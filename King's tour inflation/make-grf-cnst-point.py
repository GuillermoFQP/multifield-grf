#!/usr/bin/env python
# generate constrained Gaussian random field realization in 2D

#######################################################################

import numpy as np
from math import sqrt, pi
from scipy.linalg import solve
from numpy.random import default_rng
from scipy.fft import rfft2, irfft2, rfftfreq, fftfreq

from pathlib import Path
from argparse import ArgumentParser
import warnings

parser = ArgumentParser(description="Polynomial GRF valley with Fourier point constraints.")
parser.add_argument("output", help="New output FITS filename; existing files are not overwritten")
parser.add_argument("tour", type=int, help="One-based catalogue tour index")
parser.add_argument("--seed", type=int, default=None)
parser.add_argument("--points", type=int, default=64)
parser.add_argument("--base-grid", type=int, choices=[256,512,1024], default=256)
parser.add_argument("--amplitude", type=float, default=1.0, help="Initial V/sigma")
parser.add_argument("--mu", type=float, default=0.5)
parser.add_argument("--catalog", type=Path, default=Path("/Users/guillermo/Documents/Multifield inflation/random_potentials/King's tour inflation/tour/hex.json"))
args = parser.parse_args()
if args.tour < 1 or args.points < 8 or args.amplitude <= 0 or args.mu <= 0:
    parser.error("tour, amplitude, and mu must be positive; points must be at least 8")
if args.seed is not None and args.seed < 0:
    parser.error("seed must be nonnegative")
filename = args.output
tour_index = args.tour - 1
#if Path(filename).exists():
#    raise FileExistsError(f"Refusing to overwrite {filename}")

# number of grid points
nx = args.base_grid
ny = args.base_grid

# uniform evaluation grid
lx = 1.0
ly = 1.0
dx = 2.0*lx/nx
dy = 2.0*ly/ny
x  = np.linspace(-lx+dx/2, lx-dx/2, nx)
y  = np.linspace(-ly+dy/2, ly-dy/2, ny)

# wavenumber grids using helper functions
kx = (2.0*pi) * rfftfreq(nx, dx)
ky = (2.0*pi) * fftfreq(ny, dy)

# 2D grid iterators
X,Y   = np.meshgrid(x,y)
Kx,Ky = np.meshgrid(kx,ky)

# Laplacian operator
K2 = Kx*Kx + Ky*Ky

# partial derivatives
D = [1.0, 1.0j*Kx, 1.0j*Ky, -Kx*Kx, -Kx*Ky, -Ky*Ky]

# integrate over Fourier space
def integrate(F):
	return (2.0*np.sum(F)-np.sum(F[:,0])-np.sum(F[:,-1]))/(nx*ny)

# upscale data in Fourier space
def upgrade(F, os=8):
	ny,nn = F.shape; scaled = ny*os,(nn-1)*os+1
	# rfftfreq convention for Nyquist is negative!
	G = np.zeros(scaled, dtype=F.dtype)
	G[:ny//2, :nn] = F[:ny//2 ] * os**2
	G[-ny//2:,:nn] = F[-ny//2:] * os**2
	return G

#######################################################################

# spectrum parameters
gamma  = 1.0
cutoff = 30.0

# random generator instance
seed = args.seed if args.seed is not None else int(np.random.SeedSequence().generate_state(1)[0])
rng = default_rng(seed)

# random field power spectrum (normalized)
with np.errstate(divide='ignore'):
	P = np.exp(-K2/(2.0*cutoff**2))/K2**gamma
P[0,0] = 0.0
P[:, -1] = 0.0
P[ny//2, :] = 0.0
P /= integrate(P)

# parameters of the derivative covariance matrix
sigma,tau,nu = np.sqrt(irfft2(np.stack([P, P*Kx**2, P*Kx**4]))[:,0,0])

# Gaussian random field realization (and its derivatives)
F = 1.0*np.sqrt(P)*rfft2(rng.normal(size=[ny,nx]))
field = irfft2(np.stack([d*F for d in D]))

#######################################################################

# Impose exact band-limited point constraints at p; subtract the matching Fourier measurements.
def constraint(p, value=None, u=[1,0], v=None, du=None, dv=None, duu=None, duv=None, dvv=None, n=5, realization=None, append=None):
    # list of constraints to extend
    items=[] if append is None else append
    
    if v is None: v=[-u[1],u[0]]
    
    phase   = np.exp(-1j*(Kx*(p[0]-x[0])+Ky*(p[1]-y[0])))
    ku      = u[0]*Kx+u[1]*Ky; kv=v[0]*Kx+v[1]*Ky
    ops     = [phase, -1j*ku*phase, -1j*kv*phase, -ku*ku*phase, -ku*kv*phase, -kv*kv*phase]
    scales  = [sigma,tau,tau,nu,nu,nu]
    targets = [value,du,dv,duu,duv,dvv]
    for op,scale,target in zip(ops,scales,targets):
        if target is not None:
            actual=0. if realization is None else integrate(op.conj()*F).real
            items.append((op/scale,(target-actual)/scale))
    return items

# Compute mean field configuration resolving all constraints
def mean_field(items):
    global cov_cond, cnst_res
    H       = np.array([h for h, _ in items])
    values  = np.array([v for _, v in items])
    
    weights = 2 * P / (nx * ny)
    weights[:, 0] *= .5
    weights[:,-1] *= .5
    
    B  = (H*np.sqrt(weights)).reshape(len(items),-1)
    Br = np.ascontiguousarray(B.real); Bi=np.ascontiguousarray(B.imag)
    
    del B
    
    cov = Br@Br.T + Bi@Bi.T
    
    del Br,Bi
    
    eigen = np.linalg.eigvalsh(cov)
    
    if eigen[0] <= 0:
        raise RuntimeError("Constraint covariance is numerically singular; reduce point count or revise geometry")
    
    cov_cond = float(eigen[-1]/eigen[0])
    
    if cov_cond > 1e12:
        warnings.warn("Constraint covariance is ill-conditioned; inspect residuals and convergence")
    
    coeff = solve(cov,values)
    cnst_res = float(np.max(np.abs(cov@coeff-values)))
    
    print('condition', cov_cond, 'residual', cnst_res, flush=True)
    
    if not np.isfinite(cnst_res) or cnst_res > 1e-7:
        raise RuntimeError("Constraint solve residual exceeds tolerance")
    
    return np.einsum('i,ijk->jk', coeff, H) * P

#######################################################################

import json

'''
# king's tour trajectory on square board
desc = "King's tour trajectory on 4x4 square board"
with open('tour/king.json', 'r') as catalog:
	tour = json.load(catalog)[0]
nodes = (np.array(tour)-1.5) @ np.array([[1,0],[0,-1]])/2.0
'''

# king's tour trajectory on hexagonal board
desc = f"King's tour trajectory on hexagonal board (index {tour_index})"
with args.catalog.open('r') as catalog:
	tours = json.load(catalog)
	if tour_index >= len(tours): raise ValueError('Tour index exceeds catalogue length')
	tour = tours[tour_index]

a = 1.0
b = 1.0
basis = np.array([a * np.array([1.0, 0.0]), b * np.array([-0.5, np.sqrt(0.75)])])
nodes = np.asarray(tour, dtype=float) @ basis / 6.0

'''
# Hilbert curve of order 3
desc = "Hilbert curve of order 3 (16x16 board)"
tour = [[0,0], [0,1], [1,1], [1,0], [2,0], [3,0], [3,1], [2,1], [2,2], [3,2], [3,3], [2,3], [1,3], [1,2], [0,2], [0,3], [0,4], [1,4], [1,5], [0,5], [0,6],
        [0,7], [1,7], [1,6], [2,6], [2,7], [3,7], [3,6], [3,5], [2,5], [2,4], [3,4], [4,4], [5,4], [5,5], [4,5], [4,6], [4,7], [5,7], [5,6], [6,6], [6,7],
        [7,7], [7,6], [7,5], [6,5], [7,4], [7,3], [7,2], [6,2], [6,3], [5,3], [4,3], [4,2], [5,2], [5,1], [4,1], [4,0], [5,0], [6,0], [6,1], [7,1], [7,0]]
nodes = (np.array(tour)-3.5) @ np.array([[1,0],[0,-1]])/4.0
'''

#######################################################################

from scipy.interpolate import BSpline, CubicSpline
from scipy.integrate import cumulative_simpson as cumint

# B-spline representation
degree = 3; pts = len(tour)

knots = np.zeros(pts+degree+1); knots[-degree:] = 1.0
knots[degree:-degree] = np.linspace(0.0, 1.0, pts-degree+1)

spline_x = BSpline(knots, nodes[:,0], degree); spline_dx = spline_x.derivative()
spline_y = BSpline(knots, nodes[:,1], degree); spline_dy = spline_y.derivative()

# curve length parametrization
t = np.linspace(0.0, 1.0, 1024)
l = cumint(np.sqrt(spline_dx(t)**2 + spline_dy(t)**2), x=t, initial=0)

spline_t = CubicSpline(l/l[-1], t)

# position and unit tangent vector along the curve
def curve(s):
	t = spline_t(s)
	px = spline_x(t); ux = spline_dx(t)
	py = spline_y(t); uy = spline_dy(t)
	u = np.sqrt(ux*ux+uy*uy)
	return (np.stack([px,py], axis=-1), np.stack([ux/u,uy/u], axis=-1))

#######################################################################

# Locations along the existing spline
n = args.points
s = np.linspace(0.0, 1.0, n)

# Potential profile parameters
L     = l[-1]
V0    = args.amplitude * sigma  # Retain your current starting height
power = 2
mu    = args.mu                 # Smaller mu -> longer inflation at fixed L
chi   = L * (1.0 - s)
M     = mu**power
den   = chi**power + M
A     = V0 * (L**power + M) / L**power

# Potential and its derivatives with respect to physical arc length $\ell = Ls$
V   =  A * chi**power / den
Vu  = -A * power * M * chi**(power - 1) / den**2
Vuu = (A * power * M * chi**(power - 2) * ((power - 1)*M - (power + 1)*chi**power) / den**3)

# Transverse mass
m2 = 1.0 * nu # 2.0 * sqrt(25.0 / 27.0) * nu
## Smooth random transverse curvature along physical arc length.
#ell = L * s
#dell = L / (n - 1)
## Chosen smoothing length, resolved by several constraint points.
## This is a modeling choice, not the exact GRF correlation length.
#corr_length = max(3.0 * dell, 1.0 / cutoff)
#distance = ell[:, None] - ell[None, :]
#Cmass = np.exp(-0.5 * (distance / corr_length)**2)
## Correlated unit-variance Gaussian samples.
#chol = np.linalg.cholesky(Cmass + 1e-10 * np.eye(n))
#z = chol @ rng.normal(size=n)
## Positive, smoothly varying curvature.
## Pointwise ensemble mean is approximately mean_m2.
#floor_m2 = 0.05 * nu
#mean_m2 = 0.95 * nu
#log_scatter = 0.25
#m2 = floor_m2 + (mean_m2 - floor_m2) * np.exp(log_scatter * z - 0.5 * log_scatter**2)

# Constraints imposed on trajectory
p, u      = curve(s)
tc        = spline_t(s)
cx1,cy1   = spline_dx(tc),spline_dy(tc)
cx2,cy2   = spline_dx.derivative()(tc),spline_dy.derivative()(tc)
curvature = (cx1*cy2-cy1*cx2)/(cx1*cx1+cy1*cy1)**1.5
path      = []

for i in range(n):
    constraint(p[i], V[i], u=u[i], du=Vu[i], dv=0.0, duu=Vuu[i], duv=curvature[i]*Vu[i], dvv=m2, realization=field, append=path)

F += mean_field(path)

#######################################################################

# Upscale for plotting and analysis (in Fourier domain)
os = 4; nx *= os; ny *= os; dx /= os; dy /= os
field = irfft2(np.stack([upgrade(d*F,os) for d in D]))

#######################################################################

from astropy.io import fits

hdr = fits.Header()

# save parameter values for current realization
hdr['lx'] = lx; hdr['ly'] = ly
hdr['nx'] = nx; hdr['ny'] = ny
if seed is not None: hdr['seed'] = seed
if 'os' in locals() and os != 1: hdr['scaled'] = os
hdr['gamma'] = gamma; hdr['cutoff'] = cutoff
hdr['V0'] = V[0]; hdr['m'] = sqrt(m2)
hdr['tourix'] = tour_index
hdr['profile'] = 'POLY2'
hdr['mu'] = mu
hdr['nconpt'] = n
hdr['constr'] = 'FOURIER'
hdr['mixder'] = 'CURVATURE'
hdr['covcond'] = cov_cond
hdr['conres'] = cnst_res
hdr['COMMENT'] = "Constrained Gaussian random potential with spectrum"
hdr['COMMENT'] = f"  P(k) = exp[-k^2/{2.0*cutoff**2}]/k^{2*gamma}"
hdr['COMMENT'] = f"sampled on [{nx},{ny}] grid with |x| < {lx} and |y| < {ly}"

hdu = fits.PrimaryHDU(header=hdr); data = fits.HDUList([hdu])

# save potential and its derivatives (as 32-bit floats)
for i,name in enumerate(["V", "V,x", "V,y", "V,xx", "V,xy", "V,yy"]):
	data.append(fits.ImageHDU(data=field[i].astype(np.float32), name=name))

# save nodes of king's tour
cx = fits.Column('X', 'E', array=nodes[:,0])
cy = fits.Column('Y', 'E', array=nodes[:,1])
hdu = fits.BinTableHDU.from_columns([cx,cy])
hdu.name = 'TOUR'; data.append(hdu)

hdr = hdu.header
hdr['COMMENT'] = desc

# save potential valley spline
t = spline_t(np.linspace(0.0, 1.0, 1024))
cx = fits.Column('X', 'E', array=spline_x(t))
cy = fits.Column('Y', 'E', array=spline_y(t))
hdu = fits.BinTableHDU.from_columns([cx,cy])
hdu.name = 'SPLINE'; data.append(hdu)

hdr = hdu.header
hdr['degree'] = degree; hdr['length'] = l[-1]
hdr['COMMENT'] = f"Potential valley ({len(path)} constraints on {n} points)"

# write data to disk
data.writeto(filename, overwrite=True, checksum=False)
