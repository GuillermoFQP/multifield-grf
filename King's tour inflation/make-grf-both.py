#!/usr/bin/env python
# generate matching unconstrained and constrained Gaussian random field realizations in 2D
# using exactly the same underlying Gaussian random noise

#######################################################################

import json
import numpy as np
from math import sqrt, pi
from scipy.linalg import solve
from numpy.random import default_rng
from scipy.fft import rfft2, irfft2, rfftfreq, fftfreq
from scipy.interpolate import BSpline, CubicSpline
from scipy.integrate import cumulative_simpson as cumint
from astropy.io import fits

from sys import argv

# output filenames
filename_uncnst = 'smpout-uncnst.fits'
filename_cnst   = 'smpout-cnst.fits'

# command-line arguments:
#   argv[1] -> tour number in the trajectory catalogue, using the same
#              1-based convention as make-grf.py (default 1 -> index 0)
tour_index = 0 if len(argv) < 2 else int(argv[1]) - 1

#######################################################################

# number of grid points at the native resolution
nx = 256; ny = 256

fscale = 1.8

# uniform evaluation grid
lx = 2.0 * fscale; dx = 2.0*lx/nx; x = np.linspace(-lx+dx/2, lx-dx/2, nx)
ly = 2.0 * fscale; dy = 2.0*ly/ny; y = np.linspace(-ly+dy/2, ly-dy/2, ny)

# wavenumber grids using helper functions
kx = (2.0*pi) * rfftfreq(nx, dx)
ky = (2.0*pi) * fftfreq(ny, dy)

# 2D grid iterators
X,Y = np.meshgrid(x,y)
Kx,Ky = np.meshgrid(kx,ky)

# Laplacian operator
K2 = Kx*Kx + Ky*Ky

# partial derivatives
D = [1.0, 1.0j*Kx, 1.0j*Ky, -Kx*Kx, -Kx*Ky, -Ky*Ky]

# integrate over Fourier space
def integrate(F):
    return (2.0*np.sum(F)-np.sum(F[:,0]))/(nx*ny)

# upscale data in Fourier space
def upgrade(F, os=8):
    ny0,nn = F.shape; scaled = ny0*os,(nn-1)*os+1
    # rfftfreq convention for Nyquist is negative!
    G = np.zeros(scaled, dtype=F.dtype)
    G[:ny0//2, :nn] = F[:ny0//2 ] * os**2
    G[-ny0//2:,:nn] = F[-ny0//2:] * os**2
    return G

# return potential and its first/second derivatives, optionally upscaled
def realization_from_fourier(F, os=1):
    if os == 1:
        return irfft2(np.stack([d*F for d in D]))
    return irfft2(np.stack([upgrade(d*F, os) for d in D]))

#######################################################################

# spectrum parameters
gamma = 1.0; cutoff = 20.0 / fscale

# random generator instance
seed = None; rng = default_rng(seed)

# random field power spectrum (normalized)
with np.errstate(divide='ignore'):
    P = np.exp(-K2/(2.0*cutoff**2))/K2**gamma
P[0,0] = 0.0; P /= integrate(P)

# parameters of the derivative covariance matrix
sigma,tau,nu = np.sqrt(irfft2(np.stack([P, P*Kx**2, P*Kx**4]))[:,0,0])

#######################################################################
# Generate the Gaussian noise ONCE.
#
# F_uncnst is never modified.  The constrained realization below is made
# from a copy of this same Fourier field, guaranteeing that the two output
# potentials differ only by the deterministic constraint correction.
#######################################################################

noise = rng.normal(size=[ny,nx])
F_uncnst = np.sqrt(P)*rfft2(noise)

# Native-resolution unconstrained realization.  This exact field is used
# when evaluating the values/derivatives that enter the constraint residuals.
field_uncnst_native = realization_from_fourier(F_uncnst)

#######################################################################

# 2nd order Savitzky-Golay kernel sampling field at pixel p
def kernel(p, n=5):
    i = np.argmin((x-p[0])**2) - n//2
    j = np.argmin((y-p[1])**2) - n//2
    ddx = (X[j:j+n,i:i+n]).flatten() - p[0]
    ddy = (Y[j:j+n,i:i+n]).flatten() - p[1]
    B = np.array([np.ones(n*n),ddx,ddy,ddx*ddx,ddx*ddy,ddy*ddy])
    C = np.linalg.inv(B @ B.T) @ B
    return (i,j,C[0].reshape([n,n]))

# add a local constraint on random realization at pixel p
def constraint(p, value=None, u=[1,0], v=None, du=None, dv=None,
               duu=None, duv=None, dvv=None, n=5,
               realization=None, append=None):
    f = np.zeros([ny,nx]); i,j,w = kernel(p,n); f[j:j+n,i:i+n] = w; H = rfft2(f)

    # list of constraints to extend
    constraints = [] if append is None else append

    # directions to project on
    if v is None: v = [-u[1],u[0]]
    Ku = u[0]*Kx+u[1]*Ky; Kv = v[0]*Kx+v[1]*Ky

    # values from random realization
    f = np.zeros(6) if realization is None else [np.sum(w*realization[k,j:j+n,i:i+n]) for k in range(6)]
    fu = u[0]*f[1] + u[1]*f[2]; fv = v[0]*f[1] + v[1]*f[2]
    fuu = u[0]*u[0]*f[3] + 2.0*u[0]*u[1]*f[4] + u[1]*u[1]*f[5]
    fvv = v[0]*v[0]*f[3] + 2.0*v[0]*v[1]*f[4] + v[1]*v[1]*f[5]
    fuv = u[0]*v[0]*f[3] + (u[0]*v[1]+u[1]*v[0])*f[4] + u[1]*v[1]*f[5]

    # value and derivative constraints (preconditioned to be unitless)
    if value is not None: constraints.append((H/sigma,(value-f[0])/sigma))
    if du    is not None: constraints.append((-1.0j*Ku*H/tau,(du-fu)/tau))
    if dv    is not None: constraints.append((-1.0j*Kv*H/tau,(dv-fv)/tau))
    if duu   is not None: constraints.append((-Ku*Ku*H/nu,(duu-fuu)/nu))
    if duv   is not None: constraints.append((-Ku*Kv*H/nu,(duv-fuv)/nu))
    if dvv   is not None: constraints.append((-Kv*Kv*H/nu,(dvv-fvv)/nu))

    return constraints

# compute mean field configuration resolving all constraints
def mean_field(constraints):
    H = np.array([H for H,_ in constraints]); v = np.array([v for _,v in constraints])
    xi = np.array([[integrate(Hi.conj()*Hj*P).real for Hj in H] for Hi in H])
    return np.einsum('i,ijk', solve(xi, v), H)*P

#######################################################################

# king's tour trajectory on hexagonal board
desc = f"King's tour trajectory on hexagonal board (index {tour_index})"
with open('tour/hex.json', 'r') as catalog:
    tour = json.load(catalog)[tour_index]
nodes = np.array(tour) @ np.array([[1,0],[-0.5,sqrt(0.75)]]) * fscale / 2.5

#######################################################################

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
    tt = spline_t(s)
    px = spline_x(tt); ux = spline_dx(tt)
    py = spline_y(tt); uy = spline_dy(tt)
    unorm = np.sqrt(ux*ux+uy*uy)
    return (np.stack([px,py], axis=-1), np.stack([ux/unorm,uy/unorm], axis=-1))

#######################################################################

# nodes to impose constraints on
n = 64; s = np.linspace(0.0, 1.0, n)

# potential value and transverse mass
V = 3.0*sigma*(1.0-s); m2 = 2.0*sqrt(25.0/27.0)*nu

# constraints imposed on trajectory
p,u = curve(s); path = []

for i in range(n):
    constraint(p[i], V[i], u=u[i], du=-V[0]/l[-1], dv=0.0,
               duu=0.0, duv=0.0, dvv=m2,
               realization=field_uncnst_native, append=path)

# Apply the constraint correction to a COPY of the unconstrained Fourier field.
F_cnst = F_uncnst.copy()
F_cnst += mean_field(path)

#######################################################################

# upscale both realizations for plotting and analysis in exactly the same way
os = 4
field_uncnst = realization_from_fourier(F_uncnst, os=os)
field_cnst   = realization_from_fourier(F_cnst,   os=os)

nx_out = nx*os; ny_out = ny*os
dx_out = dx/os; dy_out = dy/os

#######################################################################

names = ["V", "V,x", "V,y", "V,xx", "V,xy", "V,yy"]

# write the unconstrained realization
def write_unconstrained(filename):
    hdr = fits.Header()
    hdr['lx'] = lx; hdr['ly'] = ly
    hdr['nx'] = nx_out; hdr['ny'] = ny_out
    if seed is not None: hdr['seed'] = seed
    if os != 1: hdr['scaled'] = os
    hdr['gamma'] = gamma; hdr['cutoff'] = cutoff
    hdr['COMMENT'] = "Unconstrained Gaussian random potential with spectrum"
    hdr['COMMENT'] = f"  P(k) = exp[-k^2/{2.0*cutoff**2}]/k^{2*gamma}"
    hdr['COMMENT'] = f"sampled on [{nx_out},{ny_out}] grid with |x| < {lx} and |y| < {ly}"
    hdr['COMMENT'] = "Shares the same underlying Gaussian realization as smpout-cnst.fits"

    hdu = fits.PrimaryHDU(header=hdr); data = fits.HDUList([hdu])
    for i,name in enumerate(names):
        data.append(fits.ImageHDU(data=field_uncnst[i].astype(np.float32), name=name))

    data.writeto(filename, overwrite=True, checksum=False)

# write the corresponding constrained realization
def write_constrained(filename):
    hdr = fits.Header()
    hdr['lx'] = lx; hdr['ly'] = ly
    hdr['nx'] = nx_out; hdr['ny'] = ny_out
    if seed is not None: hdr['seed'] = seed
    if os != 1: hdr['scaled'] = os
    hdr['gamma'] = gamma; hdr['cutoff'] = cutoff
    hdr['V0'] = V[0]; hdr['m'] = sqrt(m2)
    hdr['tourix'] = tour_index
    hdr['COMMENT'] = "Constrained Gaussian random potential with spectrum"
    hdr['COMMENT'] = f"  P(k) = exp[-k^2/{2.0*cutoff**2}]/k^{2*gamma}"
    hdr['COMMENT'] = f"sampled on [{nx_out},{ny_out}] grid with |x| < {lx} and |y| < {ly}"
    hdr['COMMENT'] = "Shares the same underlying Gaussian realization as smpout-uncnst.fits"

    hdu = fits.PrimaryHDU(header=hdr); data = fits.HDUList([hdu])

    # save potential and its derivatives (as 32-bit floats)
    for i,name in enumerate(names):
        data.append(fits.ImageHDU(data=field_cnst[i].astype(np.float32), name=name))

    # save nodes of king's tour
    cx = fits.Column('X', 'E', array=nodes[:,0])
    cy = fits.Column('Y', 'E', array=nodes[:,1])
    hdu = fits.BinTableHDU.from_columns([cx,cy])
    hdu.name = 'TOUR'; data.append(hdu)

    hdr_tour = hdu.header
    hdr_tour['COMMENT'] = desc

    # save potential valley spline
    tt = spline_t(np.linspace(0.0, 1.0, 1024))
    cx = fits.Column('X', 'E', array=spline_x(tt))
    cy = fits.Column('Y', 'E', array=spline_y(tt))
    hdu = fits.BinTableHDU.from_columns([cx,cy])
    hdu.name = 'SPLINE'; data.append(hdu)

    hdr_spline = hdu.header
    hdr_spline['degree'] = degree; hdr_spline['length'] = l[-1]
    hdr_spline['COMMENT'] = f"Potential valley ({len(path)} constraints on {n} points)"

    data.writeto(filename, overwrite=True, checksum=False)

#######################################################################

write_unconstrained(filename_uncnst)
write_constrained(filename_cnst)

print(f"Wrote {filename_uncnst}")
print(f"Wrote {filename_cnst}")
