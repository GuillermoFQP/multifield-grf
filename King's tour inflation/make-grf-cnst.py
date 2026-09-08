#!/usr/bin/env python
# generate constrained Gaussian random field realization in 2D

#######################################################################

import numpy as np
from math import sqrt, pi
from scipy.linalg import solve
from numpy.random import default_rng
from scipy.fft import rfft2, irfft2, rfftfreq, fftfreq

from sys import argv

# command-line arguments:
#   argv[1] -> output filename (default 'smpout.fits')
#   argv[2] -> tour index i into the trajectory catalogue (default 0)
filename  = 'smpout.fits' if len(argv) < 2 else argv[1]
tour_index = 0            if len(argv) < 3 else int(argv[2])-1

#######################################################################

# number of grid points
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
	ny,nn = F.shape; scaled = ny*os,(nn-1)*os+1
	# rfftfreq convention for Nyquist is negative!
	G = np.zeros(scaled, dtype=F.dtype)
	G[:ny//2, :nn] = F[:ny//2 ] * os**2
	G[-ny//2:,:nn] = F[-ny//2:] * os**2
	return G

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
#sigma,tau,nu = np.sqrt(irfft2([P, P*Kx**2, P*Kx**4])[:,0,0])
sigma,tau,nu = np.sqrt(irfft2(np.stack([P, P*Kx**2, P*Kx**4]))[:,0,0])

# Gaussian random field realization (and its derivatives)
F = np.sqrt(P)*rfft2(rng.normal(size=[ny,nx]))
#field = irfft2([d*F for d in D])
field = irfft2(np.stack([d*F for d in D]))

#######################################################################

# 2nd order Savitzky-Golay kernel sampling field at pixel p
def kernel(p, n=5):
	i = np.argmin((x-p[0])**2) - n//2
	j = np.argmin((y-p[1])**2) - n//2
	dx = (X[j:j+n,i:i+n]).flatten() - p[0]
	dy = (Y[j:j+n,i:i+n]).flatten() - p[1]
	B = np.array([np.ones(n*n),dx,dy,dx*dx,dx*dy,dy*dy])
	C = np.linalg.inv(B @ B.T) @ B
	return (i,j,C[0].reshape([n,n]))

# add a local constraint on random realization at pixel p
def constraint(p, value=None, u=[1,0], v=None, du=None, dv=None, duu=None, duv=None, dvv=None, n=5, realization=None, append=None):
	f = np.zeros([ny,nx]); i,j,w = kernel(p,n); f[j:j+n,i:i+n] = w; F = rfft2(f)
	
	# list of constraints to extend
	list = [] if append is None else append
	
	# directions to project on
	if v is None: v = [-u[1],u[0]]
	Ku = u[0]*Kx+u[1]*Ky; Kv = v[0]*Kx+v[1]*Ky
	
	# values from random realization (nearest neighbour)
	f = np.zeros(6) if realization is None else [np.sum(w*realization[k,j:j+n,i:i+n]) for k in range(6)]
	fu = u[0]*f[1] + u[1]*f[2]; fv = v[0]*f[1] + v[1]*f[2]
	fuu = u[0]*u[0]*f[3] + 2.0*u[0]*u[1]*f[4] + u[1]*u[1]*f[5]
	fvv = v[0]*v[0]*f[3] + 2.0*v[0]*v[1]*f[4] + v[1]*v[1]*f[5]
	fuv = u[0]*v[0]*f[3] + (u[0]*v[1]+u[1]*v[0])*f[4] + u[1]*v[1]*f[5]
	
	# value and derivative constraints (preconditioned to be unitless)
	if not (value is None): list.append((F/sigma,(value-f[0])/sigma))
	if not (du is None): list.append((-1.0j*Ku*F/tau,(du-fu)/tau))
	if not (dv is None): list.append((-1.0j*Kv*F/tau,(dv-fv)/tau))
	if not (duu is None): list.append((-Ku*Ku*F/nu,(duu-fuu)/nu))
	if not (duv is None): list.append((-Ku*Kv*F/nu,(duv-fuv)/nu))
	if not (dvv is None): list.append((-Kv*Kv*F/nu,(dvv-fvv)/nu))
	
	return list

# compute mean field configuration resolving all constraints
def mean_field(list):
	H = np.array([H for H,_ in list]); v = np.array([v for _,v in list])
	xi = np.array([[integrate(Hi.conj()*Hj*P).real for Hj in H] for Hi in H])
	return np.einsum('i,ijk', solve(xi, v), H)*P

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
with open('tour/hex.json', 'r') as catalog:
	tour = json.load(catalog)[tour_index]
nodes = np.array(tour) @ np.array([[1,0],[-0.5,sqrt(0.75)]]) * fscale / 2.5

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

# nodes to impose constraints on
n = 64; s = np.linspace(0.0, 1.0, n)

# potential value and transverse mass
V = 3.0*sigma*(1.0-s); m2 = 2.0*sqrt(25.0/27.0)*nu

# constraints imposed on trajectory
p,u = curve(s); path = []

for i in range(n):
	constraint(p[i], V[i], u=u[i], du=-V[0]/l[-1], dv=0.0, duu=0.0, duv=0.0, dvv=m2, realization=field, append=path)

# constrained random field realization
F += mean_field(path); #field = irfft2([d*F for d in D])

#######################################################################

# upscale for plotting and analysis (in Fourier domain)
os = 4; nx *= os; ny *= os; dx /= os; dy /= os
#field = irfft2([upgrade(d*F,os) for d in D])
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
