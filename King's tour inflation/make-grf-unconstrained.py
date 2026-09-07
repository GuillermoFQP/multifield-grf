#!/usr/bin/env python
# generate unconstrained Gaussian random field realization in 2D

#######################################################################

import numpy as np
from math import pi
from numpy.random import default_rng
from scipy.fft import rfft2, irfft2, rfftfreq, fftfreq

from sys import argv

# command-line arguments:
#   argv[1] -> output filename (default 'smpout.fits')
filename = 'smpout.fits' if len(argv) < 2 else argv[1]

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

# Gaussian random field realization (and its derivatives)
F = np.sqrt(P)*rfft2(rng.normal(size=[ny,nx]))
field = irfft2(np.stack([d*F for d in D]))

#######################################################################

# upscale for plotting and analysis (in Fourier domain)
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
hdr['COMMENT'] = "Unconstrained Gaussian random potential with spectrum"
hdr['COMMENT'] = f"  P(k) = exp[-k^2/{2.0*cutoff**2}]/k^{2*gamma}"
hdr['COMMENT'] = f"sampled on [{nx},{ny}] grid with |x| < {lx} and |y| < {ly}"

hdu = fits.PrimaryHDU(header=hdr); data = fits.HDUList([hdu])

# save potential and its derivatives (as 32-bit floats)
for i,name in enumerate(["V", "V,x", "V,y", "V,xx", "V,xy", "V,yy"]):
	data.append(fits.ImageHDU(data=field[i].astype(np.float32), name=name))

# write data to disk
data.writeto(filename, overwrite=True, checksum=False)
