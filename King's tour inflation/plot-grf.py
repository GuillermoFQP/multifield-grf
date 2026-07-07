#!/usr/bin/env python
# plot constrained Gaussian random field realization in 2D

#######################################################################

import numpy as np
from sys import argv
from astropy.io import fits

#######################################################################

# read data from disk
filename = 'smpout.fits' if len(argv) < 2 else argv[1]
hdu = fits.open(filename); hdr = hdu[0].header

# number of grid points
nx = hdr['nx']; ny = hdr['ny']

# uniform evaluation grid
lx = hdr['lx']; dx = 2.0*lx/nx; x = np.linspace(-lx+dx/2, lx-dx/2, nx)
ly = hdr['ly']; dy = 2.0*ly/ny; y = np.linspace(-ly+dy/2, ly-dy/2, ny)

# 2D grid iterators
X,Y = np.meshgrid(x,y)

# potential and its derivatives
field = np.stack([hdu[i].data.astype(np.float64) for i in range(1,7)], axis=0)
nodes = np.stack([hdu[7].data[i].astype(np.float64) for i in ['X','Y']], axis=0)
curve = np.stack([hdu[8].data[i].astype(np.float64) for i in ['X','Y']], axis=0)

# close FITS file handle
hdu.close()

#######################################################################

'''
from scipy.interpolate import RegularGridInterpolator

# nodes to check constraints on
n = 1024; s = np.linspace(0.0, 1.0, n)

# linear interpolation LUT for potential
lut = RegularGridInterpolator((y,x), field[0])

# potential along the constraint trajectory
v = lut(curve[[1,0]].T)
'''

#######################################################################

from scipy.ndimage import maximum_filter, minimum_filter

f,fx,fy,fxx,fxy,fyy = field

# field skeleton
S = fx*fy*(fxx-fyy) + fxy*(fy*fy-fx*fx)

# field extrema
peak = np.argwhere(f == maximum_filter(f, size=5, mode='wrap'))
trof = np.argwhere(f == minimum_filter(f, size=5, mode='wrap'))

#######################################################################

'''
import rlic

# line integral convolution visualizing gradient
texture = rng.choice([0.0,1.0], size=[ny,nx])
kernel = np.sin(np.linspace(0,np.pi,65)); kernel /= np.sum(kernel)
lic = rlic.convolve(texture, fx, fy, kernel=kernel, uv_mode='velocity', boundaries='periodic', iterations=4)

# normalize the distribution (optionally approximately equalize)
lic -= np.sum(lic)/np.prod(lic.shape)
lic /= np.sqrt(np.sum(lic**2)/np.prod(lic.shape))
#lic = np.tanh(np.sqrt(2.0/np.pi)*lic)
'''

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(48/5,27/5), frameon=False)
ax = fig.gca(); ax.patch.set_alpha(0.0)
ax.set_aspect('equal', adjustable='box')

fscale = 1.8

#plt.plot(s,v)

z = np.max(np.abs(f))
plt.imshow(f, extent=[-lx,lx,-ly,ly], origin='lower', vmin=-z, vmax=z, cmap='seismic', interpolation='none') # try RdBu_r?
#plt.imshow(lic, extent=[-lx,lx,-ly,ly], origin='lower', cmap='gray', interpolation='none')
plt.colorbar()

#plt.contour(X, Y, S, levels=[0.0], colors='white', linewidths=0.5, alpha=0.5)

plt.plot(x[peak[:,1]], y[peak[:,0]], '.', color='tab:red')
plt.plot(x[trof[:,1]], y[trof[:,0]], '.', color='tab:blue')

#plt.plot(nodes[0], nodes[1], 'o-')
plt.plot(curve[0], curve[1], '-', color='tab:orange', linewidth=3)

#ax.add_patch(plt.Rectangle((-1.0, -1.0), 2.0, 2.0, color='black', fill=False, linewidth=2, alpha=0.2))
ax.add_patch(plt.Circle((0, 0), 1.0 * 1.8, color='black', fill=False, linewidth=2, alpha=0.2))

# plot bounds and layout
#plt.xlim([-lx,lx]); plt.ylim([-ly,ly])
#plt.xlim([-1,1]); plt.ylim([-1,1])
plt.tight_layout()

# show in interactive console
plt.show()

#######################################################################

import pyvista as pv

pv.set_plot_theme("document")

Z = f/5.0; z = np.max(np.abs(Z))
grid = pv.StructuredGrid(X,Y,Z)
grid["V"] = grid.points[:,2]

plt = pv.Plotter(window_size=[2000, 2000], off_screen=False)
plt.add_mesh(grid, scalars='V', clim=[-z,z], cmap='RdBu_r', smooth_shading=True)
plt.view_xy(); plt.camera.elevation = -30; plt.camera.reset_clipping_range()
plt.remove_scalar_bar()

#plt.screenshot("high_res_plot.png", transparent_background=True, window_size=[4000,4000])

# show in interactive console
plt.show()
