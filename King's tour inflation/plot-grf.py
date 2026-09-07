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

# Thicker axes
plt.rc('axes', linewidth=1)

# Use LaTeX for text rendering with default LaTeX math font (Computer Modern)
plt.rc('text', usetex=True)
plt.rc('font', size=10, family='serif', serif=['Computer Modern'])

# Scaling parameter for grid
fscale = 1.8

# Plot parameters
aspect    = 1.25
width     = 7.6
ax_left   = 0.07
ax_bottom = 0.18
ax_width  = 0.80
ax_height = 0.80

# Axes limits
phi1_min, phi1_max = -lx, lx
phi2_min, phi2_max = -ly, ly

fig = plt.figure(figsize=(width/2.54, width/aspect/2.54), frameon=False)
ax = fig.add_axes([ax_left, ax_bottom, ax_width, ax_height])
ax.set_aspect('equal', adjustable='box')
ax.set_box_aspect(1)

#plt.plot(s,v)

z = np.max(np.abs(f))
im = ax.imshow(f, extent=[-lx,lx,-ly,ly], origin='lower', vmin=-z, vmax=z, cmap="RdBu_r", interpolation='none') # try RdBu_r?
#plt.imshow(lic, extent=[-lx,lx,-ly,ly], origin='lower', cmap='gray', interpolation='none')

ax.set_xlabel(r'$\varphi^1 / M_{\mathrm{Pl}}$', fontsize=10)
ax.set_ylabel(r'$\varphi^2 / M_{\mathrm{Pl}}$', fontsize=10)
ax.tick_params(axis='both', which='major', labelsize=8, length=2.5, width=1)
ticks = np.linspace(int(phi1_min), int(phi1_max), 7)
ax.set_xticks(ticks)
ax.set_yticks(ticks)

cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.08, shrink=1, extend="neither")
cbar.set_ticks(np.linspace(-4, 4, 9))
cbar.set_label(r'$V(\varphi) / 10^{-10} M_{\mathrm{Pl}}^{4}$', fontsize=10)
cbar.ax.tick_params(length=2.5, labelsize=8)

#plt.contour(X, Y, S, levels=[0.0], colors='gray', linewidths=0.5, alpha=0.5)

#plt.plot(x[peak[:,1]], y[peak[:,0]], '.', color='tab:red')
#plt.plot(x[trof[:,1]], y[trof[:,0]], '.', color='tab:blue')

#plt.plot(nodes[0], nodes[1], 'o-')
#plt.plot(curve[0], curve[1], '-', color='tab:gray', linewidth=1.0, alpha=0.25)

#ax.add_patch(plt.Rectangle((-1.0, -1.0), 2.0, 2.0, color='black', fill=False, linewidth=2, alpha=0.2))
#ax.add_patch(plt.Circle((0, 0), 1.0 * fscale, color='black', fill=False, linewidth=2, alpha=0.2))

# plot bounds and layout
#plt.xlim([-lx,lx]); plt.ylim([-ly,ly])
#plt.xlim([-1,1]); plt.ylim([-1,1])
#plt.tight_layout()

plt.savefig(f'grf-const.pdf')
plt.savefig(f'grf-cons.png', dpi=500, transparent=False)

# show in interactive console
#plt.show()

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
#plt.show()
