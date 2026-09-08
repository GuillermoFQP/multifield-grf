#!/usr/bin/env python
# plot residual between constrained and unconstrained Gaussian random fields in 2D

#######################################################################

import numpy as np
from sys import argv
from astropy.io import fits

#######################################################################

# input files:
#   argv[1] -> constrained field   (default 'smpout-cnst.fits')
#   argv[2] -> unconstrained field (default 'smpout-uncnst.fits')
filename_cnst = 'smpout-cnst.fits' if len(argv) < 2 else argv[1]
filename_uncnst = 'smpout-uncnst.fits' if len(argv) < 3 else argv[2]

# read constrained data from disk
hdu_cnst = fits.open(filename_cnst); hdr = hdu_cnst[0].header

# number of grid points
nx = hdr['nx']; ny = hdr['ny']

# uniform evaluation grid
lx = hdr['lx']; dx = 2.0*lx/nx; x = np.linspace(-lx+dx/2, lx-dx/2, nx)
ly = hdr['ly']; dy = 2.0*ly/ny; y = np.linspace(-ly+dy/2, ly-dy/2, ny)

# 2D grid iterators
X,Y = np.meshgrid(x,y)

# constrained potential and its derivatives
field_cnst = np.stack([hdu_cnst[i].data.astype(np.float64) for i in range(1,7)], axis=0)
nodes = np.stack([hdu_cnst[7].data[i].astype(np.float64) for i in ['X','Y']], axis=0)
curve = np.stack([hdu_cnst[8].data[i].astype(np.float64) for i in ['X','Y']], axis=0)

# read unconstrained data from disk
hdu_uncnst = fits.open(filename_uncnst); hdr_uncnst = hdu_uncnst[0].header

# make sure both files describe the same grid
for key in ['nx', 'ny', 'lx', 'ly']:
    if not np.isclose(hdr[key], hdr_uncnst[key]):
        raise ValueError(f"Incompatible FITS grids: header value {key} differs between files")

# unconstrained potential and its derivatives
field_uncnst = np.stack([hdu_uncnst[i].data.astype(np.float64) for i in range(1,7)], axis=0)

# residual field: constrained minus unconstrained
field_diff = field_cnst - field_uncnst

# close FITS file handles
hdu_cnst.close(); hdu_uncnst.close()

#######################################################################

from scipy.ndimage import maximum_filter, minimum_filter

f1,f1x,f1y,f1xx,f1xy,f1yy = field_uncnst
f2,f2x,f2y,f2xx,f2xy,f2yy = field_cnst
f3,f3x,f3y,f3xx,f3xy,f3yy = field_diff

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
aspect    = 2.9
width     = 15.2
ax_left   = 0.077
ax_bottom = 0.19
ax_width  = 0.26
ax_height = 0.80
ax_gap    = 0.015

# Axes limits
phi1_min, phi1_max = -lx, lx
phi2_min, phi2_max = -ly, ly

fig = plt.figure(figsize=(width/2.54, width/aspect/2.54), frameon=False)
ax1 = fig.add_axes([ax_left, ax_bottom, ax_width, ax_height])
ax2 = fig.add_axes([ax_left + ax_width + ax_gap, ax_bottom, ax_width, ax_height])
ax3 = fig.add_axes([ax_left + 2*(ax_width + ax_gap), ax_bottom, ax_width, ax_height])

for ax in (ax1, ax2, ax3):
	ax.set_aspect('equal', adjustable='box')
	ax.set_box_aspect(1)

z = np.max(np.concatenate((np.abs(f1), np.abs(f2), np.abs(f3))))

im1 = ax1.imshow(f1, extent=[-lx,lx,-ly,ly], origin='lower', vmin=-z, vmax=z, cmap="RdBu_r", interpolation='none')
im2 = ax2.imshow(f2, extent=[-lx,lx,-ly,ly], origin='lower', vmin=-z, vmax=z, cmap="RdBu_r", interpolation='none')
im3 = ax3.imshow(f3, extent=[-lx,lx,-ly,ly], origin='lower', vmin=-z, vmax=z, cmap="RdBu_r", interpolation='none')

ax2.set_xlabel(r'$\phi^1 / M_{\mathrm{Pl}}$', fontsize=10)
ax1.set_ylabel(r'$\phi^2 / M_{\mathrm{Pl}}$', fontsize=10)
ticks = np.linspace(int(phi1_min), int(phi1_max), 7)

for ax in (ax1, ax2, ax3):
	ax.tick_params(axis='both', which='major', labelsize=8, length=2.5, width=1)
	ax.set_xticks(ticks)
	ax.set_yticks(ticks)

ax2.tick_params(axis='y', labelleft=False)
ax3.tick_params(axis='y', labelleft=False)

cbar_gap    = 0.015
cbar_width  = 0.018
cbar_height = 0.75
cbar_bottom = 0.215
cax = fig.add_axes([ax_left + 3*ax_width + 2*ax_gap + cbar_gap, cbar_bottom, cbar_width, cbar_height])
cbar = fig.colorbar(im1, cax=cax, fraction=0.08, extend="neither")
cbar.set_ticks(np.linspace(-int(z), int(z), int(z)+1))
cbar.set_label(r'$M_{\mathrm{Pl}}^{-4}$', fontsize=10)
cbar.ax.tick_params(length=2.5, labelsize=8)

ax1.text(0.95, 0.95, r'$V_{0}$', transform=ax1.transAxes, ha='right', va='top', fontsize=10)
ax2.text(0.95, 0.95, r'$V_{\rm c}$', transform=ax2.transAxes, ha='right', va='top', fontsize=10)
ax3.text(0.95, 0.95, r'$\Delta V$', transform=ax3.transAxes, ha='right', va='top', fontsize=10)

plt.savefig(f'grf-cnst.pdf')
plt.savefig(f'grf-cnst.png', dpi=500, transparent=False)

# show in interactive console
#plt.show()

#######################################################################

'''
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
'''
