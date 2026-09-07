#!/usr/bin/env python
# compute number of e-folds inflating in random potential

#######################################################################

from math import *
import numpy as np

#######################################################################

from astropy.io import fits
from scipy.interpolate import RegularGridInterpolator

# scaling parameters
scale = 2.0; stretch = 2.0

# read tabulated potential data
hdu = fits.open('../smpout-cnst.fits')
hdr = hdu[0].header

# number of grid points
nx = hdr['nx']; ny = hdr['ny']

# uniform evaluation grid
lx = hdr['lx']*stretch; dx = 2.0*lx/nx; x = np.linspace(-lx+dx/2, lx-dx/2, nx)
ly = hdr['ly']*stretch; dy = 2.0*ly/ny; y = np.linspace(-ly+dy/2, ly-dy/2, ny)

# potential and its derivatives
field = np.stack([hdu[i].data.astype(np.float64) for i in range(1,4)], axis=-1)
field *= scale; field[:,:,1:] /= stretch
lut = RegularGridInterpolator((y,x), field, method='cubic')

# initial conditions
start = np.array(hdu[7].data[0])*stretch

# close FITS file handle
hdu.close()

#######################################################################

# dynamical system to be integrated
def f(t,state):
	x,y,v,w = state; V,Vx,Vy = lut((y,x)); H = sqrt((v*v+w*w)/6.0 + V/3.0)
	return np.array([v,w,-3.0*H*v-Vx,-3.0*H*w-Vy])/H

# ICs with slow-roll terminal velocity
def terminal(x,y,field=None):
	V,Vx,Vy = lut((y,x)) if field is None else field
	H0 = np.sqrt(V + sqrt(V*V + (2.0/3.0)*(Vx*Vx+Vy*Vy)))/sqrt(6.0)
	return np.array([x,y,-Vx/H0/3.0,-Vy/H0/3.0])

#######################################################################

from scipy.integrate import solve_ivp

# stopping condition (inflation ends)
def stop(t,state):
	x,y,v,w = state; V = lut((y,x))[0]
	return v*v + w*w - V

# stopping condition attributes
stop.terminal = True
stop.direction = 1.0

# ODE integration method parameters (DOP853 is fastest)
method = {'method': 'LSODA', 'atol': 1e-6, 'rtol': 1e-6}

#######################################################################


# number of grid points
nx = 256//2; ny = 256//2; nmax = 100.0

# uniform evaluation grid
lx /= 2.0; dx = 2.0*lx/nx; x = np.linspace(-lx+dx/2, lx-dx/2, nx)
ly /= 2.0; dy = 2.0*ly/ny; y = np.linspace(-ly+dy/2, ly-dy/2, ny)

# negative value flags incomplete evolution
N = -np.ones([ny,nx])

# number of e-folds until the end of inflation
for j in range(ny):
	for i in range(nx):
		field = lut((y[j],x[i]))
		if field[0] <= 0.0: continue
		ic = terminal(x[i],y[j],field)
		if stop(0,ic) > 0.0: continue
		print(ic)
		soln = solve_ivp(f, [0,nmax], ic, events=stop, **method)
		N[j,i] = soln.t[-1] if soln.status > 0 else np.inf


#######################################################################

import matplotlib.pyplot as plt
import matplotlib.colors as clr
import matplotlib.cm as cm

fig = plt.figure(); ax = fig.gca()

# field potential
#plt.imshow(field[:,:,0], extent=[-lx,lx,-ly,ly], origin='lower', cmap='seismic', norm=clr.CenteredNorm(), interpolation='none')

# number of e-folds
cmap = cm.YlOrBr; cmap.set_under('white'); cmap.set_bad('black')
plt.imshow(N, extent=[-lx,lx,-ly,ly], origin='lower', cmap=cmap, vmin=0.0, interpolation='none')
#plt.contourf(X, Y, N, levels=256, cmap='YlOrBr')
plt.colorbar()

# designer trajectory
soln = solve_ivp(f, [0,100], terminal(*start), events=stop, **method)
plt.plot(soln.y[0], soln.y[1], color="tab:blue", linewidth=1)

# clip the plot limits
plt.xlim([-lx,lx]); plt.ylim([-ly,ly])

# show in interactive console
plt.show()
