#!/usr/bin/env python
# phase diagram of chaotic inflation

#######################################################################

from math import *
import numpy as np
from scipy.optimize import root

#######################################################################

# chaotic inflation potential
def V(phi): return phi**2/2.0
def DV(phi): return phi
def H2(phi,pi): return pi*pi/6.0 + V(phi)/3.0

# dynamical system to be integrated
def f(t,state):
	phi,pi,N = state; H = sqrt(H2(phi,pi))
	return np.array([pi,-3.0*H*pi-DV(phi),H])

# slow-roll terminal velocity
def u(phi): return root(lambda x: f(0,[phi,*x,0.0])[1], 0.0).x[0]

# initial conditions
state = np.array([20.0,u(20.0),0.0])

#######################################################################

from scipy.integrate import solve_ivp

# stopping condition (inflation ends)
def stop(t,state):
	phi,pi,N = state 
	return pi*pi/2.0-H2(phi,pi)

# stopping condition attributes
stop.terminal = True
stop.direction = 1.0

# ODE integration method parameters (DOP853 is fastest)
method = {'method': 'DOP853', 'atol': 1e-13, 'rtol': 1e-11}
soln = solve_ivp(f, [0,30], state, events=stop, **method)

#######################################################################

# number of grid points
nx = 128; ny = 32

# uniform evaluation grid
lx = 20.0; x = np.linspace(-lx, lx, nx)
ly = 2.00; y = np.linspace(-ly, ly, ny)

# 2D grid iterators
X,Y = np.meshgrid(x,y)
N = np.zeros([ny,nx])

# number of e-folds until the end of inflation
for j in range(ny):
	for i in range(nx):
		ic = [x[i],y[j],0.0]
		if stop(0,ic) > 0.0: continue
		N[j,i] = solve_ivp(f, [0,30], ic, events=stop, **method).y[2,-1]

#######################################################################

import matplotlib.pyplot as plt

fig = plt.figure(); ax = fig.gca()

# number of e-folds
#plt.imshow(N, extent=[-lx,lx,-ly,ly], origin='lower', cmap='YlOrBr', interpolation='none')
plt.contourf(X, Y, N, levels=256, cmap='YlOrBr')
plt.colorbar()

# sample trajectories
for i in range(-4,5):
	soln = solve_ivp(f, [0,30], [20*i/4.5,2.0,0.0], **method)
	plt.plot(soln.y[0], soln.y[1], color="tab:blue")
	plt.plot(-soln.y[0], -soln.y[1], color="tab:blue")

# attractor trajectory
soln = solve_ivp(f, [0,50], state, **method)
plt.plot(soln.y[0], soln.y[1], color="tab:red", linewidth=3)
plt.plot(-soln.y[0], -soln.y[1], color="tab:red", linewidth=3)

# clip the plot limits
plt.xlim([-lx,lx]); plt.ylim([-ly,ly])

# show in interactive console
plt.show()
