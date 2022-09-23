#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computational Methods                                    September 23, 2022
ICI 70320
Prof. M Diaz-Maldonado

Synopsis:
Visualization of the temperature field that satisfies the Poisson equation.
Plots the (temperature) field g(t, x, y = 0.5) at the middle of the domain
with respect to the `x' coordinate and different times `t'. The plots show
that the numeric solution is in excellent agreement with the analytic one.

Copyright (c) 2022 Misael Diaz-Maldonado
This file is released under the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

References:
[0] A Gilat and V Subramaniam, Numerical Methods for Engineers and
    Scientists: An Introduction with Applications using MATLAB
[1] R Johansson, Numerical Python: Scientific Computing and Data
    Science Applications with NumPy, SciPy, and Matplotlib, 2nd edition
"""

from numpy import pi
from numpy import loadtxt
from numpy import zeros_like
from numpy import cosh, sinh, exp, sin, sqrt
from matplotlib import pyplot as plt


def fieldSteadyState(x, y):
    """
    Synopsis:
    Returns the steady-state temperature field f(x, y) that solves the
    Poisson equation with a constant heat source term.

    Inputs:
    y       (scalar) y coordinate
    x       (1st-rank numpy array) x coordinates

    Output:
    f       (1st-rank numpy array) the temperature field f(x, y)
    """

    # defines the number of elements in the Fourier series
    N = 16
    # preallocates the field
    f = zeros_like(x)
    for n in range(1, N + 1):
        ln = n * pi                     # lambda(n)
        An = 2 * (1 - (-1)**n) / ln**3  # particular solution coefficient
        Bn = (cosh(ln) - 1) / sinh(ln)  # 'geometric' parameter

        # updates the computation of the Fourier series for f(x, y)
        f += An * ( 1 - cosh(ln * y) + Bn * sinh(ln * y) ) * sin(ln * x)

    return f


def fieldTransient(t, x, y):
    """
    Synopsis:
    Returns the transient temperature field f(t, x, y) that solves the
    Poisson equation with a constant heat source term.

    Inputs:
    t       (scalar) time
    y       (scalar) y coordinate
    x       (1st-rank numpy array) x coordinates

    Output:
    f       (1st-rank numpy array) the temperature field f(t, x, y)
    """

    # defines a sufficient number of elements in the Fourier series
    N = 512
    # preallocates the field
    f = zeros_like(x)
    for n in range(1, N + 1):
        for m in range(1, N + 1):
            ln, lm = n * pi, m * pi
            lnm = ln**2 + lm**2
            Anm = 2 * ( (1 - (-1)**n) / ln ) * ( (1 - (-1)**m) / lm )
            Bnm = ( 1 / lnm + (1 - 1 / lnm) * exp(-lnm * t) )
            f += 2 * Anm * Bnm * sin(ln * x) * sin(lm * y)
    return f

# loads the numeric solution of the PDE
x, g = loadtxt('steady_2d_transport_Jacobi.dat').T

# plots the steady state solution g(t, x, y = 0.5)
plt.close('all')
plt.ion()
fig, ax = plt.subplots()
y = 0.5
f = fieldSteadyState(x, y)
ax.plot(x, f, linestyle='--', color='black', label='steady-state')
ax.plot(x[::32], g[::32], linestyle='', marker='*', markersize=12,
        color='red', label='Jacobi')
ax.set_xlabel('x')
ax.set_ylabel('g(t, x, y = 0.5)')
ax.set_title('steady temperature profile')
ax.legend()


# loads the transient data g(t, x, y = 0.5) for visualization
transient = loadtxt('transient_2d_transport_Jacobi.dat')
# unpacks the time (1st-rank array) and the field g(t,x,y) (2nd-rank array)
t, g_trans = (transient[:, 0], transient[:, 1:])


# creates a figure for ploting the field g(t, x, y=0.5) at different times
fig, ax = plt.subplots()
# defines colors for the g(t, x, y = 0.5) profiles
colors = ['blue', 'black', 'orange', 'red']
for i in range(4):
    # selects profiles sufficiently far apart with respect time
    n = 32 * i
    # plots the analytic field g(t, x, y = 0.5)
    ax.plot(x, fieldTransient(t[n], x, y), color=colors[i], linestyle='-')
    # plots the numeric field g(t, x, y = 0.5)
    ax.plot(x[::32], g_trans[n, ::32], color=colors[i], linestyle='',
            marker="*", markersize=12, label=f't = {t[n]}')

# plots the steady-state solution to provide a reference
ax.plot(x, fieldSteadyState(x, y), color='black', linestyle='--',
        label='steady-state')
ax.set_xlabel('x')
ax.set_ylabel('g(t, x, y = 0.5)')
ax.set_title('transient temperature profiles')
ax.legend()
