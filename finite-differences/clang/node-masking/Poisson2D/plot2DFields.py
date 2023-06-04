#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transient Heat Conduction                               June 04, 2022

source: plot2DFields.py
author: @misael-diaz

Synopsis:
Visualization of the steady-state temperature field that satisfies the Poisson equation.
Plots the (temperature) field g(t -> inf, x, y = 0.5) at the middle of the domain with
respect to the `x' coordinate. The plot show that the numeric solution is in excellent
agreement with the analytic one.

Copyright (c) 2023 Misael Diaz-Maldonado
This file is released under the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

References:
[0] A Gilat and V Subramaniam, Numerical Methods for Engineers and
    Scientists: An Introduction with Applications using MATLAB
[1] R Johansson, Numerical Python: Scientific Computing and Data
    Science Applications with NumPy, SciPy, and Matplotlib, 2nd edition
"""

from numpy import loadtxt
from matplotlib import pyplot as plt

# loads the analytic and numeric solutions of the PDE
x, f = loadtxt('analytic.txt').T
x, g = loadtxt('numeric.txt').T

# plots the steady state solution g(t -> infinity, x, y = 0.5)
plt.close('all')
plt.ion()
fig, ax = plt.subplots()
ax.plot(x, f, linestyle='--', color='black', label='analytic')
ax.plot(x[::4], g[::4], marker='o', markersize=12, linestyle='',
        color='red', label='numeric')
ax.set_xlabel('x')
ax.set_ylabel('g(t -> inf, x, y = 0.5)')
ax.set_title('steady temperature profile')
ax.legend()
