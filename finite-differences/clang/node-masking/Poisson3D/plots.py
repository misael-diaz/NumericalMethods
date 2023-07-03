#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transient Heat Conduction                               July 02, 2023

source: plots.py
author: @misael-diaz

Synopsis:
Visualization of the transient temperature fields that satisfy the 3D Poisson equation.

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

from numpy import arange
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.close('all')
plt.ion()
fig, ax = plt.subplots()
ids = arange(3)
colors = ['blue', 'black', 'red']
markers = ['*', 's', 'o']
# plots the transient fields f(t, x, y = 0.5, z = 0.5) with respect time and position
for i, color, marker in zip(ids, colors, markers):
    x, f = loadtxt(f'analytic{i}.txt').T
    x, g = loadtxt(f'numeric{i}.txt').T
    ax.plot(x, f, linestyle='-', color=color, label=f'analytic{i}')
    ax.plot(x[::4], g[::4], marker=marker, linestyle='', color=color, label=f'numeric{i}')

x, f = loadtxt('analytic.txt').T
# plots the analytic, steady-state, field g(t -> infinity, x, y = 0.5, z = 0.5)
ax.plot(x, f, linestyle='--', color='black', label='steady-state')
ax.set_xlabel('x')
ax.set_ylabel('g(t, x, y = 0.5, z = 0.5)')
ax.set_title('transient temperature profiles')
ax.legend()
