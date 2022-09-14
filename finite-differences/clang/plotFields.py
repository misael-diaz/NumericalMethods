#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computational Methods                                    September 14, 2022
ICI 70320
Prof. M Diaz-Maldonado

Synopsis:
Plots the (temperature) field f(t, x) with respect to position `x' and
constant time `t'.

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

from numpy import loadtxt
from matplotlib import pyplot as plt

# loads the PDE solutions
_, _, g_Jacobi      = loadtxt('pdesolJacobi.dat').T
x, f, g_GaussSeidel = loadtxt('pdesolGauss-Seidel.dat').T

plt.close('all')
plt.ion()
fig, ax = plt.subplots()
# plots the (temperature) field f(t, x) with respect to position
ax.plot(x, f, color='black', label='analytic')
ax.plot(x[::16], g_Jacobi[::16],      linestyle='', marker='s',
        markersize=12, color='blue', label='Jacobi')
ax.plot(x[::32], g_GaussSeidel[::32], linestyle='', marker='*',
        markersize=12, color='red',  label='Gauss-Seidel')
ax.set_xlabel('x')
ax.set_ylabel('f(t, x)')
ax.legend()
