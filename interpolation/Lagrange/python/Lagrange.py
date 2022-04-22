#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computational Solutions                                      April 21, 2022
IST 4360
Prof. M Diaz Maldonado

Synopsis:
Uses a Lagrange 3rd-degree polynomial to interpolate a dataset (xi, yi).


Copyright (c) 2022 Misael Diaz-Maldonado
This file is released under the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.


References:
[0] A Gilat and V Subramanian, Numerical Methods for Engineers and
    Scientists: An Introduction with Applications using MATLAB
[1] R Johansson, Numerical Python: Scientific Computing and Data
    Science Applications with NumPy, SciPy, and Matplotlib, 2nd edition
"""


from numpy import abs
from numpy import array
from numpy import linspace
from numpy.linalg import det
from numpy.linalg import cond
from numpy.linalg import solve
import matplotlib as mpl
mpl.use('qt5agg')
from matplotlib import pyplot as plt


""" Dataset """
# xi, yi coordinates
x1, y1 = -1.2, 18.8
x2, y2 =  0.2,  5.0
x3, y3 =  2.0, 16.0
x4, y4 =  3.5, 15.0


""" Lagrange 3rd-degree interpolating polynomial """
# defines the 3rd degree polynomial as a lambda (or anonymous) function
f = lambda x: (
    y1*(x - x2)*(x - x3)*(x - x4) / ( (x1 - x2)*(x1 - x3)*(x1 - x4) ) +
    y2*(x - x1)*(x - x3)*(x - x4) / ( (x2 - x1)*(x2 - x3)*(x2 - x4) ) +
    y3*(x - x1)*(x - x2)*(x - x4) / ( (x3 - x1)*(x3 - x2)*(x3 - x4) ) +
    y4*(x - x1)*(x - x2)*(x - x3) / ( (x4 - x1)*(x4 - x2)*(x4 - x3) )
)


""" post-processing """
# computes the differences between the data and the approximation returned
# by the interpolating function f(x)
diffs = (
    f'\nDifferences:\n'
    f'{abs( y1 - f(x1) ):.12f}\n'
    f'{abs( y2 - f(x2) ):.12f}\n'
    f'{abs( y3 - f(x3) ):.12f}\n'
    f'{abs( y4 - f(x4) ):.12f}\n'
)
# should display an array of zeros to indicate success
print(diffs)


# as above but accumulates the differences
xi = array([x1, x2, x3, x4])
yi = array([y1, y2, y3, y4])
diffs = (
    f'\nCumulative Difference:\n'
    f'{abs( yi - f(xi) ).sum():.12f}\n'
)
# should display a zero (scalar) to indicate success
print(diffs)


""" plots the dataset and the interpolating polynomial """
plt.close('all')
plt.ion()
fig, ax = plt.subplots()
# plots the x, y coordinates
ax.plot(xi, yi,  color="red",   marker="*", markersize=12, linestyle="")
ax.set_xlim([ -2,  5])
ax.set_ylim([-30, 50])
# plots the interpolating polynomial
x = linspace(-2, 5, 256)
ax.plot(x, f(x), color="black", linestyle="--")
ax.set_xlabel('x')
ax.set_ylabel('f(x)')
