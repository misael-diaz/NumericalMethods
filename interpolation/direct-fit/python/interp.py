#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computational Solutions                                      April 21, 2022
IST 4360
Prof. M Diaz Maldonado

Synopsis:
Obtains the coefficients of the 3rd degree interpolating polynomial for a
give dataset (xi, yi):

                    f(x) = a x^3 + b x^2 + c * x + d

where f(x) is the interpolating polynomial, and (a, b, c, d) are the
coefficients.


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


""" applies the direct-fit method """
# constructs the system of linear equations
A = array([
    [x1**3, x1**2, x1, 1],
    [x2**3, x2**2, x2, 1],
    [x3**3, x3**2, x3, 1],
    [x4**3, x4**2, x4, 1]
])

b = array([y1, y2, y3, y4])

# solves for the coefficients of the 3rd degree polynomial
x = solve(A, b)
# unpacks the coefficient array into four constants
a, b, c, d = x
# defines the 3rd degree polynomial as a lambda (or anonymous) function
f = lambda x: (a * x**3 + b * x**2 + c * x + d)


""" post-processing """
# prints the determinant condition number of the system of linear equations
# Note: if the system is ill-conditioned the interpolation is unreliable
print(f'Determinant:      { det(A):.2f}')   # ill-conditioned if singular
print(f'Condition number: {cond(A):.2f}')   # ill-conditioned if too large


# computes the differences between the data and the approximation returned
# by the interpolating function f(x)
diffs = (
    f'\nDifferences:\n'
    f'{y1 - f(x1):.12f}\n'
    f'{y2 - f(x2):.12f}\n'
    f'{y3 - f(x3):.12f}\n'
    f'{y4 - f(x4):.12f}\n'
)
# should display an array of zeros to indicate success
print(diffs)


# as above but accumulates the differences
xi = array([x1, x2, x3, x4])
yi = array([y1, y2, y3, y4])
diffs = (
    f'\nCumulative Difference:\n'
    f'{( yi - f(xi) ).sum():.12f}\n'
)
# should display a zero (scalar) to indicate success
print(diffs)


""" plots the dataset and the interpolating polynomial """
plt.close('all')
plt.ion()
fig, ax = plt.subplots()
# plots the x, y coordinates
ax.plot(xi, yi,  color="red",   marker="*", markersize=12, linestyle="")
# plots the interpolating polynomial
x = linspace(-2, 5, 256)
ax.plot(x, f(x), color="black", linestyle="--")
ax.set_xlabel('x')
ax.set_ylabel('f(x)')
fig.savefig("plots/interpolating-polynomial.png", dpi=300)
