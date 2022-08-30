#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computational Solutions of Engineering Problems         August 30, 2022
IST 4360
Prof. M Diaz-Maldondo

Synopsis:
Tests the Gaussian elimination method.

Copyright (c) 2022 Misael Diaz-Maldonado
This file is released under the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

References:
[0] A Gilat and V Subramaniam, Numerical Methods for Engineers and
    Scientists: An Introduction with Applications using MATLAB
[1] H Wendland, Numerical Linear Algebra
[2] R Johansson, Numerical Python: Scientific Computing and Data
    Science Applications with NumPy, SciPy, and Matplotlib, 2nd edition
"""


from numpy import array
from numpy import matmul
from numpy.linalg import norm
from linsys import Gaussian as solve


''' defines the system of linear equations '''
A = array([
    [2, 1, 1, 0],
    [4, 3, 3, 1],
    [8, 7, 9, 5],
    [6, 7, 9, 8],
])

b = array([1, 1, 2, 4])


''' solves the system of linear equations by Gaussian elimination '''
x = solve(A, b)


''' post-processing '''
# prints the elements of the solution vector on the console
print('Solution:')
for i in range(x.size):
    print(f'x[{i+1}]: {x[i]:8.4f}')

# computes and displays the numeric error
err = norm( matmul(A, x) - b )
print(f'Numeric Error: {err:.15f}')
