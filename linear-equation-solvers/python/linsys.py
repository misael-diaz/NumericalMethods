#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computational Solutions of Engineering Problems         August 30, 2022
IST 4360
Prof. M Diaz-Maldondo

                        Linear Equation Solver Module

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


from numpy import zeros
from numpy import newaxis
from numpy import hstack


def backSubs(U):
    """
    Synopsis:
    Possible implementation of the back substitution technique.

    Inputs:
    U       (2nk-rank numpy array) N x (N + 1) upper triangular matrix

    Output:
    x       (1st-rank numpy array) N solution vector
    """

    # gets the number of rows and columns
    N, M = rows, cols = U.shape

    # preallocates the solution vector
    x = zeros(N)
    # applies the back substitution technique
    for i in range(N):
        k = (N - 1) - i
        x[k] = U[k, M - 1]
        for j in range(k + 1, N):
            x[k] -= (U[k, j] * x[j])
        x[k] /= U[k, k]

    return x


def Gaussian(A, b):
    """
    Synopsis:
    Possible implementation of the Gaussian elimination method.

    Inputs:
    A       (2nd-rank numpy array) N x N square matrix
    b       (1st-rank numpy array) N coefficient vector

    Output:
    x       (1st-rank numpy array) N solution vector

    Comments:
    Note that this implementation does not mitigate error propagation by
    selecting the optimal pivot. It uses as a pivoting element whatever
    value that happens to be on the diagonal and for this reason it will
    **fail** sometimes. As a student of this course you are encouraged to
    improve this implementation by using the optimal pivot on each step.
    """

    # copies matrix A into the work matrix W
    W = A.astype(float)
    # appends the coefficient vector to obtain the work augmented matrix
    W = hstack([W, b[:, newaxis]])
    # gets the number of rows and columns
    N, M = rows, cols = W.shape

    # reduces the matrix W into its upper triangular form
    for p in range(N - 1):
        for i in range(p + 1, N):
            m = -W[i, p] / W[p, p]
            for j in range(p, M):
                W[i, j] += m * W[p, j]

    # obtains the solution vector `x' via back substitution
    return backSubs(W)
