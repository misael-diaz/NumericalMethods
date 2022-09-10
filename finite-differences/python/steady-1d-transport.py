#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computational Methods                                    September 09, 2022
ICI 70320
Prof. M Diaz-Maldondo

Synopsis:
Solves the one-dimensional transport problem with a uniform source term by
finite differences. The non-dimensional field variable (temperature) is
subjected to the following boundary conditions: g(x_l) = g(x_u) = 0, where
`x_l' and `x_u' are the lower and upper limits, respectively, on the x-axis
which delimit the system domain.

The field variable is obtained by solving the system of discretized
equations iteratively via the Jacobi method and the Gauss-Seidel method.
The code reports the numeric error and plots the field variable.

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

from numpy import ones
from numpy import zeros
from numpy import linspace
from numpy.linalg import norm
from matplotlib import pyplot as plt


""" default solver parameters """


# defines default tolerance and maximum number of iterations
tol_default = 1.0e-6
max_iter_default = int('0x0008FFFF', base=16)   # about 500K iters


""" defines tailored iterative linear solvers """


def Jacobi(N, H, x_lims, opts = (max_iter_default, tol_default) ):
    """
    Synopsis:
    Possible implementation of the Jacobi method for solving the
    one-dimensional heat transfer problem with a constant (heat) source
    in the domain x = [x_l, x_u]. The (temperature) field variable is
    subject to the following boundary conditions: g(x_l) = g(x_u) = 0.

    Inputs:
    N       (scalar) number of discretization intervals
    H       (scalar) non-dimensional (heat) source term
    x_lims  (tuple) storing the x-axis lower and upper limits (x_l, x_u)
    opts    (tuple) optional, solver parameters (iterations, tolerance)

    Output:
    g       (1st-rank numpy array) the (temperature) field variable

    Note:
    The implementation has been tailored to solve this particular problem.
    """

    # unpacks the maximum number of iterations and the tolerance
    max_iters, tol = opts

    # unpacks the x-axis lower and upper limits
    x_l, x_u = x_lims

    # computes the step-size
    dx = (x_u - x_l) / N

    # defines the right-hand side of the discretized (energy) equations
    b = -( (dx)**2 * H ) * ones(N + 1)

    # initializes the field (temperature) variable
    g = zeros(N + 1)
    """ applies the boundary conditions explicitly for the sake clarity """
    g[1], g[-1] = (0, 0)    # g(x_l) = g(x_u) = 0

    # defines the initial guess
    g0 = g.copy()

    """ solves the system of equations iteratively """
    for j in range(max_iters):

        g[1] = -0.5 * (b[1] - g0[2])

        for i in range(2, N - 1):
            g[i] = -0.5 * (b[i] - g0[i - 1] - g0[i + 1])

        g[N - 1] = -0.5 * (b[N - 1] - g0[N - 2])

        if ( norm(g - g0) < tol ):
            print(f"Jacobi(): solution found after {j+1} iterations")
            break

        g0 = g.copy()

    return g


def GaussSeidel(N, H, x_lims, opts = (max_iter_default, tol_default) ):
    """
    Synopsis:
    Possible implementation of the Gauss-Seidel method for solving the
    one-dimensional heat transfer problem with a constant (heat) source
    in the domain x = [x_l, x_u]. The (temperature) field variable is
    subject to the following boundary conditions: g(x_l) = g(x_u) = 0.

    Inputs:
    N       (scalar) number of discretization intervals
    H       (scalar) non-dimensional (heat) source term
    x_lims  (tuple) storing the x-axis lower and upper limits (x_l, x_u)
    opts    (tuple) optional, solver parameters (iterations, tolerance)

    Output:
    g       (1st-rank numpy array) the (temperature) field variable

    Note:
    The implementation has been tailored to solve this particular problem.
    """

    # unpacks the maximum number of iterations and the tolerance
    max_iters, tol = opts

    # unpacks the x-axis lower and upper limits
    x_l, x_u = x_lims

    # computes the step-size
    dx = (x_u - x_l) / N

    # defines the right-hand side of the discretized (energy) equations
    b = -( (dx)**2 * H ) * ones(N + 1)

    # initializes the field (temperature) variable
    g = zeros(N + 1)
    """ applies the boundary conditions explicitly for the sake clarity """
    g[1], g[-1] = (0, 0)    # g(x_l) = g(x_u) = 0

    # defines the initial guess
    g0 = g.copy()

    """ solves the system of equations iteratively """
    for j in range(max_iters):

        g[1] = -0.5 * (b[1] - g0[2])

        for i in range(2, N - 1):
            g[i] = -0.5 * (b[i] - g[i - 1] - g0[i + 1])

        g[N - 1] = -0.5 * (b[N - 1] - g[N - 2])

        if ( norm(g - g0) < tol ):
            print(f"Gauss-Seidel(): solution found after {j+1} iterations")
            break

        g0 = g.copy()

    return g


""" problem definition """


# defines the number of discretization intervals as a power of two
n = 8
N = 2**n

# defines the x-axis lower and upper limits (or bounds), x = [x_l, x_u]
x_lims = (x_l, x_u) = (-1, 1)
x = linspace(x_l, x_u, N + 1)

# defines the (heat) source term
H = 1


""" numeric solution """


# solves the discretized (energy) equations by applying the Jacobi method
g_Jacobi = Jacobi(N, H, x_lims)
# applies the Gauss-Seidel method to solve the discretized equations
g_GaussSeidel = GaussSeidel(N, H, x_lims)


""" post-processing """
# compares the numeric solutions against the analytic one


# computes the analytic solution
analytic = 0.5 * H * (1 - x**2)

# computes the average error
err_Jacobi      = norm(analytic - g_Jacobi) / analytic.size
err_GaussSeidel = norm(analytic - g_GaussSeidel) / analytic.size
errors = (
    f'Average Error:\n'
    f'Jacobi:       {err_Jacobi:.4e}\n'
    f'Gauss-Seidel: {err_GaussSeidel:.4e}\n'
)
print(errors)


""" visualization """


plt.close('all')
plt.ion()
fig, ax = plt.subplots()
ax.plot(x, g_Jacobi,      color='black', label='Jacobi')
ax.plot(x, g_GaussSeidel, color='black', label='Gauss-Seidel')
ax.plot(x, analytic,      color='red',   label='analytic')
ax.legend()
