#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computational Methods                                    September 27, 2022
ICI 70320
Prof. M Diaz-Maldonado

Synopsis:
Plots the temperature field that satisfies the Poisson equation subject to
homogeneous boundary conditions and a constant heat source.

References:
[0] R. Johansson, Numerical Python: Scientific Computing and Data
    Science Applications with NumPy, SciPy, and Matplotlib, 2nd edition
[1] A Gilat and V Subramaniam, Numerical Methods for Engineers and
    Scientists: An Introduction with Applications using MATLAB
[2] T. Bergman, A. Lavine, F. Incropera, and D. Dewitt, Fundamentals
    of Heat and Mass Transfer, 7th edition
"""


# imports supporting methods from numerical python and matplotlib
from numpy import linspace
from numpy import meshgrid
from numpy import zeros_like
from numpy import pi, sin, cosh, sinh
import matplotlib as mpl
from matplotlib import pyplot as plt


''' function definitions '''
def fieldSteadyState(x, y):
    """
    Synopsis:
    Returns the steady-state temperature field f(x, y) that solves the
    Poisson equation with a constant heat source term.

    Inputs:
    x       (2nd-rank numpy array) x coordinate matrix
    y       (2nd-rank numpy array) y coordinate matrix

    Output:
    f       (2nd-rank numpy array) the temperature field f(x, y)
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


''' main '''
# defines the number of elements of the x, y vector coordinates
numel = 257
# defines the system boundaries
x_l, x_u = 0, 1
y_l, y_u = 0, 1
# generates the mesh x, y coordinates
x, y = meshgrid(linspace(x_l, x_u, numel), linspace(y_l, y_u, numel))
# computes the steady-state (temperature) field f(x, y)
f = fieldSteadyState(x, y)


# closes existing figures and enables interactive plotting
plt.close("all")
plt.ion()
# creates figure and axes objects for plotting
fig, ax = plt.subplots()
# defines suitable normalization parameters for the colormap
norm = mpl.colors.Normalize f.min(), f.max() )
# selects the colormap
cmap = mpl.cm.get_cmap('RdBu_r')
# generates the contour plot of the steady-state field f(x, y)
p = ax.pcolor(x, y, f, norm = norm, cmap = cmap)
# inserts a colorbar to the figure
fig.colorbar(p)
