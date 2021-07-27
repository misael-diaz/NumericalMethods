#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Applied Numerical Analysis                                August 25, 2020
ME 2020 FA21
Prof. M. Diaz-Maldonado
Revised: July 27, 2021


Synopsis:
Implemements (some) numerical integration techniques.


Copyright (c) 2021 Misael Diaz-Maldonado
This file is released under the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.


References:
[0] A Gilat and V Subramanian, Numerical Methods for Engineers and
    Scientists: An Introduction with Applications using MATLAB
[1] R Johansson, Numerical Python: Scientific Computing and Data
    Science Applications with NumPy, SciPy, and Matplotlib, 2nd edition
"""


from numpy import linspace


def lsum(a, b, N, f):       # Left Riemann Sum
    """
    Synopsis:
    Integrates the function f(x) in the interval [a, b] using N intervals.
    """
    dx = (b - a) / N
    x = linspace(a, b, N + 1)
    y = f(x[0:N])
    I = dx * y.sum()
    return I


def rsum(a, b, N, f):       # Right Riemann Sum
    """
    Synopsis:
    Integrates the function f(x) in the interval [a, b] using N intervals.
    """
    dx = (b - a) / N
    x = linspace(a, b, N + 1)
    y = f(x[1:N+1])
    I = dx * y.sum()
    return I


def trap(a, b, N, f):       # Trapezoid Method
    """
    Synopsis:
    Integrates the function f(x) in the interval [a, b] using N intervals.
    """
    return ( 0.5 * ( lsum(a, b, N, f) + rsum(a, b, N, f) ) )
