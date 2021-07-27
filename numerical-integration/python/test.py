#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Applied Numerical Analysis                                    July 27, 2021
ME 2020 FA21
Prof. M. Diaz-Maldonado


Synopsis:
Tests (some) numerical integration techniques.


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

import numpy as np
from numpy import exp
from nint import lsum
from nint import rsum
from nint import trap


N = 255                         # intervals
a, b = (0, 1)                   # integration limits [a, b]
f = lambda x: exp(x)            # f(x)
EI = f(b) - f(a)                # exact integral for f(x) = exp(x)


""" numerical integrations """
NI_lsum = lsum(a, b, N, f)      # Left  Riemann Sum
NI_rsum = rsum(a, b, N, f)      # Right Riemann Sum
NI_trap = trap(a, b, N, f)      # Trapezoid Method


# error (percentage) <alternative elementwise computation>
NI = np.array([NI_lsum, NI_rsum, NI_trap])
err = np.abs(NI - EI) / EI * 100


# uses python's f-string to tabulate results
table = (
    f"Numerical Method \t   Result \t  % Error\n"
    f"----------------------------------------------------------\n"
    f"Left Riemann     \t {NI_lsum:{9}.{7}}  \t {err[0]:6.2e} % \n"
    f"Right Riemann    \t {NI_rsum:{9}.{7}}  \t {err[1]:6.2e} % \n"
    f"Trapezoid Method \t {NI_trap:{9}.{7}}  \t {err[2]:6.2e} % \n"
)


print(table)
