#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
Applied Numerical Analysis                                    July 26, 2021
ME 2020 FA21
Prof. M Diaz-Maldonado


Synopsis:
Tests bracketing nonlinear solvers.


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


from numpy import sqrt
from numpy import log10
from nlsolvers import bisect, regfal, shifter, nls_opts

# copies nonlinear solver options into a dictionary
options = nls_opts
# overrides defaults
options['tol']      = 1.0e-12
options['max_iter'] = 256
"""
options['verbose']  = False
"""

bounds = (1.0e-2, 9.0e-2)   # bracketing interval [a, b]


""" nonlinear function f(x) """
f = lambda x: (
    1.0/ sqrt(x) + 2.0 * log10(0.024651/3.7 + 2.51/(9655526.5 * sqrt(x) ) )
)


""" solves for the enclosed root with the specified method """
(x, it) = bisect (bounds, f, opts = options)   # Bisection
(x, it) = regfal (bounds, f, opts = options)   # Regula Falsi
(x, it) = shifter(bounds, f, opts = options)   # Shifter Method
