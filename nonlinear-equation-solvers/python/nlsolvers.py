#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
Applied Numerical Analysis                                    July 17, 2019
ME 2020 FA21
Prof. M Diaz-Maldonado
Revised: July 26, 2021


Synopsis:
Possible implementation of bracketing nonlinear solvers in python.


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

from numpy import abs

# sets default values for the tolerance and maximum number of iterations
TOL, MAX_ITER = (1.0e-8, 100)


def check_bounds(bounds):
    """ Synopsis: Fixes bounds if supplied in wrong order. """

    lb, ub = bounds # unpacks lower and upper bounds, respectively

    if (lb > ub):
        up = lb
        lb, ub = (ub, up)

    return (lb, ub)


def report(it):
    """ Synopsis: Reports if the method has been successful. """
    if (it != MAX_ITER):
        print(f"solution found in {it} iterations")
    else:
        print(f"method failed to find the root you may try " + 
              f"a narrower interval")
    return


def bisect(bounds, f):
    """ Synopsis: Possible implementation of the Bisection method. """

    a, b = check_bounds(bounds)

    if f(a) * f(b) > 0:
        # complains if there are no roots in bracketing interval
        print(f"no roots exist in [{a}, {b}]")
        return 0
    

    n, x = (1, (a + b) / 2)
    while ( n != MAX_ITER and abs( f(x) ) > TOL ):
        
        # updates bracketing interval [a, b]
        if f(a) * f(x) < 0:
            b = x
        else:
            a = x
            
        x = 0.5 * (a + b)
        n += 1


    report(n)
    return x


def regfal(bounds, f):
    """ Synopsis: Possible implementation of the Regula Falsi method. """

    lb, ub = check_bounds(bounds)

    if f(lb) * f(ub) > 0:
        print(f"no roots exist in [{lb}, {ub}]")
        return 0
    

    n, x = ( 1, ( lb * f(ub) - ub * f(lb) ) / ( f(ub) - f(lb) ) )
    while ( n != MAX_ITER and abs( f(x) ) > TOL ):
        
        if f(lb) * f(x) < 0:
            ub = x
        else:
            lb = x
            
        x = ( lb * f(ub) - ub * f(lb) ) / ( f(ub) - f(lb) )
        n += 1


    report(n)
    return x




"""
TODO:
    [ ] use user-supplied values for the tolerance and maximum iterations.
    [ ] raise an exception when the bracketing interval contains no root.
"""
