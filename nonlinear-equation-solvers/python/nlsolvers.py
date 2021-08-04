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
[2] docs.python.org/3/tutorial/errors.html#tut-userexceptions
"""

from numpy import abs


# sets default values for the tolerance and maximum number of iterations
TOL, MAX_ITER = (1.0e-8, 100)

nls_opts = {'tol': TOL, 'max_iter': MAX_ITER}

def check_bounds(bounds):
    """ Synopsis: Fixes bounds if supplied in wrong order. """

    lb, ub = bounds # unpacks lower and upper bounds, respectively

    if (lb > ub):
        up = lb
        lb, ub = (ub, up)

    return (lb, ub)


def check_bracket(name, bounds, f):
    """ Synopsis: Complains if there's no root in the given interval. """
    lb, ub = bounds
    errMSG = name + ": " + f"No root exist in given interval [{lb}, {ub}]"
    if ( f(lb) * f(ub) > 0. ):
        raise RuntimeError(errMSG)
    return


def report(it, name):
    """ Synopsis: Reports if the method has been successful. """
    if (it != MAX_ITER):
        print(name + " Method:")
        print(f"solution found in {it} iterations")
    else:
        print(f"method failed to find the root you may try " + 
              f"a narrower interval")
    return

def optset(**kwargs):
    """
    Synopsis:
    Uses the tolerance and maximum number of iterations options provided
    by the user if any.
    """
    opts = kwargs.get('opts', None)

    if (opts == None):
        tol      = TOL
        max_iter = MAX_ITER
    else:
        tol = opts['tol']
        max_iter = opts['max_iter']

    return (tol, max_iter)


def bisect(bounds, f, **kwargs):
    """ Synopsis: Possible implementation of the Bisection method. """

    optional = kwargs.get('opts', None)
    (tol, max_iter) = optset(opts = optional)

    name = "Bisection"
    check_bracket(name, bounds, f)
    a, b = check_bounds(bounds)


    n, x = (1, (a + b) / 2)
    while ( n != max_iter and abs( f(x) ) > tol ):
        
        # updates bracketing interval [a, b]
        if f(a) * f(x) < 0:
            b = x
        else:
            a = x
            
        x = 0.5 * (a + b)
        n += 1


    report(n, name)
    return x


def regfal(bounds, f, **kwargs):
    """ Synopsis: Possible implementation of the Regula Falsi method. """

    optional = kwargs.get('opts', None)
    (tol, max_iter) = optset(opts = optional)

    name = "Regula Falsi"
    check_bracket(name, bounds, f)
    lb, ub = check_bounds(bounds)


    n, x = ( 1, ( lb * f(ub) - ub * f(lb) ) / ( f(ub) - f(lb) ) )
    while ( n != max_iter and abs( f(x) ) > tol ):
        
        if f(lb) * f(x) < 0:
            ub = x
        else:
            lb = x
            
        x = ( lb * f(ub) - ub * f(lb) ) / ( f(ub) - f(lb) )
        n += 1


    report(n, name)
    return x




"""
TODO:
    [x] use user-supplied values for the tolerance and maximum iterations.
    [x] raise an exception when the bracketing interval contains no root.
"""
