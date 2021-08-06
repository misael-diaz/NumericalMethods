#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Applied Numerical Analysis                                November 07, 2020
ME 2020
Prof. M. Diaz-Maldonado
Revised: August 06, 2021


Synopsis:
Implements (some) numerical, Ordinary Differential Equation ODE solvers.


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


from numpy import array
from numpy import empty
from numpy import linspace


def Euler(N, trange, yi, f):
        """
        Synopsis:
        Solves a first-order ordinary differential equation with Euler's
        explicit method.

        inputs:
        N       number of time-steps for numerical integration
        trange  t[ime] range, the initial, ti, and final times, tf
        yi      initial value of the state variable, yi = y(t = ti)
        f       f(t, y), RHS of the initial-value problem: dy/dt = f(t, y)

        output:
        odesol  second-rank numpy array, [t, y]
        """

        y      = empty(N + 1)    # preallocates array for speed
        ti, tf = trange
        t,  dt = (linspace(ti, tf, N+1), (tf - ti) / N)

        y[0] = yi
        for i in range(N):
                y[i + 1] = y[i] + dt * f(t[i], y[i])

        odesol = array([t, y])   # packs into a second-rank numpy array
        return odesol


def EulerRK2(N, trange, yi, f):
        """
        Synopsis:
        Solves a first-order ordinary differential equation with
        a second-order Runge-Kutta method based on Euler's method.
        
        References:
        Equations (10.63) and (10.64) p. 406 
        A Gilat and V Subramaniam, Numerical Methods for Engineers and
        Scientists, 3rd edition
        """

        y      = empty(N + 1)
        ti, tf = trange
        t,  dt = (linspace(ti, tf, N+1), (tf - ti) / N)

        y[0] = yi
        for i in range(N):
                K1 = f(t[i], y[i])
                K2 = f(t[i] + dt, y[i] + K1 * dt)
                y[i+1] = y[i] + 0.5 * dt * (K1 + K2)

        odesol = array([t, y])
        return odesol
