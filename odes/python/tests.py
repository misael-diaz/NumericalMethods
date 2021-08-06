#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Applied Numerical Analysis                                  August 06, 2021
ME 2020 FA21
Prof. M. Diaz-Maldonado


Synopsis:
Solves the first-order Ordinary Differential Equation ODE via Euler's and
Runge-Kutta Methods:
                            dy/dt = f(t, y) = -k * y,
where k is a constant, t is the time, y is the state variable, and f(t, y)
is a function equal to the right-hand side RHS of the ODE.


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


from ode import Euler
from ode import EulerRK2 as RK2
import numpy as np
import matplotlib as mpl
mpl.use("qt5agg")
from matplotlib import pyplot as plt


# defines the right-hand side RHS of the ODE: dy/dt = f(t, y) as a lambda
k = 1.0
f = lambda t, y: (-k * y)


n = 255                         # number of integration time-steps
ti, tf = (0., 5.)               # initial time ti and final time tf
trange = (ti, tf)               # time range tuple
yi     = 1.0                    # initial value: yi = y(t = ti)
odesol = np.empty([2, n+1, 2])  # preallocates array for speed


""" solves the ODE with Euler's and second-order Runge-Kutta Methods """
odesol[:, :, 0] = Euler(n, trange, yi, f)
odesol[:, :, 1] = RK2  (n, trange, yi, f)
# unpacks the numerical solutions
t, y_Euler, y_RK2 = (odesol[0, :, 0], odesol[1, :, 0], odesol[1, :, 1])
y = np.exp(-k * t)  # analytic solution


"""plots and compares the numerical solution against the analytic one"""
plt.close("all")            # closes all existing figures
plt.ion()                   # enables interactive plotting
fig, ax = plt.subplots()    # effectively plots on the same figure

ax.plot(t, y, color="black", linewidth=2.0, label="analytic solution")
ax.plot(t[::8], y_Euler[::8], color="orange", marker="o", linestyle="",
        label="Euler's Method")
ax.plot(t[::16], y_RK2[::16], color="red",    marker="s", linestyle="",
        label="second-order Runge-Kutta Method")

ax.grid()                   # shows grid
ax.legend()                 # displays the legend
ax.set_xlabel("time, t")
ax.set_ylabel("dynamic response, y(t) = exp(-kt)")
ax.set_title("Natural Response of a First Order Dynamic System")

# exports figure as a Portable Network Graphic PNG with 300 DPI resolution
fig.savefig("figure.png", dpi=300)
