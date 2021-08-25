#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Applied Numerical Analysis                                  August 25, 2021
ME 2020 FA21
Prof. M. Diaz-Maldonado


Synopsis:
Solves for the transient response of a first-order system when subjected
to an impulse-input.


Copyright (c) 2021 Misael Diaz-Maldonado
This file is released under the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.


References:
[0] A Gilat and V Subramanian, Numerical Methods for Engineers and
    Scientists: An Introduction with Applications using MATLAB
[1] R Johansson, Numerical Python: Scientific Computing and Data
    Science Applications with NumPy, SciPy, and Matplotlib, 2nd edition
[2] CA Kluever, Dynamic Systems: Modeling, Simulation, and Control
"""


from ode import iEuler
from ode import EulerRK2 as RK2
import numpy as np
from numpy import exp
from scipy.integrate import solve_ivp
import matplotlib as mpl
mpl.use("qt5agg")
from matplotlib import pyplot as plt


# initial value, and rate and forcing constants, respectively
yi, k, b = (0.0, 1.0, 1.0)
# defines the right-hand side RHS of the ODE: dy/dt = f(t, y) as a lambda
odefun = lambda t, y: (b - k * y)
# defines lambdas for the step and impulse-response for validation
step    = lambda t: (yi - b / k) * exp(-k * t) + b / k 
impulse = lambda t: -k * (yi - b / k) * exp(-k * t)


n = 255                         # number of integration time-steps
ti, tf = (0.0, 5.0)             # initial time ti and final time tf
tspan  = (ti, tf)               # time span tuple
odesol = np.empty([2, n+1, 2])  # preallocates array for speed


""" solves the ODE with Euler's and second-order Runge-Kutta Methods """
odesol[:, :, 0] = iEuler(n, tspan, yi, odefun)
odesol[:, :, 1] = RK2   (n, tspan, yi, odefun)
# uses scipy's implementation of the fourth-order Runge-Kutta method
odesol_scipy_IVP_Solver = solve_ivp(odefun, tspan, [yi], method="RK45")
# unpacks the numerical solutions
t, y_Euler, y_RK2 = (odesol[0, :, 0], odesol[1, :, 0], odesol[1, :, 1])
t_scipy, y_scipy  = (odesol_scipy_IVP_Solver.t, odesol_scipy_IVP_Solver.y)
y_scipy = y_scipy[0, :]

"""plots and compares the numerical solutions against the analytic ones"""

# step response
plt.close("all")            # closes all existing figures
plt.ion()                   # enables interactive plotting
fig, ax = plt.subplots()    # effectively plots on the same figure

ax.plot(t, step(t), color="black", linewidth=2.0,
        label="analytic solution")
ax.plot(t_scipy, y_scipy, color="blue", marker="v", linestyle="",
        label="scipy's IVP Solver")
ax.plot(t[::8], y_Euler[::8], color="orange", marker="o", linestyle="",
        label="implicit Euler's Method")
ax.plot(t[::16], y_RK2[::16], color="red",    marker="s", linestyle="",
        label="second-order Runge-Kutta Method")

ax.grid()                   # shows grid
ax.legend()                 # displays the legend
ax.set_xlabel("time, t")
ax.set_ylabel("dynamic response, y(t)")
ax.set_title("Step Response of a First Order Dynamic System")

# exports figure as a Portable Network Graphic PNG with 300 DPI resolution
fig.savefig("output/impulse/step.png", dpi=300)



# impulse response
fig, ax = plt.subplots()

ax.plot(t, impulse(t), color="black", linewidth=2.0,
        label="analytic solution")
ax.plot(t_scipy, odefun(t_scipy, y_scipy), color="blue", marker="v",
        linestyle="", label="scipy's IVP Solver")
ax.plot(t[::8], odefun(t[::8], y_Euler[::8]), color="orange", marker="o",
        linestyle="", label="implicit Euler's Method")
ax.plot(t[::16], odefun(t[::8], y_RK2[::16]), color="red",    marker="s",
        linestyle="", label="second-order Runge-Kutta Method")

ax.grid()
ax.legend()
ax.set_xlabel("time, t")
ax.set_ylabel("dynamic response, y(t)")
ax.set_title("Impulse Response of a First Order Dynamic System")
fig.savefig("output/impulse/impulse.png", dpi=300)
