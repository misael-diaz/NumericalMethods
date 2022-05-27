#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Applied Numerical Analysis                                  August 25, 2021
ME 2020 FA21
Prof. M. Diaz-Maldonado


Synopsis:
Solves for the transient response of a nonlinear Ordinary Differential
Equation ODE that models the multimode heat transfer from a slab. The
steady-state temperature of the slab is found by applying a nonlinear
equation solver.


Copyright (c) 2021 Misael Diaz-Maldonado
This file is released under the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.


References:
[0] A Gilat and V Subramanian, Numerical Methods for Engineers and
    Scientists: An Introduction with Applications using MATLAB
[1] R Johansson, Numerical Python: Scientific Computing and Data
    Science Applications with NumPy, SciPy, and Matplotlib, 2nd edition
[2] TL Bergman, AS Lavine, FP Incropera, DP DeWitt, Fundamentals of Heat
    and Mass Transfer, 8th edition.
"""


from ode import iEuler
from ode import EulerRK2 as RK2
from nlsolvers import shifter as fzero
import numpy as np
from numpy import exp
from scipy.integrate import solve_ivp
import matplotlib as mpl
mpl.use("qt5agg")
from matplotlib import pyplot as plt


"""
defines problem and solver parameters

Physical Quantity                       Units
---------------------------------------------------
convection heat transfer coefficient    W / (m^2 K)
surroundings temperature                K
slab density                            kg / m^3
slab heat capacity                      J / (kg K)
slab thickness                          m
slab initial temperature                K
irradiation                             W / m^2
Stefan-Boltzmann constant               W / (m^2 K^4)
absorptivity
emissivity
"""
h = 20.0
T_sur = 293.15
rho, cp, tau, Ti = (2700.0, 900.0, 4.0e-3, 298.15)
G = 900.0
sigma = 5.67e-8
alpha, eps = (0.80, 0.25)

# initial value, and rate and forcing constants, respectively
yi = Ti
# defines the right-hand side RHS of the ODE: dy/dt = f(t, y) as a lambda
odefun = lambda t, T: (
    ( alpha * G - h * (T - T_sur) - eps * sigma * (T**4 - T_sur**4) ) /
    ( rho * cp * tau )
)


n = 255                         # number of integration time-steps
ti, tf = (0.0, 4.0e3)           # initial time ti and final time tf, sec
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

"""plots the numerical solutions for comparison"""
plt.close("all")
plt.ion()
fig, ax = plt.subplots()

ax.plot(t_scipy, y_scipy, color="blue", marker="v", linestyle="",
        label="scipy's IVP Solver")
ax.plot(t[::8], y_Euler[::8], color="orange", marker="o", linestyle="",
        label="implicit Euler's Method")
ax.plot(t[::16], y_RK2[::16], color="red",    marker="s", linestyle="",
        label="second-order Runge-Kutta Method")

ax.grid()
ax.legend()
ax.set_xlabel("time, t")
ax.set_ylabel("dynamic response, temperature, T(t)")
ax.set_title("Multimode Heat Transfer")

# exports figure as a Portable Network Graphic PNG with 300 DPI resolution
fig.savefig("output/heat-transfer/multimode.png", dpi=300)


""" solves for the steady-state temperature of the slab """
objf = lambda T: odefun(.0, T)  # at steady-state time t is irrelevant
ss = fzero( (320.0, 330.0), objf )
T_ss = ss[0]


"""
Final Remarks:
The slab temperature rises until a state of thermal equilibrium is reached:
the rate of energy dissipation from the slab into the surroundings matches
the input rate from irradiation. Note that the slab approaches the
equilibrium state asymptotically so that in theory an infinite amount of
time is needed. The sought steady-state temperature has been found by
supplying `odefun' to a nonlinear equation solver.
"""
