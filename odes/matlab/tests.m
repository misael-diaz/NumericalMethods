% Applied Numerical Analysis                                August 06, 2021
% ME 2020 FA21
% Prof. M. Diaz-Maldonado
%
%
% Synopsis:
% Solves the first-order Ordinary Differential Equation ODE via Euler's and
% Runge-Kutta Methods:
%                             dy/dt = f(t, y) = -k * y,
% where k is a constant, t is the time, y is the state variable, 
% and f(t, y) is a function equal to the right-hand side RHS of the ODE.
%
%
% Copyright (c) 2021 Misael Diaz-Maldonado
% This file is released under the GNU General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%
% References:
% [0] A Gilat and V Subramanian, Numerical Methods for Engineers and
%     Scientists: An Introduction with Applications using MATLAB
% [1] TO ADD

clear
close all
clc

% defines the right-hand side RHS of the ODE: dy/dt = f(t, y) as a lambda
k = 1;
f = @(t, y) (-k * y);


n = 255;             % number of integration time-steps
ti = 0.0; tf = 5.0;  % initial time ti and final time tf
trange = [ti, tf];   % time range tuple
yi     = 1.0;        % initial value: yi = y(t = ti)
odesol = zeros([n + 1, 2, 2]);


% solves the ODE with Euler's and second-order Runge-Kutta Methods
odesol(:, :, 1) = Euler   (n, trange, yi, f);
odesol(:, :, 2) = EulerRK2(n, trange, yi, f);
% unpacks the numerical solutions
t = odesol(:, 1, 1); y_Euler = odesol(:, 2, 1); y_RK2 = odesol(:, 2, 2);
y = exp(-k .* t);	% analytic solution


% plots and compares the numerical solution against the analytic one
figure(1)
plot(t, y, "color", "black", "linewidth", 2, ...
     "displayname", "analytic solution");
hold on;
plot(t(1:8:end), y_Euler(1:8:end), "color", "red", "marker", "o", ...
     "linestyle", "none", "displayname", "Euler's Method")
plot(t(1:16:end), y_RK2(1:16:end), "color", "red", "marker", "s", ...
     "linestyle", "none", "displayname", "second-order Runge-Kutta Method")

grid on                  % shows grid
legend()                 % displays the legend
xlabel("time, t")
ylabel("dynamic response, y(t) = exp(-kt)")
title("Natural Response of a First Order Dynamic System")

% exports figure as a Portable Network Graphic PNG with 300 DPI resolution
print("figure.png", "-dpng", "-r300")
