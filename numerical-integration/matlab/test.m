% Applied Numerical Analysis                                  July 27, 2021
% ME 2020 FA21
% Prof. M Diaz-Maldonado
%
%
% Synopsis:
% Tests (some) numerical integration methods.
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
%    Scientists: An Introduction with Applications using MATLAB
% [1] A Gilat, MATLAB: An Introduction with Applications, 6th edition
%

clear
close all
clc

format long g


N = 255;                	% intervals
a = 0;	b = 1;			% integration limits [a, b]
f = @(x) exp(x);        	% f(x)
EI = f(b) - f(a);       	% exact integral for f(x) = exp(x)


% numerical integrations
NI_lsum = lsum(a, b, N, f);	% Left  Riemann Sum
NI_rsum = rsum(a, b, N, f);	% Right Riemann Sum
NI_trap = trap(a, b, N, f);     % Trapezoid Method


% error (percentage)
err_lsum = abs(NI_lsum - EI) / EI * 100;
err_rsum = abs(NI_rsum - EI) / EI * 100;
err_trap = abs(NI_trap - EI) / EI * 100;


fprintf("Numerical Method \t   Result \t  %% Error\n")
fprintf("---------------------------------------------------\n")
fprintf("Left  Riemann    \t %9.6f  \t %6.2e %% \n", NI_lsum, err_lsum)
fprintf("Right Riemann    \t %9.6f  \t %6.2e %% \n", NI_rsum, err_rsum)
fprintf("Trapezoid Method \t %9.6f  \t %6.2e %% \n", NI_trap, err_trap)
