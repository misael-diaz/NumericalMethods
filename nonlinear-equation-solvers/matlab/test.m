% Applied Numerical Analysis                                 July 26, 2021
% ME 2020 FA21
% Prof. M Diaz-Maldonado
%
% source: test.m
%
% Synopsis:
% Tests bracketing methods.
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


a = 1.0e-2;	b = 9.0e-2;	% bracketing interval [a, b]


% defines the nonlinear function f(x) as a lambda function
f = @(x) 1.0 / sqrt(x) + 2.0 * log10(0.024651/3.7 + ...
         2.51/(9655526.5 * sqrt(x) ) );


% solves for the root of f(x) numerically
%[ configuration struct (optional) ]%
opt = struct('tol', 1.0e-12, 'max_iter', 256, 'verbose', 0);
x = bisect (a, b, f, opt)	% Bisection
x = regfal (a, b, f, opt)	% Regula Falsi (or False Position)
x = shifter(a, b, f, opt)	% Shifter Method


% tests throwing exceptions (interval [a, b] does not contain a root)
% a = 6.0e-2;	b = 9.0e-2;
% x = bisect(a, b, f)	passed
% x = regfal(a, b, f)	passed
% x = shifter(a, b, f)	passed
