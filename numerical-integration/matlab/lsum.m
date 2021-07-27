% Applied Numerical Analysis                               January 10, 2019
% ME 2020 FA21
% Prof. M Diaz-Maldonado
% Revised: July 27, 2021
%
%
% Synopsis:
% Possible implementation of the Left Riemann Sum in MATLAB.
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

function I = lsum(a, b, N, f)
% Integrates the function f(x) in the interval [a, b] using N intervals.
dx = (b - a) / N;
x = linspace(a, b, N + 1);
I = dx * sum( f( x(1:N) ) );
return
