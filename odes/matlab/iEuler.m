% Applied Numerical Analysis                                August 10, 2021
% ME 2020
% Prof. M. Diaz-Maldonado
%
% source: iEuler.m
%
% Synopsis:
% Implements Euler's implicit method for the numerical solution of 
% Ordinary Differential Equations ODEs of first order.
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
% [1] A Gilat, MATLAB: An Introduction with Applications, 6th edition



function odesol = iEuler(N, trange, yi, f) % implicit Euler's method

    y  = zeros(N + 1, 1);
    ti = trange(1);	tf = trange(2);
    dt = (tf - ti) / N;
    t  = linspace(ti, tf, N + 1)';

    y(1) = yi;
    for i = 1:N
	K1 = f( t(i), y(i) );
	K2 = f( t(i) + dt, y(i) + dt * K1 );
	y_lb = y(i) + dt * K1;
	y_ub = y(i) + dt * K2;
	objf = @(yn) yn - y(i) - dt * f( t(i+1), yn );
        y(i + 1) = shifter (y_lb, y_ub, objf);
    end

    odesol = [t, y];   % packs solution into a second-rank array
end
