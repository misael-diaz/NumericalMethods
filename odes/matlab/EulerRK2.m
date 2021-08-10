% Applied Numerical Analysis                              November 07, 2020
% ME 2020
% Prof. M. Diaz-Maldonado
% Revised: August 06, 2021
%
%
% Synopsis:
% Implements Euler's variant of a second-order Runge-Kutta method for the
% numerical solution of first-order Ordinary Differential Equations ODEs.
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



function odesol = EulerRK2(N, trange, yi, f)
%   Synopsis: Runge-Kutta Method of second order.

    y  = zeros(N + 1, 1);
    ti = trange(1);	tf = trange(2);
    dt = (tf - ti) / N;
    t  = linspace(ti, tf, N + 1)';

    y(1) = yi;
    for i = 1:N
	K1 = f( t(i), y(i) );
	K2 = f( t(i) + dt, y(i) + K1 * dt);
        y(i + 1) = y(i) + 0.5 * dt * (K1 + K2);
    end

    odesol = [t, y];
end
