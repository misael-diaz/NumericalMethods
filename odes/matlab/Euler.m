% Applied Numerical Analysis                              November 07, 2020
% ME 2020
% Prof. M. Diaz-Maldonado
% Revised: August 06, 2021
%
%
% Synopsis:
% Implements Euler's explicit method for the numerical solution of 
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



function odesol = Euler(N, trange, yi, f)
%   Synopsis:
%   Solves a first-order ordinary differential equation with Euler's
%   explicit method.
%
%   inputs:
%   N       number of time-steps for numerical integration
%   trange  t[ime] range, the initial, ti, and final times, tf
%   yi      initial value of the state variable, yi = y(t = ti)
%   f       f(t, y), RHS of the initial-value problem: dy/dt = f(t, y)
%
%   output:
%   odesol  second-rank numpy array, [t, y]

    y  = zeros(N + 1, 1);	% preallocates array for speed
    ti = trange(1);	tf = trange(2);
    dt = (tf - ti) / N;
    t  = linspace(ti, tf, N + 1)';

    y(1) = yi;
    for i = 1:N
        y(i + 1) = y(i) + dt * f( t(i), y(i) );
    end

    odesol = [t, y];   % packs solution into a second-rank array
end
