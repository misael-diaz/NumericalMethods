% Computational Solutions                                    April 21, 2022
% IST 4360
% Prof. M Diaz Maldonado
%
% Synopsis:
% Uses a Lagrange 3rd-degree polynomial to interpolate a dataset (xi, yi).
%
%
% Copyright (c) 2022 Misael Diaz-Maldonado
% This file is released under the GNU General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%
% References:
% [0] A Gilat and V Subramanian, Numerical Methods for Engineers and
%     Scientists: An Introduction with Applications using MATLAB
% [1] A Gilat, MATLAB: An Introduction with Applications, 6th edition


clear
close all
clc

format long g


% Dataset:
% xi, yi coordinates
x1 = -1.2; y1 = 18.8;
x2 =  0.2; y2 =  5.0;
x3 =  2.0; y3 = 16.0;
x4 =  3.5; y4 = 15.0;


% Lagrange 3rd-degree interpolating polynomial:
% Note: the anonymous function has been defined for a scalar input
f = @(x) (								...
    y1*(x - x2)*(x - x3)*(x - x4) / ( (x1 - x2)*(x1 - x3)*(x1 - x4) ) + ...
    y2*(x - x1)*(x - x3)*(x - x4) / ( (x2 - x1)*(x2 - x3)*(x2 - x4) ) + ...
    y3*(x - x1)*(x - x2)*(x - x4) / ( (x3 - x1)*(x3 - x2)*(x3 - x4) ) + ...
    y4*(x - x1)*(x - x2)*(x - x3) / ( (x4 - x1)*(x4 - x2)*(x4 - x3) )   ...
);


% Post-processing:

% computes the differences between the data and the approximation returned
% by the interpolating function f(x)
x = [x1; x2; x3; x4];
y = [y1; y2; y3; y4];

% should display an array of zeros to indicate success
diffs = 0;
fprintf("\nDifferences:\n")
for i = 1:length(x)
	diffs = diffs + abs( y(i) - f( x(i) ) );
	fprintf( "%.12f\n", abs( y(i) - f( x(i) ) ) )
end

% as above but accumulates the differences
fprintf("\nCumulative Difference:\n")
fprintf("%.12f\n", diffs)

% Plots the dataset and the interpolating polynomial:
figure(1)

plot(x, y, 'r+')		% plots the xi, yi coordinates
hold on
% performs an elementwise computation of the interpolating polynomial
x = linspace(-2, 5, 256);
y = zeros( size(x) );
for i = 1:length(x)
	y(i) = f( x(i) );	% Note: the function expects a scalar input
end
plot(x, y, '--k')		% plots interpolating polynomial
xlabel('x')
ylabel('f(x)')
