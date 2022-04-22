% Computational Solutions                                    April 21, 2022
% IST 4360
% Prof. M Diaz Maldonado
%
% Synopsis:
% Obtains the coefficients of the 3rd degree interpolating polynomial for a
% give dataset (xi, yi):
%
%                     f(x) = a x^3 + b x^2 + c * x + d
%
% where f(x) is the interpolating polynomial, and (a, b, c, d) are the
% coefficients.
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


% Applies the direct-fit method:
% constructs the system of linear equations
A = [x1^3, x1^2, x1, 1;
     x2^3, x2^2, x2, 1;
     x3^3, x3^2, x3, 1;
     x4^3, x4^2, x4, 1];

b = [y1; y2; y3; y4];

% solves for the coefficients of the 3rd degree interpolating polynomial
x = A\b;

% references the coefficients of the 3rd degree interpolating polynomial
a = x(1);	b = x(2);	c = x(3);	d = x(4);
% defines the 3rd degree polynomial as a lambda (or anonymous) function
f = @(x) (a * x.^3 + b * x.^2 + c * x + d);


% Post-processing:
% prints the determinant condition number of the system of linear equations
% Note: if the system is ill-conditioned the interpolation is unreliable
fprintf( "Determinant:      %.2f\n", det(A) )
fprintf( "Condition number: %.2f\n", cond(A) )


% computes the differences between the data and the approximation returned
% by the interpolating function f(x)
x = [x1; x2; x3; x4];
y = [y1; y2; y3; y4];

% should display an array of zeros to indicate success
fprintf("\nDifferences:\n")
for i = 1:length(x)
	fprintf( "%.12f\n", y(i) - f( x(i) ) )
end

% as above but accumulates the differences
fprintf("\nCumulative Difference:\n")
% should display a zero (scalar) to indicate success
fprintf( "%.12f\n", sum( y - f(x) ) )


% Plots the dataset and the interpolating polynomial:
figure(1)

plot(x, y, 'r+')		% plots the xi, yi coordinates
hold on
x = linspace(-2, 5, 256);
plot(x, f(x), '--k')		% plots interpolating polynomial
xlabel('x')
ylabel('f(x)')
