% Computational Solutions of Engineering Problems         August 30, 2022
% IST 4360
% Prof. M Diaz-Maldondo
%
% Synopsis:
% Tests the Gaussian elimination method.
%
% Copyright (c) 2022 Misael Diaz-Maldonado
% This file is released under the GNU General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% References:
% [0] A Gilat and V Subramaniam, Numerical Methods for Engineers and
%     Scientists: An Introduction with Applications using MATLAB
% [1] H Wendland, Numerical Linear Algebra
% [2] A Gilat, MATLAB: An Introduction with Applications


clear
close all
clc


% defines the system of linear equations
A = [
    [2, 1, 1, 0],
    [4, 3, 3, 1],
    [8, 7, 9, 5],
    [6, 7, 9, 8],
];

b = [1, 1, 2, 4]';


% solves the system of linear equations by Gaussian elimination
x = Gaussian(A, b);


% post-processing

% prints the elements of the solution vector on the console
fprintf('Solution:\n');
for i = 1:numel(x)
    fprintf('x(%d): %8.4f\n', i, x(i));
end

% computes and displays the numeric error
err = norm(A * x - b);
fprintf('Numeric Error: %.15f\n', err);
