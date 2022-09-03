% Computational Solutions of Engineering Problems         August 30, 2022
% IST 4360
% Prof. M Diaz-Maldondo
%
%                        Back Substitution Method
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

function x = backSubs(U)
% Synopsis:
% Possible implementation of the back substitution technique.
%
% Input:
% U	(2nd-rank array) N x (N + 1) upper triangular matrix
%
% Output:
% x	(1st-rank array) N x 1 solution vector

% gets the number of rows and columns
N = size(U, 1);
M = size(U, 2);

% preallocates the solution vector
x = zeros([N, 1]);
% applies the back substitution technique
for i = N:-1:1
    x(i) = U(i, M);
    for j = (i + 1):N
        x(i) = x(i) - U(i, j) * x(j);
    end
    x(i) = x(i) / U(i, i);
end

end
