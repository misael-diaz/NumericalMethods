% Computational Solutions of Engineering Problems         August 30, 2022
% IST 4360
% Prof. M Diaz-Maldondo
%
%                        Gaussian Elimination
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


function x = Gaussian(A, b)
% Synopsis:
% Possible implementation of the Gaussian elimination method.

% Inputs:
% A       (2nd-rank array) N x N square matrix
% b       (1st-rank array) N coefficient vector
%
% Output:
% x       (1st-rank array) N solution vector
%
% Comments:
% Note that this implementation does not mitigate error propagation by
% selecting the optimal pivot. It uses as a pivoting element whatever
% value that happens to be on the diagonal and for this reason it will
% **fail** sometimes. As a student of this course you are encouraged to
% improve this implementation by using the optimal pivot on each step.

% copies matrix A into the work matrix W
W = A;
% appends the coefficient vector to obtain the work augmented matrix
if (size(b, 2) > 1)
    % transposes into column vector if given a row vector
    b = b';
end
W = [W, b];
% gets the number of rows and columns, respectively
N = size(W, 1);
M = size(W, 2);

% reduces the matrix W into its upper triangular form
for p = 1:(N - 1)
    for i = (p + 1):N
        m = -W(i, p) / W(p, p);
        for j = p:M
            W(i, j) = W(i, j) +  m * W(p, j);
        end
    end
end

% obtains the solution vector `x' via back substitution
x = backSubs(W);

end
