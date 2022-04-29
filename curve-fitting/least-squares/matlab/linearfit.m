% Computational Solutions                                    April 28, 2022
% IST 4360
% Prof. M. Diaz-Maldonado
%
%
% Synopsis:
% Applies the Least Squares Method to obtain the ``best fit'' of a dataset.
%
%
% Copyright (c) 2022 Misael Diaz-Maldonado
% This file is released under the GNU General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%
% References:
% [0] B Munson, D Young, T Okiishi, and W Huebsch, Fundamentals of
%     Fluid Mechanics, 8th edition
% [1] A Gilat and V Subramanian, Numerical Methods for Engineers and
%     Scientists: An Introduction with Applications using MATLAB
% [2] R Johansson, Numerical Python: Scientific Computing and Data
%     Science Applications with NumPy, SciPy, and Matplotlib, 2nd edition


clear
close all
clc


% Dataset:
xi = linspace(1, 6, 6);
yi = [13.1, 26.0, 39.5, 52.7, 64.9, 78.6];

% applies the least squares method to obtain the ``best fit'', y = f(x)
[slope, intercept, f] = LinearFit(xi, yi);


% Plots the dataset and the ``best fit'' for visualization:
figure(1)
x = linspace(min(xi), max(xi), 256);

% plots the dataset
plot(xi, yi, 'r*')
hold on
% plots the ``best fit'' y = f(x)
plot(x, f(x), '--k')
xlabel('x')
ylabel('y')
title('Least Squares')


% Displays info about the fitting parameters and ``goodness of fit'':
fprintf('\nLinear Fit via Least Squares:\n')
fprintf('Slope: %.12f\n', slope)
fprintf('Intercept: %.12f\n', intercept)
fprintf("Pearson's correlation coefficient, R^2: %.12f\n", Pearson(xi, yi))



% Defines Methods:
function [s, i, f] = LinearFit(xi, yi)
    % Obtains the linear fit of the dataset xi, yi via Least Squares

    % solves the system of linear equations by Crammer's rule
    d = det([
        [sum(xi.^2), sum(xi)];
        [sum(xi), numel(xi)]
    ]);


    D1 = det([
        [sum(xi .* yi), sum(xi)];
        [      sum(yi), numel(xi)]
    ]);


    D2 = det([
        [sum(xi.^2), sum(xi .* yi)];
        [   sum(xi), sum(yi)]
    ]);

    % calculates the slope and intercept, respectively
    s = (D1 / d);
    i = (D2 / d);

    % defines the ``best fit'' y = f(x)
    f = @(x) (i + s * x);
end


function R = Pearson(xi, yi)
    % Returns the Pearson's coefficient of correlation, R^2

    % computes normalized deviations
    Dx = (mean(xi) - xi) / std(xi);
    Dy = (mean(yi) - yi) / std(yi);

    % computes Peasons's correlation coefficient
    [R, ~, ~] = LinearFit(Dx, Dy);
    R = R^2;
end
