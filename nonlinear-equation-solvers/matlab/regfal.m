% Applied Numerical Analysis                                 July 26, 2021
% ME 2020 FA21
% Prof. M Diaz-Maldonado
%
%
% Synopsis:
% Possible implementation of the Regula Falsi method in MATLAB.
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


function x = regfal(a, b, f, opt)
    % Synopsis:
    % Regula Falsi Method. Finds the root of the function f(x) enclosed by
    % the interval [a, b].

    name = 'Regula Falsi';
    % sets default values for the tolerance and maximum number of iterations
    TOL      = 1.0e-8;
    MAX_ITER = 100;
    VERBOSE  = 0;

    optset
    check_bounds
    check_bracket

    n = 1;
    x = ( a * f(b) - b * f(a) ) / ( f(b) - f(a) );
    while ( n ~= MAX_ITER && abs( f(x) ) > TOL )

        if ( f(a) * f(x) < 0 )
            b = x;
        else
            a = x;
        end

        x = ( a * f(b) - b * f(a) ) / ( f(b) - f(a) );
        n = n + 1;
    end

    report
    return


    %```nested functions:```%
    function optset
        % Synopsis: Uses configuration struct if provided by user.
        if ( exist('opt', 'var') )
            if ( isfield(opt, 'tol') )
                TOL      = opt.tol;
            end

            if ( isfield(opt, 'max_iter') )
                MAX_ITER = opt.max_iter;
            end

            if ( isfield(opt, 'verbose') )
                VERBOSE  = opt.verbose;
            end
	end
    end


    function check_bounds
        % Synopsis: Ensures the lower bound is less than the upper bound.
        if (a > b)
            up = a;
            a  = b;
            b  = up;
        end
    end


    function check_bracket
        % Synopsis: Complains if there's no root in the given interval.
        if ( f(a) * f(b) > 0 )
            errID  = 'NonlinearSolver:BracketingException';
            errMSG = 'No root exists in the given interval [a, b]';
            errMSG = [name, ' Method: ', errMSG];
            except = MException(errID, errMSG);
            throw(except);
        end
    end


    function report
        % Synopsis: Reports to the user if the method has been successful.
        if ( n ~= MAX_ITER )
	    if (VERBOSE)
                fprintf('%s Method:\n', name)
                fprintf('>> Solution found in %d iterations\n', n)
	    end
        else
            errID  = 'NonlinearSolver:ConvergenceException';
            errMSG = 'requires additional iterations for convergence';
            errMSG = [name, ' method ', errMSG];
            except = MException(errID, errMSG);
            throw(except);
        end
    end


end


% TODO:
% [x] Add code for user-defined tolerance and maximum number of iterations.
