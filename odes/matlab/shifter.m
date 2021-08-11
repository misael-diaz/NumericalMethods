% Applied Numerical Analysis                                 July 28, 2021
% ME 2020 FA21
% Prof. M Diaz-Maldonado
%
%
% Synopsis:
% Implements a hybrid from the bisection and the regula falsi methods.
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


function x = shifter(a, b, f, opt)
    % Synopsis:
    % Shifter Method. Finds the root of the function f(x) enclosed by
    % the interval [a, b].

    name = 'Shifter';
    TOL      = 1.0e-12;
    MAX_ITER = 256;
    VERBOSE  = 0;

    optset
    check_bounds
    check_bracket

    n = 1;
    shift
    while ( n ~= MAX_ITER && abs( f(x) ) > TOL )

    % updates the bracketing interval
        if ( f(a) * f(x) < 0 )
            b = x;
        else
            a = x;
        end

	shift
        n = n + 1;
    end

    report
    return


    %```nested functions:```%
    function shift
        % Synopsis: Shifts to the step (presumably) closest to the root.
        x1 = 0.5 * (a + b);
        x2 = ( a * f(b) - b * f(a) ) / ( f(b) - f(a) );
        if ( abs(f(x1)) < abs(f(x2)) )
            x = x1;
        else
            x = x2;
        end
    end

    function optset
        % Synopsis: Uses configuration struct if provided by user.
        if ( exist('opt', 'var') )
            TOL      = opt.tol;
            MAX_ITER = max([1, opt.max_iter]);
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
            errID  = 'NonlinearSolver:IterationException';
            errMSG = 'requires more iterations for convergence';
            errMSG = [name, ' Method: ', errMSG];
            except = MException(errID, errMSG);
	    throw(except);
        end
    end


end


% TODO:
% [x] Add code for user-defined tolerance and maximum number of iterations.
