/*
 * source: nlsolvers.c
 * author: misael-diaz
 * date:   2021/07/26
 *
 * Synopsis:
 * Implements (some) nonlinear equation solvers.
 * 
 *
 * Copyright (c) 2021 Misael Diaz-Maldonado
 *
 * This file is released under the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 *
 * References:
 * [0] A Gilat and V Subramanian, Numerical Methods for Engineers and
 *     Scientists, 3rd edition.
 * [1] A Koenig and B Moo, Accelerated C++ Practical Programming by
 *     Example.
 *
 */

#include <stdio.h>
#include "nlsolvers.h"

// implementations
double bisect ( double lb, double ub, double f(const double) )
{	// Bisection Method

	int n = 0 ;
	double xm, fm ;

        do fm = bisector (&lb, &ub, &xm, f) ;
        while (++n != MAX_ITER && fm > TOL) ;

	report (n) ;
	return xm ;

}


double regfal ( double lb, double ub, double f(const double) )
{	// Regula Falsi Method

	int n = 0 ;
	double xn, fn ;

        do fn = interp (&lb, &ub, &xn, f) ;
        while (++n != MAX_ITER && fn > TOL) ;

	report (n) ;
	return xn ;

}

void report (const int n) {
	// reports if the method has been successful
	if (n != MAX_ITER) {
		printf("solution found in %d iterations\n", n) ;
	} else {
		printf("method failed to find the root, ") ; 
		printf("you may try a narrower interval\n") ;
	}
}


double bisector ( double *lb, double *ub, double *xm, 
		  double f(const double) )
{

/*
 * Synopsis:
 * Approximates the root of the nonlinear equation f(x) with the middle
 * value, xm = (lb + ub) / 2, where `lb' and `ub' are the lower
 * and upper bounds, respectively, of the bracketing interval [lb, ub].
 * As a side effect it returns |f(xm)|, where |*| = abs(*), and `*' is
 * a wildcard.
 *
 * input/output:
 * lb		lower bound, intent(inout)
 * ub           upper bound, intent(inout)
 * xm           middle value, intent(out)
 *
 * NOTE:
 * f(x) is a function pointer (last argument) bound to the user-defined
 * function representing the nonlinear equation to be solved.
 *
 */
	
        *xm = 0.5 * (*lb + *ub) ;
	double fm = f(*xm) ;

	/* chooses the (root) bracketing interval */
	if (f(*lb) * fm < 0.) 	
		*ub = *xm ;
	else
		*lb = *xm ;


	return (fm < 0.) ? fm = -fm: fm ;	// implements abs(x)
}


double interp ( double *lb, double *ub, double *xn,
		double f(const double) )
{	// like bisector but uses interpolation to approximate the root
        *xn = ( *lb * f(*ub) - *ub * f(*lb) ) / ( f(*ub) - f(*lb) ) ;
	double fn = f(*xn) ;

	if (f(*lb) * fn < 0.)
		*ub = *xn ;
	else
		*lb = *xn ;

	return (fn < 0.) ? fn = -fn: fn ;
}


/*
 * TODO:
 * [ ] Implement guards against "empty" ranges (lower > upper bound)
 * [ ] Implement guards against applying the method on an interval
 *     that does not enclose a root.
 * [ ] Define default values for the tolerance and maximum number of
 *     iterations. The user may wish to override these parameters so these
 *     must be included in the argument lists. Use the constants as default
 *     values for these parameters.
 * [ ] report function should display the name of the numerical method
 *
 */