/*
 * Applied Numerical Analysis                                 July 26, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: nlsolvers.c
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
 * [2] www.geeksforgeeks.org/error-handling-c-programs
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "nlsolvers.h"

// implementations
double bisect_varg ( double lb, double ub, double f(const double),
		     opt_args opt )
{	// shadow function which handles initialization of optional args
	double tol = opt.tol? opt.tol: TOL ;
	int max_iter = opt.max_iter? opt.max_iter: MAX_ITER ;
	return bisect_base (lb, ub, f, tol, max_iter) ;
}


double regfal_varg ( double lb, double ub, double f(const double),
		     opt_args opt )
{	// shadow function
	double tol = opt.tol? opt.tol: TOL ;
	int max_iter = opt.max_iter? opt.max_iter: MAX_ITER ;
	return regfal_base (lb, ub, f, tol, max_iter) ;
}


double shifter_varg ( double lb, double ub, double f(const double),
		      opt_args opt )
{	// shadow function
	double tol = opt.tol? opt.tol: TOL ;
	int max_iter = opt.max_iter? opt.max_iter: MAX_ITER ;
	return shifter_base (lb, ub, f, tol, max_iter) ;
}


double bisect_base ( double lb, double ub, double f(const double),
		     double tol, int max_iter )
{	// Bisection Method

	char nm[] = "Bisection" ;
	check_bracket (lb, ub, f, nm) ;
	check_bounds  (&lb, &ub) ;

	int n = 0 ;
	double xm, fm ;

        do fm = bisector (&lb, &ub, &xm, f) ;
        while (++n != max_iter && fm > tol) ;

	report (n, nm) ;
	return xm ;

}


double regfal_base ( double lb, double ub, double f(const double),
	             double tol, int max_iter )
{	// Regula Falsi Method

	char nm[] = "Regula-Falsi" ;
	check_bracket (lb, ub, f, nm) ;
	check_bounds  (&lb, &ub) ;

	int n = 0 ;
	double xn, fn ;

        do fn = interp (&lb, &ub, &xn, f) ;
        while (++n != max_iter && fn > tol) ;

	report (n, nm) ;
	return xn ;

}


double shifter_base ( double lb, double ub, double f(const double),
	         double tol, int max_iter )
{	// Shifter Method

	char nm[] = "Shifter" ;
	check_bracket (lb, ub, f, nm) ;
	check_bounds  (&lb, &ub) ;

	int n = 0 ;
	double xn, fn ;

        do fn = shift (&lb, &ub, &xn, f) ;
        while (++n != max_iter && fn > tol) ;

	report (n, nm) ;
	return xn ;

}


void report (const int n, char nm[]) {
	// reports if the method has been successful
	if (n != MAX_ITER) {
		printf("%s Method:\n", nm) ;
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

	return absval(fm) ;
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

	return absval(fn) ;
}


double shift ( double *lb, double *ub, double *xn,
	       double f(const double) )
{	// uses hybrid approach to approximate the root
	double fn ;
	double xb = 0.5 * (*lb + *ub) ;
	double xf = ( *lb * f(*ub) - *ub * f(*lb) ) / ( f(*ub) - f(*lb) ) ;

	// shifts towards the step (presumably) closer to the root
	if ( absval(f(xb)) < absval(f(xf)) )
		*xn = xb ;
	else
		*xn = xf ;

	fn = f(*xn) ;

	if (f(*lb) * fn < 0.)
		*ub = *xn ;
	else
		*lb = *xn ;

	return absval(fn) ;
}


void check_bracket ( double lb, double ub,
		     double f(const double), char nm[] )
{
	// complains if there's no root in given interval and aborts
	if ( f(lb) * f(ub) > 0. ) {
		fprintf(stderr, "\n%s Method:\n", nm) ;
		fprintf(stderr, "No root in given interval ... \n") ;
		fprintf(stderr, "Try again, aborting execution ... \n\n") ;
		exit(EXIT_FAILURE) ;
	}
}


void check_bounds ( double *lb, double *ub ) {
	// ensures that the lower bound is less than the upper bound.
	double up = *lb ;
	if ( *lb > *ub ) {
		*lb = *ub ;
		*ub =  up ;
	}
}


/*
 * TODO:
 * [x] Implement guards against "empty" ranges (lower > upper bound)
 * [x] Implement guards against applying the method on an interval
 *     that does not enclose a root.
 * [ ] Define default values for the tolerance and maximum number of
 *     iterations. The user may wish to override these parameters so these
 *     must be included in the argument lists. Use the constants as default
 *     values for these parameters. Consider implementing a configuration
 *     struct. Possible implementation follows from stackoverflow inquiry
 *     (bk's answer):
 *     
 *     1. stackoverflow.com/questions/1472138/c-default-arguments
 *     2. modelingwithdata.org/arch/00000022.htm
 *
 *     The latter provides a detailed explanation to the proposed solution
 *     (from same stackoverflow user) which should not be overlooked!
 *
 * [x] report function should display the name of the numerical method
 *
 */
