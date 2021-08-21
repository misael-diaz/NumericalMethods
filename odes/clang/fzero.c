/*
 * Applied Numerical Analysis                               August 13, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: fzero.c
 *
 * Synopsis:
 * Implements the nonlinear equation solvers fzero.
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

#include "fzero.h"

// implementations
static void check_bracket ( double lb, double ub, double f(double, void*),
                     char nm[], void* vprms )
{
	// complains if there's no root in given interval and aborts
	if ( !signbit( f(lb, vprms) * f(ub, vprms) ) ) {
		fprintf(stderr, "\n%s Method:\n", nm) ;
		fprintf(stderr, "No root in given interval ... \n") ;
		fprintf(stderr, "Try again, aborting execution ... \n\n") ;
		exit(EXIT_FAILURE) ;
	}
}


static void check_bounds ( double *lb, double *ub ) {
	// ensures that the lower bound is less than the upper bound
	double up = *lb ;
	if ( *lb > *ub ) {
		*lb = *ub ;
		*ub =  up ;
	}
}


static void report (int max_iter, int n, char nm[], bool verbose) {
	// reports if the method has been successful
	if (n != max_iter) {

		if (verbose) {
			printf("%s Method:\n", nm) ;
			printf("solution found in %d iterations\n", n) ;
		}

	} else {
		fprintf(stderr, "%s method needs additional ", nm) ;
		fprintf(stderr, "iterations for convergence\n") ;
		exit(EXIT_FAILURE) ;
	}
}


static double shift ( double *lb, double *ub, double *xn,
	       double f(double, void*), void *vprms )
{
	double fn ;
	double xb = 0.5 * (*lb + *ub) ;
	double xf = ( *lb * f(*ub, vprms) - *ub * f(*lb, vprms) ) /
		    ( f(*ub, vprms) - f(*lb, vprms) ) ;

	if ( absval( f(xb, vprms) ) < absval( f(xf, vprms) ) )
		*xn = xb ;
	else
		*xn = xf ;

	fn = f(*xn, vprms) ;

	if ( signbit( f(*lb, vprms) * fn) )
		*ub = *xn ;
	else
		*lb = *xn ;

	return absval(fn) ;
}


// method interfaces
double fzero ( double lb, double ub, double f(double, void*), void *vprms )
{	// Shifter Method from nonlinear solvers using the default options

	char nm[] = "fzero" ;
	check_bracket (lb, ub, f, nm, vprms) ;
	check_bounds  (&lb, &ub) ;

	int n = 0 ;
	double xn, fn ;

        do fn = shift (&lb, &ub, &xn, f, vprms) ;
        while (++n != MAX_ITER && fn > TOL) ;

	report (MAX_ITER, n, nm, VERBOSE) ;
	return xn ;

}
