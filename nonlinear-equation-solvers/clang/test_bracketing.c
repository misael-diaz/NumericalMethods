/*
 * Applied Numerical Analysis
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: test_bracketing.c
 *
 * Synopsis:
 * Tests the implementation of the bisection and regula falsi methods.
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

#include <math.h>
#include <stdio.h>
#include "nlsolvers.h"

double f (double) ;	/* nonlinear function, f(x) */

int main() {
	// Solves for the positive root of the nonlinear function f(x)

	// defines solver tolerance and maximum number of iterations
	nls_conf opts = {.tol = 1.0e-12, .max_iter = 256} ;

	double lb, ub ;			/* bracketing inteval [lb, ub] */
	double x1, x2, x3 ;		// roots

	lb = 1.0e-2 ;	ub = 9.0e-2 ;

	// solves for the root of f(x) with the specified method
	x1 = bisect (lb, ub, f, opts) ;	// Bisection
	x2 = regfal (lb, ub, f, opts) ;	// Regula Falsi
	x3 = shifter(lb, ub, f, opts) ;	// Shifter


	// displays roots (methods converge to the root)
	printf("x: %24.15f %24.15f %24.15f\n", x1, x2, x3) ;
	printf("f(x): %24.6e %24.6e %24.6e\n", f(x1), f(x2), f(x3) ) ;


/*      tests non-enclosing interval        */
//	lb = 6.0e-2 ;	ub = 7.0e-2 ;
//	x1 = bisect (lb, ub, f) ;	passed
//	x2 = regfal (lb, ub, f) ;	passed
	return 0 ;
}


double f (double x) {
	// defines the nonlinear function f(x)
	return ( 1.0 / sqrt(x) + 2.0 * log10(0.024651/3.7 +
	         2.51/(9655526.5 * sqrt(x) ) ) ) ;
}
