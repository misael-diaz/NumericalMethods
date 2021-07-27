/*
 * source: test_bracketing.c
 * author: misael-diaz
 * date:   2021/07/26
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

double f (const double) ;	/* nonlinear function, f(x) */

int main() {
	// Solves for the positive root of the nonlinear function f(x)
	double lb, ub ;
	double x1, x2 ;

	lb = 2.0e-2 ;	ub = 7.0e-2 ;	/* bracketing inteval [lb, ub] */


	x1 = bisect (lb, ub, f) ;	// Bisection
	x2 = regfal (lb, ub, f) ;	// Regula Falsi


	// displays roots (both methods converge to the root)
	printf("x: %12.6f %12.6f\n", x1, x2) ;
	printf("f(x): %24.6e %24.6e\n", f(x1), f(x2) ) ;

	return 0 ;
}


double f (const double x) {
	// defines the nonlinear function f(x)
	return ( 1.0 / sqrt(x) + 2.0 * log10(0.024651/3.7 +
	         2.51/(9655526.5 * sqrt(x) ) ) ) ;
}
