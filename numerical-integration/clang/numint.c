/*
 * Applied Numerical Analysis				      July 27, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: numint.c
 *
 * Synopsis:
 * Implements (some) numerical integration methods.
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

#include "numint.h"

double lsum ( double a, double b, int N, double f(const double) ) {
// Left Riemann Sum
// Integrates the function f(x) in the interval [a, b] using N intervals.

	double dx = (b - a) / ( (double) N ) ;

	double x, sum = 0.0 ;
	for (int i = 0 ; i != N ; ++i) {
		x = a + i * dx ;
		sum += f(x) ;
	}

	return (dx * sum) ;

}


double rsum ( double a, double b, int N, double f(const double) ) {
// Right Riemann Sum
// Integrates the function f(x) in the interval [a, b] using N intervals.

	double dx = (b - a) / ( (double) N ) ;

	double x, sum = 0.0 ;
	for (int i = 1 ; i != N + 1 ; ++i) {
		x = a + i * dx ;
		sum += f(x) ;
	}

	return (dx * sum) ;

}


double trap ( double a, double b, int N, double f(const double) ) {
// Trapezoidal Method
// Integrates the function f(x) in the interval [a, b] using N intervals.

	double x, sum ;
	double dx = (b - a) / ( (double) N ) ;
	
	sum = f(a) ;
	for (int i = 1 ; i != N ; ++i) {
		x = a + i * dx ;
		sum += 2.0 * f(x) ;
	}
	sum += f(b) ;

	return (0.5 * dx * sum) ;

}


/*
 *
 * Comments:
 * The trapezoidal method could have been implemented by invoking the
 * left and right Riemann sums (as was done in Python and MATLAB) but
 * decided to take a less expensive route, involving fewer calls to
 * the function f(x).
 *
 */
