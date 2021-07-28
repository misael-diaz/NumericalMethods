/*
 * Applied Numerical Analysis				      July 27, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: test.c
 *
 * Synopsis:
 * Tests some numerical integration methods.
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
#include "numint.h"
#define absval(x) (x < 0.)? -x: x


// prototypes
double f (const double) ;	/* nonlinear function, f(x) */
void error (double exact, double numeric[3], double err[3]) ;
void report (double ni[3], double err[3]) ;


int main() {
// Synopsis:
// Integrates f(x) in the interval [a, b] using N intervals via Rectangle
// and Trapezoidal Methods.

	int N = 255 ;			// intervals
	double a = 0.0, b = 1.0 ;	// integration limits [a, b]
	double ei = f(b) - f(a) ;	// exact integral for f(x) = exp(x)
	

	double ni[3], err[3] ;		// n[umeric] i[ntegral] and err[or]
	ni[0] = lsum(a, b, N, f) ;	// Left  Riemann Sum
	ni[1] = rsum(a, b, N, f) ;	// Right Riemann Sum
	ni[2] = trap(a, b, N, f) ;	// Trapezoidal Method


	error (ei, ni, err) ;
	report (ni, err) ;
	return 0 ;
}


// implementations
void report (double ni[3], double err[3]) {
	// tabulates results
	printf("\n\n") ;
	printf("Numerical Method \t  Result \t  %% Error\n") ;
	printf("-----------------------------------------------------\n") ;
	printf("Left Riemann     \t %8.6f  \t %6.2e %% \n",
	       ni[0], err[0]) ;
	printf("Right Riemann    \t %8.6f  \t %6.2e %% \n",
	       ni[1], err[1]) ;
	printf("Trapezoid Method \t %8.6f  \t %6.2e %% \n",
	       ni[2], err[2]) ;
	printf("\n\n") ;
}


void error (double exact, double numeric[3], double err[3]) {
	// obtains the error in percentange for every element of the arrays
	for (int i = 0 ; i != 3 ; ++i)
		err[i] = absval( (numeric[i] - exact) / exact * 100.0 ) ;
}


double f (const double x) {
	// defines the nonlinear function f(x) = exp(x)
	return exp(x) ;
}


/*
 *
 * Comments:
 * The MACROS absval(x) implements the absolute value function for doubles.
 *
 */
