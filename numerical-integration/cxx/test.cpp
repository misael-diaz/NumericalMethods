/*
 * Applied Numerical Analysis			              July 27, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: test.cpp
 *
 * Synopsis:
 * Tests (some) numerical integration techniques using C++ facilities.
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

#include <cmath>
#include <iomanip>
#include <iostream>
import numerical_integrators ;
#define absval(x) (x < 0.)? -x: x

using std::cout ;
using std::endl ;
using std::fixed ;
using std::scientific ;
using std::streamsize ;
using std::setprecision ;

using numint::lsum ;
using numint::rsum ;
using numint::trap ;

void error (double exact, double numeric[3], double err[3]) ;
void report (double ni[3], double err[3]) ;


int main() {
// Synopsis:
// Integrates f(x) = exp(x) in the interval [a, b] using N intervals via
// Rectangle and Trapezoid Methods.


	/* defines f(x) = exp(x) as a lambda function */
	auto f = [](const double& x) -> double { return exp(x) ; } ;
	const int N = 255 ;		// intervals
	double a = 0., b = 1. ;		// integration limits [a, b]
	double ei = f(b) - f(a) ;	// exact integral for f(x) = exp(x)


	double ni[3], err[3] ;
	ni[0] = lsum (a, b, N, f) ;	// Left Riemann
	ni[1] = rsum (a, b, N, f) ;	// Right Riemann
	ni[2] = trap (a, b, N, f) ;	// Trapezoid Method


	error (ei, ni, err) ;
	report (ni, err) ;
	return 0 ;
}


void report (double ni[3], double err[3]) {
	// tabulates results
	streamsize prec = cout.precision() ;

	cout << endl << endl ;
	cout << "Numerical Method \t  Result \t  % Error " << endl ;
	cout << "------------------------------------------------"
	     << "-----" << endl ;
	cout << "Left Riemann     \t " << fixed << setprecision(6)
	     << ni[0] << "  \t "  << scientific << setprecision(2)
	     << err[0]<< " % " << setprecision(prec) << endl ;

	cout << "Right Riemann    \t " << fixed << setprecision(6)
	     << ni[1] << "  \t "  << scientific << setprecision(2)
	     << err[1]<< " % " << setprecision(prec) << endl ;

	cout << "Trapezoid Method \t " << fixed << setprecision(6)
	     << ni[2] << "  \t "  << scientific << setprecision(2)
	     << err[2]<< " % " << setprecision(prec) << endl ;
	cout << endl << endl ;
}


void error (double exact, double numeric[3], double err[3]) {
        // obtains the error in percentange for every element of the arrays
        for (int i = 0 ; i != 3 ; ++i)
                err[i] = absval( (numeric[i] - exact) / exact * 100.0 ) ;
}


/*
 * TODO:
 * [x] partition numerical methods into a module (since c++20)
 */
