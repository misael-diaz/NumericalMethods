/*
 * source: test-bracketing.cpp
 * author: misael-diaz
 * date:   2021/07/20
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

#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
import nonlinear_solvers ;
using std::cout ;
using std::endl ;
using std::streamsize ;
using std::setprecision ;
using std::runtime_error ;
int main() {
	// Solves for the positive root of the nonlinear function f(x)
	
	
	auto f = [](const double& x) -> double {
	    return ( 1.0 / sqrt(x) + 2.0 * log10(0.024651/3.7 +
	             2.51/(9655526.5 * sqrt(x) ) ) ) ;
	} ;
	double lb, ub ;		// lower and upper bounds [lb, ub]
	double x1, x2, x3 ;	// numerical approximation of the root


	try {
		lb = 1.0e-2 ;	ub = 9.0e-2 ;
        	x1 = nlsolver::bisect (lb, ub, f) ;	// bisection
        	x2 = nlsolver::regfal (lb, ub, f) ;	// regula falsi
		x3 = nlsolver::shifter(lb, ub, f) ;	// shifter
	} catch (runtime_error& err) {
		cout << err.what() << endl ;
		return 1 ;
	}


	// displays solutions (methods converge to the root)
	streamsize prec = cout.precision() ;
	cout << "x: " << setprecision(15) << x1 << ", " << x2
	     << ", "  << x3 << setprecision(prec) << endl ;
	cout << "f(x): " << f(x1) << ", " << f(x2)
	     << ", "     << f(x3) << endl ;
	return 0 ;
}
