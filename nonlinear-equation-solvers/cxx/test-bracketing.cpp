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

#include <iostream>
#include <stdexcept>
import nonlinear_solvers ;
int main() {
	// Solves for the positive root of f(x) = (x + 1) * (x - 1)
	
	
	// lambda function, f(x) = (x + 1) * (x - 1) = x * x - 1
	auto f = [](const double& x) -> double { return (x * x - 1.0) ; } ;
	double lb, ub ;		// lower and upper bounds [lb, ub]
	double x1, x2 ;		// numerical approximation of the root


	try {
		lb = 0.2 ;	ub = 1.5 ;
        	x1 = nlsolver::bisect (lb, ub, f) ;	// bisection
		lb = 0.2 ;	ub = 1.5 ;
        	x2 = nlsolver::regfal (lb, ub, f) ;	// regula falsi
	} catch (std::runtime_error& err) {
		std::cout << err.what() << std::endl ;
		return 1 ;
	}


	// displays solutions (both methods converge to the root)
	std::cout << "x: " << x1 << ", " << x2 << std::endl ;
	std::cout << "f(x): " << f(x1) << ", " << f(x2) << std::endl ;
	return 0 ;
}
