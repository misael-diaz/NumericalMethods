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
#include <vector>
#include <numeric>
#include <iomanip>
#include <iostream>
#include <algorithm>
/* import numerical_integrators ; */
#define absval(x) (x < 0.)? -x: x

using std::cout ;
using std::endl ;
//using std::iota ;
//using std::vector ;
//using std::transform ;
using std::fixed ;
using std::scientific ;
using std::streamsize ;
using std::setprecision ;

double lsum (double, double, const int, double f (const double&) ) ;
double rsum (double, double, const int, double f (const double&) ) ;
double trap (double, double, const int, double f (const double&) ) ;

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


// implementation of numerical integration techniques
double lsum (double a, double b, const int N, double f (const double&) ) {
	// Left Riemann Sum
	// Integrates f(x) in the interval [a, b] using N intervals.

	std::vector<int> idx(N+1) ;
	std::vector<double> x(N+1) ;
	std::iota(idx.begin(), idx.end(), 0) ;
	double dx = (b - a) / ( (double) N ) ;	// step

	// transforms index into a vector x = [a, b] of N + 1 elements
	auto t = [&a, &dx](const int& i) -> double {
		return (a + i * dx) ;
	} ;

	std::transform(idx.begin(), idx.end(), x.begin(), t) ;
	// transforms x -> f(x)
	std::transform(x.begin(), x.end(), x.begin(), f) ;

	double sum = 0 ;
	sum = std::accumulate(x.begin(), x.begin() + x.size() - 1, sum) ;

	return (dx * sum) ;
}


double rsum (double a, double b, const int N, double f (const double&) ) {
	// Right Riemann Sum
	// Integrates f(x) in the interval [a, b] using N intervals.

	std::vector<int> idx(N+1) ;
	std::vector<double> x(N+1) ;
	std::iota(idx.begin(), idx.end(), 0) ;
	double dx = (b - a) / ( (double) N ) ;

	auto t = [&a, &dx](const int& i) -> double {
		return (a + i * dx) ;
	} ;

	std::transform(idx.begin(), idx.end(), x.begin(), t) ;
	std::transform(x.begin(), x.end(), x.begin(), f) ;

	double sum = 0 ;
	sum = std::accumulate(x.begin() + 1, x.begin() + x.size(), sum) ;

	return (dx * sum) ;
}


double trap (double a, double b, const int N, double f (const double&) ) {
	// Trapezoid Method
	// Integrates f(x) in the interval [a, b] using N intervals.
	
	std::vector<int> idx(N+1) ;
	std::vector<double> x(N+1) ;
	std::iota(idx.begin(), idx.end(), 0) ;
	double dx = (b - a) / ( (double) N ) ;

	auto t = [&a, &dx](const int& i) -> double {
		return (a + i * dx) ;
	} ;

	std::transform(idx.begin(), idx.end(), x.begin(), t) ;
	std::transform(x.begin(), x.end(), x.begin(), f) ;

	double s = 0 ;
	s = std::accumulate(x.begin() + 1, x.begin() + x.size() - 1, s) ;
	s = f(a) + (2.0 * s) + f(b) ;

	return (0.5 * dx * s) ;
}


/*
 * TODO:
 * [ ] partition numerical methods into a module (since c++20)
 */


/*
 * Comments:
 * Using fully qualified names inside methods in anticipation of the
 * partitioning.
 *
 */
