/*
 * Applied Numerical Analysis                                 July 27, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: module-numerical-integrators.cpp
 *
 * Synopsis:
 * Implements (some) numerical integration techniques.
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


module ;
#include <vector>
#include <numeric>
#include <algorithm>
export module numerical_integrators ;


export namespace numint {
	double lsum( double, double, const int, double f(const double&) ) ;
	double rsum( double, double, const int, double f(const double&) ) ;
	double trap( double, double, const int, double f(const double&) ) ;
}


// implementations
double numint::lsum (double a, double b, const int N,
		     double f (const double&) )
{
	// Left Riemann Sum
	// Integrates f(x) in the interval [a, b] using N intervals.

	std::vector<int> idx(N+1) ;
	std::vector<double> x(N+1) ;
	std::iota(idx.begin(), idx.end(), 0) ;
	double dx = (b - a) / ( (double) N ) ;  // step


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


double numint::rsum (double a, double b, const int N,
		     double f (const double&) )
{
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


double numint::trap (double a, double b, const int N,
		     double f (const double&) )
{
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
