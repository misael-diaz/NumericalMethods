/*
 * Applied Numerical Analysis                                August 10, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: module-ode-solvers.cpp
 *
 * Synopsis:
 * Implements (some) Ordinary Differential Equation ODE Solvers using
 * C++ facilities.
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
#include <tuple>
#include <vector>
#include <numeric>
#include <algorithm>
export module odes ;


export namespace ode {

	// Euler's Method
	std::tuple< std::vector<double>, std::vector<double> >
	Euler (const double&, const double&, const double&, 
	       const int&, double f(const double&, const double&) ) ;
	
	// second-order Runge-Kutta Method
	std::tuple< std::vector<double>, std::vector<double> >
	EulerRK2 (const double&, const double&, const double&, 
	          const int&, double f(const double&, const double&) ) ;
}

// declarations
std::vector<double>& linspace ( std::vector<double>&, const double&, 
		                const double&, const int& ) ;

// implementations
std::tuple< std::vector<double>, std::vector<double> >
ode::Euler (const double& ti, const double& tf, const double& yi, 
	    const int& N, double f (const double&, const double&) )
{	// implements Euler's Method

	std::vector<double> t ;
	t = linspace(t, ti, tf, N + 1) ;

	std::vector<double> y(N + 1) ;
	double dt = (tf - ti) / ( (double) N ) ;
	typedef std::vector<double>::size_type size ;


	y[0] = yi ;
	for (size i = 0 ; i != t.size() ; ++i) {
		y[i + 1] = y[i] + dt * f(t[i], y[i]) ;
	}


	return std::tuple (t, y) ;
}


std::tuple< std::vector<double>, std::vector<double> >
ode::EulerRK2 (const double& ti, const double& tf, const double& yi, 
	    const int& N, double f (const double&, const double&) )
{	// implements an Euler-based, second-order, Runge-Kutta Method

	double K1, K2 ;
	std::vector<double> t ;
	t = linspace(t, ti, tf, N + 1) ;

	std::vector<double> y(N + 1) ;
	double dt = (tf - ti) / ( (double) N ) ;
	typedef std::vector<double>::size_type size ;


	y[0] = yi ;
	for (size i = 0 ; i != t.size() ; ++i) {
		K1 = f(t[i], y[i]) ;
		K2 = f(t[i] + dt, y[i] + dt * K1) ;
		y[i + 1] = y[i] + 0.5 * dt * (K1 + K2) ;
	}


	return std::tuple (t, y) ;
}


std::vector<double>& linspace (std::vector<double>& t, const double& ti, 
		               const double& tf, const int& numel)
{	// implements a numpy-like linspace method

	t.clear() ;
	std::vector<int> idx (numel) ;
	std::iota(idx.begin(), idx.end(), 0) ;
	double dt = (tf - ti) / ( (double) (numel - 1) ) ;  // time-step

	// transformational lambda function
	auto ft = [&ti, &dt](const int& i) -> double {
		return (ti + ( (double) i) * dt) ;
	} ;

	// transforms index into a vector t = [ti, tf] of numel elements
	std::transform(idx.begin(), idx.end(), std::back_inserter(t), ft) ;
	return t ;
}
