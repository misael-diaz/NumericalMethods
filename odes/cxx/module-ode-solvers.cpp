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
 * [2] www.cplusplus.com/reference/tuple/tuple
 * [3] www.geeksforgeeks.org/tuples-in-c/
 *
 */


module ;
#include <tuple>
#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>
export module odes ;


export namespace ode {

	// Euler's Method
	std::tuple< std::vector<double>, std::vector<double> >&
	Euler ( std::tuple< std::vector<double>, std::vector<double> >&,
		const double&, const double&, const double&, const int&,
		std::function< double(const double&, const double&) >& ) ;
	
	// second-order Runge-Kutta Method
	std::tuple< std::vector<double>, std::vector<double> >&
	RK2 ( std::tuple< std::vector<double>, std::vector<double> >&,
	      const double&, const double&, const double&, const int&,
	      std::function< double(const double&, const double&) >& ) ;

}

// declarations
std::vector<double>& linspace ( std::vector<double>&, const double&,
		                const double&, const int& ) ;

// implementations
std::tuple< std::vector<double>, std::vector<double> >&
ode::Euler (std::tuple< std::vector<double>, std::vector<double> >& odesol,
	    const double& ti, const double& tf, const double& yi,
	    const int& N,
	    std::function< double(const double&, const double&) >& f )
{	// possible implementation of Euler's explicit method

	double dt = (tf - ti) / ( (double) N );		// time-step, dt
	std::vector<double>& t = std::get<0> (odesol);	// time, t
	std::vector<double>& y = std::get<1> (odesol);	// y(t)

	t = linspace (t, ti, tf, N + 1);

	y.clear();		// clears (existing) data
	y.reserve(N + 1);	// preallocates vector for speed
	y[0] = yi ;
	typedef std::vector<double>::size_type size ;
	for (size i = 0 ; i != (t.size() - 1) ; ++i) {
		y[i + 1] = y[i] + dt * f(t[i], y[i]);
	}

	return odesol ;
}


std::tuple< std::vector<double>, std::vector<double> >&
ode::RK2 ( std::tuple< std::vector<double>, std::vector<double> >& odesol,
	   const double& ti, const double& tf, const double& yi,
	   const int& N,
	   std::function< double(const double&, const double&) >& f )
{	// implements an Euler-based, second-order, Runge-Kutta Method

	double K1, K2 ;					// slopes
	double dt = (tf - ti) / ( (double) N );		// time-step, dt
	std::vector<double>& t = std::get<0> (odesol);	// time, t
	std::vector<double>& y = std::get<1> (odesol);	// y(t)

	t = linspace (t, ti, tf, N + 1);

	y.clear();
	y.reserve(N + 1);
	y[0] = yi ;
	typedef std::vector<double>::size_type size ;
	for (size i = 0 ; i != (t.size() - 1) ; ++i) {
		K1 = f(t[i], y[i]);
		K2 = f(t[i] + dt, y[i] + dt * K1);
		y[i + 1] = y[i] + 0.5 * dt * (K1 + K2);
	}

	return odesol ;
}


std::vector<double>& linspace (std::vector<double>& t, const double& ti,
		               const double& tf, const int& numel)
{	// implements a numpy-like linspace method
	t.clear();				// clears (existing) data
	t.reserve (numel);			// prellocates for speed
	std::vector<int> idx (numel);		// index vector [0, numel)
	std::iota (idx.begin(), idx.end(), 0);
	double dt = (tf - ti) / ( (double) (numel - 1) );
	auto f = [&](const int& i) { return (ti + ( (double) i) * dt); };
	std::transform (idx.begin(), idx.end(), std::back_inserter(t), f);
	return t ;
}
