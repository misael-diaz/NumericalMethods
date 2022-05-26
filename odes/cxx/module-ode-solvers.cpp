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
import nonlinear_solvers ;
export module odes ;


export namespace ode {

	// Euler's Method
	std::tuple< std::vector<double>, std::vector<double> >&
	Euler ( std::tuple< std::vector<double>, std::vector<double> >&,
		double, double, double, int,
		std::function< double(double, double) >& ) ;

	// implicit Euler's method
	std::tuple< std::vector<double>, std::vector<double> >&
	iEuler ( std::tuple< std::vector<double>, std::vector<double> >&,
		 double, double, double, int,
		 std::function< double(double, double) >& ) ;
	
	// second-order Runge-Kutta Method
	std::tuple< std::vector<double>, std::vector<double> >&
	RK2 ( std::tuple< std::vector<double>, std::vector<double> >&,
	      double, double, double, int,
	      std::function< double(double, double) >& ) ;

}

// declarations
std::vector<double>& linspace (std::vector<double>&, double, double, int);

// implementations
std::tuple< std::vector<double>, std::vector<double> >&
ode::Euler (std::tuple< std::vector<double>, std::vector<double> >& odesol,
	    double ti, double tf, double yi, int N,
	    std::function< double(double, double) >& f )
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
ode::iEuler(std::tuple< std::vector<double>, std::vector<double> >& odesol,
	    double ti, double tf, double yi, int N,
	    std::function< double(double, double) >& f )
{	// possible implementation of Euler's implicit method

	double K1, K2 ;					// slopes
	double y_lb, y_ub ;				// bounds
	double dt = (tf - ti) / ( (double) N );		// time-step, dt
	std::vector<double>& t = std::get<0> (odesol);	// time, t
	std::vector<double>& y = std::get<1> (odesol);	// y(t)
	typedef std::vector<double>::size_type size ;	// typedef
	size i ;					// counter

	std::function<double(const double&)> objf = [&](const double& yn) {
		// objective function supplied to nonlinear solver fzero
		return (yn - y[i] - dt * f(t[i+1], yn));
	};

	t = linspace (t, ti, tf, N + 1);

	y.clear();		// clears (existing) data
	y.reserve(N + 1);	// preallocates vector for speed
	y[0] = yi ;
	for (i = 0 ; i != (t.size() - 1) ; ++i) {
		// brackets y[i+1] from below and above
		K1 = f(t[i], y[i]);
		K2 = f(t[i] + dt, y[i] + dt * K1);
		y_lb = y[i] + dt * K1 ;
		y_ub = y[i] + dt * K2 ;
		// obtains y[i+1] iteratively
		y[i + 1] = nlsolver::fzero (y_lb, y_ub, objf);
	}

	return odesol ;
}


std::tuple< std::vector<double>, std::vector<double> >&
ode::RK2 ( std::tuple< std::vector<double>, std::vector<double> >& odesol,
	   double ti, double tf, double yi, int N,
	   std::function< double(double, double) >& f )
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


std::vector<double>& linspace (std::vector<double>& t, double start,
			       double end, int numel)
{	// implements a numpy-like linspace method
	t.clear();				// clears (existing) data
	t.reserve (numel);			// prellocates for speed

	// calculates the time-step
	double dt = (end - start) / ( (double) (numel - 1) );

	// computes the time vector
	for (int i = 0; i != numel; ++i)
	{
		double ti = (start + ( (double) i ) * dt);	// t[i]
		t.push_back (ti);
	}

	return t ;
}

/*
 * Comments on the implementation of Euler's implicit method:
 * The objective function is defined as an instance of std::function
 * to be able to supply a lambda function with captures to the nonlinear
 * solver fzero. The fzero method has been adapted to receive an instance
 * of std::function class by reference instead of a function pointer.
 *
 * Note that the declaration of the loop index `i` has been moved above the
 * definition of the lambda for capturing.
 *
 */
