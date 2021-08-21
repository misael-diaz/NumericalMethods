/*
 * Applied Numerical Analysis                                 July 20, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: module-nonlinear-solvers.cpp
 *
 * Synopsis:
 * Implements (some) nonlinear equation solvers.
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
#include <string>
#include <iostream>
#include <stdexcept>
#include <functional>
#define absval(x) ((x < 0.)? -x: x)
export module nonlinear_solvers ;

// constants
const int MAX_ITER = 256 ;
const double TOL = 1.0e-12 ;
const bool VERBOSE = false ;

export namespace nlsolver {

	struct config {	// config[uration] struct
		config(): tol(TOL), max_iter(MAX_ITER), verbose(VERBOSE) {}
		config(double t): tol(t),   max_iter(MAX_ITER),
	                          verbose(VERBOSE) {}
		config(int    i): tol(TOL), max_iter(i),
	                          verbose(VERBOSE) {}
		config(bool   v): tol(TOL), max_iter(MAX_ITER),
	                          verbose(v) {}
		config(double t, int  i): tol(t), max_iter(i),
	                                  verbose(VERBOSE) {}
		config(double t, bool v): tol(t), max_iter(MAX_ITER),
	                                  verbose(v) {}
		config(int    i, bool v): tol(TOL), max_iter(i),
	                                  verbose(v) {}
		config(double t, int i, bool v): tol(t), max_iter(i),
	                                         verbose(v) {}
		double tol ;
		int max_iter ;
		bool verbose ;
	} ;

	double fzero(double,double,std::function<double(const double&)>&,
		     config ) ;
}


// declarations (prototypes)
void report (const int& n, const std::string& nm,
             const nlsolver::config& opt) ;
void check_bounds ( double& lb, double& ub ) ;

void check_bracket ( const double&, const double&, const std::string&,
	       	     std::function<double(const double&)>& ) ;

double shift ( double&, double&, double&, 
	       std::function<double(const double&)>& ) ;



// implementations
double nlsolver::fzero   ( double lb, double ub, 
		           std::function<double(const double&)>& f,
	       	           config opt = config() )
{	// Shifter Method

	int n = 0 ;
	double xn, fn ;
	std::string nm = "Shifter" ;

	check_bounds  (lb, ub) ;
	check_bracket (lb, ub, nm, f) ;

        do fn = shift (lb, ub, xn, f) ;
        while (++n != opt.max_iter && fn > opt.tol) ;

	report (n, nm, opt) ;
	return xn ;

}

void check_bounds ( double& lb, double& ub ) {
	// ensures the (given) interval is properly bounded [lower, upper]
	double low = ub ;
	if (lb > ub) {
		ub = lb ;
		lb = low ;
	}
}


void check_bracket ( const double& lb, const double& ub, 
		     const std::string& nm, 
		     std::function<double(const double&)>& f )
{
	// complains if there's no root in given interval
	if ( f(lb) * f(ub) > 0. ) {
		std::cout << nm << " Method:" << std::endl ;
		std::cout << "no root enclosed in " 
		          << "[" << lb << ", " << ub << "]" << std::endl ;
		throw std::runtime_error("no root in given interval") ;
	}
	
}


void report (const int& n, const std::string& nm,
             const nlsolver::config& opt)
{
	// reports if the method has been successful
	if (n != opt.max_iter) {
		if (opt.verbose) {
			std::cout << nm << " Method:" << std::endl ;
			std::cout << "solution found in " << n << " "
                                  << "iterations" << std::endl ;
		}
	} else {
		std::string errMSG = nm + " method requires more "
			"iterations for convergence" ;
		throw std::runtime_error (errMSG) ;
	}
}


double shift ( double& lb, double& ub, double& xn,
	       std::function<double(const double&)>& f )
{	// like bisector but uses the step (presumably) closer to the root
	double fn ;
	double xb = 0.5 * (lb + ub) ;
        double xf = ( lb * f(ub) - ub * f(lb) ) / ( f(ub) - f(lb) ) ;

	if ( absval(f(xb)) < absval(f(xf)) )
		xn = xb ;
	else
		xn = xf ;

	fn = f(xn) ;

	if (f(lb) * fn < 0.)
		ub = xn ;
	else
		lb = xn ;

	return absval(fn) ;
}
