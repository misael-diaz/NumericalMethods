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
#define absval(x) ((x < 0.)? -x: x)
export module nonlinear_solvers ;

// constants
const int MAX_ITER = 100 ;
const double TOL = 1.0e-8 ;

export namespace nlsolver {

	struct config {	// config[uration] struct
		config(): tol(TOL), max_iter(MAX_ITER) {}
		config(double t, int i): tol(t), max_iter(i) {}
		double tol ;
		int max_iter ;
	} ;

	double bisect ( double, double, double f(const double&), config ) ;
	double regfal ( double, double, double f(const double&), config ) ;
	double shifter( double, double, double f(const double&), config ) ;
}


// declarations (prototypes)
void report (const int& max_iter, const int& n, const std::string& nm) ;
void check_bounds ( double& lb, double& ub ) ;

void check_bracket ( const double& lb, const double& ub,
		     const std::string& nm, double f(const double&) ) ;

double bisector ( double&, double&, double&, double f(const double&) ) ;
double interp   ( double&, double&, double&, double f(const double&) ) ;
double shift    ( double&, double&, double&, double f(const double&) ) ;



// implementations
double nlsolver::bisect ( double lb, double ub,
		          double f(const double&), config opt = config() )
{	// Bisection Method

	int n = 0 ;
	double xm, fm ;
	std::string nm = "Bisection" ;

	check_bounds  (lb, ub) ;
	check_bracket (lb, ub, nm, f) ;

        do fm = bisector (lb, ub, xm, f) ;
        while (++n != opt.max_iter && fm > opt.tol) ;

	report (opt.max_iter, n, nm) ;
	return xm ;

}


double nlsolver::regfal ( double lb, double ub,
			  double f(const double&), config opt = config() )
{	// Regula Falsi Method

	int n = 0 ;
	double xn, fn ;
	std::string nm = "Regula-Falsi" ;

	check_bounds  (lb, ub) ;
	check_bracket (lb, ub, nm, f) ;

        do fn = interp (lb, ub, xn, f) ;
        while (++n != opt.max_iter && fn > opt.tol) ;

	report (opt.max_iter, n, nm) ;
	return xn ;

}


double nlsolver::shifter ( double lb, double ub,
		           double f(const double&), config opt = config() )
{	// Shifter Method

	int n = 0 ;
	double xn, fn ;
	std::string nm = "Shifter" ;

	check_bounds  (lb, ub) ;
	check_bracket (lb, ub, nm, f) ;

        do fn = shift (lb, ub, xn, f) ;
        while (++n != opt.max_iter && fn > opt.tol) ;

	report (opt.max_iter, n, nm) ;
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
		     const std::string& nm, double f(const double&) )
{
	// complains if there's no root in given interval
	if ( f(lb) * f(ub) > 0. ) {
		std::cout << nm << " Method:" << std::endl ;
		std::cout << "no root enclosed in " 
		          << "[" << lb << ", " << ub << "]" << std::endl ;
		throw std::runtime_error("no root in given interval") ;
	}
	
}


void report (const int& max_iter, const int& n, const std::string& nm) {
	// reports if the method has been successful
	if (n != max_iter) {
		std::cout << nm << " Method:" << std::endl ;
		std::cout << "solution found in " << n << " "
			  << "iterations" << std::endl ;
	} else {
		std::string errMSG = nm + " method requires more "
			"iterations for convergence" ;
		throw std::runtime_error (errMSG) ;
	}
}


double bisector ( double& lb, double& ub, double& xm, 
		  double f(const double&) )
{

/*
 * Synopsis:
 * Approximates the root of the nonlinear equation f(x) with the middle
 * value, xm = (lb + ub) / 2, where `lb' and `ub' are the lower
 * and upper bounds, respectively, of the bracketing interval [lb, ub].
 * As a side effect it returns |f(xm)|, where |*| = abs(*), and `*' is
 * a wildcard.
 *
 * input/output:
 * lb		lower bound, intent(inout)
 * ub           upper bound, intent(inout)
 * xm           middle value, intent(out)
 *
 * NOTE:
 * f(x) is a function pointer (last argument) bound to the user-defined
 * function representing the nonlinear equation to be solved.
 *
 */
	
        xm = 0.5 * (lb + ub) ;
	double fm = f(xm) ;

	/* chooses the (root) bracketing interval */
	if (f(lb) * fm < 0.) 	
		ub = xm ;
	else
		lb = xm ;

	return absval(fm) ;
}


double interp ( double& lb, double& ub, double& xn,
		double f(const double&) )
{	// like bisector but uses interpolation to approximate the root
        xn = ( lb * f(ub) - ub * f(lb) ) / ( f(ub) - f(lb) ) ;
	double fn = f(xn) ;

	if (f(lb) * fn < 0.)
		ub = xn ;
	else
		lb = xn ;

	return absval(fn) ;
}


double shift ( double& lb, double& ub, double& xn,
	       double f(const double&) )
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


/*
 * TODO:
 * [x] Define default values for the tolerance and maximum number of
 *     iterations. The user may wish to override these parameters so these
 *     must be included in the argument lists. Use the constants as default
 *     values for these parameters.
 * [x] report function should display the name of the numerical method
 *
 */
