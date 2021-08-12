#ifndef GUARD_NLSOLVERS_H
#define GUARD_NLSOLVERS_H

/*
 * source: nlsolvers.h
 * author: misael-diaz
 * date:   2021/07/26
 *
 * Synopsis:
 * Header file for the nonlinear equation solvers.
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

// Math MACROS
#define absval(x) ((x < 0.)? -x: x)

// constants (defaults)
#define MAX_ITER 100
#define TOL 1.0e-8

typedef struct {
	double tol ;	// tolerance
	int max_iter ;	// maximum number of iterations
} nls_conf ;		// nonlinear-solver configuration struct

typedef struct {
	nls_conf opts ;
} opt_args ;		/* optional arguments (companion) struct */


// declarations (prototypes):

// shadow functions (handle initialization of optional arguments)
double bisect_varg (double, double, double f(const double), opt_args) ;
double regfal_varg (double, double, double f(const double), opt_args) ;
double shifter_varg(double, double, double f(const double), opt_args) ;

// base functions (where the respective method is actually implemented)
double bisect_base (double, double, double f(const double), double, int) ;
double regfal_base (double, double, double f(const double), double, int) ;
double shifter_base(double, double, double f(const double), double, int) ;

// wrapper MACROS
#define bisect(lb, ub, f, ...)\
	bisect_varg(lb, ub, f, (opt_args){__VA_ARGS__})
#define regfal(lb, ub, f, ...)\
	regfal_varg(lb, ub, f, (opt_args){__VA_ARGS__})
#define shifter(lb, ub, f, ...)\
	shifter_varg(lb, ub, f, (opt_args){__VA_ARGS__})

// bracketing approach
double bisector ( double*, double*, double*, double f(const double) ) ;
double interp   ( double*, double*, double*, double f(const double) ) ;
double shift    ( double*, double*, double*, double f(const double) ) ;

// utility functions
void report        ( const int max_iter, const int n, char nm[] ) ;
void check_bracket ( double, double, double f(const double), char nm[] ) ;
void check_bounds  ( double*, double* ) ;
#endif




/*
 * References:
 * [*] initialization of optional arguments using Variadic MACROS VA_ARGS:
 *
 *     stackoverflow.com/questions/1472138/c-default-arguments
 *     gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
 *
 *     The implementation follows the scheme proposed by the
 *     stackoverflow user bk: define a companion struct, a shadow
 *     function that initializes optional arguments for the base function
 *     (numerical method in this case), and a variadic MACRO to define a
 *     wrapper function that hides the implementation from the user.
 *
 *     The GCC article shows how to define variadic MACROS for
 *     function with named (formal) arguments as in this case.
 *
 *
 * [*] initialization of configuration struct via MACROS: (unused)
 *     (stackoverflow.com/questions/13716913/
 *      default-value-for-struct-member-in-c)
 *
 *     Example:
 *     #define init_config(x) config x = {.tol = TOL, .max_iter = MAX_ITER}
 *
 */
