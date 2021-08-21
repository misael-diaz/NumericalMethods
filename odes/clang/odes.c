/*
 * Applied Numerical Analysis				      July 27, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: odes.c
 *
 * Synopsis:
 * Implements (some) numerical integration methods.
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


#include "odes.h"


// implementations
static double* ode_allocArray (double *x, const int numel) {
	// allocates memory for a first-rank array of doubles

	x = (double*) malloc ( numel * sizeof(double) ) ;

	if (x == NULL) {
		fprintf(stderr, "\n\n") ;
		fprintf(stderr, "ODE: failed to allocate array\n") ;
		fprintf(stderr, "aborting execution ... \n") ;
		fprintf(stderr, "\n\n") ;
		exit(EXIT_FAILURE) ;
	}

	return x ;
}


static double* linspace (double *t, double ti, double tf, const int numel)
{
        // implements a numpy-like linspace method

	double dt = (tf - ti) / ( (double) (numel - 1) ) ;

	t = ode_allocArray (t, numel) ;
        for (int i = 0 ; i != numel ; ++i)
                t[i] = ti + ( (double) i ) * dt ;

        return t ;
}


// method interfaces
double** Euler ( double **odesol, double ti, double tf, double yi,
		 const int N, double f(double t, double y, double* prms),
		 void *vprms )
{	// applies Euler's method to integrate the ODE

	// unpacks parameters
	iODE_solverParams *params = vprms ;
	double *prms = params -> prms ;

	double *t = NULL, *y = NULL ;			// t, y(t)
	double dt = (tf - ti) / ( (double) N ) ;	// time-step

	t = odesol[0] = linspace (odesol[0], ti, tf, N + 1) ;
	y = odesol[1] = ode_allocArray (odesol[1], N + 1) ;

	y[0] = yi ;
	for (int i = 0 ; i != N ; ++i)
		y[i + 1] = y[i] + dt * f(t[i], y[i], prms);

	return odesol ;
}


double** iEuler ( double **odesol, double ti, double tf, double yi,
                  const int N, double f(double t, double y, double *prms),
		  void *vprms )
{	/* Applies the implicit Euler's method */

	// unpacks the parameters for the nonlinear and ode solvers
	iODE_solverParams *params = vprms ;
	double (*objf) (double, void*) = params -> objf ;
	double *prms = params -> prms ;

	double K1, K2 ;
	double y_lb, y_ub ;
	double *t = NULL, *y = NULL ;
	double dt = (tf - ti) / ( (double) N ) ;

	t = odesol[0] = linspace (odesol[0], ti, tf, N + 1) ;
	y = odesol[1] = ode_allocArray (odesol[1], N + 1) ;

	y[0] = yi ;
	for (int i = 0 ; i != N ; ++i) {

		// bounds the solution y[i + 1] from below and above
		K1 = f(t[i], y[i], prms) ;
		K2 = f(t[i] + dt, y[i] + dt * K1, prms) ;
		y_lb = y[i] + dt * K1 ;
		y_ub = y[i] + dt * K2 ;

		// packs parameters for the nonlinear solver
		prms[0] = dt ;
		prms[1] = y[i] ;
		prms[2] = t[i + 1] ;

		/* solves for y[i + 1] iteratively */
		y[i + 1] = fzero ( y_lb, y_ub, objf, vprms ) ;
	}

	return odesol ;
}


double** EulerRK2 (double **odesol, double ti, double tf, double yi,
		   const int N, double f(double t, double y, double *prms),
		   void *vprms)
{	// applies an Euler-based, second-order, Runge-Kutta method

	// unpacks parameters
	iODE_solverParams *params = vprms ;
	double *prms = params -> prms ;

	double K1, K2 ;
	double *t = NULL, *y = NULL ;
	double dt = (tf - ti) / ( (double) N ) ;

	t = odesol[0] = linspace (odesol[0], ti, tf, N + 1) ;
	y = odesol[1] = ode_allocArray (odesol[1], N + 1) ;

	y[0] = yi ;
	for (int i = 0 ; i != N ; ++i) {
		K1 = f(t[i], y[i], prms) ;
		K2 = f(t[i] + dt, y[i] + K1 * dt, prms) ;
		y[i + 1] = y[i] + 0.5 * dt * (K1 + K2) ;
	}

	return odesol ;
}
