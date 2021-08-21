/*
 * Applied Numerical Analysis				    August 09, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: test.c
 *
 * Synopsis:
 * Tests some Ordinary Differential Equation ODE solvers.
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "odes.h"

// MACROS
#define RATE 1.0

// prototypes
double fsol (double t) ;		// analytic solution
void display (const int, double**) ;
void write   (char*, const int, double**) ;

// user-defined functions needed by the implicit Euler's method
double objf (double yn, void *vprms) ;
double odefun (double t, double y, double* prms) ;

int main() {
	// Solves first-order ODEs using Euler's and Runge-Kutta Methods


	const int N = 255 ;		// number of intervals
	const int numel = N + 1 ;	// number of elements in time array
	double ti = 0.0, tf = 5.0 ;	// initial and final times
	double yi = 1.0 ;		// initial value, y = y(t = ti)


	// allocates and packs the parameters for implicit solver
	double prms[] = {0., 0., 0., RATE} ;
	iODE_solverParams *iSolverParams =
	    (iODE_solverParams*) malloc ( sizeof(iODE_solverParams) ) ;
	iSolverParams -> objf   = objf ;
	iSolverParams -> odefun = odefun ;
	iSolverParams -> prms   = prms ;


	// creates placeholders for the numerical solutions odesol = [t, y]
	double **oderet = NULL ;	// as odesol in the Python impl*
	double *odesol[][2] = { {NULL, NULL}, {NULL, NULL}, {NULL, NULL} };


	/* solves the ode numerically via the specified methods */
	oderet = Euler   (odesol[0], ti, tf, yi, N, odefun, iSolverParams);
	oderet = EulerRK2(odesol[1], ti, tf, yi, N, odefun, iSolverParams);
	oderet = iEuler  (odesol[2], ti, tf, yi, N, odefun, iSolverParams);


	// exports numerical solutions to respective data files
	char filename[] = "output/Euler.dat" ;
	write  (filename, numel, odesol[0]) ;
	strcpy (filename, "output/EuRK2.dat") ;
	write  (filename, numel, odesol[1]) ;
	strcpy (filename, "output/iEulr.dat") ;
	write  (filename, numel, odesol[2]) ;
//	write  (filename, numel, oderet) ;	// alternate way of writing
	display (numel, oderet) ;


	// frees memory buffers
	free (iSolverParams) ;
	free (odesol[0][0]) ;
	free (odesol[0][1]) ;
	free (odesol[1][0]) ;
	free (odesol[1][1]) ;
	free (odesol[2][0]) ;
	free (odesol[2][1]) ;
	return 0 ;
}


double objf (double yn, void *vprms)
{	// Synopsis:
	// Objective function for the nonlinear solver invoked by the
	// implicit Euler method.

	/* unpacks parameters needed to evaluate the objective function */
        iODE_solverParams *params = vprms ;
        double (*fp) (double, double, double*) = params -> odefun;//f(t, y)
        double *prms = params -> prms ;
        double dt = prms[0] ;
        double yi = prms[1] ;
        double tn = prms[2] ;
//      double k  = prms[3] ;   // unused parameter, needed by odefun
        return ( yn - yi - dt * fp(tn, yn, prms) ) ;
}


double odefun (double t, double y, double* prms)
{	// Synopsis: RHS of the ODE used by the implicit Euler method.

/*	unpacks parameters used by the ODE function f(t, y) 	*/
//      double dt = prms[0] ;
//      double yi = prms[1] ;
//      double tn = prms[2] ;
	double k  = prms[3] ;
	return (-k * y) ;
}


double fsol (double t) {
	// analytic solution
	double k = RATE ;
	return exp(-k * t) ;
}


void write (char filename[], const int numel, double **odesol) {
	// writes the numerical solution to a data file
	double err ;
	double *t = odesol[0] ;
	double *y = odesol[1] ;

	FILE *pFile = NULL ;
	pFile = fopen(filename, "w") ;

	if (pFile == NULL) {
		fprintf(stderr, "I/O Error: %s\n", filename) ;
		fprintf(stderr, "aborting execution ... \n") ;
		exit(EXIT_FAILURE) ;
	}

	char fmt[] = "%23.15e %23.15e %23.15e\n" ;	// format string
	for (int i = 0 ; i != numel ; ++i) {
		err = y[i] - fsol(t[i]);
		err = absval (err) ;
		fprintf(pFile, fmt, t[i], y[i], err);
	}
	fclose(pFile) ;
}


void display (const int numel, double **odesol) {
	// displays the numerical solution to stdout

	double err ;			// absolute error
	double *t = odesol[0] ;
	double *y = odesol[1] ;

	for (int i = 0 ; i != numel ; ++i) {
		err = y[i] - fsol(t[i]);
		err = absval (err) ;
		fprintf(stdout, "%23.15e \t %23.15e \t %23.15e \n",
			t[i], y[i], err) ;
	}
}


/*
 * Comments:
 * My aim was to write a code that looks similar to the code written
 * in Python. That's why the numerical solvers return an array that
 * can be used for writing if the users wish to do so.
 *
 */
