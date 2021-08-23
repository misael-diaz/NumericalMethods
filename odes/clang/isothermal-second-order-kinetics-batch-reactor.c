/*
 * Applied Numerical Analysis				    August 23, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: isothermal-second-order-kinetics-batch-reactor.c
 *
 * Synopsis:
 * Obtains the transient response of a second-order nonlinear Ordinary
 * Differential Equation ODE, which describes the depletion of a chemical
 * species by a chemical reaction of second-order kinetics carried out in a
 * isothermal batch reactor. It's assumed that the chemical reaction takes
 * place in a liquid so that the volume of the reactor can be regarded as
 * constant.
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

// MACROS for the effective reaction rate constant
#define BETA 1.0

// prototypes
double fsol    (double t);				// exact solution
double objf    (double yn, void *vprms);		// objective fun
double odefun  (double t, double y, double* prms);	// RHS ODE, f(t, y)
void   write   (char*, const int, double**);		// writes to file
void   display (const int, double**);			// writes to stdout

int main() {
	// Solves first-order ODEs using Euler's implicit Method

	const int N = 255 ;		// number of intervals
	const int numel = N + 1 ;	// number of elements in time array
	double ti = 0.0, tf = 10.0 ;	// initial and final times
	double yi = 1.0 ;		// initial value, y = y(t = 0) = 1

	// allocates and packs the parameters for implicit solver
	double prms[] = {0., 0., 0., BETA};
	iODE_solverParams *iSolverParams =
	    (iODE_solverParams*) malloc ( sizeof(iODE_solverParams) );
	iSolverParams -> objf   = objf ;
	iSolverParams -> odefun = odefun ;
	iSolverParams -> prms   = prms ;

	// creates placeholders for the numerical solution, odesol = [t, y]
	double **oderet = NULL ;
	double *odesol[][2] = { {NULL, NULL}, {NULL, NULL} };
	// solves the ODE numerically and writes results to a data file
	oderet = iEuler  (odesol[0], ti, tf, yi, N, odefun, iSolverParams);
	oderet = EulerRK2(odesol[1], ti, tf, yi, N, odefun, iSolverParams);
	// writes numerical results
	char filename[] = "output/cheme/kinetics/iEuler.dat" ;
	write  (filename, numel, odesol[0]);
	strcpy (filename, "output/cheme/kinetics/EulRK2.dat");
	write  (filename, numel, odesol[1]);
	display (numel, oderet);

	// frees memory buffers
	free (iSolverParams);
	free (odesol[0][0]);
	free (odesol[0][1]);
	free (odesol[1][0]);
	free (odesol[1][1]);
	return 0 ;
}


// implementations
double objf (double yn, void *vprms)
{	// Synopsis:
	// Objective function for the nonlinear solver invoked by the
	// implicit Euler method.

	/* unpacks parameters needed to evaluate the objective function */
        iODE_solverParams *params = vprms ;
        double (*fp) (double, double, double*) = params -> odefun ;
        double *prms = params -> prms ;
        double dt   = prms[0];
        double yi   = prms[1];
        double tn   = prms[2];
//      double beta = prms[3];
        return ( yn - yi - dt * fp(tn, yn, prms) );
}


double odefun (double t, double y, double* prms)
{	// Synopsis: Right-hand side of the ODE.

//      double dt   = prms[0];
//      double yi   = prms[1];
//      double tn   = prms[2];
	double beta = prms[3];
	return (-beta * y * y);
}


double fsol (double t) {
/* 
 * Synopsis:
 * Analytic solution of the first-order nonlinear ODE:
 *  		y' = -BETA * y**2 		y(t = 0) = 1,
 * where BETA is the effective rate constant, and `y' is the concentration
 * of the chemical species. Note that the concentration has been made
 * non-dimensional by the initial concentration (Ca ~ C{a,0}).
 *
 */
	return ( 1.0 / (1.0 + BETA * t) );
}


void write (char filename[], const int numel, double **odesol) {
	// writes the numerical solution and error to a data file
	double *t = odesol[0];
	double *y = odesol[1];

	FILE *pFile = NULL ;
	pFile = fopen(filename, "w");

	if (pFile == NULL) {
		fprintf(stderr, "I/O Error: %s\n", filename);
		fprintf(stderr, "aborting execution ... \n");
		exit(EXIT_FAILURE);
	}

	char fmt[] = "%23.15e %23.15e %23.15e\n" ;
	for (int i = 0 ; i != numel ; ++i)
		fprintf(pFile, fmt, t[i], y[i], absval((fsol(t[i])-y[i])));
	fclose(pFile);
}


void display (const int numel, double **odesol) {
	// Displays the numerical and exact solutions, and the absolute 
	// error to the standard output.

	double err ;
	double *t = odesol[0];
	double *y = odesol[1];

	char fmt[] = "%23.15e %23.15e %23.15e %23.15e\n" ;
	for (int i = 0 ; i != numel ; ++i) {
		err = absval( (fsol(t[i]) - y[i]) );
		fprintf(stdout, fmt, t[i], y[i], fsol(t[i]), err);
	}
}
