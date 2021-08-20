/*
 * Applied Numerical Analysis				    August 20, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: ramp.c
 *
 * Synopsis:
 * Obtains the transient response of a first-order Ordinary Differential
 * Equation ODE subject to a ramp input.
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

// MACROS for the rate and external forcing constants, respectively
#define RATE 1.0
#define FEXT 1.0

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
	double ti = 0.0, tf = 5.0 ;	// initial and final times
	double yi = 0.0 ;		// initial value, y = y(t = 0) = 0

	// allocates and packs the parameters for implicit solver
	double prms[] = {0., 0., 0., RATE, FEXT};
	iODE_solverParams *iSolverParams =
	    (iODE_solverParams*) malloc ( sizeof(iODE_solverParams) );
	iSolverParams -> objf   = objf ;
	iSolverParams -> odefun = odefun ;
	iSolverParams -> prms   = prms ;

	// creates placeholders for the numerical solution, odesol = [t, y]
	double **oderet = NULL ;
	double *odesol[] = {NULL, NULL};
	// solves the ODE numerically and writes results to a data file
	oderet = iEuler (odesol, ti, tf, yi, N, odefun, iSolverParams);
	char filename[] = "output/ramp/iEuler.dat" ;
	write  (filename, numel, oderet);
	display (numel, oderet);

	// frees memory buffers
	free (iSolverParams);
	free (odesol[0]);
	free (odesol[1]);
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
        double dt = prms[0];
        double yi = prms[1];
        double tn = prms[2];
//      double k  = prms[3];	// rate constant
//      double b  = prms[4];	// external forcing constant
        return ( yn - yi - dt * fp(tn, yn, prms) );
}


double odefun (double t, double y, double* prms)
{	// Synopsis: Right-hand side of the ODE with ramp input.

//      double dt = prms[0];
//      double yi = prms[1];
//      double tn = prms[2];
	double k  = prms[3];
	double b  = prms[4];
	return (b * t - k * y);
}


double fsol (double t) {
/* 
 * Synopsis:
 * Analytic solution of the first-order ODE subject to a ramp input:
 *  		y' + k * y = b * t, 		y(t = 0) = 0,
 * where k is the rate constant, and b is the external forcing constant.
 *
 */
	double k = RATE;
	double b = FEXT;
	return (b / (k * k) * (exp(-k * t) - 1.0) + b / k * t);
}


void write (char filename[], const int numel, double **odesol) {
	// writes the numerical solution to a data file
	double *t = odesol[0];
	double *y = odesol[1];

	FILE *pFile = NULL ;
	pFile = fopen(filename, "w");

	if (pFile == NULL) {
		fprintf(stderr, "I/O Error: %s\n", filename);
		fprintf(stderr, "aborting execution ... \n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0 ; i != numel ; ++i)
		fprintf(pFile, "%23.15e \t %23.15e \n", t[i], y[i]);
	fclose(pFile);
}


void display (const int numel, double **odesol) {
	// displays the numerical solution and absolute error to stdout

	double err ;
	double *t = odesol[0];
	double *y = odesol[1];

	for (int i = 0 ; i != numel ; ++i) {
		err = absval( y[i] - fsol(t[i]) );
		fprintf(stdout, "%23.15e \t %23.15e \t %23.15e \n",
			t[i], y[i], err);
	}
}
