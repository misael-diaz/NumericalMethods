/*
 * Applied Numerical Analysis				    August 25, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: impulse.c
 *
 * Synopsis:
 * Obtains the transient response of a first-order Ordinary Differential
 * Equation ODE subject to an impulse-input.
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
#define YINI 0.0

// prototypes
double fstep   (double t);				// step solution
double fimpulse(double t);				// impulse solution
double objf    (double yn, void *vprms);		// objective fun
double odefun  (double t, double y, double* prms);	// RHS ODE, f(t, y)
void   write   (char*, const int, double**);		// writes to file
void   display (const int, double**);			// writes to stdout

int main() {
	// Solves first-order ODEs using Euler's implicit Method

	const int N = 511 ;		// number of intervals
	const int numel = N + 1 ;	// number of elements in time array
	double ti = 0.0, tf = 3.0e1 ;	// initial and final times
	double yi = YINI ;		// initial value

	// allocates and packs the parameters for implicit solver
	double prms[] = {.0, .0, .0, RATE, FEXT};
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
	char filename[] = "output/impulse/iEuler.dat" ;
	write  (filename, numel, odesol[0]);
	strcpy (filename, "output/impulse/EulRK2.dat");
	write  (filename, numel, odesol[1]);
	display(numel, oderet);

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
        double dt = prms[0];
        double yi = prms[1];
        double tn = prms[2];
//      double k  = prms[3];	// rate constant
//      double b  = prms[4];	// external forcing constant
        return ( yn - yi - dt * fp(tn, yn, prms) );
}


double odefun (double t, double y, double* prms)
{

/*
 * Synopsis: 
 * Right-hand side of the ODE with step input.
 * 
 * Comments:
 * We take advantage of the linearity of the ODE to obtain the transient
 * response of the system to an impulse-input. Note that the 
 * first-derivative of the step-response y'(t) is equal to the response
 * of the (same) system when subjected to an impulse-input. So we solve
 * for the step-response y(t) and evaluate this function to obtain 
 * the impulse-response y'(t).
 *
 */

//      double dt = prms[0];
//      double yi = prms[1];
//      double tn = prms[2];
	double k  = prms[3];
	double b  = prms[4];
	return (b - k * y);
}


double fstep (double t) {
/* 
 * Synopsis:
 * Analytic solution of the first-order ODE subject to a step input:
 *  		y' + k * y = b * u(t), 		y(t = 0) = yi,
 * where k is the rate constant, and b is the external forcing constant,
 * u(t) is the unit-step, and yi is the initial-value.
 *
 */
	double k  = RATE;
	double b  = FEXT;
	double yi = YINI;
	return ( (yi - b / k) * exp(-k * t) + b / k );
}


double fimpulse (double t) {
/* 
 * Synopsis:
 * Analytic solution of the first-order ODE subject to an impulse input:
 *  		y' + k * y = b * u(t), 		y(t = 0) = yi,
 * where k is the rate constant, and b is the external forcing constant,
 * u(t) is the unit-impulse, and yi is the initial-value. 
 *
 * Comments:
 * Note that this function returns the first-derivative of the 
 * step-response y'(t); compare the mathematical expression with that
 * defined in the fstep() function above.
 *
 */
	double k  = RATE;
	double b  = FEXT;
	double yi = YINI;
	return ( -k * (yi - b / k) * exp(-k * t) );
}



void write (char filename[], const int numel, double **odesol) {
	// writes the step and impulse-response to a data file

	double Dy ;		// first-derivative y'(t), impulse response
	double *t = odesol[0];	// time, t
	double *y = odesol[1];	// step-response, y(t)
	double err_step ;	// error
	double err_impulse ;	// error

	FILE *pFile = NULL ;
	pFile = fopen(filename, "w");

	if (pFile == NULL) {
		fprintf(stderr, "I/O Error: %s\n", filename);
		fprintf(stderr, "aborting execution ... \n");
		exit(EXIT_FAILURE);
	}

	double prms[] = {.0, .0, .0, RATE, FEXT};
	char fmt[] = "%23.15e %23.15e %23.15e %23.15e\n" ;
	for (int i = 0 ; i != numel ; ++i) {
		Dy = odefun(t[i], y[i], prms);
		err_step    = absval( (fstep(t[i]) - y[i]) );
		err_impulse = absval( 
			( fimpulse(t[i]) - odefun(t[i], y[i], prms) )
		);
		// writes the time t, step y(t), impulse y'(t), and errors
		fprintf(pFile, fmt, t[i], y[i], Dy, err_step, err_impulse);
	}
	fclose(pFile);
}


void display (const int numel, double **odesol) {
	// displays the numerical solution and absolute error to stdout

	double err ;
	double *t = odesol[0];
	double *y = odesol[1];

	char fmt[] = "%23.15e %23.15e %23.15e\n" ;
	for (int i = 0 ; i != numel ; ++i) {
		err = absval( (fstep(t[i]) - y[i]) );
		fprintf(stdout, fmt, t[i], y[i], err);
	}
}
