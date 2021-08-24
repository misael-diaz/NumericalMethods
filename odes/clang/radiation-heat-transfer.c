/*
 * Applied Numerical Analysis				    August 24, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: radiation-heat-transfer.c
 *
 * Synopsis:
 * Solves a nonlinear Ordinary Differential Equation ODE, which describes
 * the radiation heat transfer from a slab having a uniform temperature
 * distribution.
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
 * [2] TL Bergman, AS Lavine, FP Incropera, DP DeWitt, Fundamentals
 *     of Heat and Mass Transfer, 8th edition.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "odes.h"

/* MACROS
 * Stefan-Boltzmann constant, sigma, 5.67e-8 W / (m^2 K^4)
 * Emmisivity,                epsil, 0.95
 * surroundings temperature,  T_sur, 0.29315 (10^3) K
 * density,                     rho, 2700.0  kg / m^3
 * heat capacity,                 c, 900.0   J / (kg K)
 * thickness,                     L, 4.0e-3  m
 * initial temperature,          T0, 1.2000  (10^3) K
 */

#define SIGMA 5.67e-8
#define EPS   0.95
#define TSUR  0.29315
#define RHO   2700.0
#define CP    900.0
#define TAU   4.0e-3
#define TEMP0 1.2

// prototypes
double ferr    (double t, double T);			// error function
double objf    (double yn, void *vprms);		// objective fun
double odefun  (double t, double y, double* prms);	// RHS ODE, f(t, y)
void   write   (char*, const int, double**);		// writes to file
void   display (const int, double**);			// writes to stdout

int main() {
	// Solves first-order ODEs using Euler's implicit Method

	const int N = 65535 ;		// number of intervals
	const int numel = N + 1 ;	// number of elements in time array
	double ti = 0.0, tf = 5.0e1 ;	// initial and final times
	double yi = TEMP0 ;		// initial value

	// allocates and packs the parameters for implicit solver
	double prms[] = {0., 0., 0., SIGMA, EPS, TSUR, RHO, CP, TAU};
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
	char filename[] = "output/heat-transfer/radiation/iEuler.dat" ;
	write  (filename, numel, odesol[0]);
	strcpy (filename, "output/heat-transfer/radiation/EulRK2.dat");
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
        return ( yn - yi - dt * fp(tn, yn, prms) );
}


double odefun (double t, double y, double* prms)
{	// Synopsis: 
	// Right-hand side of the nonlinear ODE (scaled energy equation):
	// 		T' = T^4 - T_sur^4,
	// where T is the slab temperature, and T_sur is the surroundings
	// temperature. 
	// Note: The system time has been non-dimensionalized using 1/K
	// (see below expression for K and the final remarks at the end).

//      double dt    = prms[0];
//      double yi    = prms[1];
//      double tn    = prms[2];
//	double sigma = prms[3];	// Stefan-Boltzmann constant
//	double epsil = prms[4];	// emissivity
	double T_sur = prms[5];	// surroundings temperature
//	double rho   = prms[6];	// slab density
//	double cp    = prms[7];	// slab heat capacity
//	double L     = prms[8];	// slab thickness

	double T = y ;
//	double K = 1.0e9 * epsil * sigma / (rho * cp * L) ;
	return (  -( pow(T, 4.0) - pow(T_sur, 4.0) )  ) ;
}


double ferr (double t, double T)
{
/*
 * Synopsis:
 * Returns the error given the time t and slab temperature T.
 *
 * Comments:
 * To corroborate the correctness of the numerical solution it is
 * convenient to define the error function F(t, T), which for any
 * valid t, T pair it should evaluate to (nearly) zero.
 *
 * The mathematical expression results from separation of variables
 * and subsequent integration of the fourth-degree polynomial via
 * partial fractions. The integration constant is determined by
 * satisfying the initial condition, which we do indirectly so by
 * evaluating the integral from the initial to the present time.
 *
 */
	const double T0    = TEMP0;
	const double T_sur = TSUR;
	const double C     = 1.0 / ( 4.0 * pow(TSUR, 3.0) );
//	const double K     = 1.0e9 * EPS * SIGMA / (RHO * CP * TAU);
	return ( C * (2.0 * ( atan(T / T_sur) - atan(T0 / T_sur) ) +
	         log( (T_sur + T) / (T_sur + T0) ) -
	         log( (T_sur - T) / (T_sur - T0) ) ) - t );
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
		fprintf(pFile, fmt, t[i], y[i], absval(ferr(t[i], y[i])));
	fclose(pFile);
}


void display (const int numel, double **odesol) {
	// Displays the numerical and exact solutions, and the absolute 
	// error to the standard output.

	double err ;
	double *t = odesol[0];
	double *y = odesol[1];

	char fmt[] = "%23.15e %23.15e %23.15e\n" ;
	for (int i = 0 ; i != numel ; ++i) {
		err = absval( ferr(t[i], y[i]) );
		fprintf(stdout, fmt, t[i], y[i], err);
	}
}


/*
 * Final Remarks:
 * The slab temperature was scaled by a factor of 1000 to keep
 * propagation errors from spreading too rapidly. To account for this the
 * energy equation has been rescaled. In the end all the constants end up
 * lumped into a single constant used to non-dimensionalize time.
 * The latter is desirable since the time-step can be made as small needed
 * while keeping the final integration time to a reasonable value (as so
 * the computational cost).
 *
 * For post-processing one can convert back the time to seconds and the
 * slab temperature to Kelvins or to any other suitable units.
 *
 * I did not have to use so many time-steps (65535) to solve the problem
 * numerically. (It's a small problem so I can get away with that though.)
 * It was done to check that as the number of steps grow the error
 * function diminishes.
 *
 * At first, I did not intend to rescale the energy equation so more
 * parameters than needed are supplied to odefun() function. I am
 * leaving it like that since it shows the scaling factor 1.0e9
 * that results from scaling the temperature on both sides of the
 * energy equation, and because it also shows the resulting
 * time scale (1/K). Both might be of some use or interest eventually.
 *
 */
