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

// prototypes
double f (double, double) ;	/* ODE RHS Function f(t, y) */
void write (char filename[], const int numel, double **odesol) ;

int main() {
	// Solves an ODE using Euler's and Runge-Kutta Methods.

	const int N = 255 ;		// number of intervals
	const int numel = N + 1 ;	// number of elements in time array
	double ti = 0.0, tf = 5.0 ;	// initial and final times
	double yi = 1.0 ;		// initial value, y = y(t = ti)

	// creates placeholders for the numerical solutions odesol = [t, y]
	double **oderet = NULL ;	// as odesol in the Python impl*
	double *odesol[][2] = { {NULL, NULL}, {NULL, NULL} } ;


	/* solves the ode numerically via the specified methods */
	oderet = Euler    (odesol[0], ti, tf, yi, N, f) ;
	oderet = EulerRK2 (odesol[1], ti, tf, yi, N, f) ;


	// exports numerical solutions to data files
	char filename[] = "output/Euler.dat" ;
	strcpy (filename, "output/EuRK2.dat") ;
	write  (filename, numel, odesol[0]) ;
	write  (filename, numel, odesol[1]) ;
//	write  (filename, numel, oderet) ;	// alternate way of writing


	// frees memory buffers
	free (odesol[0][0]) ;
	free (odesol[0][1]) ;
	free (odesol[1][0]) ;
	free (odesol[1][1]) ;
	return 0 ;
}


// implementations
double f (double t, double y) {
	// ODE function is the expression at the right-hand side of the ODE
	double k = 1.0 ;
	return  (-k * y) ;
}


void write (char filename[], const int numel, double **odesol) {
	// writes the numerical solution to a data file
	double *t = odesol[0] ;
	double *y = odesol[1] ;

	FILE *pFile = NULL ;
	pFile = fopen(filename, "w") ;
	for (int i = 0 ; i != numel ; ++i)
		fprintf(pFile, "%23.15e \t %23.15e \n", t[i], y[i]) ;
	fclose(pFile) ;
}


/*
 * Comments:
 * My aim was to write a code that looks similar to the code written
 * in Python. That's why the numerical solvers return an array that
 * can be used for writing if the users wish to do so as shown in the
 * inline comment above.
 *
 */
