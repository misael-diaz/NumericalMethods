#ifndef GUARD_ODES_H
#define GUARD_ODES_H
/*
 * Applied Numerical Analysis                                 July 27, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: odes.h
 *
 * Synopsis:
 * Header file for (some) numerical integration methods.
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

#include "fzero.h"


// Type declarations
typedef struct {
	double (*objf)   (double, void*) ;
	double (*odefun) (double, double, double*) ;
	double *prms ;	// parameter array
} iODE_solverParams ;	// implicit, ODE Solver Parameters


// Methods
double**  Euler   (double**, double, double, double, const int,
                   double f(double, double) ) ;
double** iEuler   (double**, double, double, double, const int,
                   double odefun(double, double, double*), void*) ;
double** EulerRK2 (double**, double, double, double, const int,
                   double f(double, double) ) ;
#endif
