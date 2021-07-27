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


// constants
#define MAX_ITER 100
#define TOL 1.0e-8


// declarations (prototypes)

double bisect   ( double, double, double f(const double) ) ;
double regfal   ( double, double, double f(const double) ) ;
void   report   (const int n) ;
double bisector ( double*, double*, double*, double f(const double) ) ;
double interp   ( double*, double*, double*, double f(const double) ) ;

#endif
