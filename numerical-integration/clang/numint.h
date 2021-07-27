#ifndef GUARD_NUMINT_H
#define GUARD_NUMINT_H
/*
 * Applied Numerical Analysis                                 July 27, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: numint.h
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
double lsum ( double, double, int, double f(const double) ) ;
double rsum ( double, double, int, double f(const double) ) ;
double trap ( double, double, int, double f(const double) ) ;
#endif
