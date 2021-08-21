#ifndef GUARD_FZERO_H
#define GUARD_FZERO_H

/*
 * source: fzero.h
 * author: misael-diaz
 * date:   2021/08/13
 *
 * Synopsis:
 * Header file for the nonlinear equation solver fzero.
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
#include <stdbool.h>
// Math MACROS
#define absval(x) (signbit(x)? -x: x)

// constants (defaults)
#define MAX_ITER 256
#define TOL      1.0e-12
#define VERBOSE  false

// typedefs
typedef struct {	/* (yet) unused in this minimalistic version */
	double tol ;	// tolerance
	int max_iter ;	// maximum number of iterations
	bool verbose ;
} nls_conf ;		// nonlinear-solver configuration struct

// method prototypes
double fzero (double, double, double f(double, void*), void*);
#endif
