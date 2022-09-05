#ifndef GUARD_MATRIX_CLASS_H
#define GUARD_MATRIX_CLASS_H
/*
 * Computational Solutions                             September 04, 2022
 * IST 4360
 * Prof. M Diaz-Maldonado
 *
 * source: Matrix.h
 *
 * Synopsis:
 * Defines the matrix namespace.
 *
 *
 * Copyright (c) 2022 Misael Diaz-Maldonado
 *
 * This file is released under the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 *
 * References:
 * [0] A Gilat and V Subramanian, Numerical Methods for Engineers and
 *     Scientists: An Introduction with Applications using MATLAB
 * [1] A Koenig and B Moo, Accelerated C++ Practical Programming by
 *     Example.
 * [2] JJ McConnell, Analysis of Algorithms, 2nd edition
 *
 */

#include "matrix.h"

typedef struct {
	// constructor
	matrix_t* (*const zeros) (size_t, size_t);
	matrix_t* (*const create) (size_t n, size_t m, double x[n][m]);
	// destructor
	matrix_t* (*const destroy) (matrix_t*);
} matrix_namespace;

#endif
