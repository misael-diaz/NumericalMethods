#ifndef GUARD_MATRIX_H
#define GUARD_MATRIX_H
/*
 * Computational Solutions                             September 03, 2022
 * IST 4360
 * Prof. M Diaz-Maldonado
 *
 * source: matrix.h
 *
 * Synopsis:
 * Defines the matrix class.
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

typedef struct {
	// data
	size_t  rows;
	size_t  cols;
	size_t  size;
	double* buffer;
	// methods
	double (*get) (void*, size_t, size_t);
} matrix_t;

#endif
