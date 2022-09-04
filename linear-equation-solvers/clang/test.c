/*
 * Computational Solutions                             September 04, 2022
 * IST 4360
 * Prof. M Diaz-Maldonado
 *
 * source: test.c
 *
 * Synopsis:
 * Tests the matrix class.
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

#include <stdio.h>
#include "Matrix.h"
#define ROWS 4
#define COLS 5

extern matrix_namespace const matrix;	// imports the matrix namespace

int main() {

	// creates a second-rank array
	double M[ROWS][COLS] = {
		{2, 1, 1, 0, 1},
		{4, 3, 3, 1, 1},
		{8, 7, 9, 5, 2},
		{6, 7, 9, 8, 4}
	};

	// creates the matrix from the second-rank array
	matrix_t *mat = matrix.create (ROWS, COLS, M);

	// prints the elements of the matrix on the console
	for (size_t i = 0; i != (mat -> rows); ++i)
	{
		for (size_t j = 0; j != (mat -> cols); ++j)
			printf("%8.2f", mat -> get(mat, i, j));
		printf("\n");
	}

	// destroys the matrix
	mat = matrix.destroy(mat);
	return 0;
}
