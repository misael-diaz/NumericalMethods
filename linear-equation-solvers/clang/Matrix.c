/*
 * Computational Solutions                             September 04, 2022
 * IST 4360
 * Prof. M Diaz-Maldonado
 *
 * source: Matrix.c
 *
 * Synopsis:
 * Defines the methods of the matrix class.
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
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <errno.h>
#include "Matrix.h"

static double get_method (void *matrix, size_t i, size_t j)
// returns a copy of the element at Matrix[i][j]
{
	// asserts that the (universal) pointer is a matrix type pointer
	matrix_t *mat = matrix;
	size_t cols = (mat -> cols);
	double *data = (mat -> buffer);
	return (data[j + i * cols]);
}


static matrix_t* destroy (matrix_t* mat)
// destroys the matrix by freeing the allocated resources
{
	if (mat -> buffer)
	{
		free (mat -> buffer);
		mat -> buffer = NULL;
	}

	free (mat);
	mat = NULL;
	return mat;
}


static matrix_t* create (size_t rows, size_t cols, double A[rows][cols])
// constructs a matrix from a second-rank array
{
	// allocates memory for the matrix
	matrix_t *mat = (matrix_t*) malloc ( sizeof(matrix_t) );
	if (mat == NULL)
	{
		char errmsg [] = "memory allocation error: %s\n";
		fprintf (stderr, errmsg, strerror(errno) );
                exit(EXIT_FAILURE);
        }

	// allocates memory for the data buffer of the matrix
	mat -> buffer = (double*) malloc ( rows * cols * sizeof(double) );
	if (mat -> buffer == NULL)
	{
		char errmsg [] = "memory allocation error: %s\n";
		fprintf (stderr, errmsg, strerror(errno) );
                exit(EXIT_FAILURE);
        }

	double *data = mat -> buffer;
	// copies the second-rank array into the data buffer of the matrix
	for (size_t i = 0; i != rows; ++i)
	{
		for (size_t j = 0; j != cols; ++j)
			data[j + i * cols] = A[i][j];
	}

	// defines the matrix shape and size
	mat -> rows = rows;
	mat -> cols = cols;
	mat -> size = (rows * cols);
	// binds methods
	mat -> get  = get_method;
	return mat;
}


// creates a namespace for the matrix constructor and destructor
matrix_namespace const matrix = {create, destroy};
