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


static matrix_t* util_alloc_matrix_t ()
// allocates memory for a matrix object
{
	matrix_t *mat = (matrix_t*) malloc ( sizeof(matrix_t) );
	if (mat == NULL)
	{
		char errmsg [] = "memory allocation error: %s\n";
		fprintf (stderr, errmsg, strerror(errno) );
                exit(EXIT_FAILURE);
        }

	return mat;
}


static double* util_alloc_array_double_t (size_t size)
// allocates memory for an array of doubles of requested size
{
	double *array = (double*) malloc ( size * sizeof(double) );
	if (array == NULL)
	{
		char errmsg [] = "memory allocation error: %s\n";
		fprintf (stderr, errmsg, strerror(errno) );
                exit(EXIT_FAILURE);
        }

	return array;
}

static double* util_init_from_const_double_t (
	size_t size, double data[size], double c
)
// initializes the first-rank array `data' from the constant `c'
{
	for (size_t i = 0; i != size; ++i)
		data[i] = c;

	return data;
}


static double* util_init_from_array_double_t (
	size_t r, size_t c, size_t sz, double array[r][c], double data[sz]
)
// initializes the first-rank array `data' from the second-rank array
{
	size_t rows = r, cols = c;
	for (size_t i = 0; i != rows; ++i)
	{
		for (size_t j = 0; j != cols; ++j)
			data[j + i * cols] = array[i][j];
	}

	return data;
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
	matrix_t *mat = util_alloc_matrix_t ();

	size_t size = (rows * cols);
	// allocates memory for the data buffer of the matrix
	mat -> buffer = util_alloc_array_double_t (size);

	double *data = mat -> buffer;
	// copies the second-rank array into the data buffer of the matrix
	data = util_init_from_array_double_t (rows, cols, size, A, data);

	// defines the matrix shape and size
	mat -> rows = rows;
	mat -> cols = cols;
	mat -> size = size;
	// binds methods
	mat -> get  = get_method;
	return mat;
}


static matrix_t* zeros (size_t rows, size_t cols)
// constructs a matrix of zeros of requested shape
{
	// allocates memory for the matrix
	matrix_t *mat = util_alloc_matrix_t ();

	size_t size = (rows * cols);
	// allocates memory for the data buffer of the matrix
	mat -> buffer = util_alloc_array_double_t (size);

	double *data = mat -> buffer;
	// fills the data buffer of the matrix with zeros
	data = util_init_from_const_double_t (size, data, 0);

	// defines the matrix shape and size
	mat -> rows = rows;
	mat -> cols = cols;
	mat -> size = size;
	// binds methods
	mat -> get  = get_method;
	return mat;
}


static matrix_t* mult (matrix_t *matA, matrix_t *matB)
// implements the matrix multiplication C = A * B
{
	// references the data buffers of each respective matrix
	double *A = (matA -> buffer);
	double *B = (matB -> buffer);

	// gets the shapes of the matrices
	size_t rows = (matA -> rows);		// rows x knum matrix
	size_t cols = (matB -> cols);		// knum x cols matrix
	size_t knum = (matA -> cols);		// rows x cols matrix

	double prod = 0;
	// initializes the output matrix
	matrix_t *matC = zeros (rows, cols);
	double *C = (matC -> buffer);
	// performs the matrix multiplication C = A * B
	for (size_t i = 0; i != rows; ++i)
	{
		for (size_t j = 0; j != cols; ++j)
		{
			C[j + i * cols] = 0;
			for (size_t k = 0; k != knum; ++k)
			/* C[i][j] += A[i][k] * B[k][j] */
			{
				prod = A[k + i * knum] * B[j + k * cols];
				C[j + i * cols] += prod;
			}
		}
	}

	return matC;
}


// creates a namespace for the matrix class
matrix_namespace const matrix = {zeros, create, destroy, mult};
