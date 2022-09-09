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

// prototypes:
void print ();
void zeros ();
void mmult ();

int main() {

	print ();
	zeros ();
	mmult ();
	return 0;
}


// tests:
void print ()
// creates a matrix from a second-rank array, prints it, and destroys it
{
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
}


void zeros ()
// creates a matrix of zeros, prints it, and destroys it
{
	// creates a rows x cols matrix of zeros
	matrix_t *mat = matrix.zeros (ROWS, COLS);

	printf("\nzeros:\n");
	// prints the elements of the matrix on the console
	for (size_t i = 0; i != (mat -> rows); ++i)
	{
		for (size_t j = 0; j != (mat -> cols); ++j)
			printf("%8.2f", mat -> get(mat, i, j));
		printf("\n");
	}

	// destroys the matrix
	mat = matrix.destroy(mat);
}


void mmult ()
// tests matrix multiplication
{
	// creates the second-rank arrays for the test
	double A[ROWS][COLS] = {
		{2, 1, 1, 0, 1},
		{4, 3, 3, 1, 1},
		{8, 7, 9, 5, 2},
		{6, 7, 9, 8, 4}
	};

	double B[COLS][ROWS] = {
		{2, 4, 8, 6},
		{1, 3, 7, 7},
		{1, 3, 9, 9},
		{0, 1, 5, 8},
		{1, 1, 2, 4}
	};

	// creates the matrices A and B from the second-rank arrays
	matrix_t *matA = matrix.create (ROWS, COLS, A);
	matrix_t *matB = matrix.create (COLS, ROWS, B);

	// computes the matrix product: C = A * B
	matrix_t *matC = matrix.mult (matA, matB);

	printf("\nC = mult(A, B):\n");
	// prints the elements of the resulting matrix on the console
	for (size_t i = 0; i != (matC -> rows); ++i)
	{
		for (size_t j = 0; j != (matC -> cols); ++j)
			printf("%8.2f", matC -> get(matC, i, j));
		printf("\n");
	}

	// destroys the matrices
	matA = matrix.destroy (matA);
	matB = matrix.destroy (matB);
	matC = matrix.destroy (matC);
}
