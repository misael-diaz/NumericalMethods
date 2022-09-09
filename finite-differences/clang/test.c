/*
 * Computational Methods                                 September 09, 2022
 * ICI 70320
 * Prof. M Diaz-Maldonado
 *
 * source: test.c
 *
 * Synopsis:
 * Tests the vector class.
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
 * [0] A Koenig and B Moo, Accelerated C++ Practical Programming by
 *     Example.
 * [1] JJ McConnell, Analysis of Algorithms, 2nd edition
 * [2] A Gilat and V Subramaniam, Numerical Methods for Engineers and
 *     Scientists: An Introduction with Applications using MATLAB
 * [3] H Wendland, Numerical Linear Algebra
 *
 */

#include <stdio.h>
#include "Vector.h"

extern vector_namespace const vector;	// imports vector namespace

// prototypes:
void test_pushback();

int main() {

	test_pushback();
	return 0;
}

// tests:
void test_pushback() {

	size_t numel = 16;
	// creates a vector with requested storage capacity
	vector_t *vec = vector.create (numel);

	// pushes data unto the back of the vector
	for (int i = 0; i != 32; ++i)
		vec -> push_back (vec, i);

	double *array = (vec -> array);
	// prints the data on the console
	for (size_t i = 0; i != 32; ++i)
		printf("%f\n", array[i]);

	// prints the vector size on the console
	printf("size: %lu \n", vec -> size(vec));

	// clears the data stored in vector
	vec -> clear(vec);
	printf("cleared vector\n");
	printf("size: %lu \n", vec -> size(vec));

	// pushes new data unto the back of the vector
	for (int i = 0; i != 32; ++i)
		vec -> push_back (vec, 32 - i);

	array = (vec -> array);
	// prints the new data on the console
	for (size_t i = 0; i != 32; ++i)
		printf("%f\n", array[i]);

	// prints the vector size on the console
	printf("size: %lu \n", vec -> size(vec));

	// frees the memory allocated for the vector
	vec = vector.destroy (vec);
}
