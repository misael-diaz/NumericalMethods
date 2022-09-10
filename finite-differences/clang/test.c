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
void test_zeros();
void test_ones();
void test_linspace();
void test_copy();

int main() {

	test_pushback();
	test_linspace();
	test_zeros();
	test_ones();
	test_copy();
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

	// stores same data in array temporary
	double data [32];
	for (int i = 0; i != 32; ++i)
		data[i] = i;

	double diff = 0;
	double *array = (vec -> array);
	// computes differences
	for (size_t i = 0; i != 32; ++i)
		diff += (data[i] - array[i]);

	printf("push-back-method-test[0]: ");
	if (diff != 0.0)
		printf("FAIL\n");	// fails if there are differences
	else
		printf("pass\n");	// passes if the are none


	printf("push-back-method-test[1]: ");
	if (vec -> size(vec) != 32)
		printf("FAIL\n");	// fails if size does not match
	else
		printf("pass\n");	// passes otherwise


	// clears the data stored in vector
	vec -> clear(vec);
	printf("push-back-method-test[2]: ");
	if (vec -> size(vec) != 0)
		printf("FAIL\n");	// fails if size does not match
	else
		printf("pass\n");	// passes otherwise


	// pushes new data unto the back of the vector
	for (int i = 0; i != 32; ++i)
		vec -> push_back (vec, 32 - i);

	// stores same data in array temporary
	for (int i = 0; i != 32; ++i)
		data[i] = (32 - i);


	diff = 0.0;
	array = (vec -> array);
	// computes differences
	for (size_t i = 0; i != 32; ++i)
		diff += (data[i] - array[i]);

	printf("push-back-method-test[3]: ");
	if (diff != 0.0)
		printf("FAIL\n");	// fails if there are differences
	else
		printf("pass\n");	// passes if the are none


	// frees the memory allocated for the vector
	vec = vector.destroy (vec);
}


void test_zeros ()
// tests the zeros method
{
	size_t size = 16;
	// creates vector of zeros
	vector_t *vec = vector.zeros(size);

	double sum = 0;
	double *array = (vec -> array);
	// accumulates the data stored in the vector
	for (size_t i = 0; i != size; ++i)
		sum += array[i];

	printf("zeros-test[0]: ");	// checks stored data
	if (sum != 0.0)
		printf("FAIL\n");
	else
		printf("pass\n");

	printf("zeros-test[1]: ");
	if (vec -> size(vec) != size)	// checks vector size
		printf("FAIL\n");
	else
		printf("pass\n");

	// frees the memory allocated for the vector
	vec = vector.destroy (vec);
}


void test_ones ()
// tests the ones method
{
	size_t size = 16;
	// creates vector of zeros
	vector_t *vec = vector.ones(size);

	// prints the vector size on the console
	printf("size: %lu \n", vec -> size(vec));

	double *array = (vec -> array);
	// prints the vector data on the console
	for (size_t i = 0; i != size; ++i)
		printf("%7.4f\n", array[i]);

	// frees the memory allocated for the vector
	vec = vector.destroy (vec);
}


void test_linspace ()
// tests the linspace method
{
	// creates vector of 33 equally spaced elements in [-1, 1]
	vector_t *vec = vector.linspace(-1, 1, 33);

	// prints the vector size on the console
	printf("size: %lu \n", vec -> size(vec));

	double *array = (vec -> array);
	// prints the vector data on the console
	for (size_t i = 0; i != 33; ++i)
		printf("%+7.4f\n", array[i]);

	// frees the memory allocated for the vector
	vec = vector.destroy (vec);
}


void test_copy ()
// tests the (shallow) copy method of the vector class
{
	// creates a vector
	size_t size = 17;
	double x_l = -1, x_u = 1;
	vector_t *vec_x = vector.linspace(x_l, x_u, size);
	// creates another vector of the same size
	vector_t *vec_y = vector.zeros(size);
	// copies the contents of vector `x' into vector `y'
	vec_y -> copy (vec_y, vec_x);

	/* tests the copy method */

	double diff = 0;
	double *x = (vec_x -> array);
	double *y = (vec_y -> array);
	// computes the differences
	for (size_t i = 0; i != size; ++i)
		diff += (x[i] - y[i]);

	printf("copy-method-test[0](): ");
	if (diff != 0.0)
		printf("FAIL\n");	// fails if there are differences
	else
		printf("pass\n");	// passes if the are none

	// checks if the data buffers of the vectors are different objects
	printf("copy-method-test[1](): ");
	if (x == y)
		printf("FAIL\n");	// fails if objects are equal
	else
		printf("pass\n");	// passes if different

	// frees the memory allocated for the vectors
	vec_x = vector.destroy (vec_x);
	vec_y = vector.destroy (vec_y);
}
