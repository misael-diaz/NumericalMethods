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

#include <math.h>
#include <stdio.h>
#include "Vector.h"
#include "isolver.h"

extern vector_namespace const vector;	// imports vector namespace

// prototypes:
void test_pushback();
void test_zeros();
void test_ones();
void test_linspace();
void test_copy();
void test_qnorm();
vector_t** Jacobi (vector_t *pdesol[6], vector_t*, const isolver_prms_t*);
vector_t** GaussSeidel (
	vector_t *pdesol[6], vector_t*, const isolver_prms_t*
);
void test_steady_1d_transport_Jacobi();
void test_steady_1d_transport_GaussSeidel();
void test_transient_1d_transport_Jacobi ();
void test_transient_1d_transport_GaussSeidel ();

int main() {

	/*
	test_pushback();
	test_linspace();
	test_zeros();
	test_ones();
	test_copy();
	test_qnorm();
	*/

	// solves the steady 1d transport problem with the Jacobi method
	test_steady_1d_transport_Jacobi();
	// applies the Gauss-Seidel method to solve the same problem
	test_steady_1d_transport_GaussSeidel();

	// solves the transient 1d transport problem with the Jacobi method
	test_transient_1d_transport_Jacobi ();
	// solves the transient problem with the Gauss-Seidel method
	test_transient_1d_transport_GaussSeidel ();
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


	double diff = 0;
	double *array = (vec -> array);
	// computes differences
	for (size_t i = 0; i != size; ++i)
		diff += (1.0 - array[i]);

	printf("ones-test[0]: ");
	if (diff != 0.0)
		printf("FAIL\n");	// fails if there are differences
	else
		printf("pass\n");	// passes if the are none


	double sum = 0;
	// checks the vector size
	for (size_t i = 0; i != size; ++i)
		sum += array[i];

	printf("ones-test[1]: ");
	if (sum != size)
		printf("FAIL\n");	// fails if size differs
	else
		printf("pass\n");	// passes otherwise


	// frees the memory allocated for the vector
	vec = vector.destroy (vec);
}


void test_linspace ()
// tests the linspace method
{
	// defines the vector size
	size_t size = 33;
	// defines the vector limits
	double x_l = -1, x_u = 1;
	// creates vector of `size' equally spaced elements in [-1, 1]
	vector_t *vec = vector.linspace(x_l, x_u, size);


	double data [size];
	double dx = (x_u - x_l) / ( (double) (size - 1) );
	// stores same data in array temporary
	for (size_t i = 0; i != size; ++i)
		data[i] = x_l + ( (double) i ) * dx;


	double diff = 0;
	double *array = (vec -> array);
	// computes differences between the stored and expected data
	for (size_t i = 0; i != 32; ++i)
		diff += (data[i] - array[i]);

	printf("linspace-test[0]: ");
	if (diff != 0.0)
		printf("FAIL\n");	// fails if there are differences
	else
		printf("pass\n");	// passes if the are none


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


void test_qnorm()
// tests the quick-norm method
{
	// creates a vector that stores the integers in the range [0, 256)
	size_t size = 256;
	double x_l = 0, x_u = 255;
	vector_t *vec = vector.linspace(x_l, x_u, size);


	// computes the norm
	double norm = vec -> qnorm (vec);
	// computes the expected result -- the sum of squares in [0, 256)
	double sum_squares = (size * (size - 1) * (2 * size - 1) / 6);

	printf("quick-norm-method-test[0]: ");
	if (norm != sum_squares)
		printf("FAIL\n");
	else
		printf("pass\n");


	// frees the memory allocated for the vector
	vec = vector.destroy (vec);
}


vector_t** Jacobi (
	vector_t *pdesol[6], vector_t *pdevec, const isolver_prms_t* prms
)
// possible tailored implementation of the Jacobi method
{
	// gets the transient parameter
	double alpha = prms -> alpha;
	// gets the tolerance and the maximum number of iterations
	double tol = prms -> tol;
	size_t iters = prms -> iters;
	// gets the verbose parameter
	bool verbose = prms -> verbose;


	// references the position and (temperature) field variables
	vector_t *vec_x = pdesol[0];
	vector_t *vec_g = pdesol[1];
	// references the (heat) source term vector
	vector_t *vec_src = pdesol[2];
	// references the guess vector
	vector_t *vec_g0  = pdesol[3];
	// references the error vector
	vector_t *vec_err = pdesol[4];
	// references the state vector
	vector_t *vec_state = pdesol[5];
	// references the right-hand side, finite-difference, vector
	vector_t *vec_b = pdevec;


	// gets the size of the vectors
	size_t size = vec_x -> size (vec_x);
	// initializes the finite-difference vector with the source vector
	vec_b -> copy (vec_b, vec_src);


	// defines iterators
	double *b   = (vec_b   -> array);
	double *g   = (vec_g   -> array);
	double *g0  = (vec_g0  -> array);
	double *err = (vec_err -> array);
	double *state = (vec_state -> array);


	// updates the finite-difference vector with the transient vector
	for (size_t i = 0; i != size; ++i)
		b[i] += (-alpha * g[i]);


	// initializes the state
	state[0] = 1.0;
	// initializes the norm of the error vector
	double norm = 0;
	// gets the number of discretization intervals
	size_t N = (size - 1);
	// updates the diagonal coefficient with the transient contribution
	double c = 1.0 / (alpha + 2.0);
	/* applies the Jacobi method to solve for the field variable */
	for (size_t i = 0; i != iters; ++i)
	{

		// updates the node next to the lower boundary, x_l
		g[1] = -c * (b[1] - g0[2]);

		// updates the intermediate nodes
		for (size_t j = 2; j != (N - 1); ++j)
			g[j] = -c * (b[j] - g0[j - 1] - g0[j + 1]);

		// updates the node next to the upper boundary, x_u
		g[N - 1] = -c * (b[N - 1] - g0[N - 2]);


		// computes the error vector
		for (size_t n = 0; n != size; ++n)
			err[n] = (g[n] - g0[n]);


		// computes the norm of the error vector
		norm = vec_err -> qnorm (vec_err);
		/* checks for convergence */
		if (norm < tol)
		{
			state[0] = 0.0;
			char msg [] = "Jacobi(): solution found after "
				"%lu iters\n";
			if (verbose) printf(msg, i + 1);
			break;
		}


		// copies the field (g0 <- g) for the next iteration
		vec_g0 -> copy (vec_g0, vec_g);
	}

	return pdesol;
}


vector_t** GaussSeidel (
	vector_t *pdesol[6], vector_t *pdevec, const isolver_prms_t* prms
)
// possible tailored implementation of the Gauss-Seidel method
{
	// gets the transient parameter
	double alpha = prms -> alpha;
	// gets the tolerance and the maximum number of iterations
	double tol = prms -> tol;
	size_t iters = prms -> iters;
	// gets the verbose parameter
	bool verbose = prms -> verbose;


	// references the position and (temperature) field variables
	vector_t *vec_x = pdesol[0];
	vector_t *vec_g = pdesol[1];
	// references the (heat) source term vector
	vector_t *vec_src = pdesol[2];
	// references the guess vector
	vector_t *vec_g0  = pdesol[3];
	// references the error vector
	vector_t *vec_err = pdesol[4];
	// references the state vector
	vector_t *vec_state = pdesol[5];
	// references the right-hand side, finite-difference, vector
	vector_t *vec_b = pdevec;


	// gets the size of the vectors
	size_t size = vec_x -> size (vec_x);
	// initializes the finite-difference vector with the source vector
	vec_b -> copy (vec_b, vec_src);


	// defines iterators
	double *b   = (vec_b   -> array);
	double *g   = (vec_g   -> array);
	double *g0  = (vec_g0  -> array);
	double *err = (vec_err -> array);
	double *state = (vec_state -> array);


	// updates the finite-difference vector with the transient vector
	for (size_t i = 0; i != size; ++i)
		b[i] += (-alpha * g[i]);


	// initializes the state
	state[0] = 1.0;
	// initializes the norm of the error vector
	double norm = 0;
	// gets the number of discretization intervals
	size_t N = (size - 1);
	// updates the diagonal coefficient with the transient contribution
	double c = 1.0 / (alpha + 2.0);
	/* uses the Gauss-Seidel method to solve for the field variable */
	for (size_t i = 0; i != iters; ++i)
	{

		// updates the node next to the lower boundary, x_l
		g[1] = -c * (b[1] - g0[2]);

		// updates the intermediate nodes
		for (size_t j = 2; j != (N - 1); ++j)
			g[j] = -c * (b[j] - g[j - 1] - g0[j + 1]);

		// updates the node next to the upper boundary, x_u
		g[N - 1] = -c * (b[N - 1] - g[N - 2]);


		// computes the error vector
		for (size_t n = 0; n != size; ++n)
			err[n] = (g[n] - g0[n]);


		// computes the norm of the error vector
		norm = vec_err -> qnorm (vec_err);
		/* checks for convergence */
		if (norm < tol)
		{
			state[0] = 0.0;
			char msg [] = "Gauss-Seidel(): solution found "
				"after %lu iters\n";
			if (verbose) printf(msg, i + 1);
			break;
		}


		// copies the field (g0 <- g) for the next iteration
		vec_g0 -> copy (vec_g0, vec_g);
	}

	return pdesol;
}


void test_steady_1d_transport_Jacobi ()
// solves a steady one-dimensional transport problem via finite-differences
{

	/* defines the solver parameters */

	// defines the steady parameter
	double alpha = 0;
	// defines the tolerance of the linear solver
	double tol = 1.0 / ( (double) (0x400000000000) );	// ~1.4e-14
	// defines the maximum number of iterations of the linear solver
	size_t iters = (0x0008FFFF);	// about 500K iterations
	// sets the solver to be verbose
	bool verbose = true;

	// initializes the iterative solver parameters
	isolver_prms_t const prms = {
		.alpha = alpha, .tol = tol,
		.iters = iters, .verbose = verbose
	};


	/* defines the finite-differences problem */


	// sets the number of discretization intervals to 256
	size_t N = (0x00000100);
	// defines the vector size
	size_t size = (N + 1);


	// defines the system domain limits along the x-axis
	double x_l = -1.0, x_u = 1.0;
	// computes the step-size
	double dx = (x_u - x_l) / ( (double) N );


	// initializes the placeholder for the solution of the PDE
	vector_t *pdesol[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
	// initializes the position vector
	pdesol[0] = vector.linspace(x_l, x_u, size);
	// initializes the (temperature) field variable
	pdesol[1] = vector.zeros(size);
	// initializes the (heat) source vector
	pdesol[2] = vector.zeros(size);
	// initializes the guess vector
	pdesol[3] = vector.zeros(size);
	// initializes the error vector
	pdesol[4] = vector.zeros(size);
	// initializes the state vector
	pdesol[5] = vector.zeros(1);


	// creates the finite-difference vector (right-hand side of SLE)
	vector_t *pdevec = vector.zeros(size);
	// creates alias for the source vector
	vector_t *vec_b  = pdevec;


	// references the (heat) source vector
	vector_t *vec_source = pdesol[2];
	// defines an iterator for the (heat) source vector
	double *source = (vec_source -> array);

	// defines the non-dimensional volumetric (heat) source term
	double H = 1;
	// defines the constant (heat) source vector
	for (size_t i = 0; i != size; ++i)
		source[i] = -(dx) * (dx) * H;


	/* numeric solution */


	// solves for the (temperature) field variable iteratively
	vector_t **ret = Jacobi (pdesol, pdevec, &prms);


	/* post-processing */


	// references the returned position vector and the field variable
	vector_t *vec_x = ret[0];
	vector_t *vec_g = ret[1];
	vector_t *vec_src = ret[2];
	vector_t *vec_g0  = ret[3];
	vector_t *vec_err = ret[4];
	vector_t *vec_state = ret[5];


	// initializes the analytic solution vector
	vector_t *vec_analytic = vector.zeros(size);


	// defines iterators
	double *x = (vec_x -> array);
	double *g = (vec_g -> array);
	double *err = (vec_err -> array);
	double *analytic = (vec_analytic -> array);


	// computes the analytic solution
	for (size_t i = 0; i != size; ++i)
		analytic[i] = 0.5 * H * (1.0 - x[i]) * (1.0 + x[i]);

	// computes the differences between the exact and numeric solutions
	for (size_t i = 0; i != size; ++i)
		err[i] = (g[i] - analytic[i]);


	// reports the average error
	double norm = vec_err -> qnorm (vec_err);
	printf("error: %.4e\n", sqrt(norm) / size );


	// frees the vectors from memory
	vec_x = vector.destroy (vec_x);
	vec_g = vector.destroy (vec_g);
	vec_b = vector.destroy (vec_b);
	vec_g0 = vector.destroy (vec_g0);
	vec_src = vector.destroy (vec_src);
	vec_err = vector.destroy (vec_err);
	vec_state = vector.destroy (vec_state);
	vec_analytic = vector.destroy (vec_analytic);
}


void test_steady_1d_transport_GaussSeidel ()
// solves a steady one-dimensional transport problem via finite-differences
{

	/* defines the solver parameters */


	// defines the steady parameter
	double alpha = 0;
	// defines the tolerance of the linear solver
	double tol = 1.0 / ( (double) (0x400000000000) );	// ~1.4e-14
	// defines the maximum number of iterations of the linear solver
	size_t iters = (0x0008FFFF);	// about 500K iterations
	// sets the solver to be verbose
	bool verbose = true;

	// initializes the iterative solver parameters
	isolver_prms_t const prms = {
		.alpha = alpha, .tol = tol,
		.iters = iters, .verbose = verbose
	};


	/* defines the finite-differences problem */


	// sets the number of discretization intervals to 256
	size_t N = (0x00000100);
	// defines the vector size
	size_t size = (N + 1);


	// defines the system domain limits along the x-axis
	double x_l = -1.0, x_u = 1.0;
	// computes the step-size
	double dx = (x_u - x_l) / ( (double) N );


	// initializes the placeholder for the solution of the PDE
	vector_t *pdesol[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
	// initializes the position vector
	pdesol[0] = vector.linspace(x_l, x_u, size);
	// initializes the (temperature) field variable
	pdesol[1] = vector.zeros(size);
	// initializes the (heat) source vector
	pdesol[2] = vector.zeros(size);
	// initializes the guess vector
	pdesol[3] = vector.zeros(size);
	// initializes the error vector
	pdesol[4] = vector.zeros(size);
	// initializes the state vector
	pdesol[5] = vector.zeros(1);


	// creates the finite-difference vector (right-hand side of SLE)
	vector_t *pdevec = vector.zeros(size);
	// creates alias for the source vector
	vector_t *vec_b  = pdevec;


	// references the (heat) source vector
	vector_t *vec_source = pdesol[2];
	// defines an iterator for the (heat) source vector
	double *source = (vec_source -> array);

	// defines the non-dimensional volumetric (heat) source term
	double H = 1;
	// defines the constant (heat) source vector
	for (size_t i = 0; i != size; ++i)
		source[i] = -(dx) * (dx) * H;


	/* numeric solution */


	// solves for the (temperature) field variable iteratively
	vector_t **ret = GaussSeidel (pdesol, pdevec, &prms);


	/* post-processing */


	// references the returned position vector and the field variable
	vector_t *vec_x = ret[0];
	vector_t *vec_g = ret[1];
	vector_t *vec_src = ret[2];
	vector_t *vec_g0  = ret[3];
	vector_t *vec_err = ret[4];
	vector_t *vec_state = ret[5];


	// initializes the analytic solution vector
	vector_t *vec_analytic = vector.zeros(size);


	// defines iterators
	double *x = (vec_x -> array);
	double *g = (vec_g -> array);
	double *err = (vec_err -> array);
	double *analytic = (vec_analytic -> array);


	// computes the analytic solution
	for (size_t i = 0; i != size; ++i)
		analytic[i] = 0.5 * H * (1.0 - x[i]) * (1.0 + x[i]);

	// computes the differences between the exact and numeric solutions
	for (size_t i = 0; i != size; ++i)
		err[i] = (g[i] - analytic[i]);


	// reports the average error
	double norm = vec_err -> qnorm (vec_err);
	printf("error: %.4e\n", sqrt(norm) / size );


	// frees the vectors from memory
	vec_x = vector.destroy (vec_x);
	vec_g = vector.destroy (vec_g);
	vec_b = vector.destroy (vec_b);
	vec_g0 = vector.destroy (vec_g0);
	vec_src = vector.destroy (vec_src);
	vec_err = vector.destroy (vec_err);
	vec_state = vector.destroy (vec_state);
	vec_analytic = vector.destroy (vec_analytic);
}


void test_transient_1d_transport_Jacobi ()
// solves a transient 1d transport problem via finite-differences
{

	/* defines the solver parameters */

	// defines the steady parameter
	double alpha = 2.0;
	// defines the tolerance of the linear solver
	double tol = 1.0 / ( (double) (0x400000000000) );	// ~1.4e-14
	// defines the maximum number of iterations of the linear solver
	size_t iters = (0x0008FFFF);	// about 500K iterations
	// sets the solver to be verbose
	bool verbose = false;

	// initializes the iterative solver parameters
	isolver_prms_t const prms = {
		.alpha = alpha, .tol = tol,
		.iters = iters, .verbose = verbose
	};


	/* defines the finite-differences problem */


	// sets the number of discretization intervals to 256
	size_t N = (0x00000100);
	// defines the vector size
	size_t size = (N + 1);


	// defines the system domain limits along the x-axis
	double x_l = -1.0, x_u = 1.0;
	// computes the step-size
	double dx = (x_u - x_l) / ( (double) N );
	// computes the time step
	double dt = (dx * dx) / alpha;


	// initializes the placeholder for the solution of the PDE
	vector_t *pdesol[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
	// initializes the position vector
	pdesol[0] = vector.linspace(x_l, x_u, size);
	// initializes the (temperature) field variable
	pdesol[1] = vector.zeros(size);
	// initializes the (heat) source vector
	pdesol[2] = vector.zeros(size);
	// initializes the guess vector
	pdesol[3] = vector.zeros(size);
	// initializes the error vector
	pdesol[4] = vector.zeros(size);
	// initializes the state vector
	pdesol[5] = vector.zeros(1);


	// creates the finite-difference vector (right-hand side of SLE)
	vector_t *pdevec = vector.zeros(size);
	// creates alias for the source vector
	vector_t *vec_b  = pdevec;


	// references the (heat) source vector
	vector_t *vec_source = pdesol[2];
	// defines an iterator for the (heat) source vector
	double *source = (vec_source -> array);

	// defines the non-dimensional volumetric (heat) source term
	double H = 1;
	// defines the constant (heat) source vector
	for (size_t i = 0; i != size; ++i)
		source[i] = -(dx) * (dx) * H;


	/* numeric solution */


	size_t steps = (0x00100000);	// sets to ~1 million time steps
	// solves for the (temperature) field variable iteratively
	for (size_t i = 0; i != steps; ++i)
	{
		vector_t **ret = Jacobi (pdesol, pdevec, &prms);
		vector_t *vec_state = ret[5];
		double *state = (vec_state -> array);
		if (state[0] != 0.0)
		{
			printf("Jacobi solver failed\n");
			printf("try again with different solver params\n");
			break;
		}
	}


	/* post-processing */


	// references the returned position vector and the field variable
	vector_t *vec_x     = pdesol[0];
	vector_t *vec_g     = pdesol[1];
	vector_t *vec_src   = pdesol[2];
	vector_t *vec_g0    = pdesol[3];
	vector_t *vec_err   = pdesol[4];
	vector_t *vec_state = pdesol[5];


	// initializes the analytic (steady-state) solution vector
	vector_t *vec_analytic = vector.zeros(size);


	// defines iterators
	double *x = (vec_x -> array);
	double *g = (vec_g -> array);
	double *err = (vec_err -> array);
	double *analytic = (vec_analytic -> array);


	// computes the analytic (steady-state) solution
	for (size_t i = 0; i != size; ++i)
		analytic[i] = 0.5 * H * (1.0 - x[i]) * (1.0 + x[i]);

	// computes the differences between the exact and numeric solutions
	for (size_t i = 0; i != size; ++i)
		err[i] = (g[i] - analytic[i]);


	// reports the average error
	double norm = vec_err -> qnorm (vec_err);
	// computes the time t
	double t = ( (double) (steps - 1) ) * dt;
	printf("Jacobi(): steady-state solution\n");
	printf("time : %.4e\n", t);
	printf("error: %.4e\n", sqrt(norm) / size );


	// frees the vectors from memory
	vec_x = vector.destroy (vec_x);
	vec_g = vector.destroy (vec_g);
	vec_b = vector.destroy (vec_b);
	vec_g0 = vector.destroy (vec_g0);
	vec_src = vector.destroy (vec_src);
	vec_err = vector.destroy (vec_err);
	vec_state = vector.destroy (vec_state);
	vec_analytic = vector.destroy (vec_analytic);
}


void test_transient_1d_transport_GaussSeidel ()
// solves a transient 1d transport problem via finite-differences
{

	/* defines the solver parameters */

	// defines the steady parameter
	double alpha = 2.0;
	// defines the tolerance of the linear solver
	double tol = 1.0 / ( (double) (0x400000000000) );	// ~1.4e-14
	// defines the maximum number of iterations of the linear solver
	size_t iters = (0x0008FFFF);	// about 500K iterations
	// sets the solver to be verbose
	bool verbose = false;

	// initializes the iterative solver parameters
	isolver_prms_t const prms = {
		.alpha = alpha, .tol = tol,
		.iters = iters, .verbose = verbose
	};


	/* defines the finite-differences problem */


	// sets the number of discretization intervals to 256
	size_t N = (0x00000100);
	// defines the vector size
	size_t size = (N + 1);


	// defines the system domain limits along the x-axis
	double x_l = -1.0, x_u = 1.0;
	// computes the step-size
	double dx = (x_u - x_l) / ( (double) N );
	// computes the time step
	double dt = (dx * dx) / alpha;


	// initializes the placeholder for the solution of the PDE
	vector_t *pdesol[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
	// initializes the position vector
	pdesol[0] = vector.linspace(x_l, x_u, size);
	// initializes the (temperature) field variable
	pdesol[1] = vector.zeros(size);
	// initializes the (heat) source vector
	pdesol[2] = vector.zeros(size);
	// initializes the guess vector
	pdesol[3] = vector.zeros(size);
	// initializes the error vector
	pdesol[4] = vector.zeros(size);
	// initializes the state vector
	pdesol[5] = vector.zeros(1);


	// creates the finite-difference vector (right-hand side of SLE)
	vector_t *pdevec = vector.zeros(size);
	// creates alias for the source vector
	vector_t *vec_b  = pdevec;


	// references the (heat) source vector
	vector_t *vec_source = pdesol[2];
	// defines an iterator for the (heat) source vector
	double *source = (vec_source -> array);

	// defines the non-dimensional volumetric (heat) source term
	double H = 1;
	// defines the constant (heat) source vector
	for (size_t i = 0; i != size; ++i)
		source[i] = -(dx) * (dx) * H;


	/* numeric solution */


	size_t steps = (0x00100000);	// sets to ~1 million time steps
	// solves for the (temperature) field variable iteratively
	for (size_t i = 0; i != steps; ++i)
	{
		vector_t **ret = GaussSeidel (pdesol, pdevec, &prms);
		vector_t *vec_state = ret[5];
		double *state = (vec_state -> array);
		if (state[0] != 0.0)
		{
			printf("Gauss-Seidel solver failed\n");
			printf("try again with different solver params\n");
			break;
		}
	}


	/* post-processing */


	// references the returned position vector and the field variable
	vector_t *vec_x     = pdesol[0];
	vector_t *vec_g     = pdesol[1];
	vector_t *vec_src   = pdesol[2];
	vector_t *vec_g0    = pdesol[3];
	vector_t *vec_err   = pdesol[4];
	vector_t *vec_state = pdesol[5];


	// initializes the analytic (steady-state) solution vector
	vector_t *vec_analytic = vector.zeros(size);


	// defines iterators
	double *x = (vec_x -> array);
	double *g = (vec_g -> array);
	double *err = (vec_err -> array);
	double *analytic = (vec_analytic -> array);


	// computes the analytic (steady-state) solution
	for (size_t i = 0; i != size; ++i)
		analytic[i] = 0.5 * H * (1.0 - x[i]) * (1.0 + x[i]);

	// computes the differences between the exact and numeric solutions
	for (size_t i = 0; i != size; ++i)
		err[i] = (g[i] - analytic[i]);


	// reports the average error
	double norm = vec_err -> qnorm (vec_err);
	// computes the time t
	double t = ( (double) (steps - 1) ) * dt;
	printf("Gauss-Seidel(): steady-state solution\n");
	printf("time : %.4e\n", t);
	printf("error: %.4e\n", sqrt(norm) / size );


	// frees the vectors from memory
	vec_x = vector.destroy (vec_x);
	vec_g = vector.destroy (vec_g);
	vec_b = vector.destroy (vec_b);
	vec_g0 = vector.destroy (vec_g0);
	vec_src = vector.destroy (vec_src);
	vec_err = vector.destroy (vec_err);
	vec_state = vector.destroy (vec_state);
	vec_analytic = vector.destroy (vec_analytic);
}
