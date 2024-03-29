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
vector_t** JacobiSolver1D (
	vector_t *pdesol[6], vector_t*, const isolver_prms_t*
);
vector_t** GaussSeidelSolver1D (
	vector_t *pdesol[6], vector_t*, const isolver_prms_t*
);
void test_steady_1d_transport_Jacobi();
void test_steady_1d_transport_GaussSeidel();
void test_transient_1d_transport_steady_solution_Jacobi ();
void test_transient_1d_transport_steady_solution_GaussSeidel ();
void test_transient_1d_transport_Jacobi ();
void test_transient_1d_transport_GaussSeidel ();
void test_transient_2d_transport_Jacobi ();

int main() {

	/*
	test_pushback();
	test_linspace();
	test_zeros();
	test_ones();
	test_copy();
	test_qnorm();

	// solves the steady 1d transport problem with the Jacobi method
	test_steady_1d_transport_Jacobi();
	// applies the Gauss-Seidel method to solve the same problem
	test_steady_1d_transport_GaussSeidel();

	// solves the transient 1d transport problem with the Jacobi method
	test_transient_1d_transport_steady_solution_Jacobi ();
	// solves the transient problem with the Gauss-Seidel method
	test_transient_1d_transport_steady_solution_GaussSeidel ();

	// solves the transient 1d transport problem with the Jacobi method
	test_transient_1d_transport_Jacobi ();
	// solves the transient problem with the Gauss-Seidel method
	test_transient_1d_transport_GaussSeidel ();

	*/

	// solves the transient 2d transport problem with the Jacobi method
	test_transient_2d_transport_Jacobi ();
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


vector_t** JacobiSolver1D (
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


vector_t** GaussSeidelSolver1D (
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


vector_t** JacobiSolver2D (
	vector_t *pdesol[7], vector_t *pdevec, const isolver_prms_t* prms
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


	// references the components of the position vector
	vector_t *vec_x = pdesol[0];
	vector_t *vec_y = pdesol[1];
	// references the (temperature) field variable
	vector_t *vec_g = pdesol[2];
	// references the (heat) source term vector
	vector_t *vec_src = pdesol[3];
	// references the guess vector
	vector_t *vec_g0  = pdesol[4];
	// references the error vector
	vector_t *vec_err = pdesol[5];
	// references the state vector
	vector_t *vec_state = pdesol[6];
	// references the right-hand side, finite-difference, vector
	vector_t *vec_b = pdevec;


	// gets the number of nodes along the x (y) axis
	size_t nodes = vec_x -> size (vec_x);
	// gets the size of the vector that stores the field
	size_t size = vec_g -> size (vec_g);
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
		b[i] += (alpha * g[i]);


	// initializes the state
	state[0] = 1.0;
	// initializes the norm of the error vector
	double norm = 0;
	// gets the number of discretization intervals
	size_t N = (nodes - 1);
	// updates the diagonal coefficient with the transient contribution
	double c = 1.0 / (alpha + 4.0);
	/* uses Jacobi method to solve for the field variable */
	for (size_t iter = 0; iter != iters; ++iter)
	{

		// initializes g(t + dt, x, y) with the right-hand side
		for (size_t j = 1; j != N; ++j)
		{
			size_t blk = nodes * j;
			for (size_t k = (blk + 1); k != (blk + N); ++k)
				g[k] = b[k];
		}


		// updates g(t + dt, x, y) from sub and super diagonals
		for (size_t j = 1; j != N; ++j)
		{
			size_t blk = nodes * j;
			for (size_t k = (blk + 1); k != (blk + N); ++k)
				g[k] += (g0[k - 1] + g0[k + 1]);
		}


		// updates g(t + dt, x, y) from the lower band
		for (size_t j = 1; j != N; ++j)
		{
			size_t blk = nodes * j;
			for (size_t k = (blk + 1); k != (blk + N); ++k)
				g[k] += g0[k - nodes];
		}


		// updates g(t + dt, x, y) from the upper band
		for (size_t j = 1; j != N; ++j)
		{
			size_t blk = nodes * j;
			for (size_t k = (blk + 1); k != (blk + N); ++k)
				g[k] += g0[k + nodes];
		}


		// completes computation by scaling g(t + dt, x, y)
		for (size_t j = 1; j != N; ++j)
		{
			size_t blk = nodes * j;
			for (size_t k = (blk + 1); k != (blk + N); ++k)
				g[k] *= c;
		}


		// computes the error vector
		for (size_t k = 0; k != size; ++k)
			err[k] = (g[k] - g0[k]);


		// computes the norm of the error vector
		norm = vec_err -> qnorm (vec_err);
		/* checks for convergence */
		if (norm < tol)
		{
			state[0] = 0.0;
			char msg [] = "Jacobi(): solution found "
				"after %lu iters\n";
			if (verbose) printf(msg, iter + 1);
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
	vector_t **ret = JacobiSolver1D (pdesol, pdevec, &prms);


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
	vector_t **ret = GaussSeidelSolver1D (pdesol, pdevec, &prms);


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


void test_transient_1d_transport_steady_solution_Jacobi ()
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
	// initializes the (temperature) field variable g(t = 0, x) = 1
	pdesol[1] = vector.ones (size);
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


	// applies the boundary conditions
	double *gi = (pdesol[1] -> array);
	gi[0] = 0.0;	// g(t, x = x_l) = 0
	gi[N] = 0.0;	// g(t, x = x_u) = 0


	/* numeric solution */


	size_t steps = (0x00100000);	// sets to ~1 million time steps
	// solves for the (temperature) field variable iteratively
	for (size_t i = 0; i != steps; ++i)
	{
		vector_t **ret = JacobiSolver1D (pdesol, pdevec, &prms);
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
	double t = steps * dt;
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


void test_transient_1d_transport_steady_solution_GaussSeidel ()
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
	// initializes the (temperature) field variable g(t = 0, x) = 1
	pdesol[1] = vector.ones (size);
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


	// applies the boundary conditions
	double *gi = (pdesol[1] -> array);
	gi[0] = 0.0;	// g(t, x = x_l) = 0
	gi[N] = 0.0;	// g(t, x = x_u) = 0


	/* numeric solution */


	size_t steps = (0x00100000);	// sets to ~1 million time steps
	// solves for the (temperature) field variable iteratively
	for (size_t i = 0; i != steps; ++i)
	{
		vector_t **ret = GaussSeidelSolver1D (pdesol, pdevec,
						     &prms);
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
	double t = steps * dt;
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


void fpdesol (double t, double H, vector_t *vec_x, vector_t *vec_f)
// computes the analytic (temperature) field f(t, x)
{
	// gets the vector sizes
	size_t size = vec_x -> size (vec_x);
	// references the iterators
	double *x = (vec_x -> array);
	double *f = (vec_f -> array);

	size_t N = 16;
	for (size_t i = 0; i != size; ++i)
	{
		f[i] = 0.0;
		// computes the homogeneous solution
		for (size_t n = 1; n != (N + 1); ++n)
		{
			double pi = M_PI;
			double lambda = .5 * ( (double) (2 * n - 1) ) * pi;
			double lambda2 = (lambda * lambda);
			// computes the nth series coefficient
			double C = (2.0 - 2.0 / lambda2) / lambda;
			double An = (n % 2 == 0)? (-C) : C;
			// defines aliases
			double Ln = lambda, Ln2 = lambda2;
			f[i] += An * cos(Ln * x[i]) * exp(-Ln2 * t);
		}

		// completes the computation by adding the steady solution
		f[i] += 0.5 * H * (1.0 - x[i]) * (1.0 + x[i]);
	}
}


void writeField (FILE *file, vector_t *vec_g, double time, size_t record)
// writes the numeric field g(t, x) at constant time `t' in the data file
{
	// gets read-only access iterators
	const double *g = vec_g -> array;	// numeric  field g(t, x)

	// writes a new line for the current record unless its is the first
	if (record != 0) fprintf(file, "\n");

	// defines the output data format
	char fmt [] = "%18.6e";
	// writes the time `t' at the beginning (of the record)
	fprintf(file, fmt, time);
	// writes the field f(t, x) on the same record
	size_t size = vec_g -> size (vec_g);
	for (size_t i = 0; i != size; ++i)
		fprintf(file, fmt, g[i]);

	// no need to write a new line for we are accounting for that above
}


void write2DField (
	FILE *f, vector_t *vec_g, double time, size_t nodes, size_t record
)
// writes the numeric field g(t, x, y=1/2) at constant time `t' to the file
{
	// writes a new line for the current record unless its is the first
	if (record != 0) fprintf(f, "\n");

	// gets the number of discretization intervals
	size_t N = (nodes - 1);
	// gets the total number of nodes prior to y = 1/2
	size_t blk = nodes * (N / 2);
	// gets read-only access iterator for the numeric field g(t, x, y)
	const double *g = vec_g -> array;

	// defines the output data format
	char fmt [] = "%18.6e";
	// writes the time `t' at the beginning (of the record)
	fprintf(f, fmt, time);
	// writes the field g(t, x, y = 1/2) on the same record
	for (size_t k = blk; k != (blk + nodes); ++k)
		fprintf(f, fmt, g[k]);

	// no need to write a new line for we are accounting for that above
}


void write (FILE *file, vector_t *vec_x, vector_t *vec_f, vector_t *vec_g)
// writes the (temperature) fields with respect to position in a data file
{
	// gets read-only access iterators
	const double *x = vec_x -> array;	// position vector x
	const double *f = vec_f -> array;	// analytic field f(t, x)
	const double *g = vec_g -> array;	// numeric  field g(t, x)

	// defines the format for tabulating the results
	char fmt [] = ("%18.6e %18.6e %18.6e\n");
	size_t size = vec_x -> size (vec_x);
	for (size_t i = 0; i != size; ++i)
		fprintf(file, fmt, x[i], f[i], g[i]);
}


void export (char *name, vector_t *vec_x, vector_t *vec_f, vector_t *vec_g)
// exports the analytic and numeric fields f(t, x) with respect to position
{
	FILE *file = fopen (name, "w");
	if (file != NULL)
	{
		// writes results to data file and closes the file stream
		write  (file, vec_x, vec_f, vec_g);
		fclose (file);
	}
	else
	{
		/*

		Reports to the user that an Input-Output IO error has
		occurred; however, we let the code continue its execution
		so that it can free the memory allocated for all the
		vectors to avert memory leaks.

		*/
		printf("IO ERROR\n");
		printf("aborts export of numeric results\n");
	}
}


void exportSteady2DField (vector_t *vec_x, vector_t *vec_g)
// exports the steady 2d (temperature) field g(x, y=1/2)
{
	// gets the number of nodes along the x (y) axis
	size_t nodes = vec_x -> size (vec_x);
	// gets the number of discretization intervals along the x (y) axis
	size_t N = (nodes - 1);
	// creates read-only access iterators
	const double *x = (vec_x -> array);
	const double *g = (vec_g -> array);
	char pdedata [] = "steady_2d_transport_Jacobi.dat";
	FILE *pdefile = fopen(pdedata, "w");
	// exports the steady-state field at the middle, g(t_ss, x, y = .5)
	if (pdefile)
	{
		size_t j = N / 2;
		for (int i = 0; i != nodes; ++i)
		{
			size_t k = nodes * j + i;
			fprintf(pdefile, "%18.6e %18.6e\n", x[i], g[k]);
		}
		fclose (pdefile);
	}
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
	// initializes the (temperature) field variable g(t = 0, x) = 1
	pdesol[1] = vector.ones (size);
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


	// applies the boundary conditions
	double *gi = (pdesol[1] -> array);
	gi[0] = 0.0;	// g(t, x = x_l) = 0
	gi[N] = 0.0;	// g(t, x = x_u) = 0


	/* numeric solution */


	bool pdestat = true;		// assumes a successful status
	size_t steps = (0x00010000);	// sets to ~60 thousand time steps
	// solves for the (temperature) field variable iteratively
	for (size_t i = 0; i != steps; ++i)
	{
		vector_t **ret = JacobiSolver1D (pdesol, pdevec, &prms);
		vector_t *vec_state = ret[5];
		double *state = (vec_state -> array);
		if (state[0] != 0.0)
		{
			pdestat = false;// sets failure status
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


	// initializes the analytic (transient) solution vector
	vector_t *vec_f = vector.zeros(size);


	// gets iterators
	double *f = (vec_f -> array);
	double *g = (vec_g -> array);
	double *err = (vec_err -> array);


	if (pdestat)
	{
		// computes the time t
		double t = steps * dt;
		// computes the analytic (transient) solution
		fpdesol (t, H, vec_x, vec_f);

		// computes the differences
		for (size_t i = 0; i != size; ++i)
			err[i] = (g[i] - f[i]);

		// reports the average error if the solver was successful
		double norm = vec_err -> qnorm (vec_err);
		printf("Jacobi(): transient solution\n");
		printf("time : %.4e\n", t);
		printf("error: %.4e\n", sqrt(norm) / size );
	}


	// exports the fields f(t, x) if the solver was successful
	char fname [] = "pdesolJacobi.dat";
	if (pdestat) export (fname, vec_x, vec_f, vec_g);


	// frees the vectors from memory
	vec_x = vector.destroy (vec_x);
	vec_f = vector.destroy (vec_f);
	vec_g = vector.destroy (vec_g);
	vec_b = vector.destroy (vec_b);
	vec_g0 = vector.destroy (vec_g0);
	vec_src = vector.destroy (vec_src);
	vec_err = vector.destroy (vec_err);
	vec_state = vector.destroy (vec_state);
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
	// initializes the (temperature) field variable g(t = 0, x) = 1
	pdesol[1] = vector.ones (size);
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


	// applies the boundary conditions
	double *gi = (pdesol[1] -> array);
	gi[0] = 0.0;	// g(t, x = x_l) = 0
	gi[N] = 0.0;	// g(t, x = x_u) = 0


	/* numeric solution */

	size_t lines = 0;		// initializes line counter
	bool pdestat = true;		// assumes a success status
	size_t steps = (0x00010000);	// sets to ~60 thousand time steps
	char pdedata [] = "transient_Gauss-Seidel.dat";
	FILE *pdefile = fopen(pdedata, "w");
	// solves for the (temperature) field variable iteratively
	for (size_t i = 0; i != steps; ++i)
	{
		// writes g(t, x) to the data file every 256 steps
		if (pdefile && i % (0x00000100) == 0)
		{
			vector_t *vec_g = pdesol[1];
			double time = ( (double) i ) * dt;
			writeField (pdefile, vec_g, time, lines);
			++lines;
		}

		vector_t **ret = GaussSeidelSolver1D (pdesol, pdevec,
						     &prms);
		vector_t *vec_state = ret[5];
		double *state = (vec_state -> array);
		if (state[0] != 0.0)
		{
			pdestat = false;// sets failure status
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


	// initializes the analytic (transient) solution vector
	vector_t *vec_f = vector.zeros(size);


	// gets iterators
	double *f = (vec_f -> array);
	double *g = (vec_g -> array);
	double *err = (vec_err -> array);


	if (pdestat)
	{
		// computes the time t
		double t = steps * dt;
		// computes the analytic (transient) solution
		fpdesol (t, H, vec_x, vec_f);

		// computes the differences
		for (size_t i = 0; i != size; ++i)
			err[i] = (g[i] - f[i]);

		// reports the average error if the solver was successful
		double norm = vec_err -> qnorm (vec_err);
		printf("Gauss-Seidel(): transient solution\n");
		printf("time : %.4e\n", t);
		printf("error: %.4e\n", sqrt(norm) / size );
	}


	// exports the fields f(t, x) if the solver was successful
	char fname [] = "pdesolGauss-Seidel.dat";
	if (pdestat) export (fname, vec_x, vec_f, vec_g);


	// frees the vectors from memory
	vec_x = vector.destroy (vec_x);
	vec_f = vector.destroy (vec_f);
	vec_g = vector.destroy (vec_g);
	vec_b = vector.destroy (vec_b);
	vec_g0 = vector.destroy (vec_g0);
	vec_src = vector.destroy (vec_src);
	vec_err = vector.destroy (vec_err);
	vec_state = vector.destroy (vec_state);

	// closes the data file that stores the transient data
	if (pdefile) fclose(pdefile);
}


void test_transient_2d_transport_Jacobi ()
// solves a transient 2d transport problem via finite-differences
{

	/* defines the solver parameters */

	// defines the transient parameter
	double alpha = 2.0;
	// defines the tolerance of the linear solver
	double tol = 1.0 / ( (double) (0x400000000000) );	// ~1.4e-14
	// defines the maximum number of iterations of the linear solver
	size_t iters = (0x0008FFFF);	// about 500K iterations
	// silences the linear solver
	bool verbose = false;

	// initializes the iterative solver parameters
	isolver_prms_t const prms = {
		.alpha = alpha, .tol = tol,
		.iters = iters, .verbose = verbose
	};


	/* defines the finite-differences problem */


	// sets the number of discretization intervals to 256
	size_t N = (0x00000100);
	// defines the number of nodes along the x (y) axis
	size_t nodes = (N + 1);
	// defines the system size (total number of nodes)
	size_t size = nodes * nodes;


	// defines the system domain limits
	double x_l = 0.0, x_u = 1.0;
	double y_l = x_l, y_u = x_u;
	// computes the step-size along the x (y) axis
	double dx = (x_u - x_l) / ( (double) N );
	double dy = dx;
	// computes the time step
	double dt = (dx * dx) / alpha;


	// initializes the placeholder for the solution of the PDE
	vector_t *pdesol[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
	// initializes the components of the position vector
	pdesol[0] = vector.linspace(x_l, x_u, nodes);
	pdesol[1] = vector.linspace(y_l, y_u, nodes);
	/*

	sets the boundary conditions on the (temperature) field variable:

	g(t, x_l, y) = g(t, x_u, y) = g(t, x, y_l) = g(t, x, y_u) = 0

	*/
	pdesol[2] = vector.zeros(size);
	// initializes the (heat) source vector
	pdesol[3] = vector.zeros(size);
	// initializes the guess vector
	pdesol[4] = vector.zeros(size);
	// initializes the error vector
	pdesol[5] = vector.zeros(size);
	// initializes the state vector
	pdesol[6] = vector.zeros(1);


	// creates the finite-difference vector (right-hand side of SLE)
	vector_t *pdevec = vector.zeros(size);
	// creates alias for the finite-difference vector
	vector_t *vec_b  = pdevec;


	// references the (heat) source vector
	vector_t *vec_source = pdesol[3];
	// defines an iterator for the (heat) source vector
	double *source = (vec_source -> array);

	// defines the non-dimensional volumetric (heat) source term
	double H = 1;
	// defines the constant (heat) source vector
	for (size_t i = 0; i != size; ++i)
		source[i] = (dx) * (dx) * H;


	// applies the initial condition g(t = 0, x, y) = 1
	double *gi = (pdesol[2] -> array);
	for (size_t j = 1; j != N; ++j)
	{
		for (size_t i = 1; i != N; ++i)
			gi[j * nodes + i] = 1.0;
	}


	/* numeric solution */

	size_t lines = 0;		// initializes line counter
	bool pdestat = true;		// assumes a success status
	size_t steps = (0x00010000);	// sets to ~60 thousand time steps
	char pdedata [] = "transient_2d_transport_Jacobi.dat";
	FILE *pdefile = fopen(pdedata, "w");
	// solves for the (temperature) field variable iteratively
	for (size_t i = 0; i != steps; ++i)
	{
		// writes g(t, x, y = 0.5) to the data file every 256 steps
		if (pdefile && i % (0x00000100) == 0)
		{
			vector_t *vec_g = pdesol[2];
			double time = ( (double) i ) * dt;
			write2DField (pdefile, vec_g, time, nodes, lines);
			++lines;
		}

		vector_t **ret = JacobiSolver2D (pdesol, pdevec, &prms);
		vector_t *vec_state = ret[6];
		double *state = (vec_state -> array);
		if (state[0] != 0.0)
		{
			pdestat = false;// sets failure status
			printf("Jacobi 2d solver failed\n");
			printf("try again with different solver params\n");
			break;
		}
	}


	/* post-processing */


	// references the position vector, the field variable, etc.
	vector_t *vec_x     = pdesol[0];
	vector_t *vec_y     = pdesol[1];
	vector_t *vec_g     = pdesol[2];
	vector_t *vec_src   = pdesol[3];
	vector_t *vec_g0    = pdesol[4];
	vector_t *vec_err   = pdesol[5];
	vector_t *vec_state = pdesol[6];


	// exports the steady 2D field g(t, x, y = 0.5)
	exportSteady2DField (vec_x, vec_g);


	// frees the vectors from memory
	vec_x = vector.destroy (vec_x);
	vec_y = vector.destroy (vec_y);
	vec_g = vector.destroy (vec_g);
	vec_b = vector.destroy (vec_b);
	vec_g0 = vector.destroy (vec_g0);
	vec_src = vector.destroy (vec_src);
	vec_err = vector.destroy (vec_err);
	vec_state = vector.destroy (vec_state);

	// closes the data file that stores the transient data
	if (pdefile) fclose(pdefile);
}
