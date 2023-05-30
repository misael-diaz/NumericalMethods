/*
 * Transient Heat Conduction                                		May 30, 2023
 *
 * source: test.c
 * author: @misael-diaz
 *
 * Synopsis:
 * Solves the transient Poission equation iteratively with the Jacobi method
 * until the steady state is reached.
 *
 * The objective is solve the problem by masking boundary nodes. Again, it is
 * important to check which loops GCC vectorizes. The code has been written to
 * keep the vectorization that was obtained without masking. I placed inline
 * comments to indicate which loops were vectorized, as in the previous code.
 *
 *
 * Copyright (c) 2023 Misael Diaz-Maldonado
 * This file is released under the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 *
 * References:
 * [0] A Koenig and B Moo, Accelerated C++ Practical Programming by Example.
 * [1] JJ McConnell, Analysis of Algorithms, 2nd edition
 * [2] A Gilat and V Subramaniam, Numerical Methods for Engineers and
 *     Scientists: An Introduction with Applications using MATLAB
 * [3] H Wendland, Numerical Linear Algebra
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define NODE 0.0
#define SIZE 1024
#define ALPHA 2.0
#define TOLERANCE 8.673617379884035e-19
#define MAX_ITERATIONS 128
#define SUCCESS_STATE 0
#define FAILURE_STATE 1


typedef struct {
  double* x;		// position array along the x-axis
  double* f;		// exact solution array, f(t, x)
  double* g;		// estimate of the solution array, g(t + dt)
  double* g0;		// previous estimate of the solution array, g(t + dt)
  double* err;		// error array
  double* rhs;		// Right Hand Side RHS array of the PDE
  double* mask;		// mask array, one if a boundary node zero otherwise
  size_t size;		// array size
  size_t state;		// solver state
} workspace_t;


void Poisson();

int main ()
{
  Poisson();
  return 0;
}


double* alloc (size_t size)
{
  double* x = malloc( size * sizeof(double) );
  return x;
}


double* dealloc (double* x)
{
  if (x == NULL)
  {
    return x;
  }

  free(x);

  return x;
}


void zeros (size_t size, double* x)	// numpy-like zeros
{
  for (size_t i = 0; i != size; ++i)
  {
    x[i] = 0.0;
  }
}


void ones (size_t size, double* x)	// numpy-like ones
{
  for (size_t i = 0; i != size; ++i)
  {
    x[i] = 1.0;
  }
}


void linspace (double* x, double x_i, double x_f, size_t size)	// numpy-like linspace
{
  double N = (size - 1);
  double dx = (x_f - x_i) / N;
  for (size_t i = 0; i != size; ++i)
  {
    x[i] = x_i + ( (double) i ) * dx;
  }
}


double norm (size_t size, const double* x)	// gets the norm of the vector
{
  double sum = 0;
  for (size_t i = 0; i != size; ++i)
  {
    sum += (x[i] * x[i]);
  }

  return sum;
}


void copy (size_t size, const double* restrict src, double* restrict dst)
{
  for (size_t i = 0; i != size; ++i)
  {
    dst[i] = src[i];
  }
}


workspace_t* create (size_t size)
{
  workspace_t *workspace = malloc( sizeof(workspace_t) );

  if (workspace == NULL)
  {
    return workspace;
  }

  workspace -> x = alloc(size);
  workspace -> f = alloc(size);
  workspace -> g = alloc(size);
  workspace -> g0 = alloc(size);
  workspace -> err = alloc(size);
  workspace -> rhs = alloc(size);
  workspace -> mask = alloc(size);

  double* x = workspace -> x;
  double* f = workspace -> f;
  double* g = workspace -> g;
  double* g0 = workspace -> g0;
  double* err = workspace -> err;
  double* rhs = workspace -> rhs;
  double* mask = workspace -> mask;

  zeros(size, x);
  zeros(size, f);
  zeros(size, g);
  zeros(size, g0);
  zeros(size, err);
  zeros(size, rhs);
  zeros(size, mask);

  workspace -> size = size;
  workspace -> state = 1;

  return workspace;
}


workspace_t* destroy (workspace_t* workspace)
{
  if (workspace == NULL)
  {
    return workspace;
  }

  dealloc(workspace -> x);
  dealloc(workspace -> f);
  dealloc(workspace -> g);
  dealloc(workspace -> g0);
  dealloc(workspace -> err);
  dealloc(workspace -> rhs);
  dealloc(workspace -> mask);

  free(workspace);

  return workspace;
}


void set_state (workspace_t* workspace, size_t state)
{
  workspace -> state = state;
}


size_t get_state (workspace_t* workspace)
{
  return (workspace -> state);
}


// void init_rhs (size_t size, double* b, double* x, double* g)
//
// Synopsis:
// Initializes the Right-Hand-Size RHS of the Finite Difference Equations FDEs.
//
// Inputs:
// size		array size (same for both `x' and `g')
// x		x-axis position array
// g		solution array g(t) at the current time
//
// Outputs:
// b		RHS array


void init_rhs(size_t size,
	      double* restrict b,
	      const double* restrict x,
	      const double* restrict g)
{
  double alpha = ALPHA;
  double dx = (x[1] - x[0]);
  // we don't need to mask this loop since we are masking the boundary nodes elsewhere
  for (size_t i = 0; i != size; ++i)
  {
    b[i] = (dx * dx + alpha * g[i]);
  }
}


// void rhs(size_t size, double* g, double* b)
//
// Synopsis:
// Updates the solution array g(t + dt) from the RHS of the FDEs.
//
// Inputs:
// size		array size (same for both `g' and `b')
// b		RHS array
//
// Outputs:
// g		estimate of the solution array g(t + dt)


void rhs (size_t size,
	  double* restrict g,
	  const double* restrict b,
	  const double* restrict mask)
{
  for (size_t i = 0; i != size; ++i)
  {
    double const m = mask[i];
    double const value = b[i];
    double const elem = (m == NODE)? value : 0.0;
    g[i] = elem;
  }
}


// void tridiag(size_t size, double* g, double* g0)
//
// Synopsis:
// Updates the solution array g(t + dt) from the tridiagonal terms of the FDEs.
//
// Inputs:
// size		array size (same for both `g' and `g0')
// g0		previous estimate of the solution array g(t + dt)
//
// Outputs:
// g		estimate of the solution array g(t + dt)


void tridiag (size_t size,
	      double* restrict g,
	      const double* restrict g0,
	      const double* restrict mask)
{
  for (size_t i = 0; i != size; ++i)
  {
    double const m = mask[i];
    double const elem = g0[i - 1];
    g[i] += (m == NODE)? elem : 0.0;

  }

  for (size_t i = 0; i != size; ++i)
  {
    double const m = mask[i];
    double const elem = g0[i + 1];
    g[i] += (m == NODE)? elem : 0.0;
  }
}


// void scale (size_t size, double* g)
//
// Synopsis:
// Scales the solution array g(t + dt) with the main diagonal coefficient.
//
// Inputs:
// size		array size
//
// Outputs:
// g		estimate of the solution array g(t + dt)


void __attribute__ ((noinline)) scale(size_t size,
				      double* restrict g,
				      const double* restrict mask)
{
  double alpha = ALPHA;
  double c = 1.0 / (alpha + 2.0);
  for (size_t i = 0; i != size; ++i)
  {
    double const m = mask[i];
    double const elem = (m == NODE)? c : 1.0;
    g[i] *= elem;
  }
}


// void exact(size_t size, double* g, double* x)
//
// Synopsis:
// Obtains the exact steady-state solution of the Poission equation.
//
// Inputs:
// size		array size (same for both `x' and `g')
// x		x-axis position array
//
// Outputs:
// g		steady-state solution array g(t -> infinity)


void exact (size_t size, double* restrict g, const double* restrict x)
{
  for (size_t i = 0; i != size; ++i)
  {
    g[i] = 0.5 * (1.0 - x[i]) * (1.0 + x[i]);
  }
}


// void error(size_t size, double* e, double* v, double* w)
//
// Synopsis:
// Computes the elementwise differences of the vectors `v' and `w'.
//
// Inputs:
// size		array size (same for `v', `w', and `err')
// v		array
// w		array
//
// Outputs:
// e		error array, stores the elementwise differences


void error (size_t size,
	    double* restrict e,
	    const double* restrict v,
	    const double* restrict w)
{
  for (size_t i = 0; i != size; ++i)
  {
    e[i] = (v[i] - w[i]);
  }
}


// void solver (workspace_t* workspace)
//
// Synopsis:
// Solves the transient Poisson FDEs with the Jacobi method:
//
// 	-g(n + 1, i - 1) + (2 + alpha) * g(n + 1, i) - g(n + 1, i + 1) = b(i),
//
// where b(i) = dx * dx + alpha * g(n, i), dx is the distance between nodes,
// alpha = dx * dx / dt, where dt is the time-step, `n' is time index,
// `i' is the position index, and g is the (temperature) field.
//
// The solver finds the solution g(t + dt) at the next time step from the current g(t).
//
// Inputs:
// workspace	data structure containing the current data (field, position, error, etc.)
//
// Outputs:
// workspace	stores the data at the next time step


void solver (workspace_t* workspace)
{

  // parameters:

  double tol = TOLERANCE;
  size_t iters = MAX_ITERATIONS;
  size_t size = workspace -> size;

  // iterators:

  double* x = workspace -> x;
  double* b = workspace -> rhs;
  double* g = workspace -> g;
  double* g0 = workspace -> g0;
  double* err = workspace -> err;
  const double* mask = workspace -> mask;

  // initializations:

  init_rhs(size, b, x, g);			// vectorized by gcc

  // Jacobi solver:

  set_state(workspace, FAILURE_STATE);
  for (size_t i = 0; i != iters; ++i)
  {
    // updates the solution array g(t + dt):
    rhs(size, g, b, mask);			// vectorized by gcc
    tridiag(size, g, g0, mask);			// vectorized by gcc
    scale(size, g, mask);			// vectorized by gcc

    // checks for convergence:
    error(size, err, g, g0);			// vectorized by gcc
    if (norm(size, err) < tol)
    {
      set_state(workspace, SUCCESS_STATE);
      //printf("Jacobi(): solution found after %lu iters\n", i + 1);
      break;
    }

    // updates the initial solution array g(t + dt) for the next iteration:
    copy(size, g, g0);				// optimized by gcc
  }
}


void pdesol (double t, workspace_t* workspace)	// computes the exact field f(t, x)
{
  size_t size = workspace -> size;
  const double *x = workspace -> x;
  double *f = workspace -> f;

  size_t N = 256;
  for (size_t i = 0; i != size; ++i)		// obtains transient contribution f(t, x)
  {
    f[i] = 0.0;
    for (size_t n = 1; n != (N + 1); ++n)
    {
      double pi = M_PI;
      double lambda = 0.5 * ( (2.0 * ( (double) n ) - 1.0) ) * pi;
      double lambda2 = (lambda * lambda);
      double C = ( 2.0 * (1.0 - 1.0 / lambda2) ) / lambda;
      double An = ( (n % 2 == 0)? -C : C );
      double Ln = lambda, Ln2 = lambda2;
      f[i] += An * cos(Ln * x[i]) * exp(-Ln2 * t);
    }

    f[i] += 0.5 * (1.0 - x[i]) * (1.0 + x[i]);	// adds steady-state contribution f_ss(x)
  }
}


// void integrator(workspace_t* workspace)
//
// Synopsis:
// Integrates the transient Poisson FDEs via time implicit scheme.
//
// Inputs:
// workspace	data structure containing the initial data (field, position, error, etc.)
//
// Outputs:
// workspace	stores the data after completing the integration


void integrator (workspace_t* workspace)
{
  int steps = 0x00400000;
  for (int i = 0; i != steps; ++i)
  {
    solver(workspace);

    size_t state = get_state(workspace);
    if (state == FAILURE_STATE)
    {
      char msg[] = "Jacobi solver failed to converge to the solution after %d iterations";
      printf(msg, MAX_ITERATIONS);
      break;
    }

    int span = (steps / 16);
    // logs error of exact f(t+dt, x) and numeric solution g(t+dt, x) every `span' steps
    if ( ( i != 0 ) && ( (i % span) == 0 ) )
    {
      double alpha = ALPHA;
      const double* x = workspace -> x;
      double const dx = (x[1] - x[0]);
      double const dt = (dx * dx) / alpha;
      double const t = ( ( (double) i ) + 1.0 ) * dt;
      pdesol(t, workspace);

      size_t const size = workspace -> size;
      const double* f = workspace -> f;
      const double* g = workspace -> g;
      double* err = workspace -> err;
      error(size, err, f, g);
      double e = sqrt( norm(size, err) );
      printf("approximation error (transient solution t = %.4e): %e \n", t, e);
    }
  }
}


// void Poisson()
//
// Synopsis:
// Solves the 1d transient Poission equation. The initial temperature field is zero
// everywhere. There is a uniform heat source throughout the domain. The temperature
// at the system boundaries x = [-1, 1] is zero.


void Poisson ()
{
  // parameters:

  size_t size = SIZE;				// number of elements in array
  size_t N = (size - 1);			// number of discretization intervals

  // memory allocations:

  workspace_t* workspace = create(size);

  // iterators:

  double* x = workspace -> x;
  double* g = workspace -> g;
  double* g0 = workspace -> g0;
  double* err = workspace -> err;
  double* mask = workspace -> mask;

  // initializations:

  zeros(size, g0);				// inits the initial solution array
  ones(size, g);				// inits the solution array g(t=0, x) = 1
  g[0] = 0.0;					// applies BC g(t, x = -1) = 0
  g[N] = 0.0;					// applies BC g(t, x = +1) = 0

  zeros(size, mask);
  mask[0] = 1.0;
  mask[N] = 1.0;

  double x_l = -1;				// defines x-axis lower bound
  double x_u =  1;				// defines x-axis upper bound
  linspace(x, x_l, x_u, size);			// inits the x-axis position array

  // solves the transient Poisson equation:

  integrator(workspace);

  // post-processing:

  // gets the exact solution at steady-state
  exact(size, g0, x);

  // reports the approximation error
  error(size, err, g, g0);
  double e = sqrt( norm(size, err) );
  printf("approximation error (steady-state solution): %e \n", e);

  // memory deallocations:

  workspace = destroy(workspace);
}
