/*
 * Transient Heat Conduction                                	June 03, 2023
 *
 * source: test.c
 * author: @misael-diaz
 *
 * Synopsis:
 * Solves the 2d transient Poisson equation iteratively with the Jacobi method
 * until the steady state is reached.
 *
 * The objective is solve the problem by masking boundary nodes, for this simplifies
 * the programming and allows for solving the Poisson equation even for systems with
 * complicated geometries.
 *
 * Vectorized code by GCC is indicated in the source. We use a bitmask to express loops
 * with simple nested if-else conditionals into loops that GCC can auto-vectorize.
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


#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define SIZE 128
#define ALPHA 2.0
#define iNODE 0xffffffffffffffff
#define TOLERANCE 8.673617379884035e-19
#define MAX_ITERATIONS 128
#define SUCCESS_STATE 0
#define FAILURE_STATE 1
#define VERBOSE false


// we use this for auto-vectorization of loops with simple nested conditionals
typedef union
{
  uint64_t bin;		// bit pattern of the floating point data
  double data;		// floating point data
} alias_t;


typedef struct {
  double* x;		// position array along the x-axis
  double* y;		// position array along the y-axis
  double* f;		// exact solution array, f(t, x, y)
  double* g;		// estimate of the solution array, g(t + dt, x, y)
  double* g0;		// previous estimate of the solution array, g(t + dt, x, y)
  double* err;		// error array
  double* rhs;		// Right Hand Side RHS array of the PDE
  double* mask;		// bitmask is zero at boundary nodes and `iNODE' at interior ones
  size_t size;		// array size
  size_t state;		// solver state
} workspace_t;


void Poisson();

int main ()
{
  Poisson();
  return 0;
}


double* alloc (size_t const size)
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


void zeros (size_t const size, double* x)	// numpy-like zeros
{
  for (size_t i = 0; i != size; ++i)
  {
    x[i] = 0.0;
  }
}


void ones (size_t const size, double* x)	// numpy-like ones
{
  for (size_t i = 0; i != size; ++i)
  {
    x[i] = 1.0;
  }
}


void add (size_t const size, double* restrict dst, const double* restrict src)
{
  for (size_t i = 0; i != size; ++i)
  {
    dst[i] += src[i];
  }
}


void mult (size_t const size, double* restrict dst, const double* restrict src)
{
  for (size_t i = 0; i != size; ++i)
  {
    dst[i] *= src[i];
  }
}


void linspace (double* x, double x_i, double x_f, size_t const size)// numpy-like linspace
{
  double const N = (size - 1);
  double const dx = (x_f - x_i) / N;
  for (size_t i = 0; i != size; ++i)
  {
    x[i] = x_i + ( (double) i ) * dx;
  }
}


double norm (size_t const size, const double* x)	// gets the norm of the vector
{
  double sum = 0;
  for (size_t i = 0; i != size; ++i)
  {
    sum += (x[i] * x[i]);
  }

  return sum;
}


void copy (size_t const size, const double* restrict src, double* restrict dst)
{
  for (size_t i = 0; i != size; ++i)
  {
    dst[i] = src[i];
  }
}


workspace_t* create (size_t const size)
{
  workspace_t *workspace = malloc( sizeof(workspace_t) );

  if (workspace == NULL)
  {
    return workspace;
  }

  size_t const size2 = (size * size);
  workspace -> x = alloc(size);
  workspace -> y = alloc(size);
  workspace -> f = alloc(size2);
  workspace -> g = alloc(size2);
  workspace -> g0 = alloc(size2);
  workspace -> err = alloc(size2);
  workspace -> rhs = alloc(size2);
  workspace -> mask = alloc(size2);

  double* x = workspace -> x;
  double* y = workspace -> y;
  double* f = workspace -> f;
  double* g = workspace -> g;
  double* g0 = workspace -> g0;
  double* err = workspace -> err;
  double* rhs = workspace -> rhs;
  double* mask = workspace -> mask;

  zeros(size, x);
  zeros(size, y);
  zeros(size2, f);
  zeros(size2, g);
  zeros(size2, g0);
  zeros(size2, err);
  zeros(size2, rhs);
  zeros(size2, mask);

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
  dealloc(workspace -> y);
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


// void init_field (size_t size, double* g)
//
// Synopsis:
// Sets the initial (temperature) field g(t = 0, x, y). Fills with ones at the interior
// nodes and with zeros at the boundary nodes.
//
// Input:
// size		number of nodes along the x [y] axis
//
// Output:
// g		2d field array, zeros at the boundaries and ones elsewhere


void init_field (size_t const size, double* g)
{
  size_t const numel = (size * size);
  ones(numel, g);				// sets g(t = 0, x, y) = 1 (everywhere)

  for (int i = 0; i != size; ++i)
  {
    int const j = 0;
    int const k = (i + size * j);
    g[k] = 0.0;					// sets g(t = 0, x, y = 0) = 0
  }

  for (int i = 0; i != size; ++i)
  {
    int const j = (size - 1);
    int const k = (i + size * j);
    g[k] = 0.0;					// sets g(t = 0, x, y = 1) = 0
  }

  for (int j = 0; j != size; ++j)
  {
    int const i = 0;
    int const k = (i + size * j);
    g[k] = 0.0;					// sets g(t = 0, x = 0, y) = 0
  }

  for (int j = 0; j != size; ++j)
  {
    int const i = (size - 1);
    int const k = (i + size * j);
    g[k] = 0.0;					// sets g(t = 0, x = 1, y) = 0
  }
}


// void init_mask(size_t size, double* mask)
//
// Synopsis:
// Initializes node mask.
// The mask is equal to zero at boundary nodes and equal to `iNODE' at interior nodes.
// This is so that boundary node data is kept constant while updating all the nodes.
//
// Input:
// size		number of nodes along the x [y] axis
//
// Output:
// mask		array with `(size * size)'  number of elements


void init_mask (size_t const size, double* mask)
{
  alias_t* m = mask;
  size_t const numel = (size * size);
  for (size_t i = 0; i != numel; ++i)
  {
    m[i].bin = iNODE;				// sets all nodes as interior nodes
  }

  for (int i = 0; i != size; ++i)
  {
    int const j = 0;
    int const k = (i + size * j);
    mask[k] = 0.0;				// sets nodes at y = 0 as boundary nodes
  }

  for (int i = 0; i != size; ++i)
  {
    int const j = (size - 1);
    int const k = (i + size * j);
    mask[k] = 0.0;				// sets nodes at y = 1 as boundary nodes
  }

  for (int j = 0; j != size; ++j)
  {
    int const i = 0;
    int const k = (i + size * j);
    mask[k] = 0.0;				// sets nodes at x = 0 as boundary nodes
  }

  for (int j = 0; j != size; ++j)
  {
    int const i = (size - 1);
    int const k = (i + size * j);
    mask[k] = 0.0;				// sets nodes at x = 1 as boundary nodes
  }
}


// void init_rhs (size_t size, double* b, double* x, double* g)
//
// Synopsis:
// Initializes the Right-Hand-Size RHS of the Finite Difference Equations FDEs.
//
// Inputs:
// size		number of nodes along the x [y] axis
// x		x-axis position array
// g		solution array g(t) at the current time
//
// Outputs:
// b		RHS array


void init_rhs(size_t const size,
	      double* restrict b,
	      const double* restrict x,
	      const double* restrict g)
{
  double const alpha = ALPHA;
  double const dx = (x[1] - x[0]);
  size_t const numel = (size * size);
  // we don't need to mask this loop since we are masking the boundary nodes elsewhere
  for (size_t i = 0; i != numel; ++i)
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
// size		number of nodes along the x [y] axis
// b		RHS array
// mask		bitmask for protecting boundary node data
//
// Outputs:
// g		estimate of the solution array g(t + dt)


void rhs (size_t const size,
	  double* restrict g,
	  const double* restrict b,
	  const double* restrict mask)
{
  alias_t* dst = g;
  const alias_t* values = b;
  const alias_t* masks = mask;
  size_t const numel = (size * size);
  for (size_t i = 0; i != numel; ++i)
  {
    dst[i].bin = (masks[i].bin & values[i].bin);
  }
}


// void tridiag(size_t size, double* g, double* g0, double* mask)
//
// Synopsis:
// Updates the solution array g(t + dt) from the tridiagonal terms of the FDEs.
//
// Inputs:
// size		number of nodes along the x [y] axis
// g0		previous estimate of the solution array g(t + dt)
// tmp		array temporary for storing intermediate computations
// mask		bitmask for protecting boundary node data
//
// Outputs:
// g		estimate of the solution array g(t + dt)


void tridiag (size_t const size,
	      double* restrict g,
	      double* restrict tmp,
	      const double* restrict g0,
	      const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* values = g0;
  const alias_t* masks = mask;
  size_t const numel = (size * size);

  zeros(numel, tmp);

  for (size_t i = 1; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values[i - 1].bin);
  }

  add(numel, g, tmp);
  /*
  for (size_t i = 0; i != numel; ++i)
  {
    g[i] += tmp[i];
  }
  */

  for (size_t i = 0; i != (numel - 1); ++i)
  {
    t[i].bin = (masks[i].bin & values[i + 1].bin);
  }

  add(numel, g, tmp);
  /*
  for (size_t i = 0; i != numel; ++i)
  {
    g[i] += tmp[i];
  }
  */
}


// void subdiag(size_t size, double* g, double* g0, double* mask)
//
// Synopsis:
// Updates the solution array g(t + dt) from the sub-diagonal terms of the FDEs.
//
// Inputs:
// size		number of nodes along the x [y] axis
// g0		previous estimate of the solution array g(t + dt)
// tmp		array temporary for storing intermediate computations
// mask		bitmask for protecting boundary node data
//
// Outputs:
// g		estimate of the solution array g(t + dt)


void subdiag (size_t const size,
	      double* restrict g,
	      double* restrict tmp,
	      const double* restrict g0,
	      const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* values = g0;
  const alias_t* masks = mask;
  size_t const numel = (size * size);

  zeros(numel, tmp);

  for (size_t i = size; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values[i - size].bin);
  }

  add(numel, g, tmp);
  /*
  for (size_t i = 0; i != numel; ++i)
  {
    g[i] += tmp[i];
  }
  */
}


// void superdiag(size_t size, double* g, double* g0, double* mask)
//
// Synopsis:
// Updates the solution array g(t + dt) from the super-diagonal terms of the FDEs.
//
// Inputs:
// size		number of nodes along the x [y] axis
// g0		previous estimate of the solution array g(t + dt)
// tmp		array temporary for storing intermediate computations
// mask		bitmask for protecting boundary node data
//
// Outputs:
// g		estimate of the solution array g(t + dt)


void superdiag (size_t const size,
		double* restrict g,
		double* restrict tmp,
		const double* restrict g0,
		const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* values = g0;
  const alias_t* masks = mask;
  size_t const numel = (size * size);

  zeros(numel, tmp);

  for (size_t i = 0; i != (numel - size); ++i)
  {
    t[i].bin = (masks[i].bin & values[i + size].bin);
  }

  add(numel, g, tmp);
  /*
  for (size_t i = 0; i != numel; ++i)
  {
    g[i] += tmp[i];
  }
  */
}


// void scale (size_t size, double* g)
//
// Synopsis:
// Scales the solution array g(t + dt) with the main diagonal coefficient.
//
// Inputs:
// size		number of nodes along the x [y] axis
// tmp		array temporary for storing intermediate computations
// mask		bitmask for protecting boundary node data
//
// Outputs:
// g		estimate of the solution array g(t + dt)


void __attribute__ ((noinline)) scale(size_t const size,
				      double* restrict g,
				      double* restrict tmp,
				      const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* masks = mask;
  double const alpha = ALPHA;
  double const c = 1.0 / (alpha + 4.0);
  alias_t const values = { .data = c };
  size_t const numel = (size * size);

  for (size_t i = 0; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values.bin);
  }

  mult(numel, g, tmp);
  /*
  for (size_t i = 0; i != numel; ++i)
  {
    g[i] *= tmp[i];
  }
  */
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


void error (size_t const size,
	    double* restrict e,
	    const double* restrict v,
	    const double* restrict w)
{
  for (size_t i = 0; i != size; ++i)
  {
    e[i] = (v[i] - w[i]);
  }
}


// void solver(workspace_t* workspace)
//
// Synopsis:
// Solves the transient 2D Poisson FDEs with the Jacobi method.
// The solver finds the solution g(t + dt) at the next time step from the current g(t).
//
// Inputs:
// workspace	data structure containing the current data (field, position, error, etc.)
//
// Outputs:
// workspace	updates the fields g0 and g with the data of the next time step


void solver (workspace_t* workspace)
{

  // parameters:

  double const tol = TOLERANCE;
  size_t const iters = MAX_ITERATIONS;
  size_t const size = workspace -> size;
  size_t const numel = (size * size);

  // iterators:

  double* x = workspace -> x;
//double* y = workspace -> y;
  double* b = workspace -> rhs;
  double* g = workspace -> g;
  double* g0 = workspace -> g0;
  double* err = workspace -> err;
  double* tmp = workspace -> f;
  const double* mask = workspace -> mask;

  // initializations:

  init_rhs(size, b, x, g);			// vectorized by gcc

  // Jacobi solver:

  set_state(workspace, FAILURE_STATE);
  for (size_t i = 0; i != iters; ++i)
  {
    // updates the solution array g(t + dt):
    rhs(size, g, b, mask);			// vectorized by gcc
    tridiag(size, g, tmp, g0, mask);		// vectorized by gcc
    subdiag(size, g, tmp, g0, mask);		// vectorized by gcc
    superdiag(size, g, tmp, g0, mask);		// vectorized by gcc
    scale(size, g, tmp, mask);			// vectorized by gcc

    // checks for convergence:
    error(numel, err, g, g0);			// vectorized by gcc
    if (norm(numel, err) < tol)
    {
      set_state(workspace, SUCCESS_STATE);
      if (VERBOSE) printf("Jacobi(): solution found after %lu iters\n", i + 1);
      break;
    }

    // updates the initial solution array g(t + dt) for the next iteration:
    copy(numel, g, g0);				// optimized by gcc
  }
}


// void exact(size_t size, double* g, double* x)
//
// Synopsis:
// Obtains the exact steady-state solution of the Poisson equation.
//
// Inputs:
// size		number of nodes along the x [y] axis
// x		x-axis position array
//
// Outputs:
// g		steady-state solution array g(t -> infinity) (numel: size * size)


void exact (size_t const size,
	    double* restrict g,
	    const double* restrict x,
	    const double* restrict y)
{
  for (size_t j = 0; j != size; ++j)
  {
    for (size_t i = 0; i != size; ++i)
    {
      size_t const N = 64;
      size_t const k = (i + size * j);
      g[k] = 0.0;
      for (size_t n = 1; n != (N + 1); n += 2)
      {
	double const pi = M_PI;
	double const ln = ( (double) n ) * pi;
	double const ln3 = (ln * ln * ln);
	double const An = (4.0 / ln3);
	double const Bn = ( (cosh(ln) - 1.0) / sinh(ln) );
	g[k] += An * ( 1.0 - cosh(ln * y[j]) + Bn * sinh(ln * y[j]) ) * sin(ln * x[i]);
      }
    }
  }
}


void init_X (size_t const size, double* restrict X, const double* restrict x)
{
  for (size_t j = 0; j != size; ++j)
  {
    for (size_t i = 0; i != size; ++i)
    {
      size_t k = (i + j * size);
      X[k] = x[i];
    }
  }
}


void init_Y (size_t const size, double* restrict Y, const double* restrict y)
{
  for (size_t i = 0; i != size; ++i)
  {
    for (size_t j = 0; j != size; ++j)
    {
      size_t k = (j + i * size);
      Y[k] = y[i];
    }
  }
}


// void pdesol (double t, workspace_t* workspace)
//
// Synopsis:
// Computes the analytic field array f(t, x, y).
//
// Input:
// t            scalar, the current time
// workspace    data structure containing the current data (field, position, error, etc.)
//
// Output:
// workspace    updates the analytic field array `f'


void pdesol (double const t, workspace_t* workspace)// computes the exact field f(t, x, y)
{
  double* F = workspace -> f;
  double* X = workspace -> rhs;
  double* Y = workspace -> err;
  const double* x = workspace -> x;
  const double* y = workspace -> y;
  size_t const size = workspace -> size;
  size_t const numel = (size * size);
  size_t const N = 64;

  init_X(size, X, x);
  init_Y(size, Y, y);

  for (size_t k = 0; k != numel; ++k)
  {
    F[k] = 0.0;
    for (size_t m = 1; m != (N + 1); m += 2)
    {
      for (size_t n = 1; n != (N + 1); n += 2)
      {
	double const pi = M_PI;
	double const ln = ( (double) n ) * pi;
	double const lm = ( (double) m ) * pi;
	double const lnm = ( (ln * ln) + (lm * lm) );
	double const An = (2.0 / ln);
	double const Am = (2.0 / lm);
	double const Anm = 2.0 * (An * Am);
	double const Bnm = ( (1.0 / lnm) + (1.0 - 1.0 / lnm) * exp(-lnm * t) );
	F[k] += ( 2.0 * Anm * Bnm * sin(ln * X[k]) * sin(lm * Y[k]) );
      }
    }
  }
}


// double RMSE(size_t numel, double* e, double* f, double* g)
//
// Synopsis:
// Computes the Root Mean Squared Error RMSE of the numeric solution.
//
// Inputs:
// numel        number of elements (or array size)
// e            error vector of size `numel' (method ignores element values on entry)
// f            analytic field array of size `numel'
// g            numeric field array of size `numel'
//
// Outputs:
// e            error array, the elementwise difference of `f' and `g' on output
// rmse         scalar, the root mean squared error


double RMSE(size_t const numel,
	    double* restrict e,
	    const double* restrict f,
	    const double* restrict g)
{
  error(numel, e, f, g);
  double const rmse = sqrt( norm(numel, e) ) / ( (double) numel );
  return rmse;
}


// void export (char* fname, size_t size, double* g, double* x)
//
// Synopsis:
// Exports the field g(t, x, y = 0.5) array to a plain text file.
//
// Input:
// fname	filename
// size		number of nodes along the x [y] axis
// g		field array (number of elements: size * size)
// x		position array
//
// Output:
// [x, g]	writes the field array g(t, x, y = 0.5) as a function of position x


void export(const char* fname,
	    size_t const size,
	    const double* restrict g,
	    const double* restrict x)
{
  FILE* file = fopen(fname, "w");
  for (size_t i = 0; i != size; ++i)
  {
    int const j = (size / 2);
    int const k = (i + size * j);
    fprintf(file, "%.12e %.12e\n", x[i], g[k]);
  }

  fclose(file);
}


// void logger (int step, workspace_t* workspace)
//
// Synopsis:
// Logs the Root Mean Squared Error RMSE of the numeric solution.
//
// Inputs:
// step         step number (or id)
//
// Output:
// rmse         logs the RMSE = sqrt( sum( (f - g)**2 ) ) / N on the console, where
//              `f' is the analytic and `g' is the numeric field array, and `N' is the
//              array size.


void logger (int const step, workspace_t* workspace)
{
  double const alpha = ALPHA;
  const double* x = workspace -> x;
  double const dx = (x[1] - x[0]);
  double const dt = (dx * dx) / alpha;
  double const t = ( ( (double) step ) + 1.0 ) * dt;
  pdesol(t, workspace);

  size_t const size = workspace -> size;
  size_t const numel = (size * size);
  const double* f = workspace -> f;
  const double* g = workspace -> g;
  double* err = workspace -> err;
  double const e = RMSE(numel, err, f, g);
  printf("approximation error (transient solution t = %.4e): %e \n", t, e);
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
  int const steps = 0x00010000;
  for (int step = 0; step != steps; ++step)
  {
    solver(workspace);

    size_t state = get_state(workspace);
    if (state == FAILURE_STATE)
    {
      char msg[] = "Jacobi solver failed to converge to the solution after %d iterations";
      printf(msg, MAX_ITERATIONS);
      break;
    }

    int const span = (steps / 16);
    // logs error of exact f(t+dt, x) and numeric solution g(t+dt, x) every `span' steps
    if ( ( step != 0 ) && ( (step % span) == 0 ) )
    {
      logger(step, workspace);
    }
  }
}


// void Poisson()
//
// Synopsis:
// Solves the 2d transient Poisson equation. The initial temperature field is one at the
// interior nodes and zero at the boundaries. There is a uniform heat source throughout
// the domain.


void Poisson ()
{
  // parameters:

  size_t const size = SIZE;			// number of elements of position arrays
  size_t const numel = (size * size);		// size of field arrays g(t, x, y)

  // memory allocations:

  workspace_t* workspace = create(size);

  // iterators:

  double* x = workspace -> x;
  double* y = workspace -> y;
  double* f = workspace -> f;
  double* g = workspace -> g;
  double* g0 = workspace -> g0;
  double* err = workspace -> err;
  double* mask = workspace -> mask;

  // initializations:

  init_field(size, g);				// sets the initial (temperature) field
  zeros(numel, g0);				// inits the initial solution array

  init_mask(size, mask);			// masks boundary nodes

  double const x_l = 0;				// defines x-axis lower bound
  double const x_u = 1;				// defines x-axis upper bound
  linspace(x, x_l, x_u, size);			// inits the x-axis position array

  double const y_l = 0;				// defines y-axis lower bound
  double const y_u = 1;				// defines y-axis upper bound
  linspace(y, y_l, y_u, size);			// inits the y-axis position array

  // solves the transient Poisson equation:

  integrator(workspace);

  // post-processing:

  // logs the analytic and numeric solutions at steady-state
  exact(size, f, x, y);
  export("analytic.txt", size, f, x);
  export("numeric.txt", size, g, x);

  // logs the approximation error
  error(numel, err, f, g);
  double const e = RMSE(numel, err, f, g);
  printf("approximation error (steady-state solution): %e \n", e);

  // memory deallocations:

  workspace = destroy(workspace);
}


// COMMENTS:
// iNODE:       interior node, a nodes that do not lie at the boundary
