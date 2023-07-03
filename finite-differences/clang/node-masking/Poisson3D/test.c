/*
 * Transient Heat Conduction                                	June 03, 2023
 *
 * source: test.c
 * author: @misael-diaz
 *
 * Synopsis:
 * Solves the 3d transient Poisson equation iteratively with the Jacobi method
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


#define SIZE 64
#define ALPHA 2.0
#define iNODE 0xffffffffffffffff
#define TOLERANCE 8.673617379884035e-19
#define MAX_ITERATIONS 128
#define SUCCESS_STATE 0
#define FAILURE_STATE 1
#define VERBOSE true


// we use this for auto-vectorization of loops with simple nested conditionals
typedef union
{
  uint64_t bin;		// bit pattern of the floating point data
  double data;		// floating point data
} alias_t;


typedef struct {
  double* x;		// position array along the x-axis
  double* y;		// position array along the y-axis
  double* z;		// position array along the z-axis
  double* f;		// exact solution array, f(t, x, y, z)
  double* g;		// estimate of the solution array, g(t + dt, x, y, z)
  double* g0;		// previous estimate of the solution array, g(t + dt, x, y, z)
  double* err;		// error array
  double* rhs;		// Right Hand Side RHS array of the discretized PDE
  double* tmp;		// array temporary
  double* mask;		// bitmask is zero at boundary nodes and `iNODE' at interior ones
  double* data;		// contiguous data buffer (it stores x, y, z, f, g, etc.)
  size_t size;		// array size of x, y, and z position arrays
  size_t state;		// solver state
} workspace_t;


void Poisson();

int main ()
{
  Poisson();
  return 0;
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


// implements array addition
void add (size_t const size, double* restrict dst, const double* restrict src)
{
  for (size_t i = 0; i != size; ++i)
  {
    dst[i] += src[i];
  }
}


// implements array multiplication
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


double norm (size_t const size, const double* x)	// gets the norm of the array `x'
{
  double sum = 0;
  for (size_t i = 0; i != size; ++i)
  {
    sum += (x[i] * x[i]);
  }

  return sum;
}


// copies source `src' into destination `dst' array (GCC optimizes it into a library call)
void copy (size_t const size, const double* restrict src, double* restrict dst)
{
  for (size_t i = 0; i != size; ++i)
  {
    dst[i] = src[i];
  }
}


// creates the workspace for the Poisson solver
workspace_t* create (size_t const size)
{
  size_t const numel = (size * size * size);

  size_t const x_size = size;
  size_t const y_size = size;
  size_t const z_size = size;
  size_t const f_size = numel;
  size_t const g_size = numel;
  size_t const g0_size = numel;
  size_t const rhs_size = numel;
  size_t const err_size = numel;
  size_t const tmp_size = numel;
  size_t const mask_size = numel;
  size_t const data_size = x_size +
			   y_size +
			   z_size +
			   f_size +
			   g_size +
			   g0_size +
			   rhs_size +
			   err_size +
			   tmp_size +
			   mask_size;

  double* data = malloc( data_size * sizeof(double) );
  if (data == NULL)
  {
    return NULL;
  }

  workspace_t *workspace = malloc( sizeof(workspace_t) );
  if (workspace == NULL)
  {
    free(data);
    data = NULL;
    return workspace;
  }

  workspace -> x = data;
  workspace -> y = workspace -> x + size;
  workspace -> z = workspace -> y + size;
  workspace -> f = workspace -> z + size;
  workspace -> g = workspace -> f + numel;
  workspace -> g0 = workspace -> g + numel;
  workspace -> err = workspace -> g0 + numel;
  workspace -> rhs = workspace -> err + numel;
  workspace -> tmp = workspace -> rhs + numel;
  workspace -> mask = workspace -> tmp + numel;
  workspace -> data = data;

  double* x = workspace -> x;
  double* y = workspace -> y;
  double* z = workspace -> z;
  double* f = workspace -> f;
  double* g = workspace -> g;
  double* g0 = workspace -> g0;
  double* err = workspace -> err;
  double* rhs = workspace -> rhs;
  double* tmp = workspace -> tmp;
  double* mask = workspace -> mask;

  zeros(size, x);
  zeros(size, y);
  zeros(size, z);
  zeros(numel, f);
  zeros(numel, g);
  zeros(numel, g0);
  zeros(numel, err);
  zeros(numel, rhs);
  zeros(numel, tmp);
  zeros(numel, mask);

  workspace -> size = size;
  workspace -> state = FAILURE_STATE;

  data = NULL;
  return workspace;
}


// frees workspace from memory
workspace_t* destroy (workspace_t* workspace)
{
  if (workspace == NULL)
  {
    return workspace;
  }

  free(workspace -> data);
  workspace -> data = NULL;

  free(workspace);
  workspace = NULL;

  return workspace;
}


// sets the state of the numeric solver
void set_state (workspace_t* workspace, size_t state)
{
  workspace -> state = state;
}


// gets the state of the numeric solver
size_t get_state (workspace_t* workspace)
{
  return (workspace -> state);
}


// void init_field (size_t size, double* g)
//
// Synopsis:
// Sets the initial (temperature) field g(t = 0, x, y, z).
// Fills with ones at the interior nodes and with zeros at the boundary nodes.
//
// Input:
// size		number of nodes along the x, y, or z axis
//
// Output:
// g		3d field array, zeros at the boundaries and ones elsewhere


void init_field (size_t const size, double* g)
{
  size_t const size2 = (size * size);
  size_t const numel = (size * size * size);

  // sets g(t = 0, x, y, z) = 1 (everywhere):
  ones(numel, g);

  // sets g(t = 0, x = 0, y, z) = 0:
  for (size_t k = 0; k != size; ++k)
  {
    for (size_t j = 0; j != size; ++j)
    {
      size_t const i = 0;
      size_t const count = (i + size * j + size2 * k);
      g[count] = 0.0;
    }
  }

  // sets g(t = 0, x = 1, y, z) = 0:
  for (size_t k = 0; k != size; ++k)
  {
    for (size_t j = 0; j != size; ++j)
    {
      size_t const i = (size - 1);
      size_t const count = (i + size * j + size2 * k);
      g[count] = 0.0;
    }
  }

  // sets g(t = 0, x, y = 0, z) = 0:
  for (size_t k = 0; k != size; ++k)
  {
    for (size_t i = 0; i != size; ++i)
    {
      size_t const j = 0;
      size_t const count = (i + size * j + size2 * k);
      g[count] = 0.0;
    }
  }

  // sets g(t = 0, x, y = 1, z) = 0:
  for (size_t k = 0; k != size; ++k)
  {
    for (size_t i = 0; i != size; ++i)
    {
      size_t const j = (size - 1);
      size_t const count = (i + size * j + size2 * k);
      g[count] = 0.0;
    }
  }

  // sets g(t = 0, x, y, z = 0) = 0:
  for (size_t j = 0; j != size; ++j)
  {
    for (size_t i = 0; i != size; ++i)
    {
      size_t const k = 0;
      size_t const count = (i + size * j + size2 * k);
      g[count] = 0.0;
    }
  }

  // sets g(t = 0, x, y, z = 1) = 0:
  for (size_t j = 0; j != size; ++j)
  {
    for (size_t i = 0; i != size; ++i)
    {
      size_t const k = (size - 1);
      size_t const count = (i + size * j + size2 * k);
      g[count] = 0.0;
    }
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
// size		number of nodes along the x, y, or z axis
//
// Output:
// mask		array with `(size * size * size)' number of elements


void init_mask (size_t const size, double* mask)
{
  alias_t* m = mask;
  size_t const size2 = (size * size);
  size_t const numel = (size * size * size);

  // sets all nodes as interior nodes:
  for (size_t i = 0; i != numel; ++i)
  {
    m[i].bin = iNODE;
  }

  // sets nodes at x = 0 as boundary nodes:
  for (size_t k = 0; k != size; ++k)
  {
    for (size_t j = 0; j != size; ++j)
    {
      size_t const i = 0;
      size_t const count = (i + size * j + size2 * k);
      mask[count] = 0.0;
    }
  }

  // sets nodes at x = 1 as boundary nodes:
  for (size_t k = 0; k != size; ++k)
  {
    for (size_t j = 0; j != size; ++j)
    {
      size_t const i = (size - 1);
      size_t const count = (i + size * j + size2 * k);
      mask[count] = 0.0;
    }
  }

  // sets nodes at y = 0 as boundary nodes:
  for (size_t k = 0; k != size; ++k)
  {
    for (size_t i = 0; i != size; ++i)
    {
      size_t const j = 0;
      size_t const count = (i + size * j + size2 * k);
      mask[count] = 0.0;
    }
  }

  // sets nodes at y = 1 as boundary nodes:
  for (size_t k = 0; k != size; ++k)
  {
    for (size_t i = 0; i != size; ++i)
    {
      size_t const j = (size - 1);
      size_t const count = (i + size * j + size2 * k);
      mask[count] = 0.0;
    }
  }

  // sets nodes at z = 0 as boundary nodes:
  for (size_t j = 0; j != size; ++j)
  {
    for (size_t i = 0; i != size; ++i)
    {
      size_t const k = 0;
      size_t const count = (i + size * j + size2 * k);
      mask[count] = 0.0;
    }
  }

  // sets nodes at z = 1 as boundary nodes:
  for (size_t i = 0; i != size; ++i)
  {
    for (size_t j = 0; j != size; ++j)
    {
      size_t const k = (size - 1);
      size_t const count = (i + size * j + size2 * k);
      mask[count] = 0.0;
    }
  }
}


// void init_rhs (size_t size, double* b, double* x, double* g)
//
// Synopsis:
// Initializes the Right-Hand-Size RHS of the Finite Difference Equations FDEs.
//
// Inputs:
// size		number of nodes along the x, y, or z axis
// x		x-axis position array
// g		solution array g(t, x, y, z) at the current time
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
  size_t const numel = (size * size * size);
  // we don't need to mask this loop since we are masking the boundary nodes elsewhere
  for (size_t i = 0; i != numel; ++i)
  {
    b[i] = (dx * dx + alpha * g[i]);
  }
}



// void rhs(size_t size, double* g, double* b)
//
// Synopsis:
// Updates the estimate of the solution array g(t + dt, x, y, z) from the RHS of the FDEs.
//
// Inputs:
// size		number of nodes along the x, y, or z axis
// b		RHS array
// mask		bitmask for protecting boundary node data
//
// Outputs:
// g		estimate of the solution array g(t + dt, x, y, z)


void rhs (size_t const size,
	  double* restrict g,
	  const double* restrict b,
	  const double* restrict mask)
{
  alias_t* dst = g;
  const alias_t* values = b;
  const alias_t* masks = mask;
  size_t const numel = (size * size * size);
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
// size		number of nodes along the x, y, or z axis
// g0		previous estimate of the solution array g(t + dt, x, y, z)
// tmp		array temporary for storing intermediate computations
// mask		bitmask for protecting boundary node data
//
// Outputs:
// g		estimate of the solution array g(t + dt, x, y, z)


void tridiag (size_t const size,
	      double* restrict g,
	      double* restrict tmp,
	      const double* restrict g0,
	      const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* values = g0;
  const alias_t* masks = mask;
  size_t const numel = (size * size * size);

  zeros(numel, tmp);

  for (size_t i = 1; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values[i - 1].bin);
  }

  add(numel, g, tmp);

  for (size_t i = 0; i != (numel - 1); ++i)
  {
    t[i].bin = (masks[i].bin & values[i + 1].bin);
  }

  add(numel, g, tmp);
}



// void subdiag(size_t size, double* g, double* g0, double* mask)
//
// Synopsis:
// Updates the solution array g(t + dt, x, y, z) from the sub-diagonal terms of the FDEs.
//
// Inputs:
// size		number of nodes along the x, y, or z axis
// g0		previous estimate of the solution array g(t + dt)
// tmp		array temporary for storing intermediate computations
// mask		bitmask for protecting boundary node data
//
// Outputs:
// g		estimate of the solution array g(t + dt, x, y, z)


void subdiag (size_t const size,
	      double* restrict g,
	      double* restrict tmp,
	      const double* restrict g0,
	      const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* values = g0;
  const alias_t* masks = mask;
  size_t const numel = (size * size * size);

  zeros(numel, tmp);

  for (size_t i = size; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values[i - size].bin);
  }

  add(numel, g, tmp);
}



// void superdiag(size_t size, double* g, double* g0, double* mask)
//
// Synopsis:
// Updates the solution array g(t + dt, x, y, z) from the super-diagonal terms of the FDEs
//
// Inputs:
// size		number of nodes along the x, y, or z axis
// g0		previous estimate of the solution array g(t + dt, x, y, z)
// tmp		array temporary for storing intermediate computations
// mask		bitmask for protecting boundary node data
//
// Outputs:
// g		estimate of the solution array g(t + dt, x, y, z)


void superdiag (size_t const size,
		double* restrict g,
		double* restrict tmp,
		const double* restrict g0,
		const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* values = g0;
  const alias_t* masks = mask;
  size_t const numel = (size * size * size);

  zeros(numel, tmp);

  for (size_t i = 0; i != (numel - size); ++i)
  {
    t[i].bin = (masks[i].bin & values[i + size].bin);
  }

  add(numel, g, tmp);
}


// as subdiag() but uses the lowest band to update g(t + dt, x, y, z)
void subband (size_t const size,
	      double* restrict g,
	      double* restrict tmp,
	      const double* restrict g0,
	      const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* values = g0;
  const alias_t* masks = mask;
  size_t const numel = (size * size * size);
  size_t const size2 = (size * size);

  zeros(numel, tmp);

  for (size_t i = size2; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values[i - size2].bin);
  }

  add(numel, g, tmp);
}



// as superdiag() but uses the highest band to update g(t + dt, x, y, z)
void superband (size_t const size,
		double* restrict g,
		double* restrict tmp,
		const double* restrict g0,
		const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* values = g0;
  const alias_t* masks = mask;
  size_t const numel = (size * size * size);
  size_t const size2 = (size * size);

  zeros(numel, tmp);

  for (size_t i = 0; i != (numel - size2); ++i)
  {
    t[i].bin = (masks[i].bin & values[i + size2].bin);
  }

  add(numel, g, tmp);
}


// void scale (size_t size, double* g)
//
// Synopsis:
// Scales the solution array g(t + dt, x, y, z) with the main diagonal coefficient.
//
// Inputs:
// size		number of nodes along the x, y, or z axis
// tmp		array temporary for storing intermediate computations
// mask		bitmask for protecting boundary node data
//
// Outputs:
// g		estimate of the solution array g(t + dt, x, y, z)


void __attribute__ ((noinline)) scale(size_t const size,
				      double* restrict g,
				      double* restrict tmp,
				      const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* masks = mask;
  double const alpha = ALPHA;
  double const c = 1.0 / (alpha + 6.0);
  alias_t const values = { .data = c };
  size_t const numel = (size * size * size);

  for (size_t i = 0; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values.bin);
  }

  mult(numel, g, tmp);
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
// Solves the transient 3D Poisson FDEs with the Jacobi method.
// The solver finds the solution g(t + dt, x, y, z) at the next time step from the current
// g(t, x, y, z).
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
  size_t const numel = (size * size * size);

  // iterators:

  double* x = workspace -> x;
//double* y = workspace -> y;
//double* z = workspace -> z;
  double* b = workspace -> rhs;
  double* g = workspace -> g;
  double* g0 = workspace -> g0;
  double* err = workspace -> err;
  double* tmp = workspace -> tmp;
  const double* mask = workspace -> mask;

  // initializations:

  init_rhs(size, b, x, g);			// vectorized by gcc

  // Jacobi solver:

  set_state(workspace, FAILURE_STATE);
  for (size_t i = 0; i != iters; ++i)
  {
    // updates the solution array g(t + dt, x, y, z):
    rhs(size, g, b, mask);			// vectorized by gcc
    tridiag(size, g, tmp, g0, mask);		// vectorized by gcc
    subdiag(size, g, tmp, g0, mask);		// vectorized by gcc
    superdiag(size, g, tmp, g0, mask);		// vectorized by gcc
    subband(size, g, tmp, g0, mask);		// vectorized by gcc
    superband(size, g, tmp, g0, mask);		// vectorized by gcc
    scale(size, g, tmp, mask);			// vectorized by gcc

    // checks for convergence:
    error(numel, err, g, g0);			// vectorized by gcc
    if (norm(numel, err) < tol)
    {
      set_state(workspace, SUCCESS_STATE);
      if (VERBOSE) printf("Jacobi(): solution found after %lu iters\n", i + 1);
      break;
    }

    // updates the initial solution array g(t + dt, x, y, z) for the next iteration:
    copy(numel, g, g0);				// optimized by gcc
  }
}


// convenience function that stores the x position array for all nodes
void init_X (size_t const size, double* restrict X, const double* restrict x)
{
  size_t const size2 = (size * size);
  for (size_t k = 0; k != size; ++k)
  {
    for (size_t j = 0; j != size; ++j)
    {
      for (size_t i = 0; i != size; ++i)
      {
	size_t const count = (i + j * size + k * size2);
	X[count] = x[i];
      }
    }
  }
}


// convenience function that stores the y position array for all nodes
void init_Y (size_t const size, double* restrict Y, const double* restrict y)
{
  size_t const size2 = (size * size);
  for (size_t k = 0; k != size; ++k)
  {
    for (size_t j = 0; j != size; ++j)
    {
      for (size_t i = 0; i != size; ++i)
      {
	size_t const count = (i + j * size + k * size2);
	Y[count] = y[j];
      }
    }
  }
}


// convenience function that stores the z position array for all nodes
void init_Z (size_t const size, double* restrict Z, const double* restrict z)
{
  size_t const size2 = (size * size);
  for (size_t k = 0; k != size; ++k)
  {
    for (size_t j = 0; j != size; ++j)
    {
      for (size_t i = 0; i != size; ++i)
      {
	size_t const count = (i + j * size + k * size2);
	Z[count] = z[k];
      }
    }
  }
}


// void pdesol (double t, workspace_t* workspace)
//
// Synopsis:
// Computes the transient analytic field array f(t, x, y, z).
// The analytic field was obtained by applying the Finite Fourier Transform FFT method.
// Uses the (above) convenience functions so that with one index we can address the
// position and the field arrays.
//
// NOTE:
// This is the bottleneck of the numerical application.
//
// Input:
// t            scalar, the current time
// workspace    data structure containing the current data (field, position, error, etc.)
//
// Output:
// workspace    updates the analytic field array `f'


void pdesol (double const t, workspace_t* workspace)
{
  double* F = workspace -> f;
  double* X = workspace -> rhs;
  double* Y = workspace -> err;
  double* Z = workspace -> tmp;
  const double* x = workspace -> x;
  const double* y = workspace -> y;
  const double* z = workspace -> z;
  size_t const size = workspace -> size;
  size_t const numel = (size * size * size);
  // NOTE: higher order N incurs in longer application runtimes
  size_t const N = 16;

  init_X(size, X, x);
  init_Y(size, Y, y);
  init_Z(size, Z, z);

  for (size_t n = 0; n != numel; ++n)
  {
    F[n] = 0.0;
    for (size_t k = 1; k != (N + 1); k += 2)
    {
      for (size_t j = 1; j != (N + 1); j += 2)
      {
	for (size_t i = 1; i != (N + 1); i += 2)
	{
	  double const Pi = M_PI;
	  double const li = ( (double) i ) * Pi;
	  double const lj = ( (double) j ) * Pi;
	  double const lk = ( (double) k ) * Pi;
	  double const lijk = ( (li * li) + (lj * lj) + (lk * lk) );
	  double const Ai = (2.0 / li);
	  double const Aj = (2.0 / lj);
	  double const Ak = (2.0 / lk);
	  double const Aijk = (Ai * Aj * Ak);
	  double const Bijk = ( (1.0 / lijk) + (1.0 - 1.0 / lijk) * exp(-lijk * t) );
	  double const A = Aijk;
	  double const B = Bijk;
	  F[n] += ( 8.0 * A * B * sin(li * X[n]) * sin(lj * Y[n]) * sin(lk * Z[n]) );
	}
      }
    }
  }
}


// void pdesol_steady_state (workspace_t* workspace)
//
// Synopsis:
// Obtains the analytic, steady-state, solution of the Poisson equation.
// This is just like pdesol() but without the transient (exponential decay) term.
//
// Inputs:
// workspace	work data structure
//
// Outputs:
// workspace	updated with the steady-state field array


void pdesol_steady_state (workspace_t* workspace)
{
  double* F = workspace -> f;
  double* X = workspace -> rhs;
  double* Y = workspace -> err;
  double* Z = workspace -> tmp;
  const double* x = workspace -> x;
  const double* y = workspace -> y;
  const double* z = workspace -> z;
  size_t const size = workspace -> size;
  size_t const numel = (size * size * size);
  size_t const N = 16;

  init_X(size, X, x);
  init_Y(size, Y, y);
  init_Z(size, Z, z);

  for (size_t n = 0; n != numel; ++n)
  {
    F[n] = 0.0;
    for (size_t k = 1; k != (N + 1); k += 2)
    {
      for (size_t j = 1; j != (N + 1); j += 2)
      {
	for (size_t i = 1; i != (N + 1); i += 2)
	{
	  double const Pi = M_PI;
	  double const li = ( (double) i ) * Pi;
	  double const lj = ( (double) j ) * Pi;
	  double const lk = ( (double) k ) * Pi;
	  double const lijk = ( (li * li) + (lj * lj) + (lk * lk) );
	  double const Ai = (2.0 / li);
	  double const Aj = (2.0 / lj);
	  double const Ak = (2.0 / lk);
	  double const Aijk = (Ai * Aj * Ak);
	  double const Bijk = (1.0 / lijk);
	  double const A = Aijk;
	  double const B = Bijk;
	  F[n] += ( 8.0 * A * B * sin(li * X[n]) * sin(lj * Y[n]) * sin(lk * Z[n]) );
	}
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
// Exports the field g(t, x, y = 0.5, z = 0.5) array to a plain text file.
//
// Input:
// fname	filename
// size		number of nodes along the x, y, or z axis
// g		field array (number of elements: size * size * size)
// x		position array
//
// Output:
// [x, g]	writes the field g(t, x, y = 0.5, z = 0.5) with respect to position x


void export(const char* fname,
	    size_t const size,
	    const double* restrict g,
	    const double* restrict x)
{
  FILE* file = fopen(fname, "w");
  for (size_t i = 0; i != size; ++i)
  {
    size_t const size2 = (size * size);
    size_t const k = (size / 2);
    size_t const j = (size / 2);
    size_t const count = (i + size * j + size2 * k);
    fprintf(file, "%.12e %.12e\n", x[i], g[count]);
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
// count	the number of times the logger has been invoked
//
// Output:
// rmse         logs the RMSE = sqrt( sum( (f - g)**2 ) ) / N on the console, where
//              `f' is the analytic and `g' is the numeric field array, and `N' is the
//              array size.


void logger (int const step, int const count, workspace_t* workspace)
{
  double const alpha = ALPHA;
  const double* x = workspace -> x;
  double const dx = (x[1] - x[0]);
  double const dt = (dx * dx) / alpha;
  double const t = ( ( (double) step ) + 1.0 ) * dt;
  pdesol(t, workspace);

  size_t const size = workspace -> size;
  size_t const numel = (size * size * size);
  const double* f = workspace -> f;
  const double* g = workspace -> g;
  double* err = workspace -> err;
  double const e = RMSE(numel, err, f, g);
  printf("approximation error (transient solution t = %.4e): %e \n", t, e);

  char filename[80];
  sprintf(filename, "analytic%d.txt", count);
  export(filename, size, f, x);

  sprintf(filename, "numeric%d.txt", count);
  export(filename, size, g, x);
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
  int count = 0;
  int const steps = 0x00000400;
  for (int step = 0; step != steps; ++step)
  {
    solver(workspace);

    size_t const state = get_state(workspace);
    if (state == FAILURE_STATE)
    {
      char const msg[] = "Jacobi solver failed to converge after %d iterations\n";
      printf(msg, MAX_ITERATIONS);
      break;
    }

    int const span = (steps / 4);
    // logs error of exact f(t+dt, x) and numeric solution g(t+dt, x) every `span' steps
    if ( ( step != 0 ) && ( (step % span) == 0 ) )
    {
      logger(step, count, workspace);
      ++count;
    }
  }
}


// void Poisson()
//
// Synopsis:
// Solves the 3d transient Poisson equation. The initial temperature field is one at the
// interior nodes and zero at the boundaries. There is a uniform heat source throughout
// the domain.


void Poisson ()
{
  // parameters:

  size_t const size = SIZE;			// number of elements of position arrays
  size_t const numel = (size * size * size);	// size of field arrays g(t, x, y, z)

  // memory allocations:

  workspace_t* workspace = create(size);

  // iterators:

  double* x = workspace -> x;
  double* y = workspace -> y;
  double* z = workspace -> z;
  double* f = workspace -> f;
  double* g = workspace -> g;
  double* g0 = workspace -> g0;
//double* rhs = workspace -> rhs;
  double* err = workspace -> err;
//double* tmp = workspace -> tmp;
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

  double const z_l = 0;				// defines z-axis lower bound
  double const z_u = 1;				// defines z-axis upper bound
  linspace(z, z_l, z_u, size);			// inits the z-axis position array

  // checks number of boundary nodes (from faces, edges, and corners, respectively):
  double sum = 0;
  for (size_t i = 0; i != numel; ++i)
  {
    sum += g[i];
  }

  double const bnode_count = (numel - sum);
  double const bnode_count_expected = 6 * (size - 2) * (size - 2) + 12 * (size - 2) + 8;
  printf("test[0]: ");
  if (bnode_count != bnode_count_expected)
  {
    printf("FAIL\n");
  }
  else
  {
    printf("PASS\n");
  }

  sum = 0;
  const alias_t* m = mask;
  for (size_t i = 0; i != numel; ++i)
  {
    sum += (m[i].bin & 1);
  }

  // checks number of interior nodes
  double const inode_count = sum;
  double const inode_count_expected = (numel - bnode_count_expected);
  printf("test[1]: ");
  if (inode_count != inode_count_expected)
  {
    printf("FAIL\n");
  }
  else
  {
    printf("PASS\n");
  }

  // solves the transient Poisson equation:

  integrator(workspace);

  // post-processing:

  // logs the analytic and numeric solutions at steady-state
  pdesol_steady_state(workspace);
  export("analytic.txt", size, f, x);
  export("numeric.txt", size, g, x);

  // logs the approximation error of the steady-state solutions
  error(numel, err, f, g);
  double const e = RMSE(numel, err, f, g);
  printf("approximation error (steady-state solution): %e \n", e);

  // memory deallocations:

  workspace = destroy(workspace);
}


// COMMENTS:
// iNODE:       interior node, a nodes that do not lie at the boundary
