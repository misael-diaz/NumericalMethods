/*
 * Transient Heat Conduction                                	June 04, 2023
 *
 * source: test.c
 * author: @misael-diaz
 *
 * Synopsis:
 * Solves the 2d transient Poisson equation iteratively.
 * The solver executes in a multi-threaded environment via POSIX threads.
 * Auto-vectorized code by GCC is indicated in the source.
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
#include <pthread.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>


#define SIZE 128
#define ALPHA 2.0
#define iNODE 0xffffffffffffffff
#define TOLERANCE 8.673617379884035e-19
#define MAX_ITERATIONS 128
#define SUCCESS_STATE 0
#define FAILURE_STATE 1
#define VERBOSE false
#define LOG true
#define MPI_NUM_THREADS 4
#define MPI_THREAD_ROOT 0


typedef union
{
  uint64_t bin;		// binary pattern of the double precision floating-point data
  double data;		// double precision floating-point data
} alias_t;


/*
typedef struct
{
  double* x;		// position array along the x-axis
  double* y;		// position array along the y-axis
  double* f;		// exact solution array, f(t, x, y)
  double* g;		// estimate of the solution array at next step, g(t + dt, x, y)
  double* g0;		// previous estimate of the solution array, g(t + dt, x, y)
  double* err;		// error array
  double* rhs;		// Right Hand Side RHS array of the PDE
  double* mask;		// bitmask is zero at boundary nodes and `iNODE' at interior ones
  size_t size;		// array size
  size_t state;		// solver state
} workspace_t;


typedef struct
{
  size_t tid;			// thread ID
  size_t beg;			// thread begin slice
  size_t end;			// thread end slice
  workspace_t* workspace;	// placeholder for the shared data
  pthread_barrier_t* barrier;	// POSIX thread sync barrier
} sharedspace_t;


*/

typedef struct
{
  double* x;			// position array along the x-axis
  double* y;			// position array along the y-axis
  double* f;			// exact solution array, f(t, x, y)
  double* g;			// next estimate of the solution array, g(t + dt, x, y)
  double* g0;			// estimate of the solution array, g(t + dt, x, y)
  double* err;			// error array
  double* rhs;			// Right Hand Side RHS array of the PDE
  double* mask;			// bitmask for protecting boundary node data
  MPI_Status* mpi_status;	// MPI Communication status
  size_t size;			// position array size
  size_t xsize;			// extended array size
  size_t state;			// solver state
  size_t mpi_thread_rank;	// MPI Thread ID
  size_t mpi_thread_pool_size;	// number of MPI threads
  size_t mpi_thread_load_size;	// workload size of MPI thread (same for each)
} mpi_workspace_t;

void Poisson();
void test_mpi_comms();

int main ()
{
  //test_mpi_comms();
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


void linspace (double* x, double x_i, double x_f, size_t const size)// numpy-like linspace
{
  double const N = (size - 1);
  double const dx = (x_f - x_i) / N;
  for (size_t i = 0; i != size; ++i)
  {
    x[i] = x_i + ( (double) i ) * dx;
  }
}


double norm (size_t const size, const double* x)        // gets the norm of the array `x'
{
  double sum = 0;
  for (size_t i = 0; i != size; ++i)
  {
    sum += (x[i] * x[i]);
  }

  return sum;
}


void norm1 (size_t const rank,
	    size_t const numel,
	    double* restrict tmp,
	    const double* restrict err)
{
  double sum = 0;
  for (size_t i = 0; i != numel; ++i)
  {
    sum += (err[i] * err[i]);
  }

  // we could simply store at the first element but it makes no difference really
  tmp[rank] = sum;
}


double norm2 (const double* err)	// NOTE: executed by the root thread
{
  double sum = 0;
  for (size_t tid = 0; tid != MPI_NUM_THREADS; ++tid)
  {
    sum += err[tid];
  }

  return sum;
}


void add (size_t const size,
	  double* restrict dst,
	  const double* restrict src)
{
  for (size_t i = 0; i != size; ++i)
  {
    dst[i] += src[i];
  }
}


void mult(size_t const numel,
	  double* restrict dst,
	  const double* restrict src)
{
  for (size_t i = 0; i != numel; ++i)
  {
    dst[i] *= src[i];
  }
}


void copy(size_t const size,
	  const double* restrict src,
	  double* restrict dst)
{
  for (size_t i = 0; i != size; ++i)
  {
    dst[i] = src[i];
  }
}


mpi_workspace_t* create(size_t const size, const int* mpi_thread_info)
{
  size_t const numel = (size * size);
  if (numel % MPI_NUM_THREADS)
  {
    printf("the number of nodes must be a multiple of the number of threads\n");
    return NULL;
  }

  size_t const load = (numel / MPI_NUM_THREADS);
  if (load <= size)
  {
    printf("too many threads to distribute the load among threads uniformly\n");
    return NULL;
  }

  mpi_workspace_t *workspace = malloc( sizeof(mpi_workspace_t) );
  if (workspace == NULL)
  {
    printf("failed to allocate the shared workspace\n");
    return NULL;
  }

  workspace -> mpi_status = (MPI_Status*) malloc( sizeof(MPI_Status) );
  if (workspace -> mpi_status == NULL)
  {
    printf("failed to allocate MPI status placeholder!\n");
    free(workspace);
    workspace = NULL;
    return NULL;
  }

  size_t const xsize = (2 * size + load);
  // assumes successful allocations for simplicity
  workspace -> x = alloc(size);
  workspace -> y = alloc(size);
  workspace -> f = alloc(load);
  workspace -> g = alloc(load);
  // we need extra space to store data from other (contiguous) threads
  workspace -> g0 = alloc(xsize);
  workspace -> err = alloc(load);
  workspace -> rhs = alloc(load);
  workspace -> mask = alloc(load);

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
  zeros(load, f);
  zeros(load, g);
  zeros(xsize, g0);
  zeros(load, err);
  zeros(load, rhs);
  zeros(load, mask);

  workspace -> size = size;
  workspace -> xsize = xsize;
  workspace -> state = FAILURE_STATE;
  const int* mpi_thread_rank = mpi_thread_info;
  const int* mpi_thread_pool_size = (mpi_thread_info + 1);
  workspace -> mpi_thread_rank = ( (size_t) *mpi_thread_rank );
  workspace -> mpi_thread_pool_size = ( (size_t) *mpi_thread_pool_size );
  workspace -> mpi_thread_load_size = load;

  return workspace;
}


mpi_workspace_t* destroy (mpi_workspace_t* workspace)
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

  free(workspace -> mpi_status);
  workspace -> mpi_status = NULL;
  free(workspace);
  workspace = NULL;
  return workspace;
}


void set_state (mpi_workspace_t* workspace, size_t state)
{
  workspace -> state = state;
}


size_t get_state (mpi_workspace_t* workspace)
{
  return (workspace -> state);
}


/*
// void init_field (size_t size, double* g)
//
// Synopsis:
// Serial Method.
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

  for (size_t i = 0; i != size; ++i)
  {
    size_t const j = 0;
    size_t const k = (i + size * j);
    g[k] = 0.0;					// sets g(t = 0, x, y = 0) = 0
  }

  for (size_t i = 0; i != size; ++i)
  {
    size_t const j = (size - 1);
    size_t const k = (i + size * j);
    g[k] = 0.0;					// sets g(t = 0, x, y = 1) = 0
  }

  for (size_t j = 0; j != size; ++j)
  {
    size_t const i = 0;
    size_t const k = (i + size * j);
    g[k] = 0.0;					// sets g(t = 0, x = 0, y) = 0
  }

  for (size_t j = 0; j != size; ++j)
  {
    size_t const i = (size - 1);
    size_t const k = (i + size * j);
    g[k] = 0.0;					// sets g(t = 0, x = 1, y) = 0
  }
}


// void init_mask(size_t size, double* mask)
//
// Synopsis:
// Serial Method.
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
  alias_t* masks = mask;
  size_t const numel = (size * size);
  for (size_t i = 0; i != numel; ++i)
  {
    masks[i].bin = iNODE;			// sets all nodes as interior nodes
  }

  for (size_t i = 0; i != size; ++i)
  {
    size_t const j = 0;
    size_t const k = (i + size * j);
    mask[k] = 0.0;				// sets nodes at y = 0 as boundary nodes
  }

  for (size_t i = 0; i != size; ++i)
  {
    size_t const j = (size - 1);
    size_t const k = (i + size * j);
    mask[k] = 0.0;				// sets nodes at y = 1 as boundary nodes
  }

  for (size_t j = 0; j != size; ++j)
  {
    size_t const i = 0;
    size_t const k = (i + size * j);
    mask[k] = 0.0;				// sets nodes at x = 0 as boundary nodes
  }

  for (size_t j = 0; j != size; ++j)
  {
    size_t const i = (size - 1);
    size_t const k = (i + size * j);
    mask[k] = 0.0;				// sets nodes at x = 1 as boundary nodes
  }
}


// void init_rhs (size_t beg, size_t end, double* b, double* x, double* g)
//
// Synopsis:
// thread-safe method.
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
  // we don't need to mask this loop since we are masking the boundary nodes elsewhere
  for (size_t i = 0; i != size; ++i)
  {
    b[i] = (dx * dx + alpha * g[i]);
  }
}


// void rhs(size_t beg, size_t end, double* g, double* b, double* mask)
//
// Synopsis:
// thread-safe
// Updates the solution array g(t + dt) from the RHS of the FDEs.
//
// Inputs:
// beg		beginning of the data slice assigned to thread
// end		ending of the data slice assigned to thread (exclusive)
// b		RHS array
// mask         bitmask for protecting boundary node data
//
// Outputs:
// g		estimate of the solution array g(t + dt)


void rhs (size_t const beg,
	  size_t const end,
	  double* restrict g,
	  const double* restrict b,
	  const double* restrict mask)
{
  alias_t* dst = g;
  const alias_t* values = b;
  const alias_t* masks = mask;
  for (size_t i = beg; i != end; ++i)
  {
    dst[i].bin = (masks[i].bin & values[i].bin);
  }
}


// void tridiag(size_t beg, size_t end, size_t size, double* g, double* tmp, double* g0,
//              double* mask)
//
// Synopsis:
// Updates the solution array g(t + dt) from the tridiagonal terms of the FDEs.
//
// Inputs:
// beg		beginning of the data slice assigned to thread
// end		ending of the data slice assigned to thread (exclusive)
// size		number of nodes along the x [y] axis
// g0		previous estimate of the solution array g(t + dt)
// tmp          array temporary for storing intermediate computations
// mask         bitmask for protecting boundary node data
//
// Outputs:
// g		estimate of the solution array g(t + dt)


void tridiag (size_t const beg,
	      size_t const end,
	      size_t const size,
	      double* restrict g,
	      double* restrict tmp,
	      const double* restrict g0,
	      const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* values = g0;
  const alias_t* masks = mask;

  zeros(beg, end, tmp);

  if (beg == 0)
  {
    for (size_t i = 1; i != end; ++i)
    {
      t[i].bin = (masks[i].bin & values[i - 1].bin);
    }
  }
  else
  {
    for (size_t i = beg; i != end; ++i)
    {
      t[i].bin = (masks[i].bin & values[i - 1].bin);
    }
  }

  add(beg, end, g, tmp);

  zeros(beg, end, tmp);

  size_t const numel = (size * size);
  if (end == numel)
  {
    for (size_t i = beg; i != (end - 1); ++i)
    {
      t[i].bin = (masks[i].bin & values[i + 1].bin);
    }
  }
  else
  {
    for (size_t i = beg; i != end; ++i)
    {
      t[i].bin = (masks[i].bin & values[i + 1].bin);
    }
  }

  add(beg, end, g, tmp);
}


// void subdiag(size_t beg, size_t end, size_t size, double* g, double* tmp, double* g0,
//              double* mask)
//
// Synopsis:
// Updates the solution array g(t + dt) from the sub-diagonal terms of the FDEs.
//
// Inputs:
// beg		beginning of the data slice assigned to thread
// end		ending of the data slice assigned to thread (exclusive)
// size		number of nodes along the x [y] axis
// g0		previous estimate of the solution array g(t + dt)
// tmp          array temporary for storing intermediate computations
// mask         bitmask for protecting boundary node data
//
// Outputs:
// g		estimate of the solution array g(t + dt)


void subdiag (size_t const beg,
	      size_t const end,
	      size_t const size,
	      double* restrict g,
	      double* restrict tmp,
	      const double* restrict g0,
	      const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* values = g0;
  const alias_t* masks = mask;

  zeros(beg, end, tmp);

  // prevents invalid reads:
  if (beg < size)
  {
    for (size_t i = size; i != end; ++i)
    {
      t[i].bin = (masks[i].bin & values[i - size].bin);
    }
  }
  else
  {
    for (size_t i = beg; i != end; ++i)
    {
      t[i].bin = (masks[i].bin & values[i - size].bin);
    }
  }

  add(beg, end, g, tmp);
}


// void superdiag(size_t beg, size_t end, size_t size, double* g, double* tmp, double* g0,
//                double* mask)
//
// Synopsis:
// thread-safe
// Updates the solution array g(t + dt) from the super-diagonal terms of the FDEs.
//
// Inputs:
// beg		beginning of the data slice assigned to thread
// end		ending of the data slice assigned to thread (exclusive)
// size		number of nodes along the x [y] axis
// g0		previous estimate of the solution array g(t + dt)
// tmp          array temporary for storing intermediate computations
// mask         bitmask for protecting boundary node data
//
// Outputs:
// g		estimate of the solution array g(t + dt)


void superdiag (size_t const beg,
		size_t const end,
		size_t const size,
		double* restrict g,
		double* restrict tmp,
		const double* restrict g0,
		const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* values = g0;
  const alias_t* masks = mask;

  zeros(beg, end, tmp);

  // prevents invalid reads:
  size_t const numel = (size * size);
  if ( (end + size) > numel )
  {
    for (size_t i = beg; i != (end - size); ++i)
    {
      t[i].bin = (masks[i].bin & values[i + size].bin);
    }
  }
  else
  {
    for (size_t i = beg; i != end; ++i)
    {
      t[i].bin = (masks[i].bin & values[i + size].bin);
    }
  }

  add(beg, end, g, tmp);
}


// void scale (size_t beg, size_t end, double* g, double* tmp, double* mask)
//
// Synopsis:
// Scales the solution array g(t + dt) with the main diagonal coefficient.
//
// Inputs:
// beg		beginning of the data slice assigned to thread
// end		ending of the data slice assigned to thread (exclusive)
// tmp          array temporary for storing intermediate computations
// mask         bitmask for protecting boundary node data
//
// Outputs:
// g		estimate of the solution array g(t + dt)


void __attribute__ ((noinline)) scale(size_t const beg,
				      size_t const end,
				      double* restrict g,
				      double* restrict tmp,
				      const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* masks = mask;
  double const alpha = ALPHA;
  double const c = 1.0 / (alpha + 4.0);
  alias_t const values = { .data = c };

  for (size_t i = beg; i != end; ++i)
  {
    t[i].bin = (masks[i].bin & values.bin);
  }

  mult(beg, end, g, tmp);
}


// void error(size_t beg, size_t end, double* e, double* v, double* w)
//
// Synopsis:
// thread-safe
// Computes the elementwise differences of the arrays `v' and `w'.
//
// Inputs:
// beg		beginning of the data slice assigned to thread
// end		ending of the data slice assigned to thread (exclusive)
// v		array
// w		array
//
// Outputs:
// e		error array, stores the elementwise differences


void error (size_t const beg,
	    size_t const end,
	    double* restrict e,
	    const double* restrict v,
	    const double* restrict w)
{
  for (size_t i = beg; i != end; ++i)
  {
    e[i] = (v[i] - w[i]);
  }
}


// void* isolver(void* v_sharedspace)
//
// Synopsis:
// Iterative Solver.
// Solves the transient 2D Poisson FDEs with the Jacobi method.
// The solver finds the solution g(t + dt) at the next time step from the current g(t).
// Multithreaded Execution.
//
// Inputs:
// sharedspace	data structure containing the current data (field, position, error, etc.)
//
// Outputs:
// sharedspace	updates the fields g0 and g with the data of the next time step


void* isolver (void* v_sharedspace)
{
  // unpacking:

  sharedspace_t* shared = v_sharedspace;
  workspace_t* workspace = shared -> workspace;
  pthread_barrier_t* barrier = shared -> barrier;
  size_t const tid = shared -> tid;
  size_t const beg = shared -> beg;
  size_t const end = shared -> end;

  // parameters:

  double const tol = TOLERANCE;
  size_t const iters = MAX_ITERATIONS;
  size_t const size = workspace -> size;

  // iterators:

  double* x = workspace -> x;
  double* b = workspace -> rhs;
  double* g = workspace -> g;
  double* g0 = workspace -> g0;
  double* tmp = workspace -> f;
  double* err = workspace -> err;
  const double* mask = workspace -> mask;

  // initializations:

  pthread_barrier_wait(barrier);

  init_rhs(beg, end, b, x, g);			// vectorized by gcc

  // Jacobi solver:

  if (tid == 0)
  {
    set_state(workspace, FAILURE_STATE);
  }

  pthread_barrier_wait(barrier);

  // implements iterative solver:

  for (size_t i = 0; i != iters; ++i)
  {
    // updates the solution array g(t + dt):
    rhs(beg, end, g, b, mask);			// vectorized by gcc
    tridiag(beg, end, size, g, tmp, g0, mask);	// vectorized by gcc
    subdiag(beg, end, size, g, tmp, g0, mask);	// vectorized by gcc
    superdiag(beg, end, size, g, tmp, g0, mask);// vectorized by gcc
    scale(beg, end, g, tmp, mask);		// vectorized by gcc

    // computes the error array:
    pthread_barrier_wait(barrier);

    error(beg, end, err, g, g0);		// vectorized by gcc
    norm1(tid, beg, end, tmp, err);

    pthread_barrier_wait(barrier);

    // checks for convergence:
    if (norm2(beg, end, tmp) < tol)
    {
      if (tid == 0)
      {
	set_state(workspace, SUCCESS_STATE);
	if (VERBOSE)
	{
	  printf("Jacobi(): solution found after %lu iters\n", i + 1);
	}
      }
      break;
    }

    // updates the initial solution array g(t + dt) for the next iteration:

    copy(beg, end, g, g0);				// optimized by gcc

    pthread_barrier_wait(barrier);
  }

  return NULL;
}


// void exact(size_t size, double* g, double* x)
//
// Synopsis:
// serial
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


void init_X(size_t const beg,
	    size_t const end,
	    size_t const size,
	    double* restrict X,
	    const double* restrict x)
{
  for (size_t k = beg; k != end; ++k)
  {
    size_t i = (k % size);
    X[k] = x[i];
  }
}


void init_Y(size_t const beg,
	    size_t const end,
	    size_t const size,
	    double* restrict Y,
	    const double* restrict y)
{
  for (size_t k = beg; k != end; ++k)
  {
    size_t i = (k / size);
    Y[k] = y[i];
  }
}


// void pdesol (size_t beg, size_t end, double t, workspace_t* workspace)
//
// Synopsis:
// serial
// Computes the analytic field array f(t, x, y).
//
// Input:
// t            scalar, the current time
// workspace    data structure containing the current data (field, position, error, etc.)
//
// Output:
// workspace    updates the analytic field array `f'


void pdesol (size_t const beg, size_t const end, double const t, workspace_t* workspace)
{
  double* F = workspace -> f;
  double* X = workspace -> rhs;
  double* Y = workspace -> err;
  const double* x = workspace -> x;
  const double* y = workspace -> y;
  size_t const size = workspace -> size;
  size_t const numel = (size * size);
  size_t const N = 64;

  init_X(beg, end, size, X, x);
  init_Y(beg, end, size, Y, y);

  for (size_t k = beg; k != end; ++k)
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
// serial
// Computes the Root Mean Squared Error RMSE of the numeric solution.
//
// Inputs:
// numel        number of elements (or array size)
// e            error array of size `numel' (method ignores element values on entry)
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
  error(0, numel, e, f, g);
  double const rmse = sqrt( norm(numel, e) ) / ( (double) numel );
  return rmse;
}



// void export (char* fname, size_t size, double* g, double* x)
//
// Synopsis:
// serial
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
    size_t const j = (size / 2);
    size_t const k = (i + size * j);
    fprintf(file, "%.12e %.12e\n", x[i], g[k]);
  }

  fclose(file);
}


// void logger (size_t tid, size_t beg, size_t end, int step, workspace_t* workspace)
//
// Synopsis:
// thread-safe (only the master thread completes the task)
// Logs the Root Mean Squared Error RMSE of the numeric solution.
//
// Inputs:
// step         step number (or id)
//
// Output:
// rmse         logs the RMSE = sqrt( sum( (f - g)**2 ) ) / N on the console, where
//              `f' is the analytic and `g' is the numeric field array, and `N' is the
//              array size.


void logger(size_t const step, size_t const count, sharedspace_t* space)
{
  size_t const tid = space -> tid;
  size_t const beg = space -> beg;
  size_t const end = space -> end;
  workspace_t* workspace = space -> workspace;
  pthread_barrier_t* barrier = space -> barrier;

  double const alpha = ALPHA;
  const double* x = workspace -> x;
  double const dx = (x[1] - x[0]);
  double const dt = (dx * dx) / alpha;
  double const t = ( ( (double) step ) + 1.0 ) * dt;

  pthread_barrier_wait(barrier);

  pdesol(beg, end, t, workspace);

  pthread_barrier_wait(barrier);

  size_t const size = workspace -> size;
  size_t const numel = (size * size);
  const double* f = workspace -> f;
  const double* g = workspace -> g;
  double* tmp = workspace -> rhs;
  double* err = workspace -> err;

  error(beg, end, err, f, g);
  norm1(tid, beg, end, tmp, err);

  pthread_barrier_wait(barrier);

  double const e = sqrt( norm2(beg, end, tmp) ) / ( (double) numel );
  if (tid == 0)
  {
    printf("approximation error (transient solution t = %.4e): %e \n", t, e);

    char fname[80];
    sprintf(fname, "analytic%lu.txt", count);
    export(fname, size, f, x);

    sprintf(fname, "numeric%lu.txt", count);
    export(fname, size, g, x);
  }
}

// void* integrator(void* sharespace)
//
// Synopsis:
// Multithreaded.
// Integrates the transient Poisson FDEs via time implicit scheme.
//
// Inputs:
// workspace	data structure containing the initial data (field, position, error, etc.)
//
// Outputs:
// workspace	stores the data after completing the integration


void* integrator (void* v_space)
{
  sharedspace_t* space = v_space;
  size_t const tid = space -> tid;
  workspace_t* workspace = space -> workspace;
  pthread_barrier_t* barrier = space -> barrier;
  size_t count = 0;
  size_t const steps = 0x00001000;
  for (size_t step = 0; step != steps; ++step)
  {

    pthread_barrier_wait(barrier);

    isolver(space);

    pthread_barrier_wait(barrier);

    size_t state = get_state(workspace);
    if (state == FAILURE_STATE)
    {
      if (tid == 0)
      {
	char const msg [] = "Jacobi solver failed to converge to the solution "
			    "after %d iterations";
	printf(msg, MAX_ITERATIONS);
      }
      break;
    }

    if (LOG)
    {
      size_t const span = (steps / 4);
      // logs error of exact f(t+dt, x) and numeric solution g(t+dt, x) every `span' steps
      if ( ( step != 0 ) && ( (step % span) == 0 ) )
      {
	logger(step, count, space);
	++count;
      }
    }
  }

  return NULL;
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

  sharedspace_t* spaces = create(size);		// allocates and distributes thread load
  sharedspace_t* space = spaces;
  workspace_t* workspace = space -> workspace;

  // iterators:

  double* x = workspace -> x;
  double* y = workspace -> y;
  double* f = workspace -> f;
  double* g = workspace -> g;
  double* g0 = workspace -> g0;
//double* err = workspace -> err;
  double* mask = workspace -> mask;

  // initializations:

  init_field(size, g);				// sets the initial (temperature) field
  zeros(0, numel, g0);				// inits the initial solution array

  init_mask(size, mask);			// masks boundary nodes

  double const x_l = 0;				// defines x-axis lower bound
  double const x_u = 1;				// defines x-axis upper bound
  linspace(x, x_l, x_u, size);			// inits the x-axis position array

  double const y_l = 0;				// defines y-axis lower bound
  double const y_u = 1;				// defines y-axis upper bound
  linspace(y, y_l, y_u, size);			// inits the y-axis position array

  // creates parallel environment:

  pthread_barrier_t* barrier = space -> barrier;
  pthread_barrier_init(barrier, NULL, NUM_THREADS);

  pthread_t threads[NUM_THREADS];
  pthread_t* thread = threads;
  for (size_t tid = 0; tid != NUM_THREADS; ++tid)
  {
    void* v_space = space;
    pthread_create(thread, NULL, integrator, v_space);	// solves Poisson FDEs in parallel
    ++space;
    ++thread;
  }

  // waits for each thread to complete its task
  for (size_t tid = 0; tid != NUM_THREADS; ++tid)
  {
    pthread_join(threads[tid], 0);
  }

  // post-processing:

  // logs the analytic solution at steady-state
  exact(size, f, x, y);
  export("analytic.txt", size, f, x);

  // deallocates resources:

  pthread_barrier_destroy(barrier);
  spaces = destroy(spaces);
}
*/


void init_field (mpi_workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const load = workspace -> mpi_thread_load_size;
  size_t const rows = (load / size);

  double* g = workspace -> g;

  ones(load, g);				// sets g(t = 0, x, y) = 1 (everywhere)

  for (size_t j = 0; j != rows; ++j)
  {
    size_t const i = 0;
    size_t const k = (i + j * size);
    g[k] = 0.0;					// sets g(t = 0, x, y = 0) = 0
  }

  for (size_t j = 0; j != rows; ++j)
  {
    size_t const i = (size - 1);
    size_t const k = (i + j * size);
    g[k] = 0.0;					// sets g(t = 0, x, y = 1) = 0
  }

  if (rank == 0)
  {
    for (size_t i = 0; i != size; ++i)
    {
      g[i] = 0.0;				// sets g(t = 0, x, y = 0) = 0
    }
  }

  if ( rank == (MPI_NUM_THREADS - 1) )
  {
    for (size_t i = 0; i != size; ++i)
    {
      size_t const j = (rows - 1);
      size_t const k = (i + j * size);
      g[k] = 0.0;				// sets g(t = 0, x, y = 1) = 0
    }
  }
}


void init_mask (mpi_workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const load = workspace -> mpi_thread_load_size;
  size_t const rows = (load / size);

  double* mask = workspace -> mask;

  alias_t* masks = mask;
  for (size_t k = 0; k != load; ++k)
  {
    masks[k].bin = iNODE;			// sets all nodes as interior nodes
  }


  for (size_t j = 0; j != rows; ++j)
  {
    size_t const i = 0;
    size_t const k = (i + j * size);
    mask[k] = 0.0;				// masks boundary nodes at y = 0
  }

  for (size_t j = 0; j != rows; ++j)
  {
    size_t const i = (size - 1);
    size_t const k = (i + j * size);
    mask[k] = 0.0;				// masks boundary nodes at y = 1
  }

  if (rank == 0)
  {
    for (size_t i = 0; i != size; ++i)
    {
      mask[i] = 0.0;				// masks boundary nodes at x = 0
    }
  }

  if ( rank == (MPI_NUM_THREADS - 1) )
  {
    for (size_t i = 0; i != size; ++i)
    {
      size_t const j = (rows - 1);
      size_t const k = (i + j * size);
      mask[k] = 0.0;				// masks boundary nodes at x = 1
    }
  }
}


void init_guess (mpi_workspace_t* workspace)
{
  size_t const xsize = workspace -> xsize;
  double* g0 = workspace -> g0;
  zeros(xsize, g0);
}


// sends data forward (that is, to a higher rank thread)
void forward_send (const mpi_workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const load = workspace -> mpi_thread_load_size;
  const double* g0 = ( (workspace -> g0) + load );
  MPI_Send(g0, size, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
}


// receives forwarded data (that is, from a lower rank thread)
void forward_recv (mpi_workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  double* g0 = workspace -> g0;
  MPI_Status* status = workspace -> mpi_status;
  MPI_Recv(g0, size, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, status);
}


// sends data backwards (that is, to a lower rank thread)
void backward_send (const mpi_workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  const double* g0 = ( (workspace -> g0) + size );
  MPI_Send(g0, size, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD);
}


// receives backwarded data (that is, from a higher rank thread)
void backward_recv (mpi_workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const load = workspace -> mpi_thread_load_size;
  double* g0 = ( (workspace -> g0) + (load + size) );
  MPI_Status* status = workspace -> mpi_status;
  MPI_Recv(g0, size, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, status);
}


// transmits data forward (from low to high rank threads)
void forwards (mpi_workspace_t* workspace)
{
  size_t const rank = workspace -> mpi_thread_rank;
  if (rank % 2 == 0)
  {
    forward_send(workspace);
  }
  else
  {
    forward_recv(workspace);
  }

  if (rank % 2 == 1)
  {
    if ( rank != (MPI_NUM_THREADS - 1) )
    {
      forward_send(workspace);
    }
  }
  else
  {
    if (rank != 0)
    {
      forward_recv(workspace);
    }
  }
}


// transmits data backward (from high to low rank threads)
void backwards (mpi_workspace_t* workspace)
{
  size_t const rank = workspace -> mpi_thread_rank;
  if (rank % 2 == 1)
  {
    backward_send(workspace);
  }
  else
  {
    backward_recv(workspace);
  }

  if (rank % 2 == 0)
  {
    if (rank != 0)
    {
      backward_send(workspace);
    }
  }
  else
  {
    if ( rank != (MPI_NUM_THREADS - 1) )
    {
      backward_recv(workspace);
    }
  }
}


void mpi_update (mpi_workspace_t* workspace)
{
  forwards(workspace);
  backwards(workspace);
}


int* mpi_allocs ()
{
  int* mpi_thread_info = malloc( 2 * sizeof(int) );
  if (mpi_thread_info == NULL)
  {
    printf("failed to allocate mpi thread info!\n");
    return NULL;
  }

  return mpi_thread_info;
}


int* mpi_sane_checks (int* mpi_thread_info)
{
  const int* mpi_thread_rank = mpi_thread_info;
  const int* mpi_thread_pool_size = (mpi_thread_info + 1);
  if (*mpi_thread_pool_size != MPI_NUM_THREADS)
  {
    if (*mpi_thread_rank == 0)
    {
      int const request = *mpi_thread_pool_size;
      int const expects = MPI_NUM_THREADS;
      printf("solver expects %d MPI threads but user requested %d \n", expects, request);
    }

    free(mpi_thread_info);
    mpi_thread_info = NULL;
    mpi_thread_pool_size = NULL;
    mpi_thread_rank = NULL;
  }

  return mpi_thread_info;
}


void mpi_testset (mpi_workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const load = workspace -> mpi_thread_load_size;
  double* g0 = ( (workspace -> g0) + size );
  for (size_t i = 0; i != size; ++i)
  {
    g0[i] = -1.0;
  }

  g0 = ( (workspace -> g0) + load );
  for (size_t i = 0; i != size; ++i)
  {
    g0[i] = -1.0;
  }
}


void assert (size_t const rank, double const value, double const expect)
{
  if (value == expect)
  {
    printf("rank: %lu test: PASS\n", rank);
  }
  else
  {
    printf("rank: %lu test: FAIL\n", rank);
  }
}


void mpi_assert (const mpi_workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const load = workspace -> mpi_thread_load_size;

  double sum = 0;
  const double* g0 = workspace -> g0;
  for (size_t i = 0; i != size; ++i)
  {
    sum += g0[i];
  }

  g0 = ( (workspace -> g0) + (load + size) );
  for (size_t i = 0; i != size; ++i)
  {
    sum += g0[i];
  }

  if ( (rank == 0) || ( rank == (MPI_NUM_THREADS - 1) ) )
  {
    double const value = -sum;
    double const expect = size;
    assert(rank, value, expect);
  }
  else
  {
    double const value = -sum;
    double const expect = (2 * size);
    assert(rank, value, expect);
  }
}


void test_mpi_comms ()
{
  // initializes mpi parallel environment:

  int* argc = NULL;
  char*** argv = NULL;
  MPI_Init(argc, argv);

  // queries mpi thread rank and pool size:

  int* mpi_thread_info = mpi_allocs();
  if (mpi_thread_info == NULL)
  {
    MPI_Finalize();
    return;
  }

  int* mpi_thread_rank = mpi_thread_info;
  int* mpi_thread_pool_size = (mpi_thread_info + 1);

  MPI_Comm_size(MPI_COMM_WORLD, mpi_thread_pool_size);
  MPI_Comm_rank(MPI_COMM_WORLD, mpi_thread_rank);

  // sane checks:

  mpi_thread_info = mpi_sane_checks(mpi_thread_info);
  if (mpi_thread_info == NULL)
  {
    MPI_Finalize();
    return;
  }

  // allocates mpi workspace:

  size_t const size = SIZE;
  mpi_workspace_t* workspace = create(size, mpi_thread_info);
  if (workspace == NULL)
  {
    free(mpi_thread_info);
    mpi_thread_info = NULL;
    mpi_thread_pool_size = NULL;
    mpi_thread_rank = NULL;
    MPI_Finalize();
    return;
  }

  // creates test dataset by assigning -1 to data to be transmitted by threads:

  mpi_testset(workspace);

  // data transfer between threads:

  forwards(workspace);
  backwards(workspace);

  // assertions:

  mpi_assert(workspace);

  // frees allocated resources and finalizes parallel environment:

  workspace = destroy(workspace);

  free(mpi_thread_info);
  mpi_thread_info = NULL;
  mpi_thread_pool_size = NULL;
  mpi_thread_rank = NULL;
  MPI_Finalize();
}


void init_rhs(size_t const size,
	      size_t const numel,
	      double* restrict b,
	      const double* restrict x,
	      const double* restrict g)
{
  double const alpha = ALPHA;
  double const dx = (x[1] - x[0]);
  // we don't need to mask this loop since we are masking the boundary nodes elsewhere
  for (size_t i = 0; i != numel; ++i)
  {
    b[i] = (dx * dx + alpha * g[i]);
  }
}


void rhs (size_t const size,
	  size_t const numel,
	  double* restrict g,
	  const double* restrict b,
	  const double* restrict mask)
{
  alias_t* dst = g;
  const alias_t* values = b;
  const alias_t* masks = mask;
  for (size_t i = 0; i != numel; ++i)
  {
    dst[i].bin = (masks[i].bin & values[i].bin);
  }
}


void tridiag (size_t const size,
	      size_t const numel,
	      double* restrict g,
	      double* restrict tmp,
	      const double* restrict g0,
	      const double* restrict mask)
{
  alias_t* t = tmp;
  // NOTE: there's no need to check for invalid reads because `g0' has extra space
  const alias_t* values = g0;
  const alias_t* masks = mask;

  size_t const offset = size;
  for (size_t i = 0; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values[i - 1 + offset].bin);
  }

  add(numel, g, tmp);

  for (size_t i = 0; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values[i + 1 + offset].bin);
  }

  add(numel, g, tmp);
}


void subdiag (size_t const rank,
	      size_t const size,
	      size_t const numel,
	      double* restrict g,
	      double* restrict tmp,
	      const double* restrict g0,
	      const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* values = (g0 + size);
  const alias_t* masks = mask;

  zeros(numel, tmp);

  if (rank == MPI_THREAD_ROOT)
  {
    for (size_t i = size; i != numel; ++i)
    {
      t[i].bin = (masks[i].bin & values[i - size].bin);
    }
  }
  else
  {
    for (size_t i = 0; i != numel; ++i)
    {
      t[i].bin = (masks[i].bin & values[i - size].bin);
    }
  }

  add(numel, g, tmp);
}


void superdiag (size_t const rank,
		size_t const size,
		size_t const numel,
		double* restrict g,
		double* restrict tmp,
		const double* restrict g0,
		const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* values = (g0 + size);
  const alias_t* masks = mask;

  zeros(numel, tmp);

  if ( rank == (MPI_NUM_THREADS - 1) )
  {
    for (size_t i = 0; i != (numel - size); ++i)
    {
      t[i].bin = (masks[i].bin & values[i + size].bin);
    }
  }
  else
  {
    for (size_t i = 0; i != numel; ++i)
    {
      t[i].bin = (masks[i].bin & values[i + size].bin);
    }
  }

  add(numel, g, tmp);
}


void __attribute__ ((noinline)) scale(size_t const size,
				      size_t const numel,
				      double* restrict g,
				      double* restrict tmp,
				      const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* masks = mask;
  double const alpha = ALPHA;
  double const c = 1.0 / (alpha + 4.0);
  alias_t const values = { .data = c };

  for (size_t i = 0; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values.bin);
  }

  mult(numel, g, tmp);
}



void error (size_t const numel,
	    double* restrict e,
	    const double* restrict v,
	    const double* restrict w)
{
  for (size_t i = 0; i != numel; ++i)
  {
    e[i] = (v[i] - w[i]);
  }
}


double mpi_norm (mpi_workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const numel = workspace -> mpi_thread_load_size;

  double* tmp = workspace -> f;
  double* err = workspace -> err;
  const double* g = workspace -> g;
  const double* g0 = ( (workspace -> g0) + size);

  // performs the partial computation of the norm of the error array
  error(numel, err, g, g0);
  norm1(rank, numel, tmp, err);

  if (rank == MPI_THREAD_ROOT)
  {
    // gets the partial computation of the norm of the error array from other threads:
    for (size_t tid = 1; tid != MPI_NUM_THREADS; ++tid)
    {
      int const src = tid;
      int const tag = tid;
      MPI_Count const count = 1;
      double* data = (tmp + tid);
      MPI_Status* status = workspace -> mpi_status;
      MPI_Recv(data, count, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, status);
    }

    // completes the computation of the norm:
    double const e = norm2(tmp);
    // we could simply store the norm in the first element but we can afford to do this
    for (size_t tid = 0; tid != MPI_NUM_THREADS; ++tid)
    {
      tmp[tid] = e;
    }
  }
  else
  {
    MPI_Count const count = 1;
    int const dst = MPI_THREAD_ROOT;
    int const tag = rank;
    const double* data = (tmp + rank);
    // sends the partial computation of the norm to the root thread
    MPI_Send(data, count, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
  }

  // broadcasts the computed norm to the threads:

  double* data = tmp;
  int const root = MPI_THREAD_ROOT;
  MPI_Count const count = MPI_NUM_THREADS;
  // we could really just send the first element and be done with it
  MPI_Bcast(data, count, MPI_DOUBLE, root, MPI_COMM_WORLD);

  double e = tmp[rank];
  return e;
}


void isolver (mpi_workspace_t* workspace)
{
  // parameters:

  double const tol = TOLERANCE;
  size_t const iters = MAX_ITERATIONS;
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const numel = workspace -> mpi_thread_load_size;

  // iterators:

  double* x = workspace -> x;
  double* b = workspace -> rhs;
  double* g = workspace -> g;
  double* g0 = workspace -> g0;
  double* tmp = workspace -> f;
  const double* mask = workspace -> mask;

  // initializations:

  init_rhs(size, numel, b, x, g);		// vectorized by gcc

  // implements iterative solver:

  set_state(workspace, FAILURE_STATE);
  for (size_t i = 0; i != iters; ++i)
  {
    // updates the solution array g(t + dt):
    rhs(size, numel, g, b, mask);			// vectorized by gcc
    tridiag(size, numel, g, tmp, g0, mask);		// vectorized by gcc
    subdiag(rank, size, numel, g, tmp, g0, mask);	// vectorized by gcc
    superdiag(rank, size, numel, g, tmp, g0, mask);	// vectorized by gcc
    scale(size, numel, g, tmp, mask);			// vectorized by gcc

    // checks for convergence:
    if (mpi_norm(workspace) < tol)
    {
      if (rank == 0)
      {
	if (VERBOSE)
	{
	  printf("Jacobi(): solution found after %lu iters\n", i + 1);
	}
      }
      set_state(workspace, SUCCESS_STATE);
      break;
    }

    // updates the initial solution array g(t + dt) for the next iteration:

    const double* src = g;
    double* dst = (g0 + size);
    copy(numel, src, dst);				// optimized by gcc

    // transmits the ``missing'' elements of the initial solution array to other threads:

    mpi_update(workspace);
  }
}


void init_X(size_t const size,
	    size_t const numel,
	    double* restrict X,
	    const double* restrict x)
{
  for (size_t k = 0; k != numel; ++k)
  {
    size_t const i = (k % size);
    X[k] = x[i];
  }
}


void init_Y(size_t const rank,
	    size_t const size,
	    size_t const numel,
	    double* restrict Y,
	    const double* restrict y)
{
  for (size_t k = 0; k != numel; ++k)
  {
    size_t const count = (rank * numel);
    size_t const j = (k + count) / size;
    Y[k] = y[j];
  }
}


void pdesol (double const t, mpi_workspace_t* workspace)
{
  double* F = workspace -> f;
  double* X = workspace -> rhs;
  double* Y = workspace -> err;
  const double* x = workspace -> x;
  const double* y = workspace -> y;
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const numel = workspace -> mpi_thread_load_size;
  size_t const N = 64;

  init_X(size, numel, X, x);
  init_Y(rank, size, numel, Y, y);

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


void export(const char* fname,
	    size_t const size,
	    const double* restrict g,
	    const double* restrict x)
{
  FILE* file = fopen(fname, "w");
  for (size_t i = 0; i != size; ++i)
  {
    fprintf(file, "%.12e %.12e\n", x[i], g[i]);
  }

  fclose(file);
}


void logger (size_t const step, size_t const count, mpi_workspace_t* workspace)
{
  double const alpha = ALPHA;
  const double* x = workspace -> x;
  double const dx = (x[1] - x[0]);
  double const dt = (dx * dx) / alpha;
  double const t = ( ( (double) step ) + 1.0 ) * dt;

  pdesol(t, workspace);

  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const numel = workspace -> mpi_thread_load_size;
  const double* f = workspace -> f;
  const double* g = workspace -> g;
  double* tmp = workspace -> rhs;
  double* err = workspace -> err;

  error(numel, err, f, g);
  norm1(rank, numel, tmp, err);

  double const sumSquaredError = tmp[rank];
  double const rmse = sqrt(sumSquaredError) / ( (double) numel );

  if (rank == MPI_THREAD_ROOT)
  {
    printf("approximation error (transient solution t = %.4e): %e \n", t, rmse);
  }

  // exports the analytic and numeric fields g(t, x, y = 0.5) for visualization:
  if (rank == 2)
  {
    char fname[80];
    sprintf(fname, "numeric%lu.txt", count);
    export(fname, size, g, x);

    sprintf(fname, "analytic%lu.txt", count);
    export(fname, size, f, x);
  }
}


void integrator (mpi_workspace_t* workspace)
{
  size_t const rank = workspace -> mpi_thread_rank;

  size_t count = 0;
  size_t const steps = 0x00001000;
  for (size_t step = 0; step != steps; ++step)
  {
    isolver(workspace);

    size_t state = get_state(workspace);
    if (state == FAILURE_STATE)
    {
      if (rank == 0)
      {
	char const msg [] = "Jacobi solver failed to converge to the solution "
			    "after %d iterations";
	printf(msg, MAX_ITERATIONS);
      }
      break;
    }

    if (LOG)
    {
      size_t const span = (steps / 4);
      // logs error of exact f(t+dt, x) and numeric solution g(t+dt, x) every `span' steps
      if ( ( step != 0 ) && ( (step % span) == 0 ) )
      {
	logger(step, count, workspace);
	++count;
      }
    }
  }
}


void Poisson ()
{
  // initializes mpi parallel environment:

  int* argc = NULL;
  char*** argv = NULL;
  MPI_Init(argc, argv);

  // queries mpi thread rank and pool size:

  int* mpi_thread_info = mpi_allocs();
  if (mpi_thread_info == NULL)
  {
    MPI_Finalize();
    return;
  }

  int* mpi_thread_rank = mpi_thread_info;
  int* mpi_thread_pool_size = (mpi_thread_info + 1);

  MPI_Comm_size(MPI_COMM_WORLD, mpi_thread_pool_size);
  MPI_Comm_rank(MPI_COMM_WORLD, mpi_thread_rank);

  // sane checks:

  mpi_thread_info = mpi_sane_checks(mpi_thread_info);
  if (mpi_thread_info == NULL)
  {
    MPI_Finalize();
    return;
  }

  // allocates mpi workspace:

  size_t const size = SIZE;
  mpi_workspace_t* workspace = create(size, mpi_thread_info);
  if (workspace == NULL)
  {
    free(mpi_thread_info);
    mpi_thread_info = NULL;
    mpi_thread_pool_size = NULL;
    mpi_thread_rank = NULL;
    MPI_Finalize();
    return;
  }

  // initializations:

  init_mask(workspace);
  init_field(workspace);

  double const x_l = 0;                         // defines x-axis lower bound
  double const x_u = 1;                         // defines x-axis upper bound
  double* x = workspace -> x;
  linspace(x, x_l, x_u, size);                  // inits the x-axis position array

  double const y_l = 0;                         // defines y-axis lower bound
  double const y_u = 1;                         // defines y-axis upper bound
  double* y = workspace -> y;
  linspace(y, y_l, y_u, size);                  // inits the y-axis position array

  // integrates the Poisson equation:

  integrator(workspace);

  // TODO post-processing

  /*
  size_t const rank = workspace -> mpi_thread_rank;
  if (rank == 2)
  {
    const double* g = workspace -> g;

    for (size_t i = 0; i != size; ++i)
    {
      printf("%f\n", g[i]);
    }
  }
  */

  // frees allocated resources and finalizes parallel environment:

  workspace = destroy(workspace);

  free(mpi_thread_info);
  mpi_thread_info = NULL;
  mpi_thread_pool_size = NULL;
  mpi_thread_rank = NULL;
  MPI_Finalize();
}

// COMMENTS:
// iNODE:       interior node, a nodes that do not lie at the boundary


// Memcheck:
//
// To perform the memory check you can use the following command:
//
// valgrind -s mpiexec.mpich -n 4 ./test-vectorization
//
// possible output:
//
// ==9618== Memcheck, a memory error detector
// ==9618== Copyright (C) 2002-2017, and GNU GPL'd, by Julian Seward et al.
// ==9618== Using Valgrind-3.15.0 and LibVEX; rerun with -h for copyright info
// ==9618== Command: mpiexec.mpich -n 4 ./test-vectorization
// ==9618==
// rank: 0 test: PASS
// rank: 1 test: PASS
// rank: 2 test: PASS
// rank: 3 test: PASS
// ==9618==
// ==9618== HEAP SUMMARY:
// ==9618==     in use at exit: 0 bytes in 0 blocks
// ==9618==   total heap usage: 1,153 allocs, 1,153 frees, 2,372,192 bytes allocated
// ==9618==
// ==9618== All heap blocks were freed -- no leaks are possible
// ==9618==
// ==9618== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)


// TODO:
// [x] partition forward and backward methods into send_forward() recv_forward() and
//     send_backward() and recv_backward()
// [x] implement reduction computation for the vector norm
// [ ] manage allocation errors in create() method (without assumptions)
