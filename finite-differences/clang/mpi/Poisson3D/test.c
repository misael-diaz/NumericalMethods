/*
 * Transient Heat Conduction                                	July 04, 2023
 *
 * source: test.c
 * author: @misael-diaz
 *
 * Synopsis:
 * Solves the transient 3D Poisson equation iteratively.
 * Parallelizes the numeric solver with MPI.
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
 * [4] W Deen, Analysis of Transport Phenomena, 2nd edition
 *
 */


#include <stdint.h>
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
#define MPI_NUM_THREADS 16
#define MPI_THREAD_ROOT 0


typedef union
{
  uint64_t bin;		// binary pattern of the double precision floating-point data
  double data;		// double precision floating-point data
} alias_t;


typedef struct
{
  double* x;			// position array along the x-axis
  double* y;			// position array along the y-axis
  double* z;			// position array along the z-axis
  double* f;			// exact solution array, f(t, x, y, z)
  double* g;			// next estimate of the solution array, g(t + dt, x, y, z)
  double* rhs;			// Right Hand Side RHS array of the PDE
  double* err;			// error array
  double* tmp;			// array temporary
  double* mask;			// bitmask for protecting boundary node data
  double* g0;			// estimate of the solution array, g(t + dt, x, y, z)
  double* data;			// contiguous data container
  MPI_Status* mpi_status;	// MPI Communication status
  size_t mpi_thread_pool_size;	// number of MPI threads
  size_t mpi_thread_load_size;	// workload size of MPI thread (same for each)
  size_t mpi_thread_rank;	// MPI Thread ID
  size_t xsize;			// extended array size
  size_t state;			// solver state
  size_t size;			// position array size (same for x, y, z arrays)
} workspace_t;

void Poisson();
void test_mpi_comms();

int main ()
{
//test_mpi_comms();	OK!
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


// computes in parallel the norm of an array
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


// root (or master) thread performs the reduction operation to obtain the norm of an array
double norm2 (const double* err)
{
  double sum = 0;
  for (size_t tid = 0; tid != MPI_NUM_THREADS; ++tid)
  {
    sum += err[tid];
  }

  return sum;
}


// implements array addition
void add (size_t const size,
	  double* restrict dst,
	  const double* restrict src)
{
  for (size_t i = 0; i != size; ++i)
  {
    dst[i] += src[i];
  }
}


// implements array multiplication
void mult(size_t const numel,
	  double* restrict dst,
	  const double* restrict src)
{
  for (size_t i = 0; i != numel; ++i)
  {
    dst[i] *= src[i];
  }
}

// copies source `src' into destination `dst' array (GCC optimizes it into a library call)
void copy(size_t const size,
	  const double* restrict src,
	  double* restrict dst)
{
  for (size_t i = 0; i != size; ++i)
  {
    dst[i] = src[i];
  }
}


// creates the workspace for the Poisson solver
workspace_t* create (size_t const size, const int* mpi_thread_info)
{
  if (MPI_NUM_THREADS % 2)
  {
    // keep guard if you want to export the field at the center g(t, x, y ~ 0.5, z ~ 0.5)
    printf("the number of mpi threads must be even\n");
    return NULL;
  }

  size_t const numel = (size * size * size);
  if (numel % MPI_NUM_THREADS)
  {
    printf("the number of nodes must be a multiple of the number of threads\n");
    return NULL;
  }

  size_t const size2 = (size * size);
  size_t const load = (numel / MPI_NUM_THREADS);
  if (load < size2)
  {
    printf("too many threads to support the domain division scheme\n");
    return NULL;
  }

  size_t const x_size = size;
  size_t const y_size = size;
  size_t const z_size = size;
  size_t const f_size = load;
  size_t const g_size = load;
  size_t const rhs_size = load;
  size_t const err_size = load;
  size_t const tmp_size = load;
  size_t const mask_size = load;
  // allots extended size to store the data transferred by other threads
  size_t const xsize = (load + 2 * size2);
  size_t const g0_size = xsize;

  size_t const data_size = x_size +
			   y_size +
			   z_size +
			   f_size +
			   g_size +
			   rhs_size +
			   err_size +
			   tmp_size +
			   mask_size +
			   g0_size;

  double* data = malloc( data_size * sizeof(double) );
  if (data == NULL)
  {
    printf("failed to allocate the data placeholder\n");
    return NULL;
  }

  workspace_t *workspace = malloc( sizeof(workspace_t) );
  if (workspace == NULL)
  {
    free(data);
    data = NULL;
    printf("failed to allocate workspace data structure\n");
    return NULL;
  }

  workspace -> mpi_status = (MPI_Status*) malloc( sizeof(MPI_Status) );
  if (workspace -> mpi_status == NULL)
  {
    free(data);
    free(workspace);
    data = NULL;
    workspace = NULL;
    printf("failed to allocate MPI status placeholder!\n");
    return NULL;
  }

  // assumes successful allocations for simplicity
  workspace -> x = data;
  workspace -> y = workspace -> x + size;
  workspace -> z = workspace -> y + size;
  workspace -> f = workspace -> z + size;
  workspace -> g = workspace -> f + load;
  workspace -> rhs = workspace -> g + load;
  workspace -> err = workspace -> rhs + load;
  workspace -> tmp = workspace -> err + load;
  workspace -> mask = workspace -> tmp + load;
  workspace -> g0 = workspace -> mask + load;
  workspace -> data = data;

  double* x = workspace -> x;
  double* y = workspace -> y;
  double* z = workspace -> z;
  double* f = workspace -> f;
  double* g = workspace -> g;
  double* rhs = workspace -> rhs;
  double* err = workspace -> err;
  double* tmp = workspace -> tmp;
  double* mask = workspace -> mask;
  double* g0 = workspace -> g0;

  zeros(size, x);
  zeros(size, y);
  zeros(size, z);
  zeros(load, f);
  zeros(load, g);
  zeros(load, rhs);
  zeros(load, err);
  zeros(load, tmp);
  zeros(load, mask);
  zeros(xsize, g0);

  workspace -> size = size;
  workspace -> xsize = xsize;
  workspace -> state = FAILURE_STATE;
  const int* mpi_thread_rank = mpi_thread_info;
  const int* mpi_thread_pool_size = (mpi_thread_info + 1);
  workspace -> mpi_thread_rank = ( (size_t) *mpi_thread_rank );
  workspace -> mpi_thread_pool_size = ( (size_t) *mpi_thread_pool_size );
  workspace -> mpi_thread_load_size = load;

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
  free(workspace -> mpi_status);

  workspace -> data = NULL;
  workspace -> mpi_status = NULL;

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


// void init_field (workspace_t* workspace)
//
// Synopsis:
// Sets the initial (temperature) field g(t = 0, x, y, z).
// Fills with ones at the interior nodes and with zeros at the boundary nodes.
//
// Input:
// workspace	data structure containing the current data (field, position, error, etc.)
//
// Output:
// g            updated data structure with the initial field g(t = 0, x, y, z)


void init_field (workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const load = workspace -> mpi_thread_load_size;
  size_t const size2 = (size * size);
  size_t const levels = (load / size2);
  size_t const lvls = levels;

  double* g = workspace -> g;

  // sets g(t = 0, x, y, z) = 1 (everywhere)
  ones(load, g);

  // sets g(t = 0, x = 0, y, z) = 0:
  for (size_t k = 0; k != lvls; ++k)
  {
    for (size_t j = 0; j != size; ++j)
    {
      size_t const i = 0;
      size_t const count = (i + j * size + k * size2);
      g[count] = 0.0;
    }
  }

  // sets g(t = 0, x = 1, y, z) = 0:
  for (size_t k = 0; k != lvls; ++k)
  {
    for (size_t j = 0; j != size; ++j)
    {
      size_t const i = (size - 1);
      size_t const count = (i + j * size + k * size2);
      g[count] = 0.0;
    }
  }

  // sets g(t = 0, x, y = 0, z) = 0:
  for (size_t k = 0; k != lvls; ++k)
  {
    for (size_t i = 0; i != size; ++i)
    {
      size_t const j = 0;
      size_t const count = (i + j * size + k * size2);
      g[count] = 0.0;
    }
  }

  // sets g(t = 0, x, y = 1, z) = 0:
  for (size_t k = 0; k != lvls; ++k)
  {
    for (size_t i = 0; i != size; ++i)
    {
      size_t const j = (size - 1);
      size_t const count = (i + j * size + k * size2);
      g[count] = 0.0;
    }
  }

  // sets g(t = 0, x, y, z = 0) = 0:
  if (rank == 0)
  {
    for (size_t i = 0; i != size2; ++i)
    {
      g[i] = 0.0;
    }
  }

  // sets g(t = 0, x, y, z = 1) = 0:
  if ( rank == (MPI_NUM_THREADS - 1) )
  {
    g += (load - size2);
    for (size_t i = 0; i != size2; ++i)
    {
      g[i] = 0.0;
    }
  }
}


// void init_mask(size_t size, double* mask)
//
// Synopsis:
// Initializes the node masking.
// The mask is equal to zero at boundary nodes and equal to `iNODE' at interior nodes.
//
// NOTES:
// This is used to protect the boundary nodes data while updating all the nodes; this
// also simplifies the implementation of the parallelized code for the threads work is
// hardcoded into the mask.
//
// Input:
// workspace	data structure containing the current data (field, position, error, etc.)
//
// Output:
// g            updated data structure with the mask


void init_mask (workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const load = workspace -> mpi_thread_load_size;
  size_t const size2 = (size * size);
  size_t const levels = (load / size2);
  size_t const lvls = levels;

  double* mask = workspace -> mask;

  alias_t* masks = mask;
  // sets all nodes as interior nodes
  for (size_t i = 0; i != load; ++i)
  {
    masks[i].bin = iNODE;
  }

  // masks boundary nodes at x = 0:
  for (size_t k = 0; k != lvls; ++k)
  {
    for (size_t j = 0; j != size; ++j)
    {
      size_t const i = 0;
      size_t const count = (i + j * size + k * size2);
      mask[count] = 0.0;
    }
  }

  // masks boundary nodes at x = 1:
  for (size_t k = 0; k != lvls; ++k)
  {
    for (size_t j = 0; j != size; ++j)
    {
      size_t const i = (size - 1);
      size_t const count = (i + j * size + k * size2);
      mask[count] = 0.0;
    }
  }

  // masks boundary nodes at y = 0:
  for (size_t k = 0; k != lvls; ++k)
  {
    for (size_t i = 0; i != size; ++i)
    {
      size_t const j = 0;
      size_t const count = (i + j * size + k * size2);
      mask[count] = 0.0;
    }
  }

  // masks boundary nodes at y = 1:
  for (size_t k = 0; k != lvls; ++k)
  {
    for (size_t i = 0; i != size; ++i)
    {
      size_t const j = (size - 1);
      size_t const count = (i + j * size + k * size2);
      mask[count] = 0.0;
    }
  }

  // masks boundary nodes at z = 0:
  if (rank == 0)
  {
    for (size_t i = 0; i != size2; ++i)
    {
      mask[i] = 0.0;
    }
  }

  // masks boundary nodes at z = 1:
  if ( rank == (MPI_NUM_THREADS - 1) )
  {
    mask += (load - size2);
    for (size_t i = 0; i != size2; ++i)
    {
      mask[i] = 0.0;
    }
  }
}


// initializes the guess of the solution field g(t + dt, x, y, z)
void init_guess (workspace_t* workspace)
{
  size_t const xsize = workspace -> xsize;
  double* g0 = workspace -> g0;
  zeros(xsize, g0);
}


// sends data forward (that is, to a higher rank thread)
void forward_send (const workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const load = workspace -> mpi_thread_load_size;
  size_t const size2 = (size * size);
  const double* g0 = ( (workspace -> g0) + load );
  MPI_Send(g0, size2, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
}


// receives forwarded data (that is, from a lower rank thread)
void forward_recv (workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const size2 = (size * size);
  double* g0 = workspace -> g0;
  MPI_Status* status = workspace -> mpi_status;
  MPI_Recv(g0, size2, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, status);
}


// sends data backwards (that is, to a lower rank thread)
void backward_send (const workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const size2 = (size * size);
  const double* g0 = ( (workspace -> g0) + size2 );
  MPI_Send(g0, size2, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD);
}


// receives backwarded data (that is, from a higher rank thread)
void backward_recv (workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const load = workspace -> mpi_thread_load_size;
  size_t const size2 = (size * size);
  double* g0 = ( (workspace -> g0) + (size2 + load) );
  MPI_Status* status = workspace -> mpi_status;
  MPI_Recv(g0, size2, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, status);
}


// sends data forward (to higher rank threads)
void forwards (workspace_t* workspace)
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


// sends data backward (to lower rank threads)
void backwards (workspace_t* workspace)
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


// updates the field with the data from other threads for the next solver iteration
void mpi_update (workspace_t* workspace)
{
  forwards(workspace);
  backwards(workspace);
}


// performs MPI related memory allocations
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


// performs MPI sane checks
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


// creates a testset to check that the data sends/receives are working as expected
void mpi_testset (workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const size2 = (size * size);
  // we fill with negative ones the data to be transferred to other threads
  double* g0 = ( (workspace -> g0) + size2 );
  for (size_t i = 0; i != size2; ++i)
  {
    g0[i] = -1.0;
  }

  size_t const load = workspace -> mpi_thread_load_size;
  g0 = ( (workspace -> g0) + load );
  for (size_t i = 0; i != size2; ++i)
  {
    g0[i] = -1.0;
  }
}


// user-defined assertion that does not abort the numerical application
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


// this is were we verify that the data transfers are fine (we are just counting nodes)
void mpi_assert (const workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const load = workspace -> mpi_thread_load_size;

  double sum = 0;
  size_t const size2 = (size * size);
  const double* g0 = workspace -> g0;
  for (size_t i = 0; i != size2; ++i)
  {
    sum += g0[i];
  }

  g0 = ( (workspace -> g0) + (size2 + load) );
  for (size_t i = 0; i != size2; ++i)
  {
    sum += g0[i];
  }

  if ( (rank == 0) || ( rank == (MPI_NUM_THREADS - 1) ) )
  {
    double const value = -sum;
    double const expect = size2;
    assert(rank, value, expect);
  }
  else
  {
    double const value = -sum;
    double const expect = (2 * size2);
    assert(rank, value, expect);
  }
}


// driver code to check the data transfers
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
  workspace_t* workspace = create(size, mpi_thread_info);
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


// void init_rhs (size_t size, double* b, double* x, double* g)
//
// Synopsis:
// Initializes the Right-Hand-Size RHS of the Finite Difference Equations FDEs.
//
// NOTES:
// We don't need to use node masking here because we are masking the boundary nodes
// elsewhere and because the boundary conditions are zero.
//
// Inputs:
// size         number of nodes along the x, y, or z axis
// numel	number of elements (or nodes) in the partition assigned to the thread
// x            x-axis position array
// g            field array g(t, x, y, z) at the current time
//
// Outputs:
// b            RHS array


void init_rhs(size_t const size,
	      size_t const numel,
	      double* restrict b,
	      const double* restrict x,
	      const double* restrict g)
{
  double const alpha = ALPHA;
  double const dx = (x[1] - x[0]);
  for (size_t i = 0; i != numel; ++i)
  {
    b[i] = (dx * dx + alpha * g[i]);
  }
}


// void rhs(size_t size, size_t numel, double* g, double* b, double* mask)
//
// Synopsis:
// Updates the estimate of the solution array g(t + dt, x, y, z) from the RHS of the FDEs.
//
// Inputs:
// size         number of nodes along the x, y, or z axis
// numel	number of elements (or nodes) in the partition assigned to the thread
// b            (read-only) RHS array
// mask         (read-only) bitmask for protecting boundary node data
//
// Outputs:
// g            estimate of the solution array g(t + dt, x, y, z)


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


// void tridiag(size_t size, size_t numel, double* g, double* g0, double* mask)
//
// Synopsis:
// Updates the solution array g(t + dt, x, y, z) from the tridiagonal terms of the FDEs.
//
// NOTES:
// We do not need to guard against seemingly invalid reads because `g0' has extra space.
//
// Inputs:
// size         number of nodes along the x, y, or z axis
// numel	number of elements (or nodes) in the partition assigned to the thread
// tmp          array temporary for storing intermediate computations
// g0           (read-only) previous estimate of the solution array g(t + dt, x, y, z)
// mask         (read-only) bitmask for protecting boundary node data
//
// Outputs:
// g            estimate of the solution array g(t + dt, x, y, z)


void tridiag (size_t const size,
	      size_t const numel,
	      double* restrict g,
	      double* restrict tmp,
	      const double* restrict g0,
	      const double* restrict mask)
{
  alias_t* t = tmp;
  size_t const size2 = (size * size);
  size_t const offset = size2;
  const alias_t* values = (g0 + offset);
  const alias_t* masks = mask;

  for (size_t i = 0; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values[i - 1].bin);
  }

  add(numel, g, tmp);

  for (size_t i = 0; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values[i + 1].bin);
  }

  add(numel, g, tmp);
}


// void subdiag(size_t size, size_t numel, double* g, double* g0, double* mask)
//
// Synopsis:
// Updates the solution array g(t + dt, x, y, z) from the sub-diagonal terms of the FDEs.
//
// Inputs:
// size         number of nodes along the x, y, or z axis
// numel	number of elements (or nodes) in the partition assigned to the thread
// tmp          array temporary for storing intermediate computations
// g0           (read-only) previous estimate of the solution array g(t + dt, x, y, z)
// mask         (read-only) bitmask for protecting boundary node data
//
// Outputs:
// g            estimate of the solution array g(t + dt, x, y, z)


void subdiag (size_t const size,
	      size_t const numel,
	      double* restrict g,
	      double* restrict tmp,
	      const double* restrict g0,
	      const double* restrict mask)
{
  alias_t* t = tmp;
  size_t const size2 = (size * size);
  size_t const offset = size2;
  const alias_t* values = (g0 + offset);
  const alias_t* masks = mask;

  for (size_t i = 0; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values[i - size].bin);
  }

  add(numel, g, tmp);
}


// void superdiag(size_t size, size_t numel, double* g, double* g0, double* mask)
//
// Synopsis:
// Updates the solution array g(t + dt, x, y, z) from the super-diagonal terms of the FDEs
//
// Inputs:
// size         number of nodes along the x, y, or z axis
// numel        number of elements (or nodes) in the partition assigned to the thread
// tmp          array temporary for storing intermediate computations
// g0           (read-only) previous estimate of the solution array g(t + dt, x, y, z)
// mask         (read-only) bitmask for protecting boundary node data
//
// Outputs:
// g            estimate of the solution array g(t + dt, x, y, z)


void superdiag (size_t const size,
		size_t const numel,
		double* restrict g,
		double* restrict tmp,
		const double* restrict g0,
		const double* restrict mask)
{
  alias_t* t = tmp;
  size_t const size2 = (size * size);
  size_t const offset = size2;
  const alias_t* values = (g0 + offset);
  const alias_t* masks = mask;

  for (size_t i = 0; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values[i + size].bin);
  }

  add(numel, g, tmp);
}


// as subdiag() but uses the lowest band to update g(t + dt, x, y, z)
void subband (size_t const size,
	      size_t const numel,
	      double* restrict g,
	      double* restrict tmp,
	      const double* restrict g0,
	      const double* restrict mask)
{
  alias_t* t = tmp;
  size_t const size2 = (size * size);
  size_t const offset = size2;
  const alias_t* values = (g0 + offset);
  const alias_t* masks = mask;

  for (size_t i = 0; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values[i - size2].bin);
  }

  add(numel, g, tmp);
}


// as superdiag() but uses the highest band to update g(t + dt, x, y, z)
void superband (size_t const size,
		size_t const numel,
		double* restrict g,
		double* restrict tmp,
		const double* restrict g0,
		const double* restrict mask)
{
  alias_t* t = tmp;
  size_t const size2 = (size * size);
  size_t const offset = size2;
  const alias_t* values = (g0 + offset);
  const alias_t* masks = mask;

  for (size_t i = 0; i != numel; ++i)
  {
    t[i].bin = (masks[i].bin & values[i + size2].bin);
  }

  add(numel, g, tmp);
}


// void scale (size_t size, size_t numel, double* g)
//
// Synopsis:
// Scales the solution array g(t + dt, x, y, z) with the main diagonal coefficient.
//
// Inputs:
// size         number of nodes along the x, y, or z axis
// numel        number of elements (or nodes) in the partition assigned to the thread
// tmp          array temporary for storing intermediate computations
// mask         (read-only) bitmask for protecting boundary node data
//
// Outputs:
// g            estimate of the solution array g(t + dt, x, y, z)


void scale (size_t const size,
	    size_t const numel,
	    double* restrict g,
	    double* restrict tmp,
	    const double* restrict mask)
{
  alias_t* t = tmp;
  const alias_t* masks = mask;
  double const alpha = ALPHA;
  double const c = 1.0 / (alpha + 6.0);
  alias_t const values = { .data = c };

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
// size         array size (same for `v', `w', and `err')
// v            read-only array
// w            read-only array
//
// Outputs:
// e            error array, stores the elementwise differences


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


// implements the parallelized computation of the norm of the solution array estimates
double mpi_norm (workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const numel = workspace -> mpi_thread_load_size;
  size_t const size2 = (size * size);
  size_t const offset = size2;

  double* tmp = workspace -> f;
  double* err = workspace -> err;
  const double* g = workspace -> g;
  const double* g0 = (workspace -> g0 + offset);

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


// void isolver(workspace_t* workspace)
//
// Synopsis:
// Iterative Solver.
// Solves the transient 3D Poisson FDEs with the Jacobi method.
// The solver finds the solution g(t + dt, x, y, z) at the next time step from the current
// g(t, x, y, z).
//
// Inputs:
// workspace    data structure containing the current data (field, position, error, etc.)
//
// Outputs:
// workspace    updates the fields g0 and g with the data of the next time step


void isolver (workspace_t* workspace)
{
  // parameters:

  double const tol = TOLERANCE;
  size_t const iters = MAX_ITERATIONS;
  size_t const size = workspace -> size;
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const numel = workspace -> mpi_thread_load_size;
  size_t const size2 = (size * size);
  size_t const offset = size2;

  // iterators:

  double* x = workspace -> x;
  double* b = workspace -> rhs;
  double* g = workspace -> g;
  double* g0 = workspace -> g0;
  double* tmp = workspace -> tmp;
  const double* mask = workspace -> mask;

  // initializations:

  init_rhs(size, numel, b, x, g);			// vectorized by gcc

  // implements iterative solver:

  set_state(workspace, FAILURE_STATE);
  for (size_t i = 0; i != iters; ++i)
  {
    // updates the solution array g(t + dt):
    rhs(size, numel, g, b, mask);			// vectorized by gcc
    tridiag(size, numel, g, tmp, g0, mask);		// vectorized by gcc
    subdiag(size, numel, g, tmp, g0, mask);		// vectorized by gcc
    superdiag(size, numel, g, tmp, g0, mask);		// vectorized by gcc
    subband(size, numel, g, tmp, g0, mask);		// vectorized by gcc
    superband(size, numel, g, tmp, g0, mask);		// vectorized by gcc
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

    // updates the initial solution array g(t + dt, x, y, z) for the next iteration:

    const double* src = g;
    double* dst = (g0 + offset);
    copy(numel, src, dst);				// optimized by gcc

    // sends the data that other threads need for the next iteration:

    mpi_update(workspace);
  }
}


// convenience function that stores the x position array for the nodes assigned to thread
void init_X(size_t const rank,
	    size_t const size,
	    size_t const numel,
	    double* restrict X,
	    const double* restrict x)
{
  for (size_t n = 0; n != numel; ++n)
  {
    size_t const i = (n % size);
    X[n] = x[i];
  }
}


// convenience function that stores the y position array for the nodes assigned to thread
void init_Y(size_t const rank,
	    size_t const size,
	    size_t const numel,
	    double* restrict Y,
	    const double* restrict y)
{
  for (size_t n = 0; n != numel; ++n)
  {
    size_t const size2 = (size * size);
    size_t const count = (rank * numel);
    size_t const j = ( (n + count) % size2 ) / size;
    Y[n] = y[j];
  }
}


// convenience function that stores the z position array for the nodes assigned to thread
void init_Z(size_t const rank,
	    size_t const size,
	    size_t const numel,
	    double* restrict Z,
	    const double* restrict z)
{
  for (size_t n = 0; n != numel; ++n)
  {
    size_t const size2 = (size * size);
    size_t const count = (rank * numel);
    size_t const k = (n + count) / size2;
    Z[n] = z[k];
  }
}


// void pdesol (double t, workspace_t* workspace)
//
// Synopsis:
// Computes the transient analytic field array f(t, x, y, z).
// The analytic field was obtained by applying the Finite Fourier Transform FFT method.
// Uses the (above) convenience functions so that with one index we can address the
// position and the field arrays simultaneously.
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
  size_t const rank = workspace -> mpi_thread_rank;
  size_t const numel = workspace -> mpi_thread_load_size;
  size_t const N = 16;

  init_X(rank, size, numel, X, x);
  init_Y(rank, size, numel, Y, y);
  init_Z(rank, size, numel, Z, z);

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


// void export (char* fname, size_t size, double* g, double* x)
//
// Synopsis:
// Exports the field g(t, x, y ~ 0.5, z ~ 0.5) array to a plain text file.
//
// Input:
// fname        filename
// size         number of nodes along the x, y, or z axis
// g            field array (number of elements: size * size * size)
// x            position array
//
// Output:
// [x, g]       writes the field g(t, x, y ~ 0.5, z ~ 0.5) with respect to position x


void export(const char* fname,
	    size_t const size,
	    const double* restrict g,
	    const double* restrict x)
{
  FILE* file = fopen(fname, "w");
  if (file == NULL)
  {
    // we don't abort execution, for that would incur in a memory leak
    printf("IO Error: failed to open %s for writing\n", fname);
    return;
  }

  size_t const size2 = (size * size);
  size_t const offset = (size2 / 2);
  g += offset;
  for (size_t i = 0; i != size; ++i)
  {
    fprintf(file, "%.12e %.12e\n", x[i], g[i]);
  }

  fclose(file);
}


// void logger (int step, workspace_t* workspace)
//
// Synopsis:
// Logs the Root Mean Squared Error RMSE of the numeric solution.
// For simplicity we are just logging the RMSE computed by the root (or master) thread.
//
// Inputs:
// step         step number (or id)
// count        the number of times the logger has been invoked
//
// Output:
// rmse         logs the RMSE = sqrt( sum( (f - g)**2 ) ) / N on the console, where
//              `f' is the analytic and `g' is the numeric field array, and `N' is the
//              array size.


void logger (size_t const step, size_t const count, workspace_t* workspace)
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
  double* tmp = workspace -> tmp;
  double* err = workspace -> err;

  error(numel, err, f, g);
  norm1(rank, numel, tmp, err);

  double const sumSquaredError = tmp[rank];
  double const rmse = sqrt(sumSquaredError) / ( (double) numel );

  if (rank == MPI_THREAD_ROOT)
  {
    printf("approximation error (transient solution t = %.4e): %e \n", t, rmse);
  }

  // exports the analytic and numeric fields g(t, x, y ~ 0.5, z ~ 0.5) for visualization:
  if (rank == MPI_NUM_THREADS / 2)
  {
    char fname[80];
    sprintf(fname, "numeric%lu.txt", count);
    export(fname, size, g, x);

    sprintf(fname, "analytic%lu.txt", count);
    export(fname, size, f, x);
  }
}


// void integrator(workspace_t* workspace)
//
// Synopsis:
// Integrates the transient Poisson FDEs via time implicit scheme.
//
// Inputs:
// workspace    data structure containing the initial data (field, position, error, etc.)
//
// Outputs:
// workspace    contains the data after completing the integration


void integrator (workspace_t* workspace)
{
  size_t const rank = workspace -> mpi_thread_rank;

  size_t count = 0;
  size_t const steps = 0x00000800;
  for (size_t step = 0; step != steps; ++step)
  {
    isolver(workspace);

    size_t state = get_state(workspace);
    if (state == FAILURE_STATE)
    {
      if (rank == 0)
      {
	char const msg [] = "Jacobi solver failed to converge to the solution "
			    "after %d iterations\n";
	printf(msg, MAX_ITERATIONS);
      }
      break;
    }

    if (LOG)
    {
      size_t const span = (steps / 8);
      // logs error of exact f(t+dt, *) and numeric solution g(t+dt, *) every `span' steps
      if ( ( step != 0 ) && ( (step % span) == 0 ) )
      {
	logger(step, count, workspace);
	++count;
      }
    }
  }
}


// checks the interior nodes `inodes' count
void test_inodes_count (size_t const tnum,
			size_t const rank,
			size_t const load,
			size_t const size,
			size_t const count)
{
  size_t const size2 = (size * size);
  size_t const lvls = (load / size2);

  // counts the number of boundary nodes `bnodes' from the faces, edges, and corners:
  if ( rank == 0 || rank == (MPI_NUM_THREADS - 1) )
  {
    size_t const bnodes = (size - 2) * (size - 2) +		// bottom (top) face
			  4 * (size - 2) * (lvls - 1) +		// side faces
			  4 * (size - 2) +			// edges in xy-plane
			  4 * (lvls - 1) +			// edges in x(y)z-plane
			  4;					// corners
    size_t const inodes = (load - bnodes);
    if (count != inodes)
    {
      printf("rank: %lu test[%lu]: FAIL\n", tnum, rank);
    }
    else
    {
      printf("rank: %lu test[%lu]: PASS\n", tnum, rank);
    }
  }
  else
  {
    size_t const bnodes = 4 * (lvls - 2) * (size - 2) +		// sides
			  8 * (size - 2) +			// edges in xy-plane
			  4 * (lvls - 2) +			// edges in x(y)z-plane
			  8;					// corners
    size_t const inodes = (load - bnodes);
    if (count != inodes)
    {
      printf("rank: %lu test[%lu]: FAIL\n", rank, tnum);
    }
    else
    {
      printf("rank: %lu test[%lu]: PASS\n", rank, tnum);
    }
  }
}


// tests the system initialization
void test_init (const workspace_t* workspace)
{
  size_t const size = workspace -> size;
  size_t const load = workspace -> mpi_thread_load_size;
  size_t const rank = workspace -> mpi_thread_rank;

  double sum = 0;
  const double* g = workspace -> g;
  for (size_t i = 0; i != load; ++i)
  {
    sum += g[i];
  }

  test_inodes_count(0, rank, load, size, sum);

  sum = 0;
  const alias_t* m = workspace -> mask;
  for (size_t i = 0; i != load; ++i)
  {
    sum += (m[i].bin & 1);
  }

  test_inodes_count(1, rank, load, size, sum);
}


// void Poisson()
//
// Synopsis:
// Solves the 3d transient Poisson equation.
// The initial temperature field is one at the interior nodes and zero at the boundaries.
// There is a uniform heat source throughout the domain.


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
  workspace_t* workspace = create(size, mpi_thread_info);
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

  double const lb = 0.0;			// defines the axis lower bound
  double const ub = 1.0;			// defines the axis upper bound
  double const x_l = lb;			// defines x-axis lower bound
  double const x_u = ub;			// defines x-axis upper bound
  double* x = workspace -> x;
  linspace(x, x_l, x_u, size);                  // inits the x-axis position array

  double const y_l = lb;                        // defines y-axis lower bound
  double const y_u = ub;                        // defines y-axis upper bound
  double* y = workspace -> y;
  linspace(y, y_l, y_u, size);                  // inits the y-axis position array

  double const z_l = lb;                        // defines z-axis lower bound
  double const z_u = ub;                        // defines z-axis upper bound
  double* z = workspace -> z;
  linspace(z, z_l, z_u, size);                  // inits the z-axis position array

  // tests the system initialization by checking the number of interior nodes `inodes':

  test_init(workspace);

  // integrates the Poisson equation:

  integrator(workspace);

  // frees allocated resources and finalizes parallel environment:

  workspace = destroy(workspace);

  free(mpi_thread_info);
  mpi_thread_info = NULL;
  mpi_thread_pool_size = NULL;
  mpi_thread_rank = NULL;
  MPI_Finalize();
}
