#ifndef GUARD_ITERATIVE_LINEAR_SOLVER_H
#define GUARD_ITERATIVE_LINEAR_SOLVER_H
/*
 * Computational Methods                             	 September 10, 2022
 * ICI 70320
 * Prof. M Diaz-Maldonado
 *
 * source: vector.h
 *
 * Synopsis:
 * Header file for iterative linear solvers.
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
 *
 */

typedef struct {
double tol;		// tolerance
size_t iters;		// maximum number of iterations
} isolver_prms_t;

#endif
