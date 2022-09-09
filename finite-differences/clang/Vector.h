#ifndef GUARD_VECTOR_CLASS_H
#define GUARD_VECTOR_CLASS_H
/*
 * Computational Methods                                 September 09, 2022
 * ICI 70320
 * Prof. M Diaz-Maldonado
 *
 * source: Vector.h
 *
 * Synopsis:
 * Header file for the vector class.
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

#include "vector.h"

typedef struct {
// constructors:
vector_t* (*const create) (size_t);
vector_t* (*const linspace) (double, double, size_t);
// destructor:
vector_t* (*const destroy) (vector_t*);
} vector_namespace;

#endif
