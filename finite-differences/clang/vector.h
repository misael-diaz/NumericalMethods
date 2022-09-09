#ifndef GUARD_VECTOR_H
#define GUARD_VECTOR_H
/*
 * Computational Methods                             	 September 09, 2022
 * ICI 70320
 * Prof. M Diaz-Maldonado
 *
 * source: vector.h
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


typedef struct {
// data
double* begin;
double* avail;
double* limit;
double* array;
// methods
size_t (*size) (void*);
void (*clear) (void*);
void (*push_back) (void*, double);
} vector_t;
#endif
