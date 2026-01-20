
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#ifndef __PASO_UTIL_H__
#define __PASO_UTIL_H__

/****************************************************************************/

/*   Some utility routines: */

/****************************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004,2005 */
/*   author: l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"

namespace paso {
  
namespace util {  

/// Applies a sequence of N-1 Givens rotations (c,s) to v of length N which is
/// assumed to be small.
void applyGivensRotations(dim_t N, double* v, const double* c, const double* s);

/// returns the index to the largest entry in lambda
index_t arg_max(dim_t N, dim_t* lambda);

/// this int-comparison function is used by qsort/bsearch in various places
int comparIndex(const void* index1, const void* index2);

/// calculates the cumulative sum in array and returns the total sum
index_t cumsum(dim_t N, index_t* array);

index_t cumsum_maskedTrue(dim_t N, index_t* array, int* mask);

index_t cumsum_maskedFalse(dim_t N, index_t* array, int* mask);

/// returns the maximum value in integer array
index_t iMax(dim_t N, const index_t* array);

/// returns the inner product of global arrays x and y
double innerProduct(dim_t N, const double* x, const double* y,
                    escript::JMPI mpiInfo);

/// returns true if array contains value
bool isAny(dim_t N, const index_t* array, index_t value);

/// returns the global L2 norm of x
double l2(dim_t N, const double* x, escript::JMPI mpiInfo);

/// Performs an update of the form z = a*x+b*y  where y and x are long vectors.
/// If a=0, x is not used; if b=0, y is not used.
void linearCombination(dim_t N, double* z, double a, const double* x, double b,
                       const double* y);

/// returns the global Lsup of x
double lsup(dim_t N, const double* x, escript::JMPI mpiInfo);

/// returns the number of positive values in x
dim_t numPositives(dim_t N, const double* x, escript::JMPI mpiInfo);

/// Performs an update of the form x = a*x+b*y  where y and x are long vectors.
/// If b=0, y is not used.
void update(dim_t N, double a, double* x, double b, const double* y);

/// fills array x with zeroes
void zeroes(dim_t N, double* x);

/// out = in
inline void copy(dim_t N, double* out, const double* in)
{
    linearCombination(N, out, 1., in, 0., in);
}

/// x = a*x
inline void scale(dim_t N, double* x, double a)
{
    update(N, a, x, 0, x);
}

/// x = x+a*y
inline void AXPY(dim_t N, double* x, double a, const double* y)
{
    update(N, 1., x, a, y);
}

/// returns true if both arguments have the same sign, false otherwise
inline bool samesign(double a, double b)
{
    return (a>=0 && b>=0) || (a<=0 && b<=0);
}

} // namespace util
} // namespace paso

#endif // __PASO_UTIL_H__

