
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

/// Some utility routines

#ifndef __FINLEY_UTIL_H__
#define __FINLEY_UTIL_H__

#include "Finley.h"

#include <escript/Data.h>

namespace finley {
namespace util {

typedef std::pair<index_t,index_t> IndexPair;
typedef std::vector<IndexPair> ValueAndIndexList;

/// orders a ValueAndIndexList by value.
void sortValueAndIndex(ValueAndIndexList& array);

/// returns true if the data object is defined on reduced element types
inline bool hasReducedIntegrationOrder(const escript::Data& in)
{
    const int fs = in.getFunctionSpace().getTypeCode();
    return (fs == FINLEY_REDUCED_ELEMENTS || fs == FINLEY_REDUCED_FACE_ELEMENTS
                || fs == FINLEY_REDUCED_CONTACT_ELEMENTS_1
                || fs == FINLEY_REDUCED_CONTACT_ELEMENTS_2);
}

/// gathers values into array `out` from array `in` using `index`:
///   out(1:numData, 1:len) := in(1:numData, index(1:len))
void gather(int len, const index_t* index, int numData, const double* in,
            double* out);

/// adds array `in` into `out` using an `index`:
///   out(1:numData,index[p])+=in(1:numData,p) where
///   p={k=1...len, index[k]<upperBound}
template<typename Scalar>
void addScatter(int len, const index_t* index, int numData,
                const Scalar* in, Scalar* out, index_t upperBound);

/// multiplies two matrices: A(1:A1,1:A2) := B(1:A1,1:B2)*C(1:B2,1:A2)
void smallMatMult(int A1, int A2, double* A, int B2,
                  const std::vector<double>& B,
                  const std::vector<double>& C);

/// multiplies a set of matrices with a single matrix:
///   A(1:A1,1:A2,i)=B(1:A1,1:B2,i)*C(1:B2,1:A2) for i=1,len
template<typename Scalar>
void smallMatSetMult1(int len, int A1, int A2, Scalar* A, int B2,
                      const std::vector<Scalar>& B,
                      const std::vector<double>& C);

void invertSmallMat(int len, int dim, const double* A, double *invA,
                    double* det);

/// returns the normalized vector normal[dim,len] orthogonal to A(:,0,q) and
/// A(:,1,q) in the case of dim=3, or the vector A(:,0,q) in the case of dim=2
void normalVector(int len, int dim, int dim1, const double* A, double* Normal);

index_t getMinInt(int dim, dim_t N, const index_t* values);

index_t getMaxInt(int dim, dim_t N, const index_t* values);

/// calculates the minimum and maximum value from an integer array of length
/// N x dim
IndexPair getMinMaxInt(int dim, dim_t N, const index_t* values);

/// calculates the minimum and maximum value from an integer array of length N
/// disregarding the value `ignore`
IndexPair getFlaggedMinMaxInt(dim_t N, const index_t* values, index_t ignore);

/// extracts the positive entries in `mask` returning a contiguous vector of
/// those entries
std::vector<index_t> packMask(const std::vector<short>& mask);

void setValuesInUse(const int* values, dim_t numValues,
                    std::vector<int>& valuesInUse, escript::JMPI mpiInfo);

} // namespace util
} // namespace finley

#endif // __FINLEY_UTIL_H__

