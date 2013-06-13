
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/****************************************************************************

  Some utility routines

*****************************************************************************/

#ifndef __FINLEY_UTIL_H__
#define __FINLEY_UTIL_H__

#include "Finley.h"
#include "esysUtils/Esys_MPI.h"

#include <escript/Data.h>
#include <vector>

namespace finley {
namespace util {

typedef std::vector< std::pair<int,int> > ValueAndIndexList;

/// sortValueAndIndex is used to sort items by a value.
/// index points to the location of the original item array and can be used
/// to reorder the array
void sortValueAndIndex(ValueAndIndexList& array);

/// returns true if the data object is defined on reduced element types
inline bool hasReducedIntegrationOrder(const escript::Data& in)
{
    const int fs = in.getFunctionSpace().getTypeCode();
    return (fs == FINLEY_REDUCED_ELEMENTS || fs == FINLEY_REDUCED_FACE_ELEMENTS
                || fs == FINLEY_REDUCED_CONTACT_ELEMENTS_1
                || fs == FINLEY_REDUCED_CONTACT_ELEMENTS_2);
}

void gather(int len, const int* index, int numData, const double* in,
            double* out);

void addScatter(int len, const int* index, int numData, const double* in,
                double* out, int upperBound);

void smallMatMult(int A1, int A2, double* A, int B2, const double* B,
                  const double* C);

void smallMatSetMult1(int len, int A1, int A2, double* A, int B2,
                      const double* B, const double* C);

void invertSmallMat(int len, int dim, const double* A, double *invA,
                    double* det);

void normalVector(int len, int dim, int dim1, const double* A, double* Normal);

int getMinInt(int dim, int N, const int* values);

int getMaxInt(int dim, int N, const int* values);

std::pair<int,int> getMinMaxInt(int dim, int N, const int* values);

int getFlaggedMinInt(int dim, int N, const int* values, int ignore);

int getFlaggedMaxInt(int dim, int N, const int* values, int ignore);

int packMask(int N, const bool_t* mask, int* index);

void setValuesInUse(const int *values, const int numValues,
                    std::vector<int>& valuesInUse, Esys_MPIInfo* mpiinfo);

} // namespace util
} // namespace finley

#endif // __FINLEY_UTIL_H__

