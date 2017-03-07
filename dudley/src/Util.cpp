
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include "Util.h"

#include <escript/index.h>

#include <algorithm> // std::sort

namespace dudley {
namespace util {

using escript::DataTypes::real_t;
using escript::DataTypes::cplx_t;

/// comparison function for sortValueAndIndex
bool ValueAndIndexCompare(const std::pair<int,int> &i, const std::pair<int, int> &j)
{
    // to ensure we have a strict ordering as required by std
    if (i.first == j.first)
        return i.second < j.second;
    return i.first < j.first;
}

void sortValueAndIndex(ValueAndIndexList& array)
{
    std::sort(array.begin(), array.end(), ValueAndIndexCompare);
}

void gather(int len, const index_t* index, int numData, const double* in,
            double* out)
{
    for (int s = 0; s < len; s++) {
        for (int i = 0; i < numData; i++) {
            out[INDEX2(i, s, numData)] = in[INDEX2(i, index[s], numData)];
        }
    }
}

template<typename Scalar>
void addScatter(int len, const index_t* index, int numData,
                const Scalar* in, Scalar* out, index_t upperBound)
{
    for (int s = 0; s < len; s++) {
        for (int i = 0; i < numData; i++) {
            if (index[s] < upperBound) {
                out[INDEX2(i, index[s], numData)] += in[INDEX2(i, s, numData)];
            }
        }
    }
}

template
void addScatter<real_t>(int len, const index_t* index, int numData,
                                 const real_t* in, real_t* out, index_t upperBound);
template
void addScatter<cplx_t>(int len, const index_t* index, int numData,
                                 const cplx_t* in, cplx_t* out, index_t upperBound);

void smallMatMult(int A1, int A2, double* A, int B2, const double* B,
                  const double* C)
{
    for (int i = 0; i < A1; i++) {
        for (int j = 0; j < A2; j++) {
            double sum = 0.;
            for (int s = 0; s < B2; s++)
                sum += B[INDEX2(i,s,A1)] * C[INDEX2(s,j,B2)];
            A[INDEX2(i,j,A1)] = sum;
        }
    }
}

void smallMatSetMult1(int len, int A1, int A2, double* A, int B2,
                      const double* B, const double* C)
{
    for (int q = 0; q < len; q++) {
        for (int i = 0; i < A1; i++) {
            for (int j = 0; j < A2; j++) {
                double sum = 0.;
                for (int s = 0; s < B2; s++)
                    sum += B[INDEX3(i,s,q,A1,B2)] * C[INDEX2(s,j,B2)];
                A[INDEX3(i,j,q,A1,A2)] = sum;
            }
        }
    }
}

void normalVector(int len, int dim, int dim1, const double* A, double* Normal)
{
    int q;

    switch (dim) {
        case 1:
            for (q = 0; q < len; q++)
                Normal[q] = 1.;
            break;
        case 2:
            for (q = 0; q < len; q++) {
                const double A11 = A[INDEX3(0,0,q,2,dim1)];
                const double A21 = A[INDEX3(1,0,q,2,dim1)];
                const double length = sqrt(A11*A11+A21*A21);
                if (length <= 0) {
                    throw DudleyException("normalVector: area equals zero.");
                } else {
                    const double invlength = 1./length;
                    Normal[INDEX2(0,q,2)] =  A21*invlength;
                    Normal[INDEX2(1,q,2)] = -A11*invlength;
                }
            }
            break;
        case 3:
            for (q = 0; q < len; q++) {
                const double A11 = A[INDEX3(0,0,q,3,dim1)];
                const double A21 = A[INDEX3(1,0,q,3,dim1)];
                const double A31 = A[INDEX3(2,0,q,3,dim1)];
                const double A12 = A[INDEX3(0,1,q,3,dim1)];
                const double A22 = A[INDEX3(1,1,q,3,dim1)];
                const double A32 = A[INDEX3(2,1,q,3,dim1)];
                const double CO_A13 = A21*A32-A31*A22;
                const double CO_A23 = A31*A12-A11*A32;
                const double CO_A33 = A11*A22-A21*A12;
                const double length = sqrt(CO_A13*CO_A13 + CO_A23*CO_A23
                                           + CO_A33*CO_A33);
                if (length <= 0) {
                    throw DudleyException("normalVector: area equals zero.");
                } else {
                    const double invlength = 1./length;
                    Normal[INDEX2(0,q,3)] = CO_A13*invlength;
                    Normal[INDEX2(1,q,3)] = CO_A23*invlength;
                    Normal[INDEX2(2,q,3)] = CO_A33*invlength;
                }
            }
            break;
    }
}

IndexPair getMinMaxInt(int dim, dim_t N, const index_t* values)
{
    index_t vmin = escript::DataTypes::index_t_max();
    index_t vmax = escript::DataTypes::index_t_min();
    if (values && dim*N > 0) {
        vmin = vmax = values[0];
#pragma omp parallel
        {
            index_t vmin_local = vmin;
            index_t vmax_local = vmax;
#pragma omp for
            for (index_t j = 0; j < N; j++) {
                for (int i = 0; i < dim; i++) {
                    vmin_local = std::min(vmin_local, values[INDEX2(i,j,dim)]);
                    vmax_local = std::max(vmax_local, values[INDEX2(i,j,dim)]);
                }
            }
#pragma omp critical
            {
                vmin = std::min(vmin_local, vmin);
                vmax = std::max(vmax_local, vmax);
            }
        }
    }
    return IndexPair(vmin,vmax);
}

IndexPair getFlaggedMinMaxInt(dim_t N, const index_t* values, index_t ignore)
{
    index_t vmin = escript::DataTypes::index_t_max();
    index_t vmax = escript::DataTypes::index_t_min();
    if (values && N > 0) {
        vmin = vmax = values[0];
#pragma omp parallel
        {
            index_t vmin_local = vmin;
            index_t vmax_local = vmax;
#pragma omp for
            for (index_t i = 0; i < N; i++) {
                if (values[i] != ignore) {
                    vmin_local = std::min(vmin_local, values[i]);
                    vmax_local = std::max(vmax_local, values[i]);
                }
            }
#pragma omp critical
            {
                vmin = std::min(vmin_local, vmin);
                vmax = std::max(vmax_local, vmax);
            }
        }
    }
    return IndexPair(vmin,vmax);
}

std::vector<index_t> packMask(const std::vector<short>& mask)
{
    std::vector<index_t> index;
    for (index_t k = 0; k < mask.size(); k++) {
        if (mask[k] >= 0) {
            index.push_back(k);
        }
    }
    return index;
}

void setValuesInUse(const int* values, dim_t numValues,
                    std::vector<int>& valuesInUse, escript::JMPI mpiinfo)
{
    const int MAX_VALUE = std::numeric_limits<int>::max();
    int lastFoundValue = std::numeric_limits<int>::min();
    bool allFound = false;

    valuesInUse.clear();

    while (!allFound) {
        // find smallest value bigger than lastFoundValue
        int minFoundValue = MAX_VALUE;
#pragma omp parallel
        {
            int local_minFoundValue = minFoundValue;
#pragma omp for
            for (index_t i = 0; i < numValues; i++) {
                const int val = values[i];
                if (val > lastFoundValue && val < local_minFoundValue)
                    local_minFoundValue = val;
            }
#pragma omp critical
            {
                if (local_minFoundValue < minFoundValue)
                    minFoundValue = local_minFoundValue;
            }
        }
#ifdef ESYS_MPI
        int local_minFoundValue = minFoundValue;
        MPI_Allreduce(&local_minFoundValue, &minFoundValue, 1, MPI_INT,
                      MPI_MIN, mpiinfo->comm);
#endif

        // if we found a new value we need to add this to valuesInUse
        if (minFoundValue < MAX_VALUE) {
            valuesInUse.push_back(minFoundValue);
            lastFoundValue = minFoundValue;
        } else {
            allFound = true;
        }
    }
}

} // namespace util
} // namespace dudley

