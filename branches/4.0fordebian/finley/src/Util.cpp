
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************

  Some utility routines

*****************************************************************************/

#include "Finley.h"
#include "Util.h"
#include "esysUtils/index.h"

#include <algorithm> // std::sort
#include <limits>

namespace finley {
namespace util {

/// comparison function for sortValueAndIndex
bool ValueAndIndexCompare(const std::pair<int,int> &i, const std::pair<int, int> &j)
{
    // to ensure we have a strict ordering as required by std
    if (i.first == j.first)
        return i.second < j.second;
    return i.first < j.first;
}

/// orders a ValueAndIndexList by value
void sortValueAndIndex(ValueAndIndexList& array)
{
    std::sort(array.begin(), array.end(), ValueAndIndexCompare);
}

/// gathers values into vector out from vector in using index:
///   out(1:numData, 1:len) := in(1:numData, index(1:len))
void gather(int len, const int* index, int numData, const double* in, double* out)
{
    for (int s=0; s<len; s++) {
        for (int i=0; i<numData; i++) {
            out[INDEX2(i,s,numData)] = in[INDEX2(i,index[s],numData)];
        }
    }
}

/// adds a vector in into out using an index:
///   out(1:numData,index[p])+=in(1:numData,p) where
///   p={k=1...len, index[k]<upperBound}
void addScatter(const int len, const int* index, const int numData, const double* in, double* out, const int upperBound)
{
    for (int s=0; s<len; s++) {
        for (int i=0; i<numData; i++) {
            if (index[s] < upperBound) {
                out[INDEX2(i,index[s],numData)]+=in[INDEX2(i,s,numData)];
            }
        }
    }
}

/// multiplies two matrices: A(1:A1,1:A2) := B(1:A1,1:B2)*C(1:B2,1:A2)
void smallMatMult(int A1, int A2, double* A, int B2,
                  const std::vector<double>& B,
                  const std::vector<double>& C)
{
    for (int i=0; i<A1; i++) {
        for (int j=0; j<A2; j++) {
            double sum=0.;
            for (int s=0; s<B2; s++)
                sum+=B[INDEX2(i,s,A1)]*C[INDEX2(s,j,B2)];
            A[INDEX2(i,j,A1)]=sum;
        }
    }
}

/// multiplies a set of matrices with a single matrix:
///   A(1:A1,1:A2,i)=B(1:A1,1:B2,i)*C(1:B2,1:A2) for i=1,len
void smallMatSetMult1(int len, int A1, int A2, double* A, int B2,
                      const std::vector<double>& B,
                      const std::vector<double>& C)
{
    for (int q=0; q<len; q++) {
        for (int i=0; i<A1; i++) {
            for (int j=0; j<A2; j++) {
                double sum=0.;
                for (int s=0; s<B2; s++)
                    sum+=B[INDEX3(i,s,q,A1,B2)]*C[INDEX2(s,j,B2)];
                A[INDEX3(i,j,q,A1,A2)]=sum;
            }
        }
    }
}

/// inverts the set of dim x dim matrices A(:,:,1:len) with dim=1,2,3
/// the inverse and determinant are returned.
void invertSmallMat(int len, int dim, const double* A, double *invA, double* det)
{
    switch(dim) {
        case 1:
            for (int q=0; q<len; q++) {
                const double D=A[q];
                if (ABS(D) > 0) {
                    det[q]=D;
                    invA[q]=1./D;
                } else {
                    setError(ZERO_DIVISION_ERROR, "InvertSmallMat: Non-regular matrix");
                    break;
                }
            }
            break;

        case 2:
            for (int q=0; q<len; q++) {
                const double A11=A[INDEX3(0,0,q,2,2)];
                const double A12=A[INDEX3(0,1,q,2,2)];
                const double A21=A[INDEX3(1,0,q,2,2)];
                const double A22=A[INDEX3(1,1,q,2,2)];

                const double D = A11*A22-A12*A21;
                if (ABS(D) > 0) {
                    det[q]=D;
                    invA[INDEX3(0,0,q,2,2)]= A22/D;
                    invA[INDEX3(1,0,q,2,2)]=-A21/D;
                    invA[INDEX3(0,1,q,2,2)]=-A12/D;
                    invA[INDEX3(1,1,q,2,2)]= A11/D;
                } else {
                    setError(ZERO_DIVISION_ERROR, "InvertSmallMat: Non-regular matrix");
                    break;
                }
            }
            break;

        case 3:
            for (int q=0; q<len; q++) {
                const double A11=A[INDEX3(0,0,q,3,3)];
                const double A21=A[INDEX3(1,0,q,3,3)];
                const double A31=A[INDEX3(2,0,q,3,3)];
                const double A12=A[INDEX3(0,1,q,3,3)];
                const double A22=A[INDEX3(1,1,q,3,3)];
                const double A32=A[INDEX3(2,1,q,3,3)];
                const double A13=A[INDEX3(0,2,q,3,3)];
                const double A23=A[INDEX3(1,2,q,3,3)];
                const double A33=A[INDEX3(2,2,q,3,3)];

                const double D = A11*(A22*A33-A23*A32) + A12*(A31*A23-A21*A33) + A13*(A21*A32-A31*A22);
                if (ABS(D) > 0) {
                    det[q]=D;
                    invA[INDEX3(0,0,q,3,3)]=(A22*A33-A23*A32)/D;
                    invA[INDEX3(1,0,q,3,3)]=(A31*A23-A21*A33)/D;
                    invA[INDEX3(2,0,q,3,3)]=(A21*A32-A31*A22)/D;
                    invA[INDEX3(0,1,q,3,3)]=(A13*A32-A12*A33)/D;
                    invA[INDEX3(1,1,q,3,3)]=(A11*A33-A31*A13)/D;
                    invA[INDEX3(2,1,q,3,3)]=(A12*A31-A11*A32)/D;
                    invA[INDEX3(0,2,q,3,3)]=(A12*A23-A13*A22)/D;
                    invA[INDEX3(1,2,q,3,3)]=(A13*A21-A11*A23)/D;
                    invA[INDEX3(2,2,q,3,3)]=(A11*A22-A12*A21)/D;
                } else {
                    setError(ZERO_DIVISION_ERROR, "InvertSmallMat: Non-regular matrix");
                    break;
                }
            }
            break;

        default:
            setError(VALUE_ERROR, "InvertSmallMat: dim must be <=3");
            break;
    }
}

/// returns the normalized vector Normal[dim,len] orthogonal to A(:,0,q) and
/// A(:,1,q) in the case of dim=3, or the vector A(:,0,q) in the case of dim=2
void normalVector(int len, int dim, int dim1, const double* A, double* Normal)
{
    int q;
    double A11,A12,CO_A13,A21,A22,CO_A23,A31,A32,CO_A33,length,invlength;

    switch(dim) {
        case 1:
            for (q=0;q<len;q++) Normal[q]=1;
            break;
        case 2:
            for (q=0;q<len;q++) {
                A11=A[INDEX3(0,0,q,2,dim1)];
                A21=A[INDEX3(1,0,q,2,dim1)];
                length = sqrt(A11*A11+A21*A21);
                if (length <= 0) {
                    setError(ZERO_DIVISION_ERROR, __FILE__ ": area equals zero.");
                    return;
                } else {
                    invlength=1./length;
                    Normal[INDEX2(0,q,2)]=A21*invlength;
                    Normal[INDEX2(1,q,2)]=-A11*invlength;
                }
            }
            break;
        case 3:
            for (q=0;q<len;q++) {
                A11=A[INDEX3(0,0,q,3,dim1)];
                A21=A[INDEX3(1,0,q,3,dim1)];
                A31=A[INDEX3(2,0,q,3,dim1)];
                A12=A[INDEX3(0,1,q,3,dim1)];
                A22=A[INDEX3(1,1,q,3,dim1)];
                A32=A[INDEX3(2,1,q,3,dim1)];
                CO_A13=A21*A32-A31*A22;
                CO_A23=A31*A12-A11*A32;
                CO_A33=A11*A22-A21*A12;
                length=sqrt(CO_A13*CO_A13+CO_A23*CO_A23+CO_A33*CO_A33);
                if (length <= 0) {
                    setError(ZERO_DIVISION_ERROR, __FILE__ ": area equals zero.");
                    return;
                } else {
                    invlength=1./length;
                    Normal[INDEX2(0,q,3)]=CO_A13*invlength;
                    Normal[INDEX2(1,q,3)]=CO_A23*invlength;
                    Normal[INDEX2(2,q,3)]=CO_A33*invlength;
                }
            }
            break;
    }
}

/// calculates the minimum value from a dim X N integer array
int getMinInt(int dim, int N, const int* values)
{
    int out = std::numeric_limits<int>::max();
    if (values && dim*N > 0) {
        out=values[0];
#pragma omp parallel
        {
            int out_local=out;
#pragma omp for
            for (int j=0; j<N; j++) {
                for (int i=0; i<dim; i++)
                    out_local=std::min(out_local, values[INDEX2(i,j,dim)]);
            }
#pragma omp critical
            out=std::min(out_local, out);
        }
    }
    return out;
}

/// calculates the maximum value from a dim X N integer array
int getMaxInt(int dim, int N, const int* values)
{
    int out = std::numeric_limits<int>::min();
    if (values && dim*N > 0) {
        out=values[0];
#pragma omp parallel
        {
            int out_local=out;
#pragma omp for
            for (int j=0; j<N; j++) {
                for (int i=0; i<dim; i++)
                    out_local=std::max(out_local, values[INDEX2(i,j,dim)]);
            }
#pragma omp critical
            out=std::max(out_local, out);
        }
    }
    return out;
}

std::pair<int,int> getMinMaxInt(int dim, int N, const int* values)
{
    int vmin = std::numeric_limits<int>::max();
    int vmax = std::numeric_limits<int>::min();
    if (values && dim*N > 0) {
        vmin = vmax = values[0];
#pragma omp parallel
        {
            int vmin_local=vmin;
            int vmax_local=vmax;
#pragma omp for
            for (int j=0; j<N; j++) {
                for (int i=0; i<dim; i++) {
                    vmin_local=std::min(vmin_local, values[INDEX2(i,j,dim)]);
                    vmax_local=std::max(vmax_local, values[INDEX2(i,j,dim)]);
                }
            }
#pragma omp critical
            {
                vmin=std::min(vmin_local, vmin);
                vmax=std::max(vmax_local, vmax);
            }
        }
    }
    return std::pair<int,int>(vmin,vmax);
}

/// calculates the minimum and maximum value from an integer array of length N
/// disregarding the value 'ignore'
std::pair<int,int> getFlaggedMinMaxInt(int N, const int* values, int ignore)
{
    int vmin = std::numeric_limits<int>::max();
    int vmax = std::numeric_limits<int>::min();
    if (values && N > 0) {
        vmin = vmax = values[0];
#pragma omp parallel
        {
            int vmin_local=vmin;
            int vmax_local=vmax;
#pragma omp for
            for (int i=0; i<N; i++) {
                if (values[i] != ignore) {
                    vmin_local=std::min(vmin_local, values[i]);
                    vmax_local=std::max(vmax_local, values[i]);
                }
            }
#pragma omp critical
            {
                vmin=std::min(vmin_local, vmin);
                vmax=std::max(vmax_local, vmax);
            }
        }
    }
    return std::pair<int,int>(vmin,vmax);
}

/// determines the indices of the positive entries in mask returning the
/// length of index.
std::vector<int> packMask(const std::vector<short>& mask)
{
    std::vector<int> index;
    for (int k=0; k<mask.size(); k++) {
        if (mask[k] >= 0) {
            index.push_back(k);
        }
    }
    return index;
}

void setValuesInUse(const int *values, const int numValues,
                    std::vector<int>& valuesInUse, esysUtils::JMPI& mpiinfo)
{
    int lastFoundValue=INDEX_T_MIN;
    bool allFound=false;

    valuesInUse.clear();

    while (!allFound) {
        // find smallest value bigger than lastFoundValue
        int minFoundValue = INDEX_T_MAX;
#pragma omp parallel
        {
            int local_minFoundValue=minFoundValue;
#pragma omp for
            for (int i=0; i<numValues; i++) {
                const int val=values[i];
                if ((val>lastFoundValue) && (val<local_minFoundValue))
                    local_minFoundValue=val;
            }
#pragma omp critical
            {
                if (local_minFoundValue<minFoundValue)
                    minFoundValue=local_minFoundValue;
            }
        }
#ifdef ESYS_MPI
        int local_minFoundValue=minFoundValue;
        MPI_Allreduce(&local_minFoundValue, &minFoundValue, 1, MPI_INT, MPI_MIN, mpiinfo->comm);
#endif

        // if we found a new value we need to add this to valuesInUse
        if (minFoundValue < INDEX_T_MAX) {
            valuesInUse.push_back(minFoundValue);
            lastFoundValue=minFoundValue;
        } else {
            allFound=true;
        }
    }
}

} // namespace util
} // namespace finley

