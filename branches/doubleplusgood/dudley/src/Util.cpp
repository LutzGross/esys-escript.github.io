
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

/************************************************************************************/

/*   Some utility routines: */

/************************************************************************************/
#include "esysUtils/maths.h"
#include "Util.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "esysUtils/index.h"
#include "esysUtils/mem.h"
#include <limits.h>
#include "string.h"  /* for memcpy*/

/************************************************************************************/

/*   returns true if any of the values in the short array values is not equalt to Zero */

bool_t Dudley_Util_anyNonZeroDouble(dim_t N, double *values)
{
    dim_t q;
    for (q = 0; q < N; ++q)
	if (ABS(values[q]) > 0)
	    return TRUE;
    return FALSE;
}

/************************************************************************************/

/*   gathers double values out from in by index: */

/*        out(1:numData,1:len)=in(1:numData,index(1:len)) */

void Dudley_Util_Gather_double(dim_t len, index_t * index, dim_t numData, double *in, double *out)
{
    dim_t s, i;
    for (s = 0; s < len; s++)
    {
	for (i = 0; i < numData; i++)
	{
	    out[INDEX2(i, s, numData)] = in[INDEX2(i, index[s], numData)];
	}
    }
}

/************************************************************************************/

/*   gathers maybelong values out from in by index: */

/*        out(1:numData,1:len)=in(1:numData,index(1:len)) */

void Dudley_Util_Gather_int(dim_t len, index_t * index, dim_t numData, index_t * in, index_t * out)
{
    dim_t s, i;
    for (s = 0; s < len; s++)
    {
	for (i = 0; i < numData; i++)
	{
	    out[INDEX2(i, s, numData)] = in[INDEX2(i, index[s], numData)];
	}
    }
}

/************************************************************************************/

/*   adds a vector in into out using and index. */

/*        out(1:numData,index[p])+=in(1:numData,p) where p = {k=1...len , index[k]<upperBound}*/

void Dudley_Util_AddScatter(const dim_t len, const index_t * index, const dim_t numData, const double *in, double *out, const index_t upperBound)
{
    dim_t i, s;
    for (s = 0; s < len; s++)
    {
	for (i = 0; i < numData; i++)
	{
	    if (index[s] < upperBound)
	    {
		out[INDEX2(i, index[s], numData)] += in[INDEX2(i, s, numData)];
	    }
	}
    }
}

/*    multiplies two matrices */

/*          A(1:A1,1:A2)=B(1:A1,1:B2)*C(1:B2,1:A2) */

void Dudley_Util_SmallMatMult(dim_t A1, dim_t A2, double *A, dim_t B2, const double *B, const double *C)
{
    dim_t i, j, s;
    register double rtmp;
    for (i = 0; i < A1; i++)
    {
	for (j = 0; j < A2; j++)
	{
	    rtmp = 0;
	    for (s = 0; s < B2; s++)
	    {
		rtmp += B[INDEX2(i, s, A1)] * C[INDEX2(s, j, B2)];
	    }
	    A[INDEX2(i, j, A1)] = rtmp;
	}
    }
}

/*    multiplies a two sets of matries: */

/*        A(1:A1,1:A2,i)=B(1:A1,1:B2,i)*C(1:B2,1:A2,i) i=1,len */

void Dudley_Util_SmallMatSetMult(dim_t len, dim_t A1, dim_t A2, double *A, dim_t B2, const double *B, const double *C)
{
    dim_t q, i, j, s;
    register double rtmp;
    for (q = 0; q < len; q++)
    {
	for (i = 0; i < A1; i++)
	{
	    for (j = 0; j < A2; j++)
	    {
		rtmp = 0;
		for (s = 0; s < B2; s++)
		    rtmp += B[INDEX3(i, s, q, A1, B2)] * C[INDEX3(s, j, q, B2, A2)];
		A[INDEX3(i, j, q, A1, A2)] = rtmp;
	    }
	}
    }
}

/*    multiplies a set of matries with a single matrix: */

/*        A(1:A1,1:A2,i)=B(1:A1,1:B2,i)*C(1:B2,1:A2) i=1,len */

void Dudley_Util_SmallMatSetMult1(dim_t len, dim_t A1, dim_t A2, double *A, dim_t B2, const double *B, const double *C)
{
    dim_t q, i, j, s;
    register double rtmp;
    for (q = 0; q < len; q++)
    {
	for (i = 0; i < A1; i++)
	{
	    for (j = 0; j < A2; j++)
	    {
		rtmp = 0;
		for (s = 0; s < B2; s++)
		    rtmp += B[INDEX3(i, s, q, A1, B2)] * C[INDEX2(s, j, B2)];
		A[INDEX3(i, j, q, A1, A2)] = rtmp;
	    }
	}
    }
}

/*    inverts the set of dim x dim matrices A(:,:,1:len) with dim=1,2,3 */
/*    the determinante is returned. */

void Dudley_Util_InvertSmallMat(dim_t len, dim_t dim, double *A, double *invA, double *det)
{
    dim_t q;
    register double D, A11, A12, A13, A21, A22, A23, A31, A32, A33;

    switch (dim)
    {
    case 1:
	for (q = 0; q < len; q++)
	{
	    D = A[q];
	    if (ABS(D) > 0)
	    {
		det[q] = D;
		D = 1. / D;
		invA[q] = D;
	    }
	    else
	    {
		Dudley_setError(ZERO_DIVISION_ERROR, __FILE__ ": Non-regular matrix");
		return;
	    }
	}
	break;

    case 2:
	for (q = 0; q < len; q++)
	{
	    A11 = A[INDEX3(0, 0, q, 2, 2)];
	    A12 = A[INDEX3(0, 1, q, 2, 2)];
	    A21 = A[INDEX3(1, 0, q, 2, 2)];
	    A22 = A[INDEX3(1, 1, q, 2, 2)];

	    D = A11 * A22 - A12 * A21;
	    if (ABS(D) > 0)
	    {
		det[q] = D;
		D = 1. / D;
		invA[INDEX3(0, 0, q, 2, 2)] = A22 * D;
		invA[INDEX3(1, 0, q, 2, 2)] = -A21 * D;
		invA[INDEX3(0, 1, q, 2, 2)] = -A12 * D;
		invA[INDEX3(1, 1, q, 2, 2)] = A11 * D;
	    }
	    else
	    {
		Dudley_setError(ZERO_DIVISION_ERROR, __FILE__ ": Non-regular matrix");
		return;
	    }
	}
	break;

    case 3:
	for (q = 0; q < len; q++)
	{
	    A11 = A[INDEX3(0, 0, q, 3, 3)];
	    A21 = A[INDEX3(1, 0, q, 3, 3)];
	    A31 = A[INDEX3(2, 0, q, 3, 3)];
	    A12 = A[INDEX3(0, 1, q, 3, 3)];
	    A22 = A[INDEX3(1, 1, q, 3, 3)];
	    A32 = A[INDEX3(2, 1, q, 3, 3)];
	    A13 = A[INDEX3(0, 2, q, 3, 3)];
	    A23 = A[INDEX3(1, 2, q, 3, 3)];
	    A33 = A[INDEX3(2, 2, q, 3, 3)];

	    D = A11 * (A22 * A33 - A23 * A32) + A12 * (A31 * A23 - A21 * A33) + A13 * (A21 * A32 - A31 * A22);
	    if (ABS(D) > 0)
	    {
		det[q] = D;
		D = 1. / D;
		invA[INDEX3(0, 0, q, 3, 3)] = (A22 * A33 - A23 * A32) * D;
		invA[INDEX3(1, 0, q, 3, 3)] = (A31 * A23 - A21 * A33) * D;
		invA[INDEX3(2, 0, q, 3, 3)] = (A21 * A32 - A31 * A22) * D;
		invA[INDEX3(0, 1, q, 3, 3)] = (A13 * A32 - A12 * A33) * D;
		invA[INDEX3(1, 1, q, 3, 3)] = (A11 * A33 - A31 * A13) * D;
		invA[INDEX3(2, 1, q, 3, 3)] = (A12 * A31 - A11 * A32) * D;
		invA[INDEX3(0, 2, q, 3, 3)] = (A12 * A23 - A13 * A22) * D;
		invA[INDEX3(1, 2, q, 3, 3)] = (A13 * A21 - A11 * A23) * D;
		invA[INDEX3(2, 2, q, 3, 3)] = (A11 * A22 - A12 * A21) * D;
	    }
	    else
	    {
		Dudley_setError(ZERO_DIVISION_ERROR, __FILE__ ": Non-regular matrix");
		return;
	    }
	}
	break;

    }
    return;
}

/*    sets the derterminate of a set of dim x dim matrices A(:,:,1:len) with dim=1,2,3 */

void Dudley_Util_DetOfSmallMat(dim_t len, dim_t dim, double *A, double *det)
{
    dim_t q;
    register double A11, A12, A13, A21, A22, A23, A31, A32, A33;

    switch (dim)
    {
    case 1:
	for (q = 0; q < len; q++)
	{
	    det[q] = A[q];
	}
	break;

    case 2:
	for (q = 0; q < len; q++)
	{
	    A11 = A[INDEX3(0, 0, q, 2, 2)];
	    A12 = A[INDEX3(0, 1, q, 2, 2)];
	    A21 = A[INDEX3(1, 0, q, 2, 2)];
	    A22 = A[INDEX3(1, 1, q, 2, 2)];

	    det[q] = A11 * A22 - A12 * A21;
	}
	break;

    case 3:
	for (q = 0; q < len; q++)
	{
	    A11 = A[INDEX3(0, 0, q, 3, 3)];
	    A21 = A[INDEX3(1, 0, q, 3, 3)];
	    A31 = A[INDEX3(2, 0, q, 3, 3)];
	    A12 = A[INDEX3(0, 1, q, 3, 3)];
	    A22 = A[INDEX3(1, 1, q, 3, 3)];
	    A32 = A[INDEX3(2, 1, q, 3, 3)];
	    A13 = A[INDEX3(0, 2, q, 3, 3)];
	    A23 = A[INDEX3(1, 2, q, 3, 3)];
	    A33 = A[INDEX3(2, 2, q, 3, 3)];

	    det[q] = A11 * (A22 * A33 - A23 * A32) + A12 * (A31 * A23 - A21 * A33) + A13 * (A21 * A32 - A31 * A22);
	}
	break;

    }
    return;
}

/*    returns the normalized vector Normal[dim,len] orthogonal to A(:,0,q) and A(:,1,q) in the case of dim=3  */
/*    or the vector A(:,0,q) in the case of dim=2                                             */

void Dudley_NormalVector(dim_t len, dim_t dim, dim_t dim1, double *A, double *Normal)
{
    dim_t q;
    register double A11, A12, CO_A13, A21, A22, CO_A23, A31, A32, CO_A33, length, invlength;

    switch (dim)
    {
    case 1:
	for (q = 0; q < len; q++)
	    Normal[q] = 1;
	break;
    case 2:
	for (q = 0; q < len; q++)
	{
	    A11 = A[INDEX3(0, 0, q, 2, dim1)];
	    A21 = A[INDEX3(1, 0, q, 2, dim1)];
	    length = sqrt(A11 * A11 + A21 * A21);
	    if (!length > 0)
	    {
		Dudley_setError(ZERO_DIVISION_ERROR, __FILE__ ": area equals zero.");
		return;
	    }
	    else
	    {
		invlength = 1. / length;
		Normal[INDEX2(0, q, 2)] = A21 * invlength;
		Normal[INDEX2(1, q, 2)] = -A11 * invlength;
	    }
	}
	break;
    case 3:
	for (q = 0; q < len; q++)
	{
	    A11 = A[INDEX3(0, 0, q, 3, dim1)];
	    A21 = A[INDEX3(1, 0, q, 3, dim1)];
	    A31 = A[INDEX3(2, 0, q, 3, dim1)];
	    A12 = A[INDEX3(0, 1, q, 3, dim1)];
	    A22 = A[INDEX3(1, 1, q, 3, dim1)];
	    A32 = A[INDEX3(2, 1, q, 3, dim1)];
	    CO_A13 = A21 * A32 - A31 * A22;
	    CO_A23 = A31 * A12 - A11 * A32;
	    CO_A33 = A11 * A22 - A21 * A12;
	    length = sqrt(CO_A13 * CO_A13 + CO_A23 * CO_A23 + CO_A33 * CO_A33);
	    if (!length > 0)
	    {
		Dudley_setError(ZERO_DIVISION_ERROR, __FILE__ ": area equals zero.");
		return;
	    }
	    else
	    {
		invlength = 1. / length;
		Normal[INDEX2(0, q, 3)] = CO_A13 * invlength;
		Normal[INDEX2(1, q, 3)] = CO_A23 * invlength;
		Normal[INDEX2(2, q, 3)] = CO_A33 * invlength;
	    }

	}
	break;

    }
    return;
}

/*    return the length of the vector which is orthogonal to the vectors A(:,0,q) and A(:,1,q) in the case of dim=3 */
/*    or the vector A(:,0,q) in the case of dim=2                                                                   */

void Dudley_LengthOfNormalVector(dim_t len, dim_t dim, dim_t dim1, double *A, double *length)
{
    dim_t q;
    double A11, A12, CO_A13, A21, A22, CO_A23, A31, A32, CO_A33;

    switch (dim)
    {
    case 1:
	for (q = 0; q < len; q++)
	    length[q] = 1;
	break;
    case 2:
	for (q = 0; q < len; q++)
	{
	    A11 = A[INDEX3(0, 0, q, 2, dim1)];
	    A21 = A[INDEX3(1, 0, q, 2, dim1)];
	    length[q] = sqrt(A11 * A11 + A21 * A21);
	}
	break;
    case 3:
	for (q = 0; q < len; q++)
	{
	    A11 = A[INDEX3(0, 0, q, 3, dim1)];
	    A21 = A[INDEX3(1, 0, q, 3, dim1)];
	    A31 = A[INDEX3(2, 0, q, 3, dim1)];
	    A12 = A[INDEX3(0, 1, q, 3, dim1)];
	    A22 = A[INDEX3(1, 1, q, 3, dim1)];
	    A32 = A[INDEX3(2, 1, q, 3, dim1)];
	    CO_A13 = A21 * A32 - A31 * A22;
	    CO_A23 = A31 * A12 - A11 * A32;
	    CO_A33 = A11 * A22 - A21 * A12;
	    length[q] = sqrt(CO_A13 * CO_A13 + CO_A23 * CO_A23 + CO_A33 * CO_A33);
	}
	break;

    }
    return;
}

/* inverts the map map of length len */
/* there is no range checking! */
/* at output Map[invMap[i]]=i for i=0:lenInvMap */

void Dudley_Util_InvertMap(dim_t lenInvMap, index_t * invMap, dim_t lenMap, index_t * Map)
{
    dim_t i;
    for (i = 0; i < lenInvMap; i++)
	invMap[i] = 0;
    for (i = 0; i < lenMap; i++)
    {
	if (Map[i] >= 0)
	    invMap[Map[i]] = i;
    }
}

/* orders a Dudley_Util_ValueAndIndex array by value */
/* it is assumed that n is large */

int Dudley_Util_ValueAndIndex_compar(const void *arg1, const void *arg2)
{
    Dudley_Util_ValueAndIndex *e1, *e2;
    e1 = (Dudley_Util_ValueAndIndex *) arg1;
    e2 = (Dudley_Util_ValueAndIndex *) arg2;
    if (e1->value < e2->value)
	return -1;
    if (e1->value > e2->value)
	return 1;
    if (e1->index < e2->index)
	return -1;
    if (e1->index > e2->index)
	return 1;
    return 0;
}

void Dudley_Util_sortValueAndIndex(dim_t n, Dudley_Util_ValueAndIndex * array)
{
    /* OMP : needs parallelization ! */
    qsort(array, n, sizeof(Dudley_Util_ValueAndIndex), Dudley_Util_ValueAndIndex_compar);
}

/************************************************************************************/

/* calculates the minimum value from a dim X N integer array */

index_t Dudley_Util_getMinInt(dim_t dim, dim_t N, index_t * values)
{
    dim_t i, j;
    index_t out, out_local;
    out = INDEX_T_MAX;
    if (values != NULL && dim * N > 0)
    {
	out = values[0];
#pragma omp parallel private(out_local)
	{
	    out_local = out;
#pragma omp for private(i,j) schedule(static)
	    for (j = 0; j < N; j++)
	    {
		for (i = 0; i < dim; i++)
		    out_local = MIN(out_local, values[INDEX2(i, j, dim)]);
	    }
#pragma omp critical
	    out = MIN(out_local, out);
	}
    }
    return out;
}

/* calculates the maximum value from a dim X N integer array */

index_t Dudley_Util_getMaxInt(dim_t dim, dim_t N, index_t * values)
{
    dim_t i, j;
    index_t out, out_local;
    out = -INDEX_T_MAX;
    if (values != NULL && dim * N > 0)
    {
	out = values[0];
#pragma omp parallel private(out_local)
	{
	    out_local = out;
#pragma omp for private(i,j) schedule(static)
	    for (j = 0; j < N; j++)
	    {
		for (i = 0; i < dim; i++)
		{
		    out_local = MAX(out_local, values[INDEX2(i, j, dim)]);

		}
	    }
#pragma omp critical
	    out = MAX(out_local, out);
	}
    }
    return out;
}

/************************************************************************************/

/* calculates the minimum value from a dim X N integer array */

index_t Dudley_Util_getFlaggedMinInt(dim_t dim, dim_t N, index_t * values, index_t ignore)
{
    dim_t i, j;
    index_t out, out_local;
    out = INDEX_T_MAX;
    if (values != NULL && dim * N > 0)
    {
	out = values[0];
#pragma omp parallel private(out_local)
	{
	    out_local = out;
#pragma omp for private(i,j) schedule(static)
	    for (j = 0; j < N; j++)
	    {
		for (i = 0; i < dim; i++)
		    if (values[INDEX2(i, j, dim)] != ignore)
			out_local = MIN(out_local, values[INDEX2(i, j, dim)]);
	    }
#pragma omp critical
	    out = MIN(out_local, out);
	}
    }
    return out;
}

/* calculates the maximum value from a dim X N integer array */

index_t Dudley_Util_getFlaggedMaxInt(dim_t dim, dim_t N, index_t * values, index_t ignore)
{
    dim_t i, j;
    index_t out, out_local;
    out = -INDEX_T_MAX;
    if (values != NULL && dim * N > 0)
    {
	out = values[0];
#pragma omp parallel private(out_local)
	{
	    out_local = out;
#pragma omp for private(i,j) schedule(static)
	    for (j = 0; j < N; j++)
	    {
		for (i = 0; i < dim; i++)
		    if (values[INDEX2(i, j, dim)] != ignore)
			out_local = MAX(out_local, values[INDEX2(i, j, dim)]);
	    }
#pragma omp critical
	    out = MAX(out_local, out);
	}
    }
    return out;
}

/* set the index of the positive entries in mask. The length of index is returned. */

dim_t Dudley_Util_packMask(dim_t N, index_t * mask, index_t * index)
{
    dim_t out, k;
    out = 0;
    /*OMP */
    for (k = 0; k < N; k++)
    {
	if (mask[k] >= 0)
	{
	    index[out] = k;
	    out++;
	}
    }
    return out;
}

/* returns true if array contains value */
bool_t Dudley_Util_isAny(dim_t N, index_t * array, index_t value)
{
    bool_t out = FALSE;
    dim_t i;
#pragma omp parallel for private(i) schedule(static) reduction(||:out)
    for (i = 0; i < N; i++)
	out = out || (array[i] == value);
    return out;
}

/* calculates the cummultative sum in array and returns the total sum */
index_t Dudley_Util_cumsum(dim_t N, index_t * array)
{
    index_t out = 0, tmp;
    dim_t i;
#ifdef _OPENMP
    index_t *partial_sums = NULL, sum;
    partial_sums = new  index_t[omp_get_max_threads()];
#pragma omp parallel private(sum,i,tmp)
    {
	sum = 0;
#pragma omp for schedule(static)
	for (i = 0; i < N; ++i)
	    sum += array[i];
	partial_sums[omp_get_thread_num()] = sum;
#pragma omp barrier
#pragma omp master
	{
	    out = 0;
	    for (i = 0; i < omp_get_max_threads(); ++i)
	    {
		tmp = out;
		out += partial_sums[i];
		partial_sums[i] = tmp;
	    }
	}
#pragma omp barrier
	sum = partial_sums[omp_get_thread_num()];
#pragma omp for schedule(static)
	for (i = 0; i < N; ++i)
	{
	    tmp = sum;
	    sum += array[i];
	    array[i] = tmp;
	}
    }
    delete[] partial_sums;
#else
    for (i = 0; i < N; ++i)
    {
	tmp = out;
	out += array[i];
	array[i] = tmp;
    }
#endif
    return out;
}

void Dudley_Util_setValuesInUse(const index_t * values, const dim_t numValues, dim_t * numValuesInUse,
				index_t ** valuesInUse, Esys_MPIInfo * mpiinfo)
{
    dim_t i;
    index_t lastFoundValue = INDEX_T_MIN, minFoundValue, local_minFoundValue, *newValuesInUse = NULL;
    register index_t itmp;
    bool_t allFound = FALSE;
    dim_t nv = 0;

    while (!allFound)
    {
	/* 
	 *  find smallest value bigger than lastFoundValue 
	 */
	minFoundValue = INDEX_T_MAX;
#pragma omp parallel private(local_minFoundValue)
	{
	    local_minFoundValue = minFoundValue;
#pragma omp for private(i,itmp) schedule(static)
	    for (i = 0; i < numValues; i++)
	    {
		itmp = values[i];
		if ((itmp > lastFoundValue) && (itmp < local_minFoundValue))
		    local_minFoundValue = itmp;
	    }
#pragma omp critical
	    {
		if (local_minFoundValue < minFoundValue)
		    minFoundValue = local_minFoundValue;
	    }

	}
#ifdef ESYS_MPI
	local_minFoundValue = minFoundValue;
	MPI_Allreduce(&local_minFoundValue, &minFoundValue, 1, MPI_INT, MPI_MIN, mpiinfo->comm);
#endif
	/* if we found a new tag we need to add this too the valuesInUseList */

	if (minFoundValue < INDEX_T_MAX)
	{
	    newValuesInUse = MEMALLOC(nv + 1, index_t);
	    if (*valuesInUse != NULL)
	    {
		memcpy(newValuesInUse, *valuesInUse, sizeof(index_t) * nv);
		MEMFREE(*valuesInUse);
	    }
	    newValuesInUse[nv] = minFoundValue;
	    *valuesInUse = newValuesInUse;
	    newValuesInUse = NULL;
	    nv++;
	    lastFoundValue = minFoundValue;
	}
	else
	{
	    allFound = TRUE;
	}
    }
    *numValuesInUse = nv;
}

#ifdef ESYS_MPI
void Dudley_printDoubleArray(FILE * fid, dim_t n, double *array, char *name)
{
    index_t i;

    if (name)
	fprintf(fid, "%s [ ", name);
    else
	fprintf(fid, "[ ");
    for (i = 0; i < (n < 60 ? n : 60); i++)
	fprintf(fid, "%g ", array[i]);
    if (n >= 30)
	fprintf(fid, "... ");
    fprintf(fid, "]\n");
}

void Dudley_printIntArray(FILE * fid, dim_t n, int *array, char *name)
{
    index_t i;

    if (name)
	fprintf(fid, "%s [ ", name);
    else
	fprintf(fid, "[ ");
    for (i = 0; i < (n < 60 ? n : 60); i++)
	fprintf(fid, "%d ", array[i]);
    if (n >= 30)
	fprintf(fid, "... ");
    fprintf(fid, "]\n");
}

void Dudley_printMaskArray(FILE * fid, dim_t n, int *array, char *name)
{
    index_t i;

    if (name)
	fprintf(fid, "%s [ ", name);
    else
	fprintf(fid, "[ ");
    for (i = 0; i < (n < 60 ? n : 60); i++)
	if (array[i] != -1)
	    fprintf(fid, "%3d ", array[i]);
	else
	    fprintf(fid, "  * ");
    if (n >= 30)
	fprintf(fid, "... ");
    fprintf(fid, "]\n");
}
#endif
