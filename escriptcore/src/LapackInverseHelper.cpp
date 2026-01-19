
/*****************************************************************************
*
* Copyright (c) 2009-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "LapackInverseHelper.h"

#ifdef ESYS_HAVE_LAPACK

#ifdef ESYS_MKL_LAPACK
#include <mkl_lapack.h>
#else	// assuming lapacke
#include <lapacke.h>
#endif

#endif


using namespace escript;

#define SUCCESS 0
#define NEEDLAPACK 5
#define ERRFACTORISE 6
#define ERRINVERT 7

LapackInverseHelper::LapackInverseHelper(int N)
{
	piv=0;
	work=0;
	lwork=0;
	this->N=N;
#ifdef ESYS_HAVE_LAPACK
	piv=new int[N];
	int blocksize=64;	// this is arbitrary. For implementations that require work array 
				// maybe we should look into the Lapack ILAENV function
#ifdef ESYS_MKL_LAPACK
	int minus1=-1;
	double dummyd=0;
	int result=0;
	dgetri(&N, &dummyd, &N, &N, &dummyd, &minus1, &result);		// The only param that matters is the second -1
	if (result==0)
	{
		blocksize=static_cast<int>(dummyd);
	}
	// If there is an error, then fail silently.
	// Why: I need this to be threadsafe so I can't throw and If this call fails, 
	//      It will fail again when we try to invert the matrix 
#endif	
	lwork=N*blocksize;
	work=new double[lwork];
#endif
}

LapackInverseHelper::~LapackInverseHelper()
{
	if (piv!=0)
	{
	    delete[] piv;
	}
	if (work!=0)
	{
	    delete[] work;
	}
}

int 
LapackInverseHelper::invert(double* matrix)
{
#ifndef ESYS_HAVE_LAPACK
	return NEEDLAPACK;
#else
#ifdef ESYS_MKL_LAPACK
	int res=0;
	int size=N;
	dgetrf(&N,&N,matrix,&N,piv,&res);
	if (res!=0)
	{
	    return ERRFACTORISE;
	}
	dgetri(&N, matrix, &N, piv, work, &lwork, &res);
	if (res!=0)
	{
	    return ERRINVERT;
	}
#else		// assuming lapacke
	int res=LAPACKE_dgetrf(LAPACK_COL_MAJOR, N,N,matrix,N, piv);
	if (res!=0)
	{
	    return ERRFACTORISE;
	}
	res=LAPACKE_dgetri(LAPACK_COL_MAJOR ,N,matrix,N , piv);
	if (res!=0)
	{
	    return ERRINVERT;
	}
#endif
	return SUCCESS;
#endif
}

