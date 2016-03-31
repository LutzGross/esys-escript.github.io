
/*****************************************************************************
*
* Copyright (c) 2009-2016 by The University of Queensland
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

#include "LapackInverseHelper.h"

#ifdef USE_LAPACK

#ifdef MKL_LAPACK
#include <mkl_lapack.h>
#else	// assuming clapack
extern "C"
{
#include <clapack.h>
}
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
#ifdef USE_LAPACK
	piv=new int[N];
	int blocksize=64;	// this is arbitrary. For implementations that require work array 
				// maybe we should look into the Lapack ILAENV function
#ifdef MKL_LAPACK
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
#ifndef USE_LAPACK
	return NEEDLAPACK;
#else
#ifdef MKL_LAPACK
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
#else		// assuming clapack
	int res=clapack_dgetrf(CblasColMajor, N,N,matrix,N, piv);
	if (res!=0)
	{
	    return ERRFACTORISE;
	}
	res=clapack_dgetri(CblasColMajor ,N,matrix,N , piv);
	if (res!=0)
	{
	    return ERRINVERT;
	}
#endif
	return SUCCESS;
#endif
}

