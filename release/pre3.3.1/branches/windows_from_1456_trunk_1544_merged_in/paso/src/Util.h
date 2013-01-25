/* $Id:$ */
/*******************************************************
 *
 *       Copyright 2008 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#ifndef INC_PASO_UTIL
#define INC_PASO_UTIL

#include "Common.h"
#include "Paso_MPI.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void Paso_zeroes(const dim_t n, double* x);
void Paso_Update(const dim_t n, const double a, double* x, const double b, const double* y);
void Paso_LinearCombination(const dim_t n, double*z, const double a,const double* x, const double b, const double* y);
double Paso_InnerProduct(const dim_t n,const double* x, const double* y, Paso_MPIInfo* mpiinfo);
double Paso_l2(const dim_t n, const double* x, Paso_MPIInfo* mpiinfo);
void ApplyGivensRotations(const dim_t n,double* v,const double* c,const double* s);

#endif
