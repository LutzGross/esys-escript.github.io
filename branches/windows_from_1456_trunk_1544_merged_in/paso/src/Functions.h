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

#ifndef INC_PASO_FUNCTIONS
#define INC_PASO_FUNCTIONS


#include "Common.h"
#include "Paso_MPI.h"
typedef struct Paso_Function {
  int kind;
  dim_t local_n;
  Paso_MPIInfo *mpi_info;
} Paso_Function;

err_t Paso_FunctionDerivative(double* J0w, const double* w, Paso_Function* F, const double *f0, const double *x0, double* setoff);
err_t Paso_FunctionCall(Paso_Function * F,double* value, const double* arg);

#endif
