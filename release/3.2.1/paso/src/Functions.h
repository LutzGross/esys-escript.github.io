
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#ifndef INC_PASO_FUNCTIONS
#define INC_PASO_FUNCTIONS


#include "Common.h"
#include "esysUtils/Esys_MPI.h"
#include "performance.h"

enum Paso_FunctionType {
  LINEAR_SYSTEM,
  FCT
};

typedef enum Paso_FunctionType Paso_FunctionType;

typedef struct Paso_Function {
  Paso_FunctionType kind;
  dim_t n;
  Esys_MPIInfo *mpi_info;
  double *b;
  double *tmp;
  void *more;
} Paso_Function;

err_t Paso_FunctionDerivative(double* J0w, const double* w, Paso_Function* F, const double *f0, const double *x0, double* setoff, Paso_Performance *pp);
err_t Paso_FunctionCall(Paso_Function * F,double* value, const double* arg, Paso_Performance *pp);
void Paso_Function_free(Paso_Function * F);

#endif
