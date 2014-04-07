
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


#ifndef __PASO_FUNCTIONS_H__
#define __PASO_FUNCTIONS_H__

#include "Common.h"
#include "esysUtils/Esys_MPI.h"
#include "performance.h"

enum Paso_FunctionType {
    LINEAR_SYSTEM
};

typedef enum Paso_FunctionType Paso_FunctionType;

struct Paso_Function
{
    Paso_FunctionType kind;
    dim_t n;
    Esys_MPIInfo *mpi_info;
    double* b;
    double* tmp;
};

err_t Paso_FunctionDerivative(double* J0w, const double* w, Paso_Function* F,
                              const double* f0, const double* x0,
                              double* setoff, Paso_Performance* pp);

err_t Paso_FunctionCall(Paso_Function* F, double* value, const double* arg,
                        Paso_Performance* pp);

void Paso_Function_free(Paso_Function* F);

#endif // __PASO_FUNCTIONS_H__

