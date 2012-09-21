
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
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


#ifndef INC_PASOREACTIVE
#define INC_PASOREACTIVE

#include "Transport.h"

#define PASO_RT_EXP_LIM_MIN  sqrt(EPSILON) /* exp(h)-1 ~ h + h**2/2 for abs(h) <  PASO_RT_EXP_LIM_MIN */
#define PASO_RT_EXP_LIM_MAX  log(1./sqrt(EPSILON)) /* it is assumed that exp(h) with  h>PASO_RT_EXP_LIM_MAX is not reliable */ 

    
typedef struct Paso_ReactiveSolver {
  double A;
  double dt;
} Paso_ReactiveSolver;


PASO_DLL_API
err_t Paso_ReactiveSolver_solve(Paso_ReactiveSolver* support, Paso_TransportProblem* fctp, double* u, double* u_old,  const double* source, Paso_Options* options, Paso_Performance *pp);

PASO_DLL_API
Paso_ReactiveSolver* Paso_ReactiveSolver_alloc(Paso_TransportProblem* fctp);

PASO_DLL_API
void Paso_ReactiveSolver_free(Paso_ReactiveSolver* in);

PASO_DLL_API
double Paso_ReactiveSolver_getSafeTimeStepSize(Paso_TransportProblem* fctp);

PASO_DLL_API
void Paso_ReactiveSolver_initialize(const double dt, Paso_ReactiveSolver* rsolver, Paso_Options* options);

#endif /* #ifndef INC_PASOREACTIVE */
