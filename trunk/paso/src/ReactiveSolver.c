
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: reactive solver (D is a diagonal matrix) 
 *        
 *   - Mv_t=Dv+q   v(0)=u       
 *
 *  to return v(dt)
 *
*/

/**************************************************************/

/* Author: l.gross@uq.edu.au                                */

/**************************************************************/


#include "ReactiveSolver.h"
#include "PasoUtil.h"
#include "Solver.h"

err_t  Paso_ReactiveSolver_solve(Paso_ReactiveSolver* support, Paso_TransportProblem* fctp, double* u, double dt, double* source, Paso_Options* options, Paso_Performance *pp)
{
    return SOLVER_NO_ERROR;
}

Paso_ReactiveSolver* Paso_ReactiveSolver_alloc(Paso_TransportProblem* fctp)
{
    Paso_ReactiveSolver* out=NULL;
    out=MEMALLOC(1,Paso_ReactiveSolver);
    if (Paso_checkPtr(out)) return NULL;
    return out;
}

void Paso_ReactiveSolver_free(Paso_ReactiveSolver* in)
{
   if (in!=NULL) {
      MEMFREE(in);
   }
}
double Paso_ReactiveSolver_getSafeTimeStepSize(Paso_TransportProblem* fctp)
{
   return LARGE_POSITIVE_FLOAT;
}