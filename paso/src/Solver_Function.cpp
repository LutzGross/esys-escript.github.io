
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


#include "Common.h"
#include "Functions.h"
#include "PasoUtil.h"
#include "Solver.h"

/*
 * generate Linear System (mainly for test purposes)
 *
 */
Paso_Function_LinearSystem* Paso_Function_LinearSystem_alloc(
        paso::SystemMatrix_ptr A, double* b, paso::Options* options)
{
    Paso_Function_LinearSystem* out=NULL;
    A->setPreconditioner(options);
    if (!Esys_noError()) return NULL;
    out=new Paso_Function_LinearSystem();
    out->kind = LINEAR_SYSTEM;
    out->mpi_info = Esys_MPIInfo_getReference(A->mpi_info);
    out->n = A->getTotalNumRows();
    out->mat = A;
    out->b = b;
    out->tmp = new double[out->n];
    if (Esys_noError()) {
        return out;
    } else {
        Paso_Function_LinearSystem_free(out);
        return NULL;
    }
}

void Paso_Function_LinearSystem_free(Paso_Function_LinearSystem* F) 
{
   if (F!=NULL) {
       Esys_MPIInfo_free(F->mpi_info);
       delete[] F->tmp;
       delete F;
   }
}
/*
 * evaluates value=P*(b-Ax)
 *
 */
err_t Paso_Function_LinearSystem_call(Paso_Function_LinearSystem* F,
                                      double* value, const double* arg,
                                      Paso_Performance *pp)
{
    Paso_Copy(F->n, F->tmp, F->b); /* tmp=b */
    paso::SystemMatrix_MatrixVector_CSR_OFFSET0(PASO_ONE, F->mat, arg, -PASO_ONE, F->tmp); /* tmp=(A*arg-tmp) */
    F->mat->solvePreconditioner(value, F->tmp);  /* value=P*tmp */
    return NO_ERROR;
}

