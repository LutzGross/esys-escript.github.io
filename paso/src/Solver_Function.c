
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


#include "Common.h"
#include "Functions.h"
#include "PasoUtil.h"
#include "Solver.h"
/*
 * generate Linear System (mainly for test purposes)
 *
 */
Paso_Function * Paso_Function_LinearSystem_alloc(Paso_SystemMatrix* A, double* b, Paso_Options* options)
{
    Paso_Function * out=NULL;
    Paso_SystemMatrix_setPreconditioner(A,options);
    if (! Esys_noError()) return NULL;
    out=MEMALLOC(1,Paso_Function);
    if (! Esys_checkPtr(out)) {
        out->kind=LINEAR_SYSTEM;
        out->mpi_info=Esys_MPIInfo_getReference(A->mpi_info);
        out->n=Paso_SystemMatrix_getTotalNumRows(A);
        out->more=(void*)Paso_SystemMatrix_getReference(A);
        out->b=b;
        out->tmp=MEMALLOC(out->n, double);
        Esys_checkPtr(out->tmp);
    }
    if (Esys_noError()) {
        return out;
    } else {
        Paso_Function_LinearSystem_free(out);
        return NULL;
    }
}
void Paso_Function_LinearSystem_free(Paso_Function * F) 
{
   if (F!=NULL) {
       Esys_MPIInfo_free(F->mpi_info);
       Paso_SystemMatrix_free((Paso_SystemMatrix*)(F->more));
       MEMFREE(F->tmp);
       MEMFREE(F);
   }
}
/*
 * evaluates value=P*(b-Ax)
 *
 */
err_t Paso_Function_LinearSystem_call(Paso_Function * F,double* value, const double* arg, Paso_Performance *pp)
{
    Paso_SystemMatrix* A=(Paso_SystemMatrix*)(F->more);
    Paso_Copy(F->n,F->tmp,F->b); /* tmp=b */
    Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(PASO_ONE, A, arg,-PASO_ONE, F->tmp); /* tmp=(A*arg-tmp) */
    Paso_SystemMatrix_solvePreconditioner(A,value,F->tmp);  /* value=P*tmp */
    return NO_ERROR;
}
