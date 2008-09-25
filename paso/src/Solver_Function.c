
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


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
    Paso_Solver_setPreconditioner(A,options);
    if (! Paso_noError()) return NULL;
    out=MEMALLOC(1,Paso_Function);
    if (! Paso_checkPtr(out)) {
        out->kind=LINEAR_SYSTEM;
        out->mpi_info=Paso_MPIInfo_getReference(A->mpi_info);
        out->n=Paso_SystemMatrix_getTotalNumRows(A);
        out->more=(void*)Paso_SystemMatrix_reference(A);
        out->b=b;
        out->tmp=MEMALLOC(out->n, double);
        Paso_checkPtr(out->tmp);
    }
    if (Paso_noError()) {
        return out;
    } else {
        Paso_Function_LinearSystem_free(out);
        return NULL;
    }
}
void Paso_Function_LinearSystem_free(Paso_Function * F) 
{
   if (F!=NULL) {
       Paso_MPIInfo_free(F->mpi_info);
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
    Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(ONE, A, arg,-ONE, F->tmp); /* tmp=(A*arg-tmp) */
    Paso_Solver_solvePreconditioner(A,value,F->tmp);  /* value=P*tmp */
    return NO_ERROR;
}
