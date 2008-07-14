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

#include "Common.h"
#include "SolverFCT.h"
#include "PasoUtil.h"

/*
 * these are the calls to allocate, free and call the function dfineing the FCT eqution
 *
 */
Paso_Function * Paso_Function_FCT_alloc(Paso_MPIInfo *mpi_info)
{
    Paso_Function * out=NULL;
    out=MEMALLOC(1,Paso_Function);
    if (! Paso_checkPtr(out)) {
        out->kind=FCT;
    }
    if (Paso_noError()) {
        return out;
    } else {
        Paso_Function_FCT_free(out);
        return NULL;
    }
}
void Paso_Function_FCT_free(Paso_Function * F) 
{
   if (F!=NULL) {
       MEMFREE(F);
   }
}

err_t Paso_Function_FCT_call(Paso_Function * F,double* value, const double* arg)
{

}

