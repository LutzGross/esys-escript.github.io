
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


/****************************************************************************

 * Paso: preconditioner set up

 ****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"
#include "PasoUtil.h"
#include "Preconditioner.h"

namespace paso {

void Preconditioner_free(Preconditioner* in)
{
    if (in!=NULL) {
        Preconditioner_Smoother_free(in->jacobi);
        Preconditioner_Smoother_free(in->gs);
        Preconditioner_AMG_Root_free(in->amg);
        Solver_ILU_free(in->ilu);
        Solver_RILU_free(in->rilu);
        delete in;
    }
}

Preconditioner* Preconditioner_alloc(SystemMatrix_ptr A, Options* options)
{
    Preconditioner* prec = new Preconditioner;
    prec->type=UNKNOWN;
    prec->jacobi=NULL;
    prec->gs=NULL;
    prec->amg=NULL;
    prec->rilu=NULL;
    prec->ilu=NULL;
    
    if (options->verbose && options->use_local_preconditioner)
        printf("Paso: Applying preconditioner locally only.\n");

    switch (options->preconditioner) {
        default:
        case PASO_JACOBI:
            if (options->verbose) {
                if (options->sweeps >0 ) {
                    printf("Preconditioner: Jacobi(%d) preconditioner is used.\n",options->sweeps);
                } else {
                    printf("Preconditioner: Jacobi preconditioner is used.\n");
                }
            }
            prec->jacobi=Preconditioner_Smoother_alloc(A, TRUE, options->use_local_preconditioner, options->verbose);
            prec->type=PASO_JACOBI;
            prec->sweeps=options->sweeps;
            break;

        case PASO_GS:
            if (options->verbose)  {
                if (options->sweeps >0) {
                    printf("Preconditioner: Gauss-Seidel(%d) preconditioner is used.\n",options->sweeps);
                } else {
                    printf("Preconditioner: Gauss-Seidel preconditioner is used.\n");
                }
            }
            prec->gs = Preconditioner_Smoother_alloc(A, false, options->use_local_preconditioner, options->verbose);
            prec->type = PASO_GS;
            prec->sweeps = options->sweeps;
            break;

        case PASO_BOOMERAMG:
        case PASO_AMLI:
        case PASO_AMG:
            prec->amg = Preconditioner_AMG_Root_alloc(A, options);
            prec->type = PASO_AMG;
            break;

        case PASO_ILU0:
            if (options->verbose)
                printf("Preconditioner: ILU preconditioner is used.\n");
            prec->ilu = Solver_getILU(A->mainBlock, options->verbose);
            prec->type = PASO_ILU0;
            Esys_MPIInfo_noError(A->mpi_info);
            break;

        case PASO_RILU:
            if (options->verbose)
                printf("Preconditioner: RILU preconditioner is used.\n");
            prec->rilu = Solver_getRILU(A->mainBlock,options->verbose);
            Esys_MPIInfo_noError(A->mpi_info);
            prec->type=PASO_RILU;
            break;

        case PASO_NO_PRECONDITIONER:
            if (options->verbose)
                printf("Preconditioner: no preconditioner is applied.\n");
            prec->type=PASO_NO_PRECONDITIONER;
            break;
    }
    if (!Esys_noError()) {
        Preconditioner_free(prec);
        return NULL;
    }
    return prec;
}

/* Applies the preconditioner. */
/* Has to be called within a parallel region. */
/* Barrier synchronization is performed before the evaluation to make sure that the input vector is available */
void Preconditioner_solve(Preconditioner* prec, SystemMatrix_ptr A,
                          double* x, double* b)
{
    dim_t n=0;

    switch (prec->type) {
        default:
        case PASO_JACOBI:
            Preconditioner_Smoother_solve(A, prec->jacobi, x, b, prec->sweeps, false);
            break;   
        case PASO_GS:
            Preconditioner_Smoother_solve(A, prec->gs, x, b, prec->sweeps, false);
            break;
        case PASO_AMG:
            Preconditioner_AMG_Root_solve(A, prec->amg, x, b);
            break;
        case PASO_ILU0:
            Solver_solveILU(A->mainBlock, prec->ilu, x, b);
            break;
        case PASO_RILU:
            Solver_solveRILU(prec->rilu, x, b);
            break;
        case PASO_NO_PRECONDITIONER:
            n = MIN(A->getTotalNumCols(), A->getTotalNumRows());
            Paso_Copy(n,x,b);
            break;
    }
}

} // namespace paso

