
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/


/****************************************************************************

 * Paso: preconditioner set up

 ****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Preconditioner.h"
#include "Options.h"
#include "PasoUtil.h"
#include "SystemMatrix.h"

namespace paso {

void Preconditioner_free(Preconditioner* in)
{
    if (in!=NULL) {
        Preconditioner_Smoother_free(in->jacobi);
        Preconditioner_Smoother_free(in->gs);
        Solver_ILU_free(in->ilu);
        Solver_RILU_free(in->rilu);
        delete in;
    }
}

Preconditioner* Preconditioner_alloc(SystemMatrix_ptr<double> A, Options* options)
{
    Preconditioner* prec = new Preconditioner;
    prec->jacobi=NULL;
    prec->gs=NULL;
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
            prec->jacobi=Preconditioner_Smoother_alloc(A, true, options->use_local_preconditioner, options->verbose);
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

        case PASO_ILU0:
            if (options->verbose)
                printf("Preconditioner: ILU preconditioner is used.\n");
            prec->ilu = Solver_getILU(A->mainBlock, options->verbose);
            prec->type = PASO_ILU0;
            break;

        case PASO_RILU:
            if (options->verbose)
                printf("Preconditioner: RILU preconditioner is used.\n");
            prec->rilu = Solver_getRILU(A->mainBlock,options->verbose);
            prec->type=PASO_RILU;
            break;

        case PASO_NO_PRECONDITIONER:
            if (options->verbose)
                printf("Preconditioner: no preconditioner is applied.\n");
            prec->type=PASO_NO_PRECONDITIONER;
            break;
    }
    return prec;
}

/* Applies the preconditioner. */
/* Has to be called within a parallel region. */
/* Barrier synchronization is performed before the evaluation to make sure that the input vector is available */
void Preconditioner_solve(Preconditioner* prec, SystemMatrix_ptr<double> A,
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
        case PASO_ILU0:
            Solver_solveILU(A->mainBlock, prec->ilu, x, b);
            break;
        case PASO_RILU:
            Solver_solveRILU(prec->rilu, x, b);
            break;
        case PASO_NO_PRECONDITIONER:
            n = std::min(A->getTotalNumCols(), A->getTotalNumRows());
            util::copy(n,x,b);
            break;
    }
}

} // namespace paso

