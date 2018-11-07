
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************/

/* Paso: AMG set-ups                                          */

/****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "Options.h"
#include "Preconditioner.h"

#include <iostream>

namespace paso {

void Preconditioner_AMG_Root_free(Preconditioner_AMG_Root* in)
{
    if (in) {
        Preconditioner_AMG_free(in->amg);
        Preconditioner_LocalAMG_free(in->localamg);
        Preconditioner_Smoother_free(in->amgsubstitute);
        delete in;
    }
}

Preconditioner_AMG_Root* Preconditioner_AMG_Root_alloc(SystemMatrix_ptr A,
                                                       Options* options)
{
    Preconditioner_AMG_Root* prec=new Preconditioner_AMG_Root;
    prec->amg = NULL;
    prec->localamg = NULL;
    prec->amgsubstitute = NULL;
    
    prec->is_local = (A->mpi_info->size == 1) || options->use_local_preconditioner;
    if (prec->is_local) {
        prec->localamg = Preconditioner_LocalAMG_alloc(A->mainBlock, 1, options);
    } else {
        prec->amg = Preconditioner_AMG_alloc(A, 1, options);
    }

    if (options->verbose) {
        if (prec->localamg || prec->amg) {
            std::cout << "Preconditioner_AMG_Root:  Smoother is ";
            if (options->smoother == PASO_JACOBI) {
                std::cout << "Jacobi";
            } else {
                std::cout << "Gauss-Seidel";
            }
            std::cout << " with " << options->pre_sweeps << "/"
                << options->post_sweeps << " pre/post sweeps";
            if (options->interpolation_method == PASO_CLASSIC_INTERPOLATION) {
                std::cout << " and classical interpolation.";
            } else if (options->interpolation_method == PASO_CLASSIC_INTERPOLATION_WITH_FF_COUPLING) {
                std::cout << " and classical interpolation with enforced FF coupling.";
            } else {
                std::cout << " and direct interpolation.";
            }
            std::cout << std::endl;
        } else {
            std::cout << "Preconditioner_AMG_Root:  no coarsening constructed." << std::endl;
        }
    } // verbose?


    if (prec->localamg != NULL) {
        options->num_level=Preconditioner_LocalAMG_getMaxLevel(prec->localamg);
        options->coarse_level_sparsity=Preconditioner_LocalAMG_getCoarseLevelSparsity(prec->localamg);
        options->num_coarse_unknowns=Preconditioner_LocalAMG_getNumCoarseUnknowns(prec->localamg);
    } else if (prec->amg != NULL) {
        options->num_level=Preconditioner_AMG_getMaxLevel(prec->amg);
        options->coarse_level_sparsity=Preconditioner_AMG_getCoarseLevelSparsity(prec->amg);
        options->num_coarse_unknowns=Preconditioner_AMG_getNumCoarseUnknowns(prec->amg);
    } else {
        prec->sweeps=options->sweeps;
        prec->amgsubstitute=Preconditioner_Smoother_alloc(A, (options->smoother == PASO_JACOBI), prec->is_local, options->verbose);
        options->num_level=0;
        if (options->verbose) {
            if (options->smoother == PASO_JACOBI) {
                std::cout << "Preconditioner: Jacobi(" << prec->sweeps
                    << ") preconditioner is used." << std::endl;
            } else {
                std::cout << "Preconditioner: Gauss-Seidel("
                    << prec->sweeps << ") preconditioner is used."
                    << std::endl;
            }
        }
    }
    return prec;
}

/* Applies the preconditioner. */
/* Has to be called within a parallel region. */
/* Barrier synchronization is performed before the evaluation to make sure that the input vector is available */
void Preconditioner_AMG_Root_solve(SystemMatrix_ptr A,
                                   Preconditioner_AMG_Root* prec,
                                   double* x, double* b)
{
    if (prec->localamg != NULL) {
        Preconditioner_LocalAMG_solve(A->mainBlock, prec->localamg, x, b);
    } else if ( prec->amg != NULL) {
        Preconditioner_AMG_solve(A, prec->amg,x,b);
    } else {
        Preconditioner_Smoother_solve(A, prec->amgsubstitute, x, b,
                                      prec->sweeps, false);
    }
}

} // namespace paso

