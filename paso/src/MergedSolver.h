
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


/************************************************************************************/

/* Paso: AMG preconditioner  (local version)                  */

/************************************************************************************/

/* Author: lgao@uq.edu.au, l.gross@uq.edu.au                 */

/************************************************************************************/

#ifndef INC_PASO_MERGEDSOLVER
#define INC_PASO_MERGEDSOLVER

#include "Paso.h"
#include "SystemMatrix.h"
#include "Options.h"
#include "esysUtils/Esys_MPI.h"
#include "Paso.h"


typedef struct Paso_MergedSolver
{
    esysUtils::JMPI mpi_info;
    Paso_SparseMatrix *A;

    double* x;
    double* b;
    index_t *counts; 
    index_t *offset;
    index_t reordering;
    index_t refinements;
    index_t verbose;
    index_t sweeps;

} Paso_MergedSolver;

Paso_SparseMatrix* Paso_MergedSolver_mergeSystemMatrix(Paso_SystemMatrix* A);
Paso_MergedSolver* Paso_MergedSolver_alloc(Paso_SystemMatrix *A, Paso_Options* options);
void Paso_MergedSolver_free(Paso_MergedSolver* in);
void Paso_MergedSolver_solve(Paso_MergedSolver* ms, double* local_x, double* local_b) ;

#endif
