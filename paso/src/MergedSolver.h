
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

/* Paso: AMG preconditioner  (local version)                  */

/**************************************************************/

/* Author: lgao@uq.edu.au, l.gross@uq.edu.au                 */

/**************************************************************/

#ifndef INC_PASO_MERGEDSOLVER
#define INC_PASO_MERGEDSOLVER

#include "Paso.h"
#include "SystemMatrix.h"
#include "Options.h"
#include "esysUtils/Esys_MPI.h"
#include "Paso.h"


typedef struct Paso_MergedSolver
{
    Esys_MPIInfo *mpi_info;
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
