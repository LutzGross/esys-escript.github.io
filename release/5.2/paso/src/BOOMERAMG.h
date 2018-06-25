
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

/* Paso: interface to HYPRE, a software library of high performance
 *       preconditioners and solvers. We use the BoomerAMG part provided
 *       in this library
 */

/****************************************************************************/

/* Author: Lin Gao, l.gao@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_BOOMERAMG_H__
#define __PASO_BOOMERAMG_H__

#include "SystemMatrix.h"

#ifdef ESYS_HAVE_BOOMERAMG
#include <HYPRE_krylov.h>
#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#endif

namespace paso {

struct Preconditioner_BoomerAMG
{
#ifdef ESYS_HAVE_BOOMERAMG
    HYPRE_IJMatrix A;
    HYPRE_ParCSRMatrix parcsr_A;
    HYPRE_IJVector b;
    HYPRE_ParVector par_b;
    HYPRE_IJVector x;
    HYPRE_ParVector par_x;
    HYPRE_Solver solver;
#else
    void* n;
#endif
};

void Preconditioner_BoomerAMG_free(Preconditioner_BoomerAMG* in);
Preconditioner_BoomerAMG* Preconditioner_BoomerAMG_alloc(SystemMatrix_ptr A,
                                                         Options* options);

void Preconditioner_BoomerAMG_solve(SystemMatrix_ptr A,
                                    Preconditioner_BoomerAMG* amg,
                                    double* x, double * b);


} // namespace paso

#endif // __PASO_BOOMERAMG_H__

