
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

/* Paso: SystemMatrix: interface to HYPRE, a software library of
   high performance preconditioners and solvers. We use the 
   BoomerAMG provided in this library */

/**************************************************************/

/* Copyrights by ACcESS Australia 2011 */
/* Author: Lin Gao, l.gao@uq.edu.au */

/**************************************************************/

#ifndef INC_PASO_BOOMERAMG
#define INC_PASO_BOOMERAMG

#include "SystemMatrix.h"
#include "performance.h"

#ifdef BOOMERAMG
#include <HYPRE_krylov.h> 
#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#endif

typedef struct {
#ifdef BOOMERAMG
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
} Paso_BOOMERAMG_Handler;
#endif
