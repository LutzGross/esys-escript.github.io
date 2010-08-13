
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

/* Paso: jacobi preconditioner: */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004, 2005 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Solver.h"

/**************************************************************/

/* frees the Jacobi preconditioner */

void Paso_Solver_Jacobi_free(Paso_Solver_Jacobi * in) {
  if (in!=NULL) {
     MEMFREE(in->values);
     MEMFREE(in->pivot);
     MEMFREE(in);
  }
}
Paso_Solver_Jacobi* Paso_Solver_getJacobi(Paso_SystemMatrix * A_p)
{
   Paso_Solver_Jacobi* jac=Paso_Solver_getLocalJacobi(A_p->mainBlock);
   if (Paso_MPIInfo_noError(A_p->mpi_info)) {
         return jac;
   } else {
         return NULL;
   } 
}

/**************************************************************/

/* Jacobi precondioner set up */

Paso_Solver_Jacobi* Paso_Solver_getLocalJacobi(Paso_SparseMatrix * A_p) {
  Paso_Solver_Jacobi* out=NULL;
  dim_t n = A_p->numCols;
  dim_t n_block=A_p->row_block_size;
  dim_t block_size=A_p->block_size;
  
  /* check matrix is square */
  if (A_p->col_block_size !=A_p->row_block_size) {
    Paso_setError(TYPE_ERROR, "Paso_Solver_getJacobi: Jacobi preconditioner square block size.");
    return NULL;
  }
   
  out=MEMALLOC(1,Paso_Solver_Jacobi);
  if (! Paso_checkPtr(out)) {
     out->values = MEMALLOC( ((size_t) n) * ((size_t) block_size),double);
     out->pivot = MEMALLOC(  ((size_t) n) * ((size_t) n_block), index_t);

     
     if (! ( Paso_checkPtr(out->values) || Paso_checkPtr(out->pivot) )) {
         Paso_SparseMatrix_invMain(A_p, out->values, out->pivot ); 
     }
     
  }
  if (Paso_noError()) {
     return out;
  } else {
     Paso_Solver_Jacobi_free(out);
     return NULL;
  }
} 
/**************************************************************/

/* applies Jacobi preconditioner */

/* should be called within a parallel region                                              */
/* barrier synconization should be performed to make sure that the input vector available */

void Paso_Solver_solveJacobi(Paso_SystemMatrix * A_p, Paso_Solver_Jacobi * prec, double * x, double * b) {
   Paso_Solver_solveLocalJacobi( A_p->mainBlock, prec, x,  b) ;
   return;
}

void Paso_Solver_solveLocalJacobi(Paso_SparseMatrix * A_p, Paso_Solver_Jacobi * prec, double * x, double * b) {
     dim_t n = A_p->numCols;
     dim_t n_block=A_p->row_block_size;
     Paso_Solver_applyBlockDiagonalMatrix(n_block,n,prec->values,prec->pivot,x,b);
     return;
}

