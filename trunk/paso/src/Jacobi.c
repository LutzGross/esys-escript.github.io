
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
#include "Preconditioner.h"
#include "BlockOps.h"

/**************************************************************/

/* frees the Jacobi preconditioner */

void Paso_Preconditioner_Jacobi_free(Paso_Preconditioner_Jacobi * in) {
  if (in!=NULL) {
     MEMFREE(in->values);
     MEMFREE(in->pivot);
     MEMFREE(in);
  }
}
Paso_Preconditioner_Jacobi* Paso_Preconditioner_Jacobi_alloc(Paso_SystemMatrix * A_p)
{
   Paso_Preconditioner_Jacobi* jac=Paso_Preconditioner_LocalJacobi_alloc(A_p->mainBlock);
   if (Paso_MPIInfo_noError(A_p->mpi_info)) {
         return jac;
   } else {
         return NULL;
   } 
}

/**************************************************************/

/* Jacobi precondioner set up */

Paso_Preconditioner_Jacobi* Paso_Preconditioner_LocalJacobi_alloc(Paso_SparseMatrix * A_p) {
  Paso_Preconditioner_Jacobi* out=NULL;
  dim_t n = A_p->numCols;
  dim_t n_block=A_p->row_block_size;
  dim_t block_size=A_p->block_size;
  
  /* check matrix is square */
  if (A_p->col_block_size !=A_p->row_block_size) {
    Paso_setError(TYPE_ERROR, "Paso_Solver_getJacobi: Jacobi preconditioner square block size.");
    return NULL;
  }
   
  out=MEMALLOC(1,Paso_Preconditioner_Jacobi);
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
     Paso_Preconditioner_Jacobi_free(out);
     return NULL;
  }
} 
/**************************************************************/

/* applies Jacobi preconditioner */

/* should be called within a parallel region                                              */
/* barrier synconization should be performed to make sure that the input vector available */

void Paso_Preconditioner_Jacobi_solve(Paso_SystemMatrix * A_p, Paso_Preconditioner_Jacobi * prec, double * x, double * b) {
   Paso_Preconditioner_LocalJacobi_solve( A_p->mainBlock, prec, x,  b) ;
   return;
}

void Paso_Preconditioner_LocalJacobi_solve(Paso_SparseMatrix * A_p, Paso_Preconditioner_Jacobi * prec, double * x, double * b) {
     dim_t n = A_p->numCols;
     dim_t n_block=A_p->row_block_size;
     Paso_BlockOps_allMV(n_block,n,prec->values,prec->pivot,x,b);
     return;
}

