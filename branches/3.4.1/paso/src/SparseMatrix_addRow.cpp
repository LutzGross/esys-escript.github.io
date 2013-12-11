
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/*************************************************************************************************/

/* Paso: SparseMatrix:  adds the row entries to an array */

/*************************************************************************************************/

/* Author: l.gross@uq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"

void Paso_SparseMatrix_addRow_CSR_OFFSET0(Paso_SparseMatrix* A, double* array) {
   dim_t ir,irow,icb,irb;
   index_t iptr;
   register double fac;
   #pragma omp parallel for private(ir,irb,irow,fac,iptr,icb) schedule(static)
   for (ir=0;ir< A->pattern->numOutput;ir++) {
       for (irb=0;irb< A->row_block_size;irb++) {
	  irow=irb+A->row_block_size*ir;
	  fac=0.;
	  for (iptr=A->pattern->ptr[ir];iptr<A->pattern->ptr[ir+1]; iptr++) {
	      for (icb=0;icb< A->col_block_size;icb++) 
                 fac+=A->val[iptr*A->block_size+irb+A->row_block_size*icb];

          }
          array[irow]+=fac;
       }
   }
}
void Paso_SparseMatrix_copyBlockToMainDiagonal(Paso_SparseMatrix * A_p, const double* in)
{
   index_t ir;
   const dim_t n = A_p->pattern->numOutput;
   const dim_t nblk = A_p->block_size;
   const size_t nblk_size=sizeof(double)*nblk  ;
   const index_t* main_ptr=Paso_SparseMatrix_borrowMainDiagonalPointer(A_p);
   #pragma omp parallel for private(ir) schedule(static)
   for (ir=0;ir< n;ir++) {
       memcpy((void *)(&(A_p->val[main_ptr[ir]*nblk])), (void *)( &in[nblk*ir] ),  nblk_size);
   }

}
void Paso_SparseMatrix_copyBlockFromMainDiagonal(Paso_SparseMatrix * A_p, double* out)
{
   index_t ir;
   const dim_t n = A_p->pattern->numOutput;
   const dim_t nblk =A_p->block_size;
   const size_t nblk_size=sizeof(double)*nblk  ;
   const index_t* main_ptr=Paso_SparseMatrix_borrowMainDiagonalPointer(A_p);
   #pragma omp parallel for private(ir) schedule(static)
   for (ir=0;ir< n;ir++) {
       memcpy((void *)( &out[nblk*ir] ), (void *)(&(A_p->val[main_ptr[ir]*nblk])), nblk_size);
   }
}

void Paso_SparseMatrix_copyFromMainDiagonal(Paso_SparseMatrix * A_p, double* out)
{
   index_t ir, ib;
   const dim_t n = A_p->pattern->numOutput;
   const dim_t nblk = A_p->block_size;
   const dim_t blk = MIN(A_p->row_block_size, A_p->col_block_size);
   const index_t* main_ptr=Paso_SparseMatrix_borrowMainDiagonalPointer(A_p);
   #pragma omp parallel for private(ir,ib) schedule(static)
   for (ir=0;ir< n;ir++) {
       for (ib=0;ib<blk; ib++) {
	    out[ir*blk+ib] = A_p->val[main_ptr[ir]*nblk+ib+A_p->row_block_size*ib];
       }
   }

}

void Paso_SparseMatrix_copyToMainDiagonal(Paso_SparseMatrix * A_p, const double* in)
{
   index_t ir, ib;
   const dim_t n = A_p->pattern->numOutput;
   const dim_t nblk = A_p->block_size;
   const dim_t blk = MIN(A_p->row_block_size, A_p->col_block_size);
   const index_t* main_ptr=Paso_SparseMatrix_borrowMainDiagonalPointer(A_p);
   #pragma omp parallel for private(ir,ib) schedule(static)
   for (ir=0;ir< n;ir++) {
       for (ib=0;ib<blk; ib++) {
	     A_p->val[main_ptr[ir]*nblk+ib+A_p->row_block_size*ib] = in[ir*blk+ib];
       }
   }

}

