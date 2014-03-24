
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


/****************************************************************************/

/* Paso: SparseMatrix:  adds the row entries to an array */

/****************************************************************************/

/* Author: l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"

namespace paso {

void SparseMatrix_addRow_CSR_OFFSET0(const SparseMatrix* A, double* array)
{
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
void SparseMatrix_copyBlockToMainDiagonal(SparseMatrix* A, const double* in)
{
   index_t ir;
   const dim_t n = A->pattern->numOutput;
   const dim_t nblk = A->block_size;
   const size_t nblk_size=sizeof(double)*nblk  ;
   const index_t* main_ptr=SparseMatrix_borrowMainDiagonalPointer(A);
   #pragma omp parallel for private(ir) schedule(static)
   for (ir=0;ir< n;ir++) {
       memcpy((void *)(&(A->val[main_ptr[ir]*nblk])), (void *)( &in[nblk*ir] ),  nblk_size);
   }

}
void SparseMatrix_copyBlockFromMainDiagonal(const SparseMatrix* A, double* out)
{
   index_t ir;
   const dim_t n = A->pattern->numOutput;
   const dim_t nblk =A->block_size;
   const size_t nblk_size=sizeof(double)*nblk;
   const index_t* main_ptr=SparseMatrix_borrowMainDiagonalPointer(A);
   #pragma omp parallel for private(ir) schedule(static)
   for (ir=0;ir< n;ir++) {
       memcpy((void *)( &out[nblk*ir] ), (void *)(&(A->val[main_ptr[ir]*nblk])), nblk_size);
   }
}

void SparseMatrix_copyFromMainDiagonal(const SparseMatrix* A, double* out)
{
   index_t ir, ib;
   const dim_t n = A->pattern->numOutput;
   const dim_t nblk = A->block_size;
   const dim_t blk = MIN(A->row_block_size, A->col_block_size);
   const index_t* main_ptr=SparseMatrix_borrowMainDiagonalPointer(A);
   #pragma omp parallel for private(ir,ib) schedule(static)
   for (ir=0;ir< n;ir++) {
       for (ib=0;ib<blk; ib++) {
	    out[ir*blk+ib] = A->val[main_ptr[ir]*nblk+ib+A->row_block_size*ib];
       }
   }
}

void SparseMatrix_copyToMainDiagonal(SparseMatrix* A, const double* in)
{
   index_t ir, ib;
   const dim_t n = A->pattern->numOutput;
   const dim_t nblk = A->block_size;
   const dim_t blk = MIN(A->row_block_size, A->col_block_size);
   const index_t* main_ptr=SparseMatrix_borrowMainDiagonalPointer(A);
   #pragma omp parallel for private(ir,ib) schedule(static)
   for (ir=0;ir< n;ir++) {
       for (ib=0;ib<blk; ib++) {
	     A->val[main_ptr[ir]*nblk+ib+A->row_block_size*ib] = in[ir*blk+ib];
       }
   }
}

} // namespace paso

