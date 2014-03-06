
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


/*************************************************************************************************/

/* Paso: SparseMatrix: adds the absolute values of row entries to an array */

/*************************************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"

void Paso_SparseMatrix_addAbsRow_CSR_OFFSET0(const Paso_SparseMatrix* A, double* array) {
   dim_t ir,icb,irb;
   index_t iptr;
   #pragma omp parallel for private(ir,irb,iptr,icb) schedule(static)
   for (ir=0;ir< A->pattern->numOutput;ir++) {
       for (irb=0;irb< A->row_block_size;irb++) {
	  const dim_t irow=irb+A->row_block_size*ir;
          register double fac=0.;
	  for (iptr=A->pattern->ptr[ir];iptr<A->pattern->ptr[ir+1]; iptr++) {
	      for (icb=0;icb< A->col_block_size;icb++) {
                 const double rtmp=A->val[iptr*A->block_size+irb+A->row_block_size*icb];
                 fac+=ABS(rtmp);
              }
          }
          array[irow]+=fac;
        }
   }
}
void Paso_SparseMatrix_maxAbsRow_CSR_OFFSET0(const Paso_SparseMatrix* A, double* array) {
   dim_t ir,icb,irb;
   index_t iptr;
   #pragma omp parallel for private(ir,irb,iptr,icb) schedule(static)
   for (ir=0;ir< A->pattern->numOutput;ir++) {
       for (irb=0;irb< A->row_block_size;irb++) {
	  const dim_t irow=irb+A->row_block_size*ir;
          register double fac=0.;
	  for (iptr=A->pattern->ptr[ir];iptr<A->pattern->ptr[ir+1]; iptr++) {
	      for (icb=0;icb< A->col_block_size;icb++)  {
                 const double rtmp=A->val[iptr*A->block_size+irb+A->row_block_size*icb];
                 fac=MAX(fac,rtmp);
              }
          }
          array[irow]=MAX(array[irow], fac);
        }
   }
}
