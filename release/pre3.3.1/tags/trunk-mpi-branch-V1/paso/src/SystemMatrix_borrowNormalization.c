/* $Id$ */


/*
********************************************************************************
*               Copyright   2006 by ACcESS MNRF                                *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**************************************************************/

/* Paso: SystemMatrix                                           */

/*  returns a borrowed reference to the matrix normaliztion vector */

/*  if the vector is in valid a new vector is calculate as the inverse */
/*  of the sum of the absolute value in each row/column                */

/**************************************************************/

/* Copyrights by ACcESS Australia 2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

double* Paso_SystemMatrix_borrowNormalization(Paso_SystemMatrix* A) {
   dim_t ir,irow,ic,icb,icol,irb;
   index_t iptr,icol_failed,irow_failed;
   double fac;
   if (!A->normalizer_is_valid) {
      if ((A->type & MATRIX_FORMAT_CSC) || (A->type & MATRIX_FORMAT_SYM) || (A->type & MATRIX_FORMAT_OFFSET1)) {
        Paso_setError(TYPE_ERROR,"No normalization available for compressed sparse column, symmetric storage scheme or index offset 1.");
      } else {
          irow_failed=-1;
          #pragma omp parallel for private(ir,irb,irow,fac,iptr,icb) schedule(static)
          for (ir=0;ir< A->pattern->myNumOutput;ir++) {
	    for (irb=0;irb< A->row_block_size;irb++) {
	      irow=irb+A->row_block_size*ir;
              fac=0.;
	      for (iptr=A->pattern->ptr[ir];iptr<A->pattern->ptr[ir+1]; iptr++) {
	          for (icb=0;icb< A->col_block_size;icb++) 
                    fac+=ABS(A->val[iptr*A->block_size+irb+A->row_block_size*icb]);
              }
              if (ABS(fac)>0) {
                 A->normalizer[irow]=1./fac;
              } else {
                 A->normalizer[irow]=1.;
                 irow_failed=irow;
              }
            }
          }
          if (irow_failed>=0) 
             Paso_setError(ZERO_DIVISION_ERROR,"There is a row containing zero entries only.");
          A->normalizer_is_valid=TRUE;
      }
      Paso_MPI_noError(A->mpi_info );
   }
   return A->normalizer;
}
