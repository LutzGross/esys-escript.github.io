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

/* Paso: SystemMatrix */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004,2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Paso_MPI.h"
#include "SystemMatrix.h"
#include "TRILINOS.h"

/**************************************************************/

/* allocates a SystemMatrix of type type using the given matrix pattern 
   if type is UNKOWN CSR is used.
   if CSC or CSC_BLK1 is used pattern has to give the CSC pattern.
   if CSR or CSR_BLK1 is used pattern has to give the CSR pattern.
   Values are initialized by zero.  */

Paso_SystemMatrix* Paso_SystemMatrix_alloc(Paso_SystemMatrixType type,Paso_SystemMatrixPattern *pattern, int row_block_size, int col_block_size) {

  double time0;
  Paso_SystemMatrix*out=NULL;
  Paso_resetError();
  time0=Paso_timer();
  dim_t n_norm,i;
  out=MEMALLOC(1,Paso_SystemMatrix);
  if (! Paso_checkPtr(out)) {  
     out->pattern=NULL;  
     out->solver_package=PASO_PASO;  
     out->solver=NULL;  
     out->val=NULL;  
     out->trilinos_data=NULL;
     out->reference_counter=1;
     out->type=type;
     /* ====== compressed sparse columns === */
     if (type & MATRIX_FORMAT_CSC) {
        if (type & MATRIX_FORMAT_SYM) {
           Paso_setError(TYPE_ERROR,"Generation of matrix pattern for symmetric CSC is not implemented yet.");
           return NULL;
        } else {
           if ((type & MATRIX_FORMAT_BLK1) || row_block_size!=col_block_size || col_block_size>3) {
              if (type & MATRIX_FORMAT_OFFSET1) {
                  out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,PATTERN_FORMAT_OFFSET1,col_block_size,row_block_size);
              } else {
                  out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,PATTERN_FORMAT_DEFAULT,col_block_size,row_block_size);
              }
              out->row_block_size=1;
              out->col_block_size=1;
           } else {
              if ( (type & MATRIX_FORMAT_OFFSET1) ==(pattern->type & PATTERN_FORMAT_OFFSET1)) {
                  out->pattern=Paso_SystemMatrixPattern_reference(pattern);
              } else {
                  out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,(type & MATRIX_FORMAT_OFFSET1)? PATTERN_FORMAT_OFFSET1:  PATTERN_FORMAT_DEFAULT,1,1);
              }
              out->row_block_size=row_block_size;
              out->col_block_size=col_block_size;
           }
        }
        out->row_distribution=Paso_Distribution_getReference(out->pattern->input_distribution);
        out->col_distribution=Paso_Distribution_getReference(out->pattern->output_distribution);
     } else {
     /* ====== compressed sparse row === */
        if (type & MATRIX_FORMAT_SYM) {
           Paso_setError(TYPE_ERROR,"Generation of matrix pattern for symmetric CSR is not implemented yet.");
           return NULL;
        } else {
           if ((type & MATRIX_FORMAT_BLK1) || row_block_size!=col_block_size || col_block_size>3)  {
              if (type & MATRIX_FORMAT_OFFSET1) {
                  out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,PATTERN_FORMAT_OFFSET1,row_block_size,col_block_size);
              } else {
                  out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,PATTERN_FORMAT_DEFAULT,row_block_size,col_block_size);
              }
              out->row_block_size=1;
              out->col_block_size=1;
           } else {
              if ((type & MATRIX_FORMAT_OFFSET1)==(pattern->type & PATTERN_FORMAT_OFFSET1)) {
                  out->pattern=Paso_SystemMatrixPattern_reference(pattern);
              } else {
                  out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,(type & MATRIX_FORMAT_OFFSET1)? PATTERN_FORMAT_OFFSET1:  PATTERN_FORMAT_DEFAULT,1,1);
              }
              out->row_block_size=row_block_size;
              out->col_block_size=col_block_size;
           }
        }
        out->row_distribution=Paso_Distribution_getReference(out->pattern->output_distribution);
        out->col_distribution=Paso_Distribution_getReference(out->pattern->input_distribution);
     }
     out->logical_row_block_size=row_block_size;
     out->logical_col_block_size=col_block_size;
     out->logical_block_size=out->logical_row_block_size*out->logical_block_size;
     out->block_size=out->row_block_size*out->col_block_size;
     out->myLen=(size_t)(out->pattern->myLen)*(size_t)(out->block_size);

     out->numRows = out->row_distribution->myNumComponents;
     out->myNumRows = out->row_distribution->numComponents;
     out->numCols = out->col_distribution->myNumComponents;
     out->myNumCols = out->col_distribution->numComponents;
     out->mpi_info = Paso_MPIInfo_getReference(out->pattern->mpi_info);
     /* allocate memory for matrix entries */
     if (type & MATRIX_FORMAT_TRILINOS_CRS) {
         Paso_TRILINOS_alloc(out->trilinos_data, out->pattern,out->row_block_size,out->col_block_size);
     } else {
         if (type & MATRIX_FORMAT_CSC) {
            n_norm = out->myNumCols * out->col_block_size;
         } else {
            n_norm = out->myNumRows * out->row_block_size;
         }
         out->val=MEMALLOC(out->myLen,double);
         out->normalizer=MEMALLOC(n_norm,double);
         out->normalizer_is_valid=FALSE;
         if (! Paso_checkPtr(out->val)) Paso_SystemMatrix_setValues(out,DBLE(0));
         if (! Paso_checkPtr(out->normalizer)) {
             #pragma omp parallel for private(i) schedule(static)
             for (i=0;i<n_norm;++i) out->normalizer[i]=0.;
         }
     }
  }
  /* all done: */
  if (! Paso_noError()) {
    Paso_SystemMatrix_dealloc(out);
    return NULL;
  } else {
    #ifdef Paso_TRACE
    printf("timing: system matrix %.4e sec\n",Paso_timer()-time0);
    printf("Paso_SystemMatrix_alloc: %ld x %ld system matrix has been allocated.\n",(long)out->numRows,(long)out->numCols);
    #endif
    return out;
  }
}

/* returns a reference to Paso_SystemMatrix in */

Paso_SystemMatrix* Paso_SystemMatrix_reference(Paso_SystemMatrix* in) {
   if (in!=NULL) ++(in->reference_counter);
   return NULL;
}

/* deallocates a SystemMatrix: */

void Paso_SystemMatrix_dealloc(Paso_SystemMatrix* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<=0) {
        Paso_SystemMatrixPattern_dealloc(in->pattern);
        Paso_Distribution_dealloc(in->row_distribution);
        Paso_Distribution_dealloc(in->col_distribution);
        Paso_MPIInfo_dealloc(in->mpi_info);
        MEMFREE(in->val);
        MEMFREE(in->normalizer);
        Paso_TRILINOS_free(in->trilinos_data);
        Paso_solve_free(in); 
        MEMFREE(in);
        #ifdef Paso_TRACE
        printf("Paso_SystemMatrix_dealloc: system matrix as been deallocated.\n");
        #endif
     }
   }
}
/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:38  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.2  2005/09/07 00:59:08  gross
 * some inconsistent renaming fixed to make the linking work.
 *
 * Revision 1.1.2.1  2005/09/05 06:29:47  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
