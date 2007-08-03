/* $Id:$ */

/*
********************************************************************************
*               Copyright 2006, 2007 by ACcESS MNRF                            *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**************************************************************/

/* Paso: SparseMatrix */

/**************************************************************/

/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"
#include "TRILINOS.h"

/**************************************************************/

/* allocates a SparseMatrix of type type using the given matrix pattern 
   if type is UNKOWN CSR is used.
   if CSC or CSC_BLK1 is used pattern has to give the CSC pattern.
   if CSR or CSR_BLK1 is used pattern has to give the CSR pattern.
   Values are initialized by zero.  */

Paso_SparseMatrix* Paso_SparseMatrix_alloc(Paso_SparseMatrixType type,Paso_Pattern *pattern, int row_block_size, int col_block_size) {

  double time0;
  Paso_SparseMatrix*out=NULL;
  dim_t n_norm,i;

  Paso_resetError();
  time0=Paso_timer();
  out=MEMALLOC(1,Paso_SparseMatrix);
  if (! Paso_checkPtr(out)) {  
     out->pattern=NULL;  
     out->val=NULL;  
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
                  out->pattern=Paso_Pattern_unrollBlocks(pattern,PATTERN_FORMAT_OFFSET1,col_block_size,row_block_size);
              } else {
                  out->pattern=Paso_Pattern_unrollBlocks(pattern,PATTERN_FORMAT_DEFAULT,col_block_size,row_block_size);
              }
              out->row_block_size=1;
              out->col_block_size=1;
           } else {
              if ( (type & MATRIX_FORMAT_OFFSET1) ==(pattern->type & PATTERN_FORMAT_OFFSET1)) {
                  out->pattern=Paso_Pattern_getReference(pattern);
              } else {
                  out->pattern=Paso_Pattern_unrollBlocks(pattern,(type & MATRIX_FORMAT_OFFSET1)? PATTERN_FORMAT_OFFSET1:  PATTERN_FORMAT_DEFAULT,1,1);
              }
              out->row_block_size=row_block_size;
              out->col_block_size=col_block_size;
           }
           out->numRows = out->pattern->numInput;
           out->numCols = out->pattern->numOutput;
        }
     } else {
     /* ====== compressed sparse row === */
        if (type & MATRIX_FORMAT_SYM) {
           Paso_setError(TYPE_ERROR,"Generation of matrix pattern for symmetric CSR is not implemented yet.");
           return NULL;
        } else {
           if ((type & MATRIX_FORMAT_BLK1) || row_block_size!=col_block_size || col_block_size>3)  {
              if (type & MATRIX_FORMAT_OFFSET1) {
                  out->pattern=Paso_Pattern_unrollBlocks(pattern,PATTERN_FORMAT_OFFSET1,row_block_size,col_block_size);
              } else {
                  out->pattern=Paso_Pattern_unrollBlocks(pattern,PATTERN_FORMAT_DEFAULT,row_block_size,col_block_size);
              }
              out->row_block_size=1;
              out->col_block_size=1;
           } else {
              if ((type & MATRIX_FORMAT_OFFSET1)==(pattern->type & PATTERN_FORMAT_OFFSET1)) {
                  out->pattern=Paso_Pattern_getReference(pattern);
              } else {
                  out->pattern=Paso_Pattern_unrollBlocks(pattern,(type & MATRIX_FORMAT_OFFSET1)? PATTERN_FORMAT_OFFSET1:  PATTERN_FORMAT_DEFAULT,1,1);
              }
              out->row_block_size=row_block_size;
              out->col_block_size=col_block_size;
           }
           out->numRows = out->pattern->numOutput;
           out->numCols = out->pattern->numInput;
        }
     }
     out->logical_row_block_size=row_block_size;
     out->logical_col_block_size=col_block_size;
     out->logical_block_size=out->logical_row_block_size*out->logical_block_size;
     out->block_size=out->row_block_size*out->col_block_size;
     out->len=(size_t)(out->pattern->len)*(size_t)(out->block_size);

     out->val=MEMALLOC(out->len,double);
     if (! Paso_checkPtr(out->val)) Paso_SparseMatrix_setValues(out,DBLE(0));
  }
  /* all done: */
  if (! Paso_noError()) {
    Paso_SparseMatrix_free(out);
    return NULL;
  } else {
    #ifdef Paso_TRACE
    printf("timing: system matrix %.4e sec\n",Paso_timer()-time0);
    printf("Paso_SparseMatrix_alloc: %ld x %ld system matrix has been allocated.\n",(long)out->numRows,(long)out->numCols);
    #endif
    return out;
  }
}

/* returns a reference to Paso_SparseMatrix in */

Paso_SparseMatrix* Paso_SparseMatrix_getReference(Paso_SparseMatrix* in) {
   if (in!=NULL) ++(in->reference_counter);
   return in;
}

/* deallocates a SparseMatrix: */

void Paso_SparseMatrix_free(Paso_SparseMatrix* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<=0) {
        Paso_Pattern_free(in->pattern);
        MEMFREE(in->val);
        MEMFREE(in);
        #ifdef Paso_TRACE
        printf("Paso_SparseMatrix_free: system matrix as been deallocated.\n");
        #endif
     }
   }
}
