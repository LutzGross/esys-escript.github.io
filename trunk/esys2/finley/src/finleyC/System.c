/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrix */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "System.h"

/**************************************************************/

/* allocates a SystemMatrix of type type using the given matrix pattern 
   if type is UNKOWN CSR is used.
   if CSC or CSC_BLK1 is used pattern has to give the CSC pattern.
   if CSR or CSR_BLK1 is used pattern has to give the CSR pattern.
   Values are initialized by zero.  */

Finley_SystemMatrix* Finley_SystemMatrix_alloc(Finley_SystemMatrixType type,Finley_SystemMatrixPattern *pattern, int row_block_size, int col_block_size) {
  double time0;
  Finley_SystemMatrix*out=NULL;
  Finley_ErrorCode=NO_ERROR;
  time0=Finley_timer();
  out=MEMALLOC(1,Finley_SystemMatrix);
  if (! Finley_checkPtr(out)) {  
     out->pattern=NULL;  
     out->direct=NULL;  
     out->iterative=NULL;
     out->val=NULL;  
     out->reference_counter=1;
     /* check the matrix type */
     switch(type) {
        case CSC:
          out->type=CSC;
          if (row_block_size!=col_block_size || col_block_size>3) {
             out->row_block_size=1;
             out->col_block_size=1;
             out->pattern=Finley_SystemMatrixPattern_unrollBlocks(pattern,col_block_size,row_block_size);
          } else {
             out->pattern=Finley_SystemMatrixPattern_reference(pattern);
             out->row_block_size=row_block_size;
             out->col_block_size=col_block_size;
          }
          break;
        case CSR:
           out->type=CSR;
           if (row_block_size!=col_block_size || col_block_size>3) {
              out->row_block_size=1;
              out->col_block_size=1;
              out->pattern=Finley_SystemMatrixPattern_unrollBlocks(pattern,row_block_size,col_block_size);
           } else { 
              out->pattern=Finley_SystemMatrixPattern_reference(pattern);
              out->row_block_size=row_block_size;
              out->col_block_size=col_block_size;
          }
          break;
        case CSC_BLK1:
          out->type=CSC;
          out->row_block_size=1;
          out->col_block_size=1;
          if (row_block_size==1 && col_block_size==1) {
              out->pattern=Finley_SystemMatrixPattern_reference(pattern);
          } else {
             out->pattern=Finley_SystemMatrixPattern_unrollBlocks(pattern,col_block_size,row_block_size);
          }
          break;
        case CSR_BLK1:
          out->type=CSR;
          out->row_block_size=1;
          out->col_block_size=1;
          if (row_block_size==1 && col_block_size==1) {
              out->pattern=Finley_SystemMatrixPattern_reference(pattern);
          } else {
             out->pattern=Finley_SystemMatrixPattern_unrollBlocks(pattern,row_block_size,col_block_size);
          }
          break;
        case CSC_SYM:
        case CSC_BLK1_SYM:
          out->type=CSC_SYM;
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"convertion of matrix pattern for symmetric CSC is not implemented yet.");
          return NULL;
        default:
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"unknown matrix type identifier %d.",type);
          return NULL;
     }
     if (out->type==CSC || out->type==CSC_SYM ) {
         out->num_rows=out->pattern->n_index;
         out->num_cols=out->pattern->n_ptr;
     } else {
         out->num_rows=out->pattern->n_ptr;
         out->num_cols=out->pattern->n_index;
     } 
     out->logical_row_block_size=row_block_size;
     out->logical_col_block_size=col_block_size;
     out->logical_block_size=out->logical_row_block_size*out->logical_block_size;
     out->row_block_size=out->row_block_size;
     out->col_block_size=out->col_block_size;
     out->block_size=out->row_block_size*out->col_block_size;
     out->len=(size_t)(out->pattern->len)*(size_t)(out->block_size);
     /* allocate memory for matrix entries */
     out->val=MEMALLOC(out->len,double);
     if (! Finley_checkPtr(out->val)) {
        Finley_SystemMatrix_setValues(out,DBLE(0));
     }
  }
  /* all done: */
  if (Finley_ErrorCode!=NO_ERROR) {
    Finley_SystemMatrix_dealloc(out);
    return NULL;
  } else {
    printf("timing: system matrix %.4e sec\n",Finley_timer()-time0);
    #ifdef Finley_TRACE
    printf("Finley_SystemMatrix_alloc: %ld x %ld system matrix has been allocated.\n",(long)out->num_rows,(long)out->num_cols);
    #endif
    return out;
  }
}

/* returns a reference to Finley_SystemMatrix in */

Finley_SystemMatrix* Finley_SystemMatrix_reference(Finley_SystemMatrix* in) {
   if (in!=NULL) ++(in->reference_counter);
   return NULL;
}

/* deallocates a SystemMatrix: */

void Finley_SystemMatrix_dealloc(Finley_SystemMatrix* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<=0) {
        MEMFREE(in->val);
        Finley_SystemMatrixPattern_dealloc(in->pattern);
        Finley_SystemMatrix_solve_free(in); 
        MEMFREE(in);
        #ifdef Finley_TRACE
        printf("Finley_SystemMatrix_dealloc: system matrix as been deallocated.\n");
        #endif
     }
   }
}
