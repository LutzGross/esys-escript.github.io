/* $Id$ */

/**************************************************************/

/* Paso: SystemMatrix */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004,2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

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
             out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,col_block_size,row_block_size);
          } else {
             out->pattern=Paso_SystemMatrixPattern_reference(pattern);
             out->row_block_size=row_block_size;
             out->col_block_size=col_block_size;
          }
          break;
        case CSR:
           out->type=CSR;
           if (row_block_size!=col_block_size || col_block_size>3) {
              out->row_block_size=1;
              out->col_block_size=1;
              out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,row_block_size,col_block_size);
           } else { 
              out->pattern=Paso_SystemMatrixPattern_reference(pattern);
              out->row_block_size=row_block_size;
              out->col_block_size=col_block_size;
          }
          break;
        case CSC_BLK1:
          out->type=CSC;
          out->row_block_size=1;
          out->col_block_size=1;
          if (row_block_size==1 && col_block_size==1) {
              out->pattern=Paso_SystemMatrixPattern_reference(pattern);
          } else {
             out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,col_block_size,row_block_size);
          }
          break;
        case CSR_BLK1:
          out->type=CSR;
          out->row_block_size=1;
          out->col_block_size=1;
          if (row_block_size==1 && col_block_size==1) {
              out->pattern=Paso_SystemMatrixPattern_reference(pattern);
          } else {
             out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,row_block_size,col_block_size);
          }
          break;
        case CSC_SYM:
        case CSC_BLK1_SYM:
          out->type=CSC_SYM;
          Paso_setError(TYPE_ERROR,"convertion of matrix pattern for symmetric CSC is not implemented yet.");
          return NULL;
        default:
          Paso_setError(TYPE_ERROR,"unknown matrix type identifier.");
          return NULL;
     }
     if (out->type==CSC || out->type==CSC_SYM ) {
         out->num_rows=out->pattern->n_index;
         out->num_cols=out->pattern->n_ptr;
         n_norm = out->num_cols * out->col_block_size;

     } else {
         out->num_rows=out->pattern->n_ptr;
         out->num_cols=out->pattern->n_index;
         n_norm = out->num_rows * out->row_block_size;
     } 
     out->logical_row_block_size=row_block_size;
     out->logical_col_block_size=col_block_size;
     out->logical_block_size=out->logical_row_block_size*out->logical_block_size;
     out->block_size=out->row_block_size*out->col_block_size;
     out->len=(size_t)(out->pattern->len)*(size_t)(out->block_size);
     /* allocate memory for matrix entries */
     out->val=MEMALLOC(out->len,double);
     out->normalizer=MEMALLOC(n_norm,double);
     out->normalizer_is_valid=FALSE;
     if (! Paso_checkPtr(out->val)) {
        Paso_SystemMatrix_setValues(out,DBLE(0));
     }
     if (! Paso_checkPtr(out->normalizer)) {
         #pragma omp parallel for private(i) schedule(static)
         for (i=0;i<n_norm;++i) out->normalizer[i]=0.;
     }
  }
  /* all done: */
  if (! Paso_noError()) {
    Paso_SystemMatrix_dealloc(out);
    return NULL;
  } else {
    #ifdef Paso_TRACE
    printf("timing: system matrix %.4e sec\n",Paso_timer()-time0);
    printf("Paso_SystemMatrix_alloc: %ld x %ld system matrix has been allocated.\n",(long)out->num_rows,(long)out->num_cols);
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
        MEMFREE(in->val);
        MEMFREE(in->normalizer);
        Paso_SystemMatrixPattern_dealloc(in->pattern);
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
