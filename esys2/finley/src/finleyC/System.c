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
  dim_t n_norm,i;
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
     if (! Finley_checkPtr(out->val)) {
        Finley_SystemMatrix_setValues(out,DBLE(0));
     }
     if (! Finley_checkPtr(out->normalizer)) {
         #pragma omp parallel for private(i) schedule(static)
         for (i=0;i<n_norm;++i) out->normalizer[i]=0.;
     }
  }
  /* all done: */
  if (Finley_ErrorCode!=NO_ERROR) {
    Finley_SystemMatrix_dealloc(out);
    return NULL;
  } else {
    #ifdef Finley_TRACE
    printf("timing: system matrix %.4e sec\n",Finley_timer()-time0);
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
        MEMFREE(in->normalizer);
        Finley_SystemMatrixPattern_dealloc(in->pattern);
        Finley_SystemMatrix_solve_free(in); 
        MEMFREE(in);
        #ifdef Finley_TRACE
        printf("Finley_SystemMatrix_dealloc: system matrix as been deallocated.\n");
        #endif
     }
   }
}
/*
 * $Log$
 * Revision 1.9  2005/09/01 03:31:36  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-01
 *
 * Revision 1.8  2005/08/23 01:24:29  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-08-23
 *
 * Revision 1.7.2.3  2005/08/24 02:02:18  gross
 * timing output switched off. solver output can be swiched through getSolution(verbose=True) now.
 *
 * Revision 1.7.2.2  2005/08/19 02:44:09  gross
 * stopping criterion modified to cope with badly balanced equations
 *
 * Revision 1.7.2.1  2005/08/15 12:02:53  gross
 * memory leak fixed. it is still not clear if there is no problem anymore
 *
 * Revision 1.7  2005/07/08 04:07:56  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.6  2005/06/30 01:53:56  gross
 * a bug in coloring fixed
 *
 * Revision 1.1.1.1.2.5  2005/03/15 07:23:55  gross
 * Finley's interface to the SCSL library can deal with systems of PDEs now. tests shows that the SCSL library cannot deal with problems with more then 200000 unknowns. problem has been reported to SGI.
 *
 * Revision 1.1.1.1.2.4  2005/03/02 23:35:05  gross
 * reimplementation of the ILU in Finley. block size>1 still needs some testing
 *
 * Revision 1.1.1.1.2.3  2004/11/24 01:37:15  gross
 * some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
 *
 * Revision 1.1.1.1.2.2  2004/11/12 06:58:18  gross
 * a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
 *
 * Revision 1.1.1.1.2.1  2004/10/28 22:59:24  gross
 * finley's RecTest.py is running now: problem in SystemMatrixAdapater fixed
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
