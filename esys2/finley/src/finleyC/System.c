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
  double *val=NULL;
  Finley_SystemMatrixType out_type;
  Finley_ErrorCode=NO_ERROR;
  
  /* check the matrix type */
  switch(type) {
     case CSC:
        out_type=CSC;
        break;
     case CSR:
        out_type=CSR;
        break;
     case CSC_BLK1:
        out_type=CSC;
        if (row_block_size!=1 || col_block_size!=1) {
            Finley_ErrorCode=TYPE_ERROR;
            sprintf(Finley_ErrorMsg,"convertion of matrix pattern for block size one for logical block size > 1 is not implemented yet");
            return NULL;
        }
        break;
     case CSR_BLK1:
        out_type=CSR;
        if (row_block_size!=1 || col_block_size!=1) {
            Finley_ErrorCode=TYPE_ERROR;
            sprintf(Finley_ErrorMsg,"convertion of matrix pattern for block size one for logical block size > 1 is not implemented yet");
            return NULL;
        }
        break;
     case CSC_SYM:
     case CSC_BLK1_SYM:
        out_type=CSC_SYM;
        Finley_ErrorCode=TYPE_ERROR;
        sprintf(Finley_ErrorMsg,"convertion of matrix pattern for symmetric CSC is not implemented yet.");
        return NULL;
     default:
        Finley_ErrorCode=TYPE_ERROR;
        sprintf(Finley_ErrorMsg,"unknown matrix type identifier %d.",type);
        return NULL;
  }
  time0=Finley_timer();
  /*  allocate the return value */

  out=MEMALLOC(1,Finley_SystemMatrix);
  if (! Finley_checkPtr(out)) {

     /* is block size 1 enforced ? */
     if (type==CSC_BLK1 || type==CSR_BLK1 || type==CSR_BLK1_SYM || type==CSC_BLK1_SYM) {
        out->row_block_size=1;
        out->col_block_size=1;
     } else {
        out->row_block_size=row_block_size;
        out->col_block_size=col_block_size;
     }
     if (out_type==CSC || type==CSC_BLK1 ||  type==CSC_SYM || type==CSC_BLK1_SYM ) {
         out->num_cols=pattern->n_index;
         out->num_rows=pattern->n_ptr;
     } else {
         out->num_rows=pattern->n_index;
         out->num_cols=pattern->n_ptr;
     } 

     out->type=out_type;
     out->logical_row_block_size=row_block_size;
     out->logical_col_block_size=col_block_size;
     out->logical_block_size=out->logical_row_block_size*out->logical_block_size;
     out->block_size=out->row_block_size*out->col_block_size;
     out->pattern=Finley_SystemMatrixPattern_reference(pattern);
     out->len=(size_t)(out->pattern->len)*(size_t)(out->block_size);
     out->reference_counter=1;
     out->direct=NULL;  
     out->iterative=NULL;

   
     /* allocate memory for matrix entries */
     val=MEMALLOC(out->len,double);
     if (! Finley_checkPtr(val)) {
        out->val=val;
        Finley_SystemMatrix_setValues(out,DBLE(0));
     }
   }
  /* all done: */

  printf("timing: system matrix %.4e sec\n",Finley_timer()-time0);
  if (Finley_ErrorCode!=NO_ERROR) {
    MEMFREE(val);
    Finley_SystemMatrix_dealloc(out);
    return NULL;
  } else {
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
/*
 * $Log$
 * Revision 1.2  2004/12/14 05:39:30  jgs
 * *** empty log message ***
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
