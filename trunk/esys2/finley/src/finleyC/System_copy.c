/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrix :                           */
/*  copies a SystemMatrix values into an output array */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "System.h"


void Finley_SystemMatrix_copy(Finley_SystemMatrix* in,double* array) {
  maybelong i,j,iptr,len_index_list=0;
  maybelong block_size=in->row_block_size*in->col_block_size;
  switch(in->type) {
    case CSR:
        len_index_list=in->num_rows;
        break;
    case CSC:
        len_index_list=in->num_cols;
        break;
    default:
        Finley_ErrorCode = TYPE_ERROR;
        sprintf(Finley_ErrorMsg, "Unknown matrix type.");
        return;
  }
  #pragma omp parallel for private(i,iptr,j) schedule(static)
  for (i=0;i< len_index_list;i++) {
     for (iptr=in->ptr[i]-PTR_OFFSET;iptr<in->ptr[i+1]-PTR_OFFSET; iptr++) {
         for (j=0;j<block_size;j++) array[iptr*block_size+j]=in->val[iptr*block_size+j];
     }
  }
}
/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
