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
  maybelong i,j,iptr;
  #pragma omp parallel for private(i,iptr,j) schedule(static)
  for (i=0;i< in->pattern->n_ptr;++i) {
     for (iptr=(in->pattern->ptr[i])-PTR_OFFSET;iptr<(in->pattern->ptr[i+1])-PTR_OFFSET; ++iptr) {
         for (j=0;j<in->block_size;j++) array[iptr*(in->block_size)+j]=in->val[iptr*(in->block_size)+j];
     }
  }
}
/*
 * $Log$
 * Revision 1.4  2004/12/15 07:08:33  jgs
 * *** empty log message ***
 *
 *
 *
 */
