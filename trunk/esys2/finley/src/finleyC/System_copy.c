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
 * Revision 1.2  2004/12/14 05:39:31  jgs
 * *** empty log message ***
 *
 * Revision 1.1.1.1.2.1  2004/11/12 06:58:19  gross
 * a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
