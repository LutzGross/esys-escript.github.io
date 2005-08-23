/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrix :                           */
/*  sets the values of the system matrix to a value */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "System.h"


/**************************************************************/

void  Finley_SystemMatrix_setValues(Finley_SystemMatrix* in,double value) {
  dim_t i,j;
  index_t iptr;
  #pragma omp parallel for private(i,iptr,j) schedule(static)
  for (i=0;i< in->pattern->n_ptr;++i) {
     for (iptr=(in->pattern->ptr[i])-PTR_OFFSET;iptr<(in->pattern->ptr[i+1])-PTR_OFFSET;++iptr) {
         for (j=0;j<(in->block_size);++j) in->val[iptr*(in->block_size)+j]=value;
     }
  }
  in->normalizer_is_valid=FALSE;
}
/*
 * $Log$
 * Revision 1.6  2005/08/23 01:24:30  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-08-23
 *
 * Revision 1.5.2.1  2005/08/19 02:44:09  gross
 * stopping criterion modified to cope with badly balanced equations
 *
 * Revision 1.5  2005/07/08 04:07:59  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.4  2004/12/15 07:08:34  jgs
 * *** empty log message ***
 * Revision 1.1.1.1.2.2  2005/06/29 02:34:57  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1.2.1  2004/11/12 06:58:19  gross
 * a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
 *
 *
 *
 */
