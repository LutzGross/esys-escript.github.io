/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrixPatternPattern */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004, 2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "Common.h"
#include "Util.h"
#include "SystemPattern.h"

/**************************************************************/

/* creates SystemMatrixPattern  */

Finley_SystemMatrixPattern* Finley_SystemMatrixPattern_unrollBlocks(Finley_SystemMatrixPattern* pattern, \
                                           dim_t row_block_size,dim_t col_block_size) {
  Finley_SystemMatrixPattern*out=NULL;
  index_t *ptr=NULL,*index=NULL,iPtr;
  dim_t i,j,k,l;
  Finley_ErrorCode=NO_ERROR;
  dim_t block_size=row_block_size*col_block_size;
  dim_t new_n_ptr=(pattern->n_ptr)*row_block_size;
  dim_t new_len=(pattern->len)*block_size;

  ptr=MEMALLOC(new_n_ptr+1,index_t);
  index=MEMALLOC(new_len,index_t);


  if (! ( Finley_checkPtr(ptr) || Finley_checkPtr(index) ) )  {
     #pragma omp parallel
     {
        #pragma omp for private(i) schedule(static)
        for (i=0;i<new_n_ptr+1;++i) ptr[i]=0;

        #pragma omp master
        ptr[new_n_ptr]=new_len;

        #pragma omp for private(i,k) schedule(static) 
        for (i=0;i<pattern->n_ptr;++i) 
            for (k=0;k<row_block_size;++k) ptr[i*row_block_size+k]=(pattern->ptr[i]-PTR_OFFSET)*block_size+(pattern->ptr[i+1]-pattern->ptr[i])*col_block_size*k;
          
        #pragma omp for private(i,iPtr) schedule(static) 
        for (i=0;i<new_n_ptr;++i) 
            for (iPtr=ptr[i];iPtr<ptr[i+1];++iPtr) index[iPtr]=0;

        #pragma omp for private(i,j,iPtr,k) schedule(static) 
        for (i=0;i<pattern->n_ptr;++i) {
           for (iPtr=pattern->ptr[i];iPtr<pattern->ptr[i+1];++iPtr)  {
              for (k=0;k<row_block_size;++k) {
                 for (j=0;j<col_block_size;++j) {
                    index[ptr[i*row_block_size+k]+(iPtr-pattern->ptr[i])*col_block_size+j]=(pattern->index[iPtr]-INDEX_OFFSET)*col_block_size+j;
                 }
              }
           }
        }
     }
     /* create return value */
     out=Finley_SystemMatrixPattern_alloc(new_n_ptr,ptr,index);
  }  
  if (Finley_ErrorCode!=NO_ERROR) {
     MEMFREE(index);
     MEMFREE(ptr);
  }
  return out;
}
/*
 * $Log$
 * Revision 1.3  2005/07/08 04:07:57  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.2  2005/04/01 05:48:56  jgs
 * *** empty log message ***
 *
 * Revision 1.1.2.4  2005/07/01 07:02:13  gross
 * some bug with OPENMP fixed
 *
 * Revision 1.1.2.3  2005/06/30 01:53:56  gross
 * a bug in coloring fixed
 *
 * Revision 1.1.2.2  2005/06/29 02:34:56  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.2.1  2005/03/15 07:23:55  gross
 * Finley's interface to the SCSL library can deal with systems of PDEs now. tests shows that the SCSL library cannot deal with problems with more then 200000 unknowns. problem has been reported to SGI.
 *
 *
 *
 */
