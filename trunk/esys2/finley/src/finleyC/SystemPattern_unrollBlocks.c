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
                                           maybelong row_block_size,maybelong col_block_size) {
  Finley_SystemMatrixPattern*out=NULL;
  maybelong *ptr=NULL,*index=NULL;
  maybelong i,k,j,l;
  Finley_ErrorCode=NO_ERROR;
  maybelong block_size=row_block_size*col_block_size;
  maybelong new_n_ptr=(pattern->n_ptr)*row_block_size;
  maybelong new_len=(pattern->len)*block_size;

  ptr=MEMALLOC(new_n_ptr+1,maybelong);
  index=MEMALLOC(new_len,maybelong);


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
          
        #pragma omp for private(i,k) schedule(static)
        for (i=0;i<new_n_ptr;++i) 
            for (k=ptr[i];k<ptr[i+1];++k) index[k]=0;

        #pragma omp for private(i,j,k,l) schedule(static)
        for (i=0;i<pattern->n_ptr;++i) {
           for (l=pattern->ptr[i];l<pattern->ptr[i+1];++l)  {
              for (k=0;k<row_block_size;++k) {
                 for (j=0;j<col_block_size;++j) {
                    index[ptr[i*row_block_size+k]+(l-pattern->ptr[i])*col_block_size+j]=(pattern->index[l]-INDEX_OFFSET)*col_block_size+j;
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
 * Revision 1.2  2005/04/01 05:48:56  jgs
 * *** empty log message ***
 *
 * Revision 1.1.2.1  2005/03/15 07:23:55  gross
 * Finley's interface to the SCSL library can deal with systems of PDEs now. tests shows that the SCSL library cannot deal with problems with more then 200000 unknowns. problem has been reported to SGI.
 *
 *
 *
 */
