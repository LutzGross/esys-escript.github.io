
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: inverts the main diagonal of a SparseMatrix: */

/**************************************************************/

/* Copyrights by ACcESS Australia 2010 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "SparseMatrix.h"

Paso_SparseMatrix* Paso_SparseMatrix_getTranspose(Paso_SparseMatrix* P)
{
   
   Paso_Pattern *outpattern=NULL;
   Paso_SparseMatrix *out=NULL;
   
   const dim_t C=P->numCols;
   const dim_t F=P->numRows-C;
   const dim_t n=C+F;
   const dim_t block_size=P->row_block_size;
   dim_t i,j,k=0;
   index_t iptr,jptr;
   Paso_IndexListArray* index_list = Paso_IndexListArray_alloc(C);
   
   
   for (i=0;i<n;++i) {
      for (iptr=P->pattern->ptr[i];iptr<P->pattern->ptr[i+1]; ++iptr) {
	 j=P->pattern->index[iptr];
	 Paso_IndexListArray_insertIndex(index_list,j,i);
	 }
	 }
	 
	 outpattern=Paso_Pattern_fromIndexListArray(0,index_list,0,C+F,0);
	 out=Paso_SparseMatrix_alloc(P->type,outpattern,block_size,block_size,FALSE);
	 
	 
	 if (block_size==1) {
	    for (i=0;i<out->numRows;++i) {
	       for (iptr=out->pattern->ptr[i];iptr<out->pattern->ptr[i+1]; ++iptr) {
		  j=out->pattern->index[iptr];
		  /*This can be replaced by bsearch!!*/
		  for (jptr=P->pattern->ptr[j];jptr<P->pattern->ptr[j+1]; ++jptr) {
		     k=P->pattern->index[jptr];
		     if(k==i) {
			out->val[iptr]=P->val[jptr];
			}
			}
			}
			}
			} else if (block_size==2) {
			   for (i=0;i<out->numRows;++i) {
			      for (iptr=out->pattern->ptr[i];iptr<out->pattern->ptr[i+1]; ++iptr) {
				 j=out->pattern->index[iptr];
				 /*This can be replaced by bsearch!!*/
				 for (jptr=P->pattern->ptr[j];jptr<P->pattern->ptr[j+1]; ++jptr) {
				    k=P->pattern->index[jptr];
				    if(k==i) {
				       out->val[iptr*block_size*block_size]=P->val[jptr*block_size*block_size];
				       out->val[iptr*block_size*block_size+1]=P->val[jptr*block_size*block_size+2];
				       out->val[iptr*block_size*block_size+2]=P->val[jptr*block_size*block_size+1];
				       out->val[iptr*block_size*block_size+3]=P->val[jptr*block_size*block_size+3];
				       }
				       }
				       }
				       }
				       } else if (block_size==3) {
					  for (i=0;i<out->numRows;++i) {
					     for (iptr=out->pattern->ptr[i];iptr<out->pattern->ptr[i+1]; ++iptr) {
						j=out->pattern->index[iptr];
						/*This can be replaced by bsearch!!*/
						for (jptr=P->pattern->ptr[j];jptr<P->pattern->ptr[j+1]; ++jptr) {
						   k=P->pattern->index[jptr];
						   if(k==i) {
						      out->val[iptr*block_size*block_size]=P->val[jptr*block_size*block_size];
						      out->val[iptr*block_size*block_size+1]=P->val[jptr*block_size*block_size+3];
						      out->val[iptr*block_size*block_size+2]=P->val[jptr*block_size*block_size+6];
						      out->val[iptr*block_size*block_size+3]=P->val[jptr*block_size*block_size+1];
						      out->val[iptr*block_size*block_size+4]=P->val[jptr*block_size*block_size+4];
						      out->val[iptr*block_size*block_size+5]=P->val[jptr*block_size*block_size+7];
						      out->val[iptr*block_size*block_size+6]=P->val[jptr*block_size*block_size+2];
						      out->val[iptr*block_size*block_size+7]=P->val[jptr*block_size*block_size+5];
						      out->val[iptr*block_size*block_size+8]=P->val[jptr*block_size*block_size+8];
						      }
						      }
						      }
						      }
						      }
						      
						      /* clean up */
						      Paso_IndexListArray_free(index_list);
						      Paso_Pattern_free(outpattern);
						      return out;
						      }
						      