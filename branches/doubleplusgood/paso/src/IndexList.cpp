
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/



/************************************************************************************/

/* Paso: Index List                                           */

/************************************************************************************/
 
/* Author: Lutz Gross, l.gross@uq.edu.au */

/************************************************************************************/

#include "IndexList.h"
#include "Paso.h"

/* inserts row index row into the Paso_IndexList in if it does not exist */

void Paso_IndexList_insertIndex(Paso_IndexList* in, index_t index) {
  dim_t i;
  /* is index in in? */
  for (i=0;i<in->n;i++) {
    if (in->index[i]==index)  return;
  }
  /* index could not be found */
  if (in->n==INDEXLIST_LENGTH) {
     /* if in->index is full check the extension */
     if (in->extension==NULL) {
        in->extension=new Paso_IndexList;
        if (Esys_checkPtr(in->extension)) return;
        in->extension->n=0;
        in->extension->extension=NULL;
     }
     Paso_IndexList_insertIndex(in->extension,index);
  } else {
     /* insert index into in->index*/
     in->index[in->n]=index;
     in->n++;
  }
}

/* counts the number of row indices in the Paso_IndexList in */

dim_t Paso_IndexList_count(Paso_IndexList* in, index_t range_min,index_t range_max) {
  dim_t i;
  dim_t out=0;
  register index_t itmp;
  if (in==NULL) {
     return 0;
  } else {
    for (i=0;i<in->n;i++) {
          itmp=in->index[i];
          if ((itmp>=range_min) && (range_max>itmp)) ++out;
    }
     return out+Paso_IndexList_count(in->extension, range_min,range_max);
  }
}

/* counts the number of row indices in the Paso_IndexList in */

void Paso_IndexList_toArray(Paso_IndexList* in, index_t* array, index_t range_min,index_t range_max, index_t index_offset) {
  dim_t i, ptr;
  register index_t itmp;
  if (in!=NULL) {
    ptr=0;
    for (i=0;i<in->n;i++) {
          itmp=in->index[i];
          if ((itmp>=range_min) && (range_max>itmp)) {
             array[ptr]=itmp+index_offset;
             ptr++;
          }

    }
    Paso_IndexList_toArray(in->extension,&(array[ptr]), range_min, range_max, index_offset);
  }
}

/* deallocates the Paso_IndexList in by recursive calls */

void Paso_IndexList_free(Paso_IndexList* in) {
  if (in!=NULL) {
    Paso_IndexList_free(in->extension);
    delete in;
  }
}

Paso_IndexListArray* Paso_IndexListArray_alloc(const dim_t n)
{
   register dim_t i;
   Paso_IndexListArray* out = new Paso_IndexListArray;
   
   if (! Esys_checkPtr(out)) {
      
       out->n=n;
       out->index_list=new Paso_IndexList[out->n];
   
       if (Esys_checkPtr(out->index_list)) {
	  Paso_IndexListArray_free(out);
       } else {
	  #pragma omp parallel for private(i) schedule(static)
	  for(i=0; i<out->n; ++i) {
	     out->index_list[i].extension=NULL;
	     out->index_list[i].n=0;
	  }
       }
       
   }
   return out;
}

void Paso_IndexListArray_free(Paso_IndexListArray* in)
{
   register dim_t i;
   if (in !=NULL) {

      if (in->index_list!=NULL) {
	 #pragma omp parallel for private(i) schedule(static)
	 for(i=0; i<in->n; ++i) Paso_IndexList_free(in->index_list[i].extension);
      
         delete[] in->index_list;
      }
      delete in;

   }
}

