
/* $Id: Pattern.c 1306 2007-09-18 05:51:09Z ksteube $ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/* Paso: Pattern */

/**************************************************************/
 
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Pattern.h"

/**************************************************************/

/* allocates a Pattern  */

Paso_Pattern* Paso_Pattern_alloc(int type, dim_t input_block_size, dim_t output_block_size, dim_t numOutput, dim_t numInput, index_t* ptr, index_t* index) {
  Paso_Pattern*out=NULL;
  index_t index_offset=(type & PATTERN_FORMAT_OFFSET1 ? 1:0);
  index_t loc_min_index,loc_max_index,min_index=index_offset,max_index=index_offset-1;
  dim_t i, sum=0;
  Paso_resetError();

  if (type & PATTERN_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_Pattern_alloc: symmetric matrix pattern is not supported yet");
    return NULL;
  }
  if (ptr!=NULL && index != NULL) {
    #pragma omp parallel private(loc_min_index,loc_max_index,i)
     {
        loc_min_index=index_offset;
        loc_max_index=index_offset-1;
        if (type & PATTERN_FORMAT_OFFSET1) {
           #pragma omp for schedule(static) 
           for (i=0;i<numOutput;++i) {
               if (ptr[i]<ptr[i+1]) {
                 #ifdef USE_QSORTG
                    qsortG(&(index[ptr[i]-1]),(size_t)(ptr[i+1]-ptr[i]),sizeof(index_t),Paso_comparIndex); 
                 #else
                    qsort(&(index[ptr[i]-1]),(size_t)(ptr[i+1]-ptr[i]),sizeof(index_t),Paso_comparIndex); 
                 #endif
                 loc_min_index=MIN(loc_min_index,index[ptr[i]-1]);
                 loc_max_index=MAX(loc_max_index,index[ptr[i+1]-2]);
               }
           }
        } else {
           #pragma omp for schedule(static) 
           for (i=0;i<numOutput;++i) {
               if (ptr[i]<ptr[i+1]) {
                 #ifdef USE_QSORTG
                    qsortG(&(index[ptr[i]]),(size_t)(ptr[i+1]-ptr[i]),sizeof(index_t),Paso_comparIndex); 
                 #else
                    qsort(&(index[ptr[i]]),(size_t)(ptr[i+1]-ptr[i]),sizeof(index_t),Paso_comparIndex); 
                 #endif
                 loc_min_index=MIN(loc_min_index,index[ptr[i]]);
                 loc_max_index=MAX(loc_max_index,index[ptr[i+1]-1]);
               }
           }
        }
        #pragma omp critical
        {
           min_index=MIN(loc_min_index,min_index);
           max_index=MAX(loc_max_index,max_index);
        }
    }
    if ( (min_index<index_offset) || (max_index>=numInput+index_offset) ) {
      Paso_setError(TYPE_ERROR,"Paso_Pattern_alloc: Pattern index out of range.");
      return NULL;
    }
  }
  out=MEMALLOC(1,Paso_Pattern);
  if (! Paso_checkPtr(out)) {
      out->type=type;
      out->reference_counter=1;
      out->numOutput=numOutput;
      out->numInput=numInput;
      out->ptr=ptr;
      out->index=index;
      out->input_block_size=input_block_size;
      out->output_block_size=output_block_size;
      out->block_size=out->input_block_size * out->output_block_size;
      if (out->ptr == NULL) {
          out->len=0;
      } else {
          out->len=out->ptr[out->numOutput] - index_offset;
      }
  }
  #ifdef Paso_TRACE
  printf("Paso_Pattern_alloc: system matrix pattern as been allocated.\n");
  #endif
  return out;
}

/* returns a reference to in */

Paso_Pattern* Paso_Pattern_getReference(Paso_Pattern* in) {
     if (in!=NULL) {
        ++(in->reference_counter);
     }
     return in;
}
  
/* deallocates a Pattern: */

void Paso_Pattern_free(Paso_Pattern* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<=0) {
        MEMFREE(in->ptr);
        MEMFREE(in->index);
        MEMFREE(in);
        #ifdef Paso_TRACE
        printf("Paso_Pattern_free: pattern as been deallocated.\n");
        #endif
     }
   }
}
/* *************************************************************/

/*  some routines which help to get the matrix pattern from elements: */

/*  this routine is used by qsort called in Paso_Pattern_alloc */

int Paso_comparIndex(const void *index1,const void *index2){
   index_t Iindex1,Iindex2;
   Iindex1=*(index_t*)index1;
   Iindex2=*(index_t*)index2;
   if (Iindex1<Iindex2) {
      return -1;
   } else {
      if (Iindex1>Iindex2) {
         return 1;
      } else {
         return 0;
      }
   }
}

bool_t Paso_Pattern_isEmpty(Paso_Pattern* in) {
     if (in != NULL) {
         if ((in->ptr != NULL) && (in->index != NULL)) return FALSE;
     }
     return TRUE;
}
