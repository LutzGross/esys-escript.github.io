
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
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


/* computes the pattern coming from matrix-matrix multiplication
*
**/

Paso_Pattern* Paso_Pattern_multiply(int type, Paso_Pattern* A, Paso_Pattern* B) {
  Paso_Pattern*out=NULL;
  index_t index_offset=(type & PATTERN_FORMAT_OFFSET1 ? 1:0);
  index_t iptrA,iptrB;
  dim_t i,j,k;
  Paso_IndexList* index_list=NULL;

  index_list=TMPMEMALLOC(A->numOutput,Paso_IndexList);
  if (! Paso_checkPtr(index_list)) {
  
      #pragma omp parallel private(i)
      {
        #pragma omp for schedule(static)
        for(i=0;i<A->numOutput;++i) {
             index_list[i].extension=NULL;
             index_list[i].n=0;
        }
      }
  }
  
  for(i = 0; i < A->numOutput; i++) {
     for(iptrA = A->ptr[i]; iptrA < A->ptr[i+1]; ++iptrA) {
      j = A->index[iptrA];
      for(iptrB = B->ptr[j]; iptrB < B->ptr[j+1]; ++iptrB) {
    	k = B->index[iptrB];
        Finley_IndexList_insertIndex(&(index_list[i]),k);
     }
    }
  }
    
  out=Paso_IndexList_createPattern(0, A->numOutput,index_list,0,INDEXLIST_LENGTH,0);

  #ifdef Paso_TRACE
  printf("Paso_Pattern_multipy: new pattern has been allocated.\n");
  #endif

 /* clean up */
   if (index_list!=NULL) {
        #pragma omp parallel for private(i) 
        for(i=0;i<A->numOutput;++i) Paso_IndexList_free(index_list[i].extension);
     }
  TMPMEMFREE(index_list);
  
return out;
}



/*
 * Computes the pattern  of C = A binary operation B for CSR matrices A,B
 *
 * Note: we do not check whether A_ij(op)B_ij=0
 *
 */
Paso_Pattern* Paso_Pattern_binop(int type, Paso_Pattern* A, Paso_Pattern* B) {
  Paso_Pattern*out=NULL;
  index_t index_offset=(type & PATTERN_FORMAT_OFFSET1 ? 1:0);
  index_t iptrA,iptrB,*A_row=NULL,*B_row=NULL;
  dim_t i,j,k;

  Paso_IndexList* index_list=NULL;

 index_list=TMPMEMALLOC(A->numOutput,Paso_IndexList);
   if (! Paso_checkPtr(index_list)) {
  
      #pragma omp parallel private(i)
      {
        #pragma omp for schedule(static)
        for(i=0;i<A->numOutput;++i) {
             index_list[i].extension=NULL;
             index_list[i].n=0;
        }
      }
  }
  for(i = 0; i < B->numOutput; i++){
    iptrA = A->ptr[i],
    iptrB = B->ptr[i];
    
    while (iptrA < A->ptr[i+1] && iptrB < B->ptr[i+1]) {
        j = A->index[iptrA];
        k = B->index[iptrB];
        if (j<k) {
           Finley_IndexList_insertIndex(&(index_list[i]),j);
           iptrA++;
        } else if (j>k) {
            Finley_IndexList_insertIndex(&(index_list[i]),k);
            iptrB++;
        } else if (j==k) {
            Finley_IndexList_insertIndex(&(index_list[i]),j);
            iptrB++;
            iptrA++;
        }
    }
    while(iptrA < A->ptr[i+1]) {
        j = A->index[iptrA];
        Finley_IndexList_insertIndex(&(index_list[i]),j);
        iptrA++;
    }
    while(iptrB < B->ptr[i+1]) {
        k = B->index[iptrB];
        Finley_IndexList_insertIndex(&(index_list[i]),k);
        iptrB++;
    }
  }
 
  out=Paso_IndexList_createPattern(0, A->numOutput,index_list,0,INDEXLIST_LENGTH,0);

  #ifdef Paso_TRACE
  printf("Paso_Pattern_binop: new pattern has been allocated.\n");
  #endif

 /* clean up */
   if (index_list!=NULL) {
        #pragma omp parallel for private(i) 
        for(i=0;i<A->numOutput;++i) Paso_IndexList_free(index_list[i].extension);
     }
  TMPMEMFREE(index_list);

  return out;
}

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
        in->extension=TMPMEMALLOC(1,Paso_IndexList);
        if (Paso_checkPtr(in->extension)) return;
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

/* count the number of row indices in the Paso_IndexList in */

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
    TMPMEMFREE(in);
  }
}

/* creates a Paso_pattern from a range of indices */
Paso_Pattern* Paso_IndexList_createPattern(dim_t n0, dim_t n,Paso_IndexList* index_list,index_t range_min,index_t range_max,index_t index_offset)
{
   dim_t *ptr=NULL;
   register dim_t s,i,itmp;
   index_t *index=NULL;
   Paso_Pattern* out=NULL;

   ptr=MEMALLOC(n+1-n0,index_t);
   if (! Paso_checkPtr(ptr) ) {
       /* get the number of connections per row */
       #pragma omp parallel for schedule(static) private(i)
       for(i=n0;i<n;++i) {
              ptr[i-n0]=Paso_IndexList_count(&index_list[i],range_min,range_max);
       }
       /* accumulate ptr */
       s=0;
       for(i=n0;i<n;++i) {
               itmp=ptr[i-n0];
               ptr[i-n0]=s;
               s+=itmp;
       }
       ptr[n-n0]=s;
       /* fill index */
       index=MEMALLOC(ptr[n-n0],index_t);
       if (! Paso_checkPtr(index)) {
              #pragma omp parallel for schedule(static)
              for(i=n0;i<n;++i) {
                  Paso_IndexList_toArray(&index_list[i],&index[ptr[i-n0]],range_min,range_max,index_offset);
              }
              out=Paso_Pattern_alloc(PATTERN_FORMAT_DEFAULT,1,1,n-n0,range_max+index_offset,ptr,index);
       }
  }
  if (! Paso_noError()) {
        MEMFREE(ptr);
        MEMFREE(index);
        Paso_Pattern_free(out);
  }
  return out;
}
