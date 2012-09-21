
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
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

/* Paso: Pattern */

/************************************************************************************/
 
/* Author: Lutz Gross, l.gross@uq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "Pattern.h"

/************************************************************************************/

/* allocates a Pattern  */

Paso_Pattern* Paso_Pattern_alloc(int type, dim_t numOutput, dim_t numInput, index_t* ptr, index_t* index) {
  Paso_Pattern *out=NULL;
  index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  index_t loc_min_index,loc_max_index,min_index=index_offset,max_index=index_offset-1;
  dim_t i;
  Esys_resetError();

  if (ptr!=NULL && index != NULL) {
    #pragma omp parallel private(loc_min_index,loc_max_index,i)
     {
        loc_min_index=index_offset;
        loc_max_index=index_offset-1;
        if (type & MATRIX_FORMAT_OFFSET1) {
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
      Esys_setError(TYPE_ERROR,"Paso_Pattern_alloc: Pattern index out of range.");
      return NULL;
    }
  }
  out=MEMALLOC(1,Paso_Pattern);
  if (! Esys_checkPtr(out)) {
      out->type=type;
      out->reference_counter=1;
      out->numOutput=numOutput;
      out->numInput=numInput;
      out->ptr=ptr;
      out->index=index;
      out->main_iptr = NULL;
      out->coloring = NULL;
      out->numColors=-1;

      if (out->ptr == NULL) {
          out->len=0;
      } else {
          out->len=out->ptr[out->numOutput] - index_offset;
      }
  }
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
	MEMFREE(in->main_iptr);
	MEMFREE(in->coloring);
        MEMFREE(in);
     }
   }
}
/* ***********************************************************************************/

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


/* creates a Paso_pattern from a range of indices */
Paso_Pattern* Paso_Pattern_fromIndexListArray(dim_t n0, Paso_IndexListArray* index_list_array,index_t range_min,index_t range_max,index_t index_offset)
{
   const dim_t n=index_list_array->n;
   Paso_IndexList* index_list = index_list_array->index_list;
   dim_t *ptr=NULL;
   register dim_t s,i,itmp;
   index_t *index=NULL;
   Paso_Pattern* out=NULL;

   ptr=MEMALLOC(n+1-n0,index_t);
   if (! Esys_checkPtr(ptr) ) {
       /* get the number of connections per row */ 
       #pragma omp parallel for private(i) schedule(static)
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
       if (! Esys_checkPtr(index)) {
              #pragma omp parallel for private(i) schedule(static) 
              for(i=n0;i<n;++i) {
                  Paso_IndexList_toArray(&index_list[i],&index[ptr[i-n0]],range_min,range_max,index_offset);
              }
              out=Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT,n-n0,range_max+index_offset,ptr,index);
       }
  }
  if (! Esys_noError()) {
        MEMFREE(ptr);
        MEMFREE(index);
        Paso_Pattern_free(out);
  }
  return out;
}

index_t* Paso_Pattern_borrowMainDiagonalPointer(Paso_Pattern* A) 
{
    const dim_t n=A->numOutput;
    int fail=0;
    index_t *index,*where_p, i;
    
     if (A->main_iptr == NULL) {
         A->main_iptr=MEMALLOC(n,index_t);
         if (! Esys_checkPtr(A->main_iptr) ) {
	     #pragma omp parallel 
             {
                 /* identify the main diagonals */
                 #pragma omp for schedule(static) private(i, index, where_p)
                 for (i = 0; i < n; ++i) {
		    index=&(A->index[A->ptr[i]]);
		    where_p=bsearch(&i,
				    index,
				    (size_t) (A->ptr[i + 1]-A->ptr[i]),
			             sizeof(index_t),
			             Paso_comparIndex);
					      
		    if (where_p==NULL) {
		        fail=1;
		    } else {
		       A->main_iptr[i]=A->ptr[i]+(index_t)(where_p-index);
		    }
                 }
     
             }
	     if (fail > 0) {
	       MEMFREE(A->main_iptr);
	       A->main_iptr=NULL;
	     }

	 }
     }
     return A->main_iptr;
}
			  

dim_t Paso_Pattern_getNumColors(Paso_Pattern* A)
{
   Paso_Pattern_borrowColoringPointer(A);  /* make sure numColors is defined */
   return A->numColors;
}
index_t* Paso_Pattern_borrowColoringPointer(Paso_Pattern* A)
{
   dim_t n=A->numInput;
   /* is coloring available ? */
   if (A->coloring == NULL) {
      
      A->coloring=MEMALLOC(n,index_t);
      if ( ! Esys_checkPtr(A->coloring)) {
	 Paso_Pattern_color(A,&(A->numColors),A->coloring);
	 if (! Esys_noError()) {
	    MEMFREE(A->coloring);
	 }
      } 
   }
   return A->coloring;
}
dim_t Paso_Pattern_maxDeg(Paso_Pattern* A)
{
   dim_t deg=0, loc_deg=0, i;
   const dim_t n=A->numInput;

   #pragma omp parallel private(i, loc_deg) 
   {
         loc_deg=0;
	 #pragma omp for schedule(static)
	 for (i = 0; i < n; ++i) {
	    loc_deg=MAX(loc_deg, A->ptr[i+1]-A->ptr[i]);
	 }
	 #pragma omp critical
	 {
	    deg=MAX(deg, loc_deg);
         }
   }
   return deg;
}

