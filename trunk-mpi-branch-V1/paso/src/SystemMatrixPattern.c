/* $Id$ */

/*
********************************************************************************
*               Copyright   2006 by ACcESS MNRF                                *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**************************************************************/

/* Paso: SystemMatrixPatternPattern */

/**************************************************************/
 
/* Copyrights by ACcESS Australia 2003, 2004,2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrixPattern.h"

/**************************************************************/

/* allocates a SystemMatrixPattern  */

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_alloc(int type, 
                                                         Paso_Distribution* output_distribution,
                                                         Paso_Distribution* input_distribution,
                                                         index_t* ptr,
                                                         index_t* index,
                                                         dim_t numHops,
                                                         index_t *hop) {
  Paso_SystemMatrixPattern*out;
  Paso_MPIInfo* mpi_info=output_distribution->mpi_info;
  index_t index_offset=(type & PATTERN_FORMAT_OFFSET1 ? 1:0);
  index_t loc_min_index,loc_max_index,min_index=index_offset,max_index=index_offset-1;
  dim_t i, sum=0;
  Paso_resetError();

  if (input_distribution->mpi_info != output_distribution->mpi_info) {
    Paso_setError(VALUE_ERROR,"Paso_SystemMatrixPattern_alloc: row and column distribution must base on the same communicator.");
    return NULL;
  }

  if (type & PATTERN_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_SystemMatrixPattern_alloc: symmetric matrix pattern is not supported yet");
    return NULL;
  }
  dim_t myNumOutput=output_distribution->myNumComponents;
  #pragma omp parallel private(loc_min_index,loc_max_index,i)
   {
      loc_min_index=index_offset;
      loc_max_index=index_offset-1;
      if (type & PATTERN_FORMAT_OFFSET1) {
         #pragma omp for schedule(static) 
         for (i=0;i<myNumOutput;++i) {
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
         for (i=0;i<myNumOutput;++i) {
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
  #ifdef Paso_MPI
     loc_min_index=min_index;
     loc_max_index=max_index;
     MPI_Reduce (&loc_max_index,&max_index,1,PASO_MPI_INT,MPI_MAX,mpi_info->comm)
     MPI_Reduce (&loc_min_index,&min_index,1,PASO_MPI_INT,MPI_MIN,mpi_info->comm)
  #endif
  if (min_index<index_offset) {
    Paso_setError(TYPE_ERROR,"Paso_SystemMatrixPattern_alloc: Pattern index out of index offset range.");
    return NULL;
  }
  if (min_index<input_distribution->firstComponent) {
    Paso_setError(TYPE_ERROR,"Paso_SystemMatrixPattern_alloc: Minimum pattern index out of input distribution range.");
    return NULL;
  }
  if (input_distribution->firstComponent+input_distribution->numComponents <= max_index) {
    Paso_setError(TYPE_ERROR,"Paso_SystemMatrixPattern_alloc: Maximum pattern index out of input distribution range.");
    return NULL;
  }
  sum=0;
  for (i=0;i<numHops;++i) sum+=hop[i];
  if (! (sum == output_distribution->mpi_info->size) ) {
    Paso_setError(SYSTEM_ERROR,"Paso_SystemMatrixPattern_alloc: processor hops do not define a closed loop.");
    return NULL;
  }

  out=MEMALLOC(1,Paso_SystemMatrixPattern);
  if (Paso_checkPtr(out)) return NULL;
  out->hop=MEMALLOC((output_distribution->mpi_info->size),index_t);
  if (Paso_checkPtr(hop)) {
     MEMFREE(out);
     return NULL;
  }
 
  out->type=type;
  out->reference_counter=1;
  out->myNumOutput=myNumOutput;
  out->maxNumOutput=output_distribution->maxNumComponents;
  out->numOutput=output_distribution->numComponents;
  out->myNumInput=input_distribution->myNumComponents;
  out->maxNumInput=input_distribution->maxNumComponents;
  out->numInput=input_distribution->numComponents;
  out->myLen=ptr[myNumOutput]-ptr[0];
  out->ptr=ptr;
  out->index=index;
  out->input_distribution=Paso_Distribution_getReference(input_distribution);
  out->output_distribution=Paso_Distribution_getReference(output_distribution);
  out->mpi_info = Paso_MPIInfo_getReference(mpi_info);
  out->numHops=numHops;
  for (i=0;i<out->numHops;++i) out->hop[i]=hop[i];
  /* this has to go */
  #ifdef PASO_MPI
      out->output_node_distribution=NULL;
      out->input_node_distribution=NULL;
  #endif
  #ifdef Paso_TRACE
  printf("Paso_SystemMatrixPattern_dealloc: system matrix pattern as been allocated.\n");
  #endif
  return out;
}

/* returns a reference to in */

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_reference(Paso_SystemMatrixPattern* in) {
     if (in!=NULL) {
        ++(in->reference_counter);
     }
     return in;
}
  
/* deallocates a SystemMatrixPattern: */

void Paso_SystemMatrixPattern_dealloc(Paso_SystemMatrixPattern* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<=0) {
        Paso_Distribution_free(in->output_distribution);
        Paso_Distribution_free(in->input_distribution);
        Paso_MPIInfo_dealloc(in->mpi_info);
        MEMFREE(in->ptr);
        MEMFREE(in->index);
        MEMFREE(in->hop);
        /* this has to go */
        #ifdef PASO_MPI
           Finley_NodeDistribution_dealloc(in->output_node_distribution);
           Finley_NodeDistribution_dealloc(in->input_node_distribution);
        #endif
        MEMFREE(in);
        #ifdef Paso_TRACE
        printf("Paso_SystemMatrixPattern_dealloc: system matrix pattern as been deallocated.\n");
        #endif
     }
   }
}
/* *************************************************************/

/*  some routines which help to get the matrix pattern from elements: */

/*  this routine is used by qsort called in Paso_SystemMatrixPattern_alloc */

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
