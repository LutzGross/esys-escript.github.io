/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrix */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "IndexList.h"
#include "System.h"
#include "Mesh.h"

/**************************************************************/

/* allocates a SystemMatrix in the CSR/Harwell boeing format from a
   Finley_Mesh. Values are initialized by zero.  */

Finley_SystemMatrix* Finley_SystemMatrix_alloc(Finley_Mesh *mesh,
  Finley_SystemMatrixType type, int symmetric,int row_block_size, int reduce_row_order,
  int col_block_size, int reduce_col_order) {
  double time0,time1;
  maybelong packed_row_block_size, packed_col_block_size, unpacked_row_block_size, unpacked_col_block_size,i,block_size,j,iptr;
  maybelong *rowLabel=NULL,*colLabel=NULL,len_index_list=0;
  Finley_SystemMatrix*out;
  Finley_IndexList* index_list=NULL;
  Finley_ErrorCode=NO_ERROR;
  
  time0=Finley_timer();
  /* the matrix block row_block_size x col_block_size is fully
  incorporated into the matrix pattern.  this is required for using
  the SGI scientific library. row_offset is the offset for the row
  index typically */

  packed_row_block_size=row_block_size;
  packed_col_block_size=col_block_size;
  unpacked_row_block_size=1;
  unpacked_col_block_size=1;
 

  /*  allocate the return value */

  out=(Finley_SystemMatrix*)MEMALLOC(sizeof(Finley_SystemMatrix));
  if (Finley_checkPtr(out)) goto clean;

  if (type==UNKNOWN) {
     out->type=FINLEY_DEFAULT_MATRIX_TYPE;
  } else {
     out->type=type;
  }
  
  out->total_row_block_size=row_block_size;
  out->total_col_block_size=col_block_size;

  out->row_block_size=unpacked_row_block_size;
  out->col_block_size=unpacked_col_block_size;

  block_size=out->row_block_size*out->col_block_size;

/* ******************************** */
  if (reduce_col_order) {
     out->num_cols=packed_col_block_size*mesh->Nodes->reducedNumDegreesOfFreedom;
     colLabel=mesh->Nodes->reducedDegreeOfFreedom;
  } else {
     out->num_cols=packed_col_block_size*mesh->Nodes->numDegreesOfFreedom;
     colLabel=mesh->Nodes->degreeOfFreedom;
  }
     
  if (reduce_row_order) {
     out->num_rows=packed_row_block_size*mesh->Nodes->reducedNumDegreesOfFreedom;
     rowLabel=mesh->Nodes->reducedDegreeOfFreedom;
  } else {
     out->num_rows=packed_row_block_size*mesh->Nodes->numDegreesOfFreedom;
     rowLabel=mesh->Nodes->degreeOfFreedom;
  }
/* ******************************** */

  out->symmetric=symmetric;
  out->index=NULL;
  out->val=NULL;
  out->ptr=NULL;
  out->reference_counter=0;
  out->lenOfVal=0;
  out->solve=NULL;
  out->iterative=NULL;

  /*  initialize the (temporary) list index_list of the colums indices: */
  switch(out->type) {
  case CSR:
    len_index_list=out->num_rows;
    break;
  case CSC:
    len_index_list=out->num_cols;
    break;
  default:
    Finley_ErrorCode=TYPE_ERROR;
    sprintf(Finley_ErrorMsg,"Unknown matrix type.");
    goto clean;
  } /* switch out->type */

  index_list=(Finley_IndexList*) TMPMEMALLOC(len_index_list*sizeof(Finley_IndexList));
  if (Finley_checkPtr(index_list)) goto clean;

  #pragma omp parallel for private(i) schedule(static)
  for(i=0;i<len_index_list;i++) {
       index_list[i].extension=NULL;
       index_list[i].n=0;
  }

  /*  insert contributions from element matrices into colums index index_list: */
  Finley_IndexList_insertElements(index_list,mesh->Elements,
                                  reduce_row_order,packed_row_block_size,rowLabel,
                                  reduce_col_order,packed_col_block_size,colLabel,
                                  symmetric, out->type);
  Finley_IndexList_insertElements(index_list,mesh->FaceElements,
                                  reduce_row_order,packed_row_block_size,rowLabel,
                                  reduce_col_order,packed_col_block_size,colLabel,
                                  symmetric, out->type);
  Finley_IndexList_insertElements(index_list,mesh->ContactElements,
                                  reduce_row_order,packed_row_block_size,rowLabel,
                                  reduce_col_order,packed_col_block_size,colLabel,
                                  symmetric, out->type);
  Finley_IndexList_insertElements(index_list,mesh->Points,
                                  reduce_row_order,packed_row_block_size,rowLabel,
                                  reduce_col_order,packed_col_block_size,colLabel,
                                  symmetric, out->type);
  if (Finley_ErrorCode!=NO_ERROR) goto clean;

  /* allocate the ptr: */

  out->ptr=(maybelong*)MEMALLOC(((len_index_list)+1)*sizeof(maybelong));
  if (Finley_checkPtr(out->ptr)) {
    Finley_SystemMatrix_dealloc(out);
    goto clean;
  }
  #pragma omp parallel for private(i) schedule(static)
  for(i=0;i<len_index_list;i++) out->ptr[i]=0;
  out->ptr[len_index_list]=0;

  /* count entries in each column and store to pointer vector: */
  out->ptr[0]=PTR_OFFSET;
  /* OMP */
  for(i=0;i<len_index_list;i++)
    out->ptr[i+1]=out->ptr[i]+Finley_IndexList_count(&index_list[i]);

  /* allocate index and val: */

  out->index=(maybelong*)MEMALLOC((out->ptr[len_index_list]-PTR_OFFSET)*sizeof(maybelong));
  out->lenOfVal=(out->ptr[len_index_list]-PTR_OFFSET)*out->row_block_size*out->col_block_size;
  out->val=(double*)MEMALLOC(out->lenOfVal*sizeof(double));
  if (Finley_checkPtr(out->index) || Finley_checkPtr(out->val) ) {
    Finley_SystemMatrix_dealloc(out);
    goto clean;
  }
  #pragma omp parallel firstprivate(out,index_list,len_index_list)
  { 

#pragma omp master
{
time1=Finley_timer();
}
     #pragma omp for private(i,iptr,j) schedule(static) 
     for (i=0;i< len_index_list;i++) {
        for (iptr=out->ptr[i]-PTR_OFFSET;iptr<out->ptr[i+1]-PTR_OFFSET; iptr++) {
               out->index[iptr]=0;
               for (j=0;j<block_size;j++) out->val[iptr*block_size+j]=0.;
        }
     }
#pragma omp master
{
printf("timing: matrix pattern: instantation: %.4e sec\n",Finley_timer()-time1);
}


     /* change the list of the row indicies into an array:  */
     /* each index is ordered by increased size  */
   
     #pragma omp for private(i) schedule(static)
     for(i=0;i<len_index_list;i++) {
       Finley_IndexList_toArray(&index_list[i],&(out->index[out->ptr[i]-PTR_OFFSET]));
       qsort(&(out->index[out->ptr[i]-PTR_OFFSET]),(int)(out->ptr[i+1]-out->ptr[i]),sizeof(maybelong), Finley_comparIndex); 
     }

  }
  /* out->val is initialized to zero: */

time1=Finley_timer();
  Finley_SystemMatrix_setValues(out,DBLE(0));
printf("timing: matrix pattern: initialize: %.4e sec\n",Finley_timer()-time1);
  /* all done: */


 clean:
  if (index_list!=NULL) {
    #pragma omp parallel for private(i) 
    for(i=0;i<len_index_list;i++) Finley_IndexList_free(index_list[i].extension);
  }
  TMPMEMFREE(index_list);

  printf("timing: matrix pattern: %.4e sec\n",Finley_timer()-time0);

  if (Finley_ErrorCode!=NO_ERROR) {
    return NULL;
  } else {
    #ifdef Finley_TRACE
    printf("Finley_SystemMatrix_alloc: %ld x %ld system matrix has been allocated.\n",(long)out->num_rows,(long)out->num_cols);
    #endif
    out->reference_counter++;
    return out;
  }
}

/* deallocates a SystemMatrix: */

void Finley_SystemMatrix_dealloc(Finley_SystemMatrix* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<=0) {
        MEMFREE(in->ptr);
        MEMFREE(in->index);
        MEMFREE(in->val);
        Finley_SystemMatrix_solve_free(in);
        Finley_SystemMatrix_iterative_free(in);
        MEMFREE(in);
        #ifdef Finley_TRACE
        printf("Finley_SystemMatrix_dealloc: system matrix as been deallocated.\n");
        #endif
     }
   }
}
/* *************************************************************/

/*  some routines which help to get the matrix pattern from elements: */

/*  this routine is used by qsort called in Finley_SystemMatrix_alloc */

int Finley_comparIndex(const void *index1,const void *index2){
   maybelong Iindex1,Iindex2;
   Iindex1=*(maybelong*)index1;
   Iindex2=*(maybelong*)index2;
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
/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
