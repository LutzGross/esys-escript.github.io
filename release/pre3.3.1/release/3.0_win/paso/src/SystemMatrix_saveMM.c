
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: SystemMatrix is saved to Matrix Market format      */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004 */
/* Author: davies@access.edu.au */

/**************************************************************/

#include "SystemMatrix.h"
#include "mmio.h"

void Paso_SystemMatrix_saveMM(Paso_SystemMatrix * A_p, char * fileName_p) {
  FILE * fileHandle_p = NULL;
  dim_t N,M,i, iptr_ij;
  MM_typecode matcode;                        

  if (A_p->mpi_info->size > 1) {
       Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_saveMM: currently single processor runs are supported.\n");
       return;
  }
  if (A_p->col_block_size !=A_p->row_block_size) {
    Paso_setError(TYPE_ERROR, "Paso_SystemMatrix_saveMM: currently only square blocks are supported.");
    return;
  }
  if (A_p->row_block_size>3) {
       Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_saveMM: currently only block size 3 is supported.\n");
       return;
  }

  if (A_p->type & MATRIX_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_saveMM does not support symmetric storage scheme");
    return;
  }
  /* open the file */
  fileHandle_p = fopen(fileName_p, "w");
  if (fileHandle_p==NULL) {
    Paso_setError(IO_ERROR,"file could not be opened for writing");
    return;
  }
  if (A_p->type & MATRIX_FORMAT_CSC) {
    Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_saveMM does not support CSC yet.");
  } else {
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);

    N=Paso_SystemMatrix_getGlobalNumRows(A_p);
    M=Paso_SystemMatrix_getGlobalNumCols(A_p);
    mm_write_banner(fileHandle_p, matcode); 
    mm_write_mtx_crd_size(fileHandle_p, N*A_p->row_block_size, M*A_p->col_block_size, A_p->mainBlock->pattern->ptr[N]*A_p->block_size);
    if(A_p->row_block_size==1)
      for (i=0; i<N; i++) 
        for (iptr_ij=A_p->mainBlock->pattern->ptr[i];iptr_ij<A_p->mainBlock->pattern->ptr[i+1]; ++iptr_ij) {
          fprintf(fileHandle_p, "%d %d %25.15e\n", i+1, A_p->mainBlock->pattern->index[iptr_ij]+1, A_p->mainBlock->val[iptr_ij]);
        }

    if(A_p->row_block_size==2)
      for (i=0; i<N; i++) {
        for (iptr_ij=A_p->mainBlock->pattern->ptr[i];iptr_ij<A_p->mainBlock->pattern->ptr[i+1]; ++iptr_ij) {
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*2+1, 2*(A_p->mainBlock->pattern->index[iptr_ij])+1, A_p->mainBlock->val[iptr_ij*4]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*2+1, 2*(A_p->mainBlock->pattern->index[iptr_ij])+1+1, A_p->mainBlock->val[iptr_ij*4+2]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*2+1+1, 2*(A_p->mainBlock->pattern->index[iptr_ij])+1, A_p->mainBlock->val[iptr_ij*4+1]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*2+1+1, 2*(A_p->mainBlock->pattern->index[iptr_ij])+1+1, A_p->mainBlock->val[iptr_ij*4+3]);
        }
      }

    if(A_p->row_block_size==3)
      for (i=0; i<N; i++) {
       for (iptr_ij=A_p->mainBlock->pattern->ptr[i];iptr_ij<A_p->mainBlock->pattern->ptr[i+1]; ++iptr_ij) {
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+1, 3*(A_p->mainBlock->pattern->index[iptr_ij])+1, A_p->mainBlock->val[iptr_ij*9]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+1, 3*(A_p->mainBlock->pattern->index[iptr_ij])+1+1, A_p->mainBlock->val[iptr_ij*9+3]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+1, 3*(A_p->mainBlock->pattern->index[iptr_ij])+2+1, A_p->mainBlock->val[iptr_ij*9+6]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+1+1, 3*(A_p->mainBlock->pattern->index[iptr_ij])+1, A_p->mainBlock->val[iptr_ij*9+1]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+1+1, 3*(A_p->mainBlock->pattern->index[iptr_ij])+1+1, A_p->mainBlock->val[iptr_ij*9+4]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+1+1, 3*(A_p->mainBlock->pattern->index[iptr_ij])+2+1, A_p->mainBlock->val[iptr_ij*9+7]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+2+1, 3*(A_p->mainBlock->pattern->index[iptr_ij])+1, A_p->mainBlock->val[iptr_ij*9+2]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+2+1, 3*(A_p->mainBlock->pattern->index[iptr_ij])+1+1, A_p->mainBlock->val[iptr_ij*9+5]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+2+1, 3*(A_p->mainBlock->pattern->index[iptr_ij])+2+1, A_p->mainBlock->val[iptr_ij*9+8]);
        }
      } 
     
  }

  /* close the file */
  fclose(fileHandle_p);
  
  return;
}
