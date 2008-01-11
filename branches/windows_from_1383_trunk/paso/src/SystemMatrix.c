
/* $Id$ */

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

/* Paso: SystemMatrix */

/**************************************************************/

/* Author: gross@access.edu.au */

/**************************************************************/

#include "SystemMatrix.h"

/**************************************************************/

/* allocates a SystemMatrix of type type using the given matrix pattern 
   if type is UNKOWN CSR is used.
   if CSC or CSC_BLK1 is used pattern has to give the CSC pattern.
   if CSR or CSR_BLK1 is used pattern has to give the CSR pattern.
   Values are initialized by zero.  */

Paso_SystemMatrix* Paso_SystemMatrix_alloc(Paso_SystemMatrixType type,Paso_SystemMatrixPattern *pattern, int row_block_size, int col_block_size) {

  double time0;
  Paso_SystemMatrix*out=NULL;
  dim_t n_norm,i;
  Paso_SystemMatrixType pattern_format_out;

  Paso_resetError();
  time0=Paso_timer();
  out=MEMALLOC(1,Paso_SystemMatrix);
  if (! Paso_checkPtr(out)) {  
     out->type=type;
     out->pattern=NULL;  
     out->row_distribution=NULL;
     out->col_distribution=NULL;
     out->mpi_info=Paso_MPIInfo_getReference(pattern->mpi_info);
     out->mainBlock=NULL;
     out->coupleBlock=NULL;
     out->normalizer_is_valid=FALSE;
     out->normalizer=NULL; 
     out->solver_package=PASO_PASO;  
     out->solver=NULL;  
     out->trilinos_data=NULL;
     out->reference_counter=1;

     pattern_format_out= (type & MATRIX_FORMAT_OFFSET1)? PATTERN_FORMAT_OFFSET1:  PATTERN_FORMAT_DEFAULT;
     /* ====== compressed sparse columns === */
     if (type & MATRIX_FORMAT_CSC) {
        if (type & MATRIX_FORMAT_SYM) {
           Paso_setError(TYPE_ERROR,"Generation of matrix pattern for symmetric CSC is not implemented yet.");
        } else {
           if ((type & MATRIX_FORMAT_BLK1) || row_block_size!=col_block_size || col_block_size>3) {
              out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,pattern_format_out,col_block_size,row_block_size);
              out->row_block_size=1;
              out->col_block_size=1;
           } else {
              out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,pattern_format_out,1,1);
              out->row_block_size=row_block_size;
              out->col_block_size=col_block_size;
           }
        }
        out->row_distribution=Paso_Distribution_getReference(out->pattern->input_distribution);
        out->col_distribution=Paso_Distribution_getReference(out->pattern->output_distribution);
     } else {
     /* ====== compressed sparse row === */
        if (type & MATRIX_FORMAT_SYM) {
           Paso_setError(TYPE_ERROR,"Generation of matrix pattern for symmetric CSR is not implemented yet.");
        } else {
           if ((type & MATRIX_FORMAT_BLK1) || row_block_size!=col_block_size || col_block_size>3)  {
              out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,pattern_format_out,row_block_size,col_block_size);
              out->row_block_size=1;
              out->col_block_size=1;
           } else {
              out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,pattern_format_out,1,1);
              out->row_block_size=row_block_size;
              out->col_block_size=col_block_size;
           }
        }
        out->row_distribution=Paso_Distribution_getReference(out->pattern->output_distribution);
        out->col_distribution=Paso_Distribution_getReference(out->pattern->input_distribution);
     }
     out->logical_row_block_size=row_block_size;
     out->logical_col_block_size=col_block_size;
     out->logical_block_size=out->logical_row_block_size*out->logical_block_size;
     out->block_size=out->row_block_size*out->col_block_size;
     /* this should be bypassed if trilinos is used */
     if (type & MATRIX_FORMAT_TRILINOS_CRS) {
        #ifdef TRILINOS
        out->trilinos_data=Paso_TRILINOS_alloc();
        #endif
     } else {
        out->solver_package=PASO_PASO;  
        out->mainBlock=Paso_SparseMatrix_alloc(type,out->pattern->mainPattern,row_block_size,col_block_size);
        out->coupleBlock=Paso_SparseMatrix_alloc(type,out->pattern->couplePattern,row_block_size,col_block_size);
        /* allocate memory for matrix entries */
        if (type & MATRIX_FORMAT_CSC) {
           n_norm = out->mainBlock->numCols * out->col_block_size;
        } else {
           n_norm = out->mainBlock->numRows * out->row_block_size;
        }
        out->normalizer=MEMALLOC(n_norm,double);
        out->normalizer_is_valid=FALSE;
        if (! Paso_checkPtr(out->normalizer)) {
           #pragma omp parallel for private(i) schedule(static)
           for (i=0;i<n_norm;++i) out->normalizer[i]=0.;
       }
    }
  }
  /* all done: */
  if (! Paso_noError()) {
    Paso_SystemMatrix_free(out);
    return NULL;
  } else {
    #ifdef Paso_TRACE
    printf("timing: system matrix %.4e sec\n",Paso_timer()-time0);
    printf("Paso_SystemMatrix_alloc: system matrix has been allocated.\n");
    #endif
    return out;
  }
}

/* returns a reference to Paso_SystemMatrix in */

Paso_SystemMatrix* Paso_SystemMatrix_reference(Paso_SystemMatrix* in) {
   if (in!=NULL) ++(in->reference_counter);
   return NULL;
}

/* deallocates a SystemMatrix: */

void Paso_SystemMatrix_free(Paso_SystemMatrix* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<=0) {
        Paso_SystemMatrixPattern_free(in->pattern);
        Paso_Distribution_free(in->row_distribution);
        Paso_Distribution_free(in->col_distribution);
        Paso_MPIInfo_free(in->mpi_info);
        Paso_SparseMatrix_free(in->mainBlock);
        Paso_SparseMatrix_free(in->coupleBlock);
        MEMFREE(in->normalizer);
        Paso_solve_free(in); 
        #ifdef TRILINOS
        Paso_TRILINOS_free(in->trilinos_data);
        #endif
        MEMFREE(in);
        #ifdef Paso_TRACE
        printf("Paso_SystemMatrix_free: system matrix as been deallocated.\n");
        #endif
     }
   }
}
void Paso_SystemMatrix_allocBuffer(Paso_SystemMatrix* A) {
   if (! Paso_Coupler_bufferIsAllocated(A->pattern->coupler) ) {
       Paso_Coupler_allocBuffer(A->pattern->coupler,A->col_block_size);
   }
}
void Paso_SystemMatrix_freeBuffer(Paso_SystemMatrix* A) {
    Paso_Coupler_freeBuffer(A->pattern->coupler);
}
void  Paso_SystemMatrix_startCollect(Paso_SystemMatrix* A, double* in)
{
  Paso_Coupler_startCollect(A->pattern->coupler, in);
}
double* Paso_SystemMatrix_finishCollect(Paso_SystemMatrix* A)
{
 Paso_Coupler_finishCollect(A->pattern->coupler);
 return A->pattern->coupler->recv_buffer;
}

dim_t Paso_SystemMatrix_getTotalNumRows(Paso_SystemMatrix* A){
  return A->mainBlock->numRows * A->row_block_size;
}

dim_t Paso_SystemMatrix_getTotalNumCols(Paso_SystemMatrix* A){
  return A->mainBlock->numCols * A->col_block_size;
}
dim_t Paso_SystemMatrix_getGlobalNumRows(Paso_SystemMatrix* A) {
  if (A->type & MATRIX_FORMAT_CSC) {
      return  Paso_Distribution_getGlobalNumComponents(A->pattern->input_distribution);
  }  else {
      return  Paso_Distribution_getGlobalNumComponents(A->pattern->output_distribution);
  }
}
dim_t Paso_SystemMatrix_getGlobalNumCols(Paso_SystemMatrix* A) {
  if (A->type & MATRIX_FORMAT_CSC) {
      return  Paso_Distribution_getGlobalNumComponents(A->pattern->output_distribution);
  }  else {
      return  Paso_Distribution_getGlobalNumComponents(A->pattern->input_distribution);
  }

}
dim_t Paso_SystemMatrix_getNumOutput(Paso_SystemMatrix* A) {
   return Paso_SystemMatrixPattern_getNumOutput(A->pattern);
}

