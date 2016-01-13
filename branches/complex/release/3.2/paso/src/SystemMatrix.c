
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

/* Paso: SystemMatrix */

/**************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "SystemMatrix.h"
#include "Preconditioner.h"

/**************************************************************/

/* allocates a SystemMatrix of type type using the given matrix pattern 
   Values are initialized by zero. 
   if patternIsUnrolled and type & MATRIX_FORMAT_BLK1, it is assumed that the pattern is allready unrolled to match the requested block size
   and offsets otherwise unrolling and offset adjustment will be performed. 
*/

Paso_SystemMatrix* Paso_SystemMatrix_alloc(Paso_SystemMatrixType type,Paso_SystemMatrixPattern *pattern, int row_block_size, int col_block_size,
  const bool_t patternIsUnrolled) {

  Paso_SystemMatrix*out=NULL;
  dim_t n_norm,i;
  Paso_SystemMatrixType pattern_format_out;
  bool_t unroll=FALSE;

  pattern_format_out= (type & MATRIX_FORMAT_OFFSET1)? MATRIX_FORMAT_OFFSET1:  MATRIX_FORMAT_DEFAULT;
  Esys_resetError();
  if (patternIsUnrolled) {
     if ( ! XNOR(type & MATRIX_FORMAT_OFFSET1, pattern->type & MATRIX_FORMAT_OFFSET1) ) {
         Esys_setError(TYPE_ERROR,"Paso_SystemMatrix_alloc: requested offset and pattern offset does not match.");
         return NULL;
     }
  }
  /* do we need to apply unrolling ? */
  unroll  
        /* we don't like non-square blocks */
    =   (row_block_size!=col_block_size)
        /* or any block size bigger than 3 */
    ||  (col_block_size>3) 
        /* or if lock size one requested and the block size is not 1 */
    ||  ((type & MATRIX_FORMAT_BLK1) &&  (col_block_size>1) )
        /* or the offsets are wrong */
    ||  ((type & MATRIX_FORMAT_OFFSET1) != ( pattern->type & MATRIX_FORMAT_OFFSET1));
  
  out=MEMALLOC(1,Paso_SystemMatrix);
  if (! Esys_checkPtr(out)) {  
     out->type=type;
     out->pattern=NULL;  
     out->row_distribution=NULL;
     out->col_distribution=NULL;
     out->mpi_info=Esys_MPIInfo_getReference(pattern->mpi_info);
     out->row_coupler=NULL;
     out->col_coupler=NULL;
     out->mainBlock=NULL;
     out->row_coupleBlock=NULL;
     out->col_coupleBlock=NULL;
     out->is_balanced=FALSE;
     out->balance_vector=NULL; 
     out->solver_package=PASO_PASO;  
     out->solver_p=NULL;  
     out->trilinos_data=NULL;
     out->reference_counter=1;
     out->logical_row_block_size=row_block_size;
     out->logical_col_block_size=col_block_size;


     if (type & MATRIX_FORMAT_CSC) {
         if (unroll) {
               if (patternIsUnrolled) {
                  out->pattern=Paso_SystemMatrixPattern_getReference(pattern);
               } else {
                  out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,pattern_format_out,col_block_size,row_block_size);
               }
               out->row_block_size=1;
               out->col_block_size=1;
         } else {
              out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,pattern_format_out,1,1);
              out->row_block_size=row_block_size;
              out->col_block_size=col_block_size;
         }
         if (Esys_noError()) {
           out->row_distribution=Paso_Distribution_getReference(out->pattern->input_distribution);
           out->col_distribution=Paso_Distribution_getReference(out->pattern->output_distribution);
         }
     } else {
         if (unroll) {
              if (patternIsUnrolled) {
                  out->pattern=Paso_SystemMatrixPattern_getReference(pattern);
              } else {
                  out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,pattern_format_out,row_block_size,col_block_size);
              }
              out->row_block_size=1;
              out->col_block_size=1;
         } else {
              out->pattern=Paso_SystemMatrixPattern_unrollBlocks(pattern,pattern_format_out,1,1);
              out->row_block_size=row_block_size;
              out->col_block_size=col_block_size;
         }
         if (Esys_noError()) {
              out->row_distribution=Paso_Distribution_getReference(out->pattern->output_distribution);
              out->col_distribution=Paso_Distribution_getReference(out->pattern->input_distribution);
         }
     }
     if (Esys_noError()) {
        out->block_size=out->row_block_size*out->col_block_size;
        out->col_coupler=Paso_Coupler_alloc(pattern->col_connector,out->col_block_size);
        out->row_coupler=Paso_Coupler_alloc(pattern->row_connector,out->row_block_size);
        /* this should be bypassed if trilinos is used */
        if (type & MATRIX_FORMAT_TRILINOS_CRS) {
           #ifdef TRILINOS
           out->trilinos_data=Paso_TRILINOS_alloc();
           #endif
        } else {
           out->solver_package=PASO_PASO;  
           out->mainBlock=Paso_SparseMatrix_alloc(type,out->pattern->mainPattern,row_block_size,col_block_size,TRUE);
           out->col_coupleBlock=Paso_SparseMatrix_alloc(type,out->pattern->col_couplePattern,row_block_size,col_block_size,TRUE);
           out->row_coupleBlock=Paso_SparseMatrix_alloc(type,out->pattern->row_couplePattern,row_block_size,col_block_size,TRUE);
           if (Esys_noError()) {
              /* allocate memory for matrix entries */
              n_norm = MAX(out->mainBlock->numCols * out->col_block_size, out->mainBlock->numRows * out->row_block_size);
	      out->balance_vector=MEMALLOC(n_norm,double);
	      out->is_balanced=FALSE;
	      if (! Esys_checkPtr(out->balance_vector)) {
                 #pragma omp parallel for private(i) schedule(static)
                 for (i=0;i<n_norm;++i) out->balance_vector[i]=1.;
              }
           }
        }
     }
  }
  /* all done: */
  if (! Esys_noError()) {
    Paso_SystemMatrix_free(out);
    return NULL;
  } else {
    return out;
  }
}

/* returns a reference to Paso_SystemMatrix in */

Paso_SystemMatrix* Paso_SystemMatrix_getReference(Paso_SystemMatrix* in) {
   if (in!=NULL) ++(in->reference_counter);
   return in;
}

/* deallocates a SystemMatrix: */

void Paso_SystemMatrix_free(Paso_SystemMatrix* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<=0) {
	Paso_solve_free(in);
        Paso_SystemMatrixPattern_free(in->pattern);
        Paso_Distribution_free(in->row_distribution);
        Paso_Distribution_free(in->col_distribution);
        Esys_MPIInfo_free(in->mpi_info);
        Paso_Coupler_free(in->row_coupler);
        Paso_Coupler_free(in->col_coupler);
        Paso_SparseMatrix_free(in->mainBlock);
        Paso_SparseMatrix_free(in->col_coupleBlock);
        Paso_SparseMatrix_free(in->row_coupleBlock);
	MEMFREE(in->balance_vector);
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
void  Paso_SystemMatrix_startCollect(Paso_SystemMatrix* A,const double* in)
{
  Paso_SystemMatrix_startColCollect(A,in);
}
double* Paso_SystemMatrix_finishCollect(Paso_SystemMatrix* A)
{
 return Paso_SystemMatrix_finishColCollect(A);
}

void  Paso_SystemMatrix_startColCollect(Paso_SystemMatrix* A,const double* in)
{
  Paso_Coupler_startCollect(A->col_coupler, in);
}
double* Paso_SystemMatrix_finishColCollect(Paso_SystemMatrix* A)
{
 Paso_Coupler_finishCollect(A->col_coupler);
 return A->col_coupler->recv_buffer;
}
void  Paso_SystemMatrix_startRowCollect(Paso_SystemMatrix* A,const double* in)
{
  Paso_Coupler_startCollect(A->row_coupler, in);
}
double* Paso_SystemMatrix_finishRowCollect(Paso_SystemMatrix* A)
{
 Paso_Coupler_finishCollect(A->row_coupler);
 return A->row_coupler->recv_buffer;
}

dim_t Paso_SystemMatrix_getTotalNumRows(const Paso_SystemMatrix* A){
  return A->mainBlock->numRows * A->row_block_size;
}

dim_t Paso_SystemMatrix_getTotalNumCols(const Paso_SystemMatrix* A){
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

index_t* Paso_SystemMatrix_borrowMainDiagonalPointer(Paso_SystemMatrix * A_p) 
{
    index_t* out=NULL;
    int fail=0;
    out=Paso_SparseMatrix_borrowMainDiagonalPointer(A_p->mainBlock);
    if (out==NULL) fail=1;
    #ifdef ESYS_MPI
    {
         int fail_loc = fail;
         MPI_Allreduce(&fail_loc, &fail, 1, MPI_INT, MPI_MAX, A_p->mpi_info->comm);
    }
    #endif
    if (fail>0) Esys_setError(VALUE_ERROR, "Paso_SystemMatrix_borrowMainDiagonalPointer: no main diagonal");
    return out;
}

void Paso_SystemMatrix_makeZeroRowSums(Paso_SystemMatrix * A_p, double* left_over) 
{
   index_t ir, ib, irow;
   register double rtmp1, rtmp2;
   const dim_t n = Paso_SystemMatrixPattern_getNumOutput(A_p->pattern);
   const dim_t nblk = A_p->block_size;
   const dim_t blk = A_p->row_block_size;
   const index_t* main_ptr=Paso_SystemMatrix_borrowMainDiagonalPointer(A_p);
   
   
   Paso_SystemMatrix_rowSum(A_p, left_over); /* left_over now hold the row sum */

   #pragma omp parallel for private(ir,ib, rtmp1, rtmp2) schedule(static)
   for (ir=0;ir< n;ir++) {
       for (ib=0;ib<blk; ib++) {
	     irow=ib+blk*ir;
	     rtmp1=left_over[irow];
	     rtmp2=A_p->mainBlock->val[main_ptr[ir]*nblk+ib+blk*ib];
	     A_p->mainBlock->val[main_ptr[ir]*nblk+ib+blk*ib] = -rtmp1;
	     left_over[irow]=rtmp2+rtmp1;
       }
   }
}
void Paso_SystemMatrix_copyBlockFromMainDiagonal(Paso_SystemMatrix * A_p, double* out)
{
    Paso_SparseMatrix_copyBlockFromMainDiagonal(A_p->mainBlock, out);
    return;
}
void Paso_SystemMatrix_copyBlockToMainDiagonal(Paso_SystemMatrix * A_p, const double* in) 
{
    Paso_SparseMatrix_copyBlockToMainDiagonal(A_p->mainBlock, in);
    return;
}
void Paso_SystemMatrix_copyFromMainDiagonal(Paso_SystemMatrix * A_p, double* out)
{
    Paso_SparseMatrix_copyFromMainDiagonal(A_p->mainBlock, out);
    return;
}
void Paso_SystemMatrix_copyToMainDiagonal(Paso_SystemMatrix * A_p, const double* in) 
{
    Paso_SparseMatrix_copyToMainDiagonal(A_p->mainBlock, in);
    return;
}

void Paso_SystemMatrix_setPreconditioner(Paso_SystemMatrix* A,Paso_Options* options) {
   if (A->solver_p==NULL) {
      A->solver_p=Paso_Preconditioner_alloc(A,options);
   }
}

/* applies the preconditioner */
/* has to be called within a parallel reqion */
/* barrier synchronization is performed before the evaluation to make sure that the input vector is available */
void Paso_SystemMatrix_solvePreconditioner(Paso_SystemMatrix* A,double* x,double* b){
   Paso_Preconditioner* prec=(Paso_Preconditioner*) A->solver_p;
   Paso_Preconditioner_solve(prec, A,x,b);
}
void Paso_SystemMatrix_freePreconditioner(Paso_SystemMatrix* A) {
   Paso_Preconditioner* prec=NULL;
   if (A!=NULL) {
      prec=(Paso_Preconditioner*) A->solver_p;
      Paso_Preconditioner_free(prec);
      A->solver_p=NULL;
   }
}
