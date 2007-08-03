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

/* Paso: SparseMatrix is saved to Matrix Market format      */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004 */
/* Author: davies@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "mmio.h"
#include "SparseMatrix.h"

void Paso_SparseMatrix_saveMM_CSC(Paso_SparseMatrix * A_p, FILE * fileHandle_p) {
  int iRow, iCol, iPtr,ir,ic;
  index_t index_offset;
  MM_typecode matrixCode;
  index_offset=(A_p->type & MATRIX_FORMAT_OFFSET1 ? 1:0);

  if (A_p->val == NULL) return;
  /* set the matrix code */
  mm_initialize_typecode(&matrixCode);
  mm_set_matrix(&matrixCode);
  mm_set_coordinate(&matrixCode);
  mm_set_real(&matrixCode);
  
  mm_write_banner(fileHandle_p, matrixCode);
  mm_write_mtx_crd_size(fileHandle_p, A_p->numRows*A_p->row_block_size, A_p->numCols*A_p->col_block_size,A_p->len);

  for (iCol = 0; iCol< A_p->pattern->numOutput; iCol++) {
     for (ic = 0; ic< A_p->col_block_size; ic++) {
         for (iPtr = A_p->pattern->ptr[iCol] - index_offset;iPtr < A_p->pattern->ptr[iCol+1] - index_offset; iPtr++) {
            for (ir = 0; ir< A_p->row_block_size; ir++) {
                 fprintf(fileHandle_p, "%12d %12d %e\n",
		          (A_p->pattern->index[iPtr]-index_offset)*A_p->row_block_size+ir+1, 
		          iCol*A_p->col_block_size+ic+ 1, 
		          A_p->val[iPtr*A_p->block_size+ir+A_p->row_block_size*ic]);
               }
            }
         }
      }

  return;
}

void Paso_SparseMatrix_saveMM_CSR(Paso_SparseMatrix * A_p, FILE * fileHandle_p) {
  int iRow, iCol, iPtr,ir,ic;
  index_t index_offset;
  MM_typecode matrixCode;
  index_offset=(A_p->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
   if (A_p->val == NULL) {
       return;
   }
  /* set the matrix code */
  mm_initialize_typecode(&matrixCode);
  mm_set_matrix(&matrixCode);
  mm_set_coordinate(&matrixCode);
  mm_set_real(&matrixCode);
  
  mm_write_banner(fileHandle_p, matrixCode);
  mm_write_mtx_crd_size(fileHandle_p, A_p->numRows*A_p->row_block_size, A_p->numCols*A_p->col_block_size,A_p->len);

  for (iRow = 0; iRow< A_p->pattern->numOutput; iRow++) {
      for (ir = 0; ir< A_p->row_block_size; ir++) {
         for (iPtr = A_p->pattern->ptr[iRow] - index_offset;iPtr < A_p->pattern->ptr[iRow+1] - index_offset; iPtr++) {
            for (ic = 0; ic< A_p->col_block_size; ic++) {
                fprintf(fileHandle_p, "%12d %12d %e\n",
		          iRow*A_p->row_block_size+ir+ 1, 
		          (A_p->pattern->index[iPtr]-index_offset)*A_p->col_block_size+ic+1, 
		          A_p->val[iPtr*A_p->block_size+ir+A_p->row_block_size*ic]);
             }
         }
      }
  }
  return;
}
