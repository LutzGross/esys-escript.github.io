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

/* Paso: SystemMatrix is saved to Matrix Market format      */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004 */
/* Author: davies@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "mmio.h"
#include "SystemMatrix.h"

void Paso_SystemMatrix_saveMM(Paso_SystemMatrix * A_p, char * fileName_p) {

  int iRow, iCol, iPtr,ir,ic;
  index_t index_offset=(A_p->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  FILE * fileHandle_p = NULL;
  MM_typecode matrixCode;

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

  /* set the matrix code */
  mm_initialize_typecode(&matrixCode);
  mm_set_matrix(&matrixCode);
  mm_set_coordinate(&matrixCode);
  mm_set_real(&matrixCode);
  
  mm_write_banner(fileHandle_p, matrixCode);
  mm_write_mtx_crd_size(fileHandle_p, A_p->num_rows*A_p->row_block_size, A_p->num_cols*A_p->col_block_size,A_p->len);

  if (A_p->type & MATRIX_FORMAT_CSC) {
      for (iCol = 0; iCol< A_p->pattern->n_ptr; iCol++) {
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
  } else { 
      for (iRow = 0; iRow< A_p->pattern->n_ptr; iRow++) {
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
  }

  /* close the file */
  fclose(fileHandle_p);
  
  return;
}

/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:39  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.1  2005/09/05 06:29:48  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
