/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrix is saved to Matrix Market format      */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004 */
/* Author: davies@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "mmio.h"
#include "System.h"

void Finley_SystemMatrix_saveMM(Finley_SystemMatrix * A_p, char * fileName_p) {

  int iRow, iCol, iPtr;

  /* open the file */
  FILE * fileHandle_p = fopen(fileName_p, "w");
  if (fileHandle_p==NULL) {
    Finley_ErrorCode =IO_ERROR;
    sprintf(Finley_ErrorMsg, "file %s could not be opened for writing",fileName_p);
    return;
  }

  /* set the matrix code */
  MM_typecode matrixCode;
  mm_initialize_typecode(&matrixCode);
  mm_set_matrix(&matrixCode);
  mm_set_coordinate(&matrixCode);
  mm_set_real(&matrixCode);
  
  mm_write_banner(fileHandle_p, matrixCode);
  mm_write_mtx_crd_size(fileHandle_p, A_p->num_rows, A_p->num_cols, A_p->ptr[A_p->num_cols] - PTR_OFFSET);

  switch(A_p->type) {
  case CSR:
    #if INDEX_OFFSET==0
      /* add 1 for 1-based indexing */
      for (iRow = 0; iRow < A_p->num_rows; iRow++)
	for (iPtr = A_p->ptr[iRow] - PTR_OFFSET;
	     iPtr < A_p->ptr[iRow+1] - PTR_OFFSET; iPtr++)
	  fprintf(fileHandle_p, "%12d %12d %22.15g\n",
		  iRow + 1,  A_p->index[iPtr] + 1, A_p->val[iPtr]);
    #else 
      /* matrix market uses 1-based indexing */
      for (iRow = 0; iRow< A_p->num_rows; iRow++)
	for (iPtr = A_p->ptr[iRow] - PTR_OFFSET;
	     iPtr < A_p->ptr[iRow+1] - PTR_OFFSET; iPtr++)
	  fprintf(fileHandle_p, "%12d %12d %22.15g\n",
		  iRow + INDEX_OFFSET, A_p->index[iPtr], A_p->val[iPtr]);
    #endif
    break;
  case CSC:
    #if INDEX_OFFSET==0
      /* add 1 for 1-based indexing */
      for (iCol = 0; iCol < A_p->num_cols; iCol++)
	for (iPtr = A_p->ptr[iCol] - PTR_OFFSET;
	     iPtr < A_p->ptr[iCol+1] - PTR_OFFSET; iPtr++)
	  fprintf(fileHandle_p, "%12d %12d %22.15g\n", A_p->index[iPtr] + 1, iCol + 1, A_p->val[iPtr]);
    #else
      /* matrix market uses 1-based indexing */
      for (iCol = 0; iCol< A_p->num_cols; iCol++)
	for (iPtr = A_p->ptr[iCol] - PTR_OFFSET;
	     iPtr < A_p->ptr[iCol+1] - PTR_OFFSET; iPtr++)
	  fprintf(fileHandle_p, "%12d %12d %22.15g\n", A_p->index[iPtr],
		  iCol + INDEX_OFFSET, A_p->val[iPtr]);
    #endif
  }

  /* close the file */
  fclose(fileHandle_p);
  
  return;
}
/*
 * $Log$
 * Revision 1.3  2004/12/15 03:48:47  jgs
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
