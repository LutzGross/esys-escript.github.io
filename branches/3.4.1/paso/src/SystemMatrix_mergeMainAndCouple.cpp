
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
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


/**********************************************************************************/
/* Paso: SystemMatrix                                       */
/*							    */
/*  Merge the MainBlock and CoupleBlock in the matrix       */
/*  Input: SystemMatrix A				    */
/*  Output: 						    */
/*	p_ptr: the pointer to a vector of locations that    */
/*	       start a row.				    */
/*	p_idx: the pointer to the column indices for each   */
/*	       of the rows, ordered by rows.	            */
/*	p_val: the pointer to the data corresponding 	    */
/*	       directly to the column entries in p_idx.     */
/**********************************************************************************/

/* Copyrights by ACcESS Australia 2003 */
/* Author: Lin Gao, l.gao@uq.edu.au */

/**********************************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"

#ifdef _OPENMP
#include <omp.h>
#endif


void Paso_SystemMatrix_mergeMainAndCouple(Paso_SystemMatrix* A, index_t** p_ptr, index_t** p_idx, double** p_val){
  if (A->type & MATRIX_FORMAT_DEFAULT) {
      Paso_SystemMatrix_mergeMainAndCouple_CSR_OFFSET0(A, p_ptr, p_idx, p_val);
  } else if (A->type & MATRIX_FORMAT_CSC) {
      /* CSC part is for PASTIX */
      if (A->type & (MATRIX_FORMAT_OFFSET1 + MATRIX_FORMAT_BLK1)) {
          Paso_SystemMatrix_mergeMainAndCouple_CSC_OFFSET1(A, p_ptr, p_idx, p_val);
      } else {
          Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_mergeMainAndCouple: CSC with index 0 or block size larger than 1 is not supported.");
      }
  } else if (A->type & MATRIX_FORMAT_TRILINOS_CRS) {
      Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_mergeMainAndCouple: TRILINOS is not supported.");
  } else {
      Esys_setError(SYSTEM_ERROR,"Paso_SystemMatrix_mergeMainAndCouple: CRS is not supported.");
  }
}


void Paso_SystemMatrix_mergeMainAndCouple_CSR_OFFSET0(Paso_SystemMatrix* A, index_t** p_ptr, index_t** p_idx, double** p_val) {

  index_t i, j, i_ub, j_lb, j_ub, row, num_vals, main_num_vals;
  index_t couple_num_vals, rank, row_offset, ij_ptr=0, idx=0, idx2=0;
  index_t main_num_rows, couple_num_rows, col_offset, num_cols;
  index_t *main_ptr, *main_idx, *couple_ptr, *couple_idx, *global_id;
  double  *main_val, *couple_val, *rows=NULL;
  Paso_Coupler* coupler=NULL;

  if (A->mainBlock->col_block_size!=1 ||
      A->mainBlock->row_block_size!=1 ||
      A->col_coupleBlock->col_block_size!=1 ||
      A->col_coupleBlock->row_block_size!=1) {
      Paso_SystemMatrix_mergeMainAndCouple_CSR_OFFSET0_Block(A, p_ptr, p_idx, p_val);
      return;
  }

  if (A->mpi_info->size == 1) {
      /* initialisation */
      main_num_rows=A->mainBlock->numRows;
      main_ptr=A->mainBlock->pattern->ptr;
      main_idx=A->mainBlock->pattern->index;
      main_val=A->mainBlock->val;

     /* allocate arrays "ptr", "index" and "val" */
     num_vals = main_ptr[main_num_rows]-1;
     *p_ptr = new index_t[main_num_rows+1];
     *p_idx = new index_t[num_vals];
     *p_val = new double[num_vals];

     #pragma omp for schedule(static) private(i,ij_ptr,j_lb,j_ub)
     for (i=0; i<main_num_rows; i++) {
	j_lb = main_ptr[i];
	j_ub = main_ptr[i+1];
	(*p_ptr)[i] = j_lb;
	for (ij_ptr=j_lb; ij_ptr<j_ub; ++ij_ptr) {
	    (*p_idx)[ij_ptr] = main_idx[ij_ptr];
	    (*p_val)[ij_ptr] = main_val[ij_ptr];
	}
     }
     (*p_ptr)[main_num_rows] = main_ptr[main_num_rows];
     return;
  }

  main_num_rows=A->mainBlock->numRows;
  couple_num_rows=A->col_coupleBlock->numRows;
  rank = A->mpi_info->rank;

  if (main_num_rows != couple_num_rows) {
      Esys_setError(TYPE_ERROR,"Paso_SystemMatrix_mergeMainAndCouple_CSR_OFFSET0: number of rows do not match.");
      return;
  }

  if (A->global_id == NULL) {
    /* prepare for global coordinates in colCoupleBlock, the results are
       in coupler->recv_buffer */
    rows=new double[main_num_rows];
    row_offset = A->row_distribution->first_component[rank];
    #pragma omp parallel for private(i) schedule(static)
    for (i=0; i<main_num_rows; ++i) rows[i]=row_offset+i;
    coupler= Paso_Coupler_alloc(A->col_coupler->connector, 1);
    Paso_Coupler_startCollect(coupler, rows);
  }

  /* initialisation, including allocate arrays "ptr", "index" and "val" */
  main_ptr=A->mainBlock->pattern->ptr;
  main_idx=A->mainBlock->pattern->index;
  main_val=A->mainBlock->val;
  couple_ptr=A->col_coupleBlock->pattern->ptr;
  couple_idx=A->col_coupleBlock->pattern->index;
  couple_val=A->col_coupleBlock->val;
  col_offset = A->col_distribution->first_component[rank];
  main_num_vals = main_ptr[main_num_rows]-main_ptr[0];
  couple_num_vals = couple_ptr[couple_num_rows]-couple_ptr[0];
  num_vals = main_num_vals + couple_num_vals;
  *p_ptr = new index_t[main_num_rows+1];
  *p_idx = new index_t[num_vals];
  *p_val = new double[num_vals];
  (*p_ptr)[0] = 0;
  i = 0;
  j = 0;

  if (A->global_id == NULL) {
    Paso_Coupler_finishCollect(coupler);
    delete[] rows;
    num_cols = A->col_coupleBlock->numCols;
    global_id = new index_t[num_cols];
    #pragma omp parallel for private(i) schedule(static)
    for (i=0; i<num_cols; ++i)
        global_id[i] = coupler->recv_buffer[i];
    A->global_id = global_id;
    Paso_Coupler_free(coupler);
  }
  global_id = A->global_id;

  /* merge mainBlock and col_coupleBlock */
  for (row=1; row<=main_num_rows; row++) {
      i_ub = main_ptr[row];
      j_ub = couple_ptr[row];
      while (i < i_ub || j < j_ub) {
	  ij_ptr = i + j;
          if (j < j_ub) {
	      idx = global_id[couple_idx[j]];
	  }
	  if (i < i_ub) {
	      idx2 = main_idx[i] + col_offset;
	  }
          if (j == j_ub || (i < i_ub && idx2 < idx)){
              (*p_idx)[ij_ptr] = idx2;
              (*p_val)[ij_ptr] = main_val[i];
              i++;
          } else {
              (*p_idx)[ij_ptr] = idx;
              (*p_val)[ij_ptr] = couple_val[j];
              j++;
          }
      }
      (*p_ptr)[row] = ij_ptr+1;
  }

}


void Paso_SystemMatrix_mergeMainAndCouple_CSR_OFFSET0_Block(Paso_SystemMatrix* A, index_t** p_ptr, index_t** p_idx, double** p_val) {

  index_t i, j, i_ub, j_lb, j_ub, row, num_vals, main_num_vals, ib;
  index_t couple_num_vals, rank, row_offset, ij_ptr=0, idx=0, idx2=0;
  index_t main_num_rows, couple_num_rows, col_offset, num_cols, block_size;
  index_t *main_ptr, *main_idx, *couple_ptr, *couple_idx, *global_id;
  double  *main_val, *couple_val, *rows=NULL;
  Paso_Coupler* coupler=NULL;

  block_size = A->block_size;
  if (A->mpi_info->size == 1) {
      /* initialisation */
      main_num_rows=A->mainBlock->numRows;
      main_ptr=A->mainBlock->pattern->ptr;
      main_idx=A->mainBlock->pattern->index;
      main_val=A->mainBlock->val;

     /* allocate arrays "ptr", "index" and "val" */
     num_vals = main_ptr[main_num_rows]-1;
     *p_ptr = new index_t[main_num_rows+1];
     *p_idx = new index_t[num_vals];
     *p_val = new double[num_vals * block_size];

     #pragma omp for schedule(static) private(i,ij_ptr,j_lb,j_ub, ib)
     for (i=0; i<main_num_rows; i++) {
	j_lb = main_ptr[i];
	j_ub = main_ptr[i+1];
	(*p_ptr)[i] = j_lb;
	for (ij_ptr=j_lb; ij_ptr<j_ub; ij_ptr++) {
	    (*p_idx)[ij_ptr] = main_idx[ij_ptr];
	    for (ib=0; ib<block_size; ib++)
		(*p_val)[ij_ptr*block_size+ib] = main_val[ij_ptr*block_size+ib];
	}
     }
     (*p_ptr)[main_num_rows] = main_ptr[main_num_rows];
     return;
  }

  main_num_rows=A->mainBlock->numRows;
  couple_num_rows=A->col_coupleBlock->numRows;
  rank = A->mpi_info->rank;

  if (main_num_rows != couple_num_rows) {
      Esys_setError(TYPE_ERROR,"Paso_SystemMatrix_mergeMainAndCouple_CSR_OFFSET0: number of rows do not match.");
      return;
  }

  if (A->global_id == NULL) {
    /* prepare for global coordinates in colCoupleBlock, the results are
       in coupler->recv_buffer */
    rows=new double[main_num_rows];
    row_offset = A->row_distribution->first_component[rank];
    #pragma omp parallel for private(i) schedule(static)
    for (i=0; i<main_num_rows; ++i) rows[i]=row_offset+i;
    coupler= Paso_Coupler_alloc(A->col_coupler->connector, 1);
    Paso_Coupler_startCollect(coupler, rows);
  }

  /* initialisation, including allocate arrays "ptr", "index" and "val" */
  main_ptr=A->mainBlock->pattern->ptr;
  main_idx=A->mainBlock->pattern->index;
  main_val=A->mainBlock->val;
  couple_ptr=A->col_coupleBlock->pattern->ptr;
  couple_idx=A->col_coupleBlock->pattern->index;
  couple_val=A->col_coupleBlock->val;
  col_offset = A->col_distribution->first_component[rank];
  main_num_vals = main_ptr[main_num_rows]-main_ptr[0];
  couple_num_vals = couple_ptr[couple_num_rows]-couple_ptr[0];
  num_vals = main_num_vals + couple_num_vals;
  *p_ptr = new index_t[main_num_rows+1];
  *p_idx = new index_t[num_vals];
  *p_val = new double[num_vals*block_size];
  (*p_ptr)[0] = 0;
  i = 0;
  j = 0;

  if (A->global_id == NULL) {
    Paso_Coupler_finishCollect(coupler);
    delete[] rows;
    num_cols = A->col_coupleBlock->numCols;
    global_id = new index_t[num_cols];
    #pragma omp parallel for private(i) schedule(static)
    for (i=0; i<num_cols; ++i)
        global_id[i] = coupler->recv_buffer[i];
    A->global_id = global_id;
    Paso_Coupler_free(coupler);
  }
  global_id = A->global_id;

  /* merge mainBlock and col_coupleBlock */
  for (row=1; row<=main_num_rows; row++) {
      i_ub = main_ptr[row];
      j_ub = couple_ptr[row];
      while (i < i_ub || j < j_ub) {
	  ij_ptr = i + j;
          if (j < j_ub) {
	      idx = global_id[couple_idx[j]];
	  }
	  if (i < i_ub) {
	      idx2 = main_idx[i] + col_offset;
	  }
          if (j == j_ub || (i < i_ub && idx2 < idx)){
              (*p_idx)[ij_ptr] = idx2;
	      for (ib=0; ib<block_size; ib++)
		(*p_val)[ij_ptr*block_size+ib] = main_val[i*block_size+ib];
              i++;
          } else {
              (*p_idx)[ij_ptr] = idx;
	      for (ib=0; ib<block_size; ib++)
		(*p_val)[ij_ptr*block_size+ib] = couple_val[j*block_size+ib];
              j++;
          }
      }
      (*p_ptr)[row] = ij_ptr+1;
  }

}


void Paso_SystemMatrix_mergeMainAndCouple_CSC_OFFSET1(Paso_SystemMatrix* A, index_t** p_ptr, index_t** p_idx, double** p_val) {
/*
  index_t i, j, i_ub, j_ub, col, num_vals, main_num_vals;
  index_t couple_num_vals, idx, rank, main_offset;
  index_t main_num_cols=A->mainBlock->pattern->numOutput;
  index_t couple_num_cols=A->col_coupleBlock->pattern->numOutput;
  index_t *main_ptr=A->mainBlock->pattern->ptr;
  index_t *main_idx=A->mainBlock->pattern->index;
  double  *main_val=A->mainBlock->val;
  index_t *couple_ptr=A->col_coupleBlock->pattern->ptr;
  index_t *couple_idx=A->col_coupleBlock->pattern->index;
  double  *couple_val=A->col_coupleBlock->val;
  index_t *couple_global=NULL;
  Paso_Coupler* coupler=A->col_coupler;
*/
}


