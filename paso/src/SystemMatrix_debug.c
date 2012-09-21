
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
* 
*****************************************************
************************************************************************************

 Paso: SystemMatrix: debugging tools

************************************************************************************

Author: Lutz Gross, l.gross@uq.edu.au 

************************************************************************************/

#include "SystemMatrix.h"
#include "esysUtils/error.h"

/************************************************************************************/

/* fills the matrix with values i+f1*j where i and j are the global row
 * and column indices of the matrix entry */
void Paso_SystemMatrix_fillWithGlobalCoordinates(Paso_SystemMatrix *A, const double f1)
{
   dim_t ib, iPtr, q, i,p;
   const dim_t n=Paso_SystemMatrix_getNumRows(A);
   const dim_t m=Paso_SystemMatrix_getNumCols(A);
   const index_t me = A->mpi_info->rank;
   const dim_t block_size=A->block_size;
   const index_t row_offset = A->row_distribution->first_component[me];
   const index_t col_offset = A->col_distribution->first_component[me];
   Paso_Coupler* col_coupler=NULL, *row_coupler=NULL;
   double *cols=NULL, *rows=NULL; 
   
   cols=TMPMEMALLOC(m, double);
   rows=TMPMEMALLOC(n, double);
   col_coupler= Paso_Coupler_alloc(A->col_coupler->connector, 1);
   row_coupler = Paso_Coupler_alloc(A->col_coupler->connector, 1);
   
   #pragma omp parallel for private(i)
   for (i=0; i<n; ++i) rows[i]=row_offset+i;
   Paso_Coupler_startCollect(col_coupler, rows);
   
   #pragma omp parallel for private(i)
   for (i=0; i<m; ++i) cols[i]=col_offset+i;
   Paso_Coupler_startCollect(row_coupler, cols);
   
   

   /* main block : */
   for (q=0; q< n; ++q){
      for (iPtr =A->mainBlock->pattern->ptr[q]; iPtr<A->mainBlock->pattern->ptr[q+1]; ++iPtr) {
	 p =A->mainBlock->pattern->index[iPtr];
	 for (ib=0; ib<block_size; ib++) A->mainBlock->val[iPtr*block_size+ib]=f1*rows[q]+cols[p];
      }
   }

   Paso_Coupler_finishCollect(col_coupler);
   if ( A->col_coupleBlock != NULL ){
      for (q=0; q< A->col_coupleBlock->pattern->numOutput; ++q){
	 for (iPtr =A->col_coupleBlock->pattern->ptr[q]; iPtr<A->col_coupleBlock->pattern->ptr[q+1]; ++iPtr) {
	    p =A->col_coupleBlock->pattern->index[iPtr];
	    for (ib=0; ib<block_size; ib++) A->col_coupleBlock->val[iPtr*block_size+ib]=f1*rows[q]+col_coupler->recv_buffer[p];
	 }
      }
   }
   
   Paso_Coupler_finishCollect(row_coupler);
   if ( A->row_coupleBlock != NULL ){
      for (p=0; p< A->row_coupleBlock->pattern->numOutput; ++p){
	 for (iPtr =A->row_coupleBlock->pattern->ptr[p]; iPtr<A->row_coupleBlock->pattern->ptr[p+1]; ++iPtr) {
	    q =A->row_coupleBlock->pattern->index[iPtr];
	    for (ib=0; ib<block_size; ib++) A->row_coupleBlock->val[iPtr*block_size+ib]=row_coupler->recv_buffer[p]*f1+cols[q];
	 }
      }      
   }
   
   TMPMEMFREE(cols);
   TMPMEMFREE(rows);
   Paso_Coupler_free(row_coupler);
   Paso_Coupler_free(col_coupler);
   return; 		      
}

void Paso_SystemMatrix_print(Paso_SystemMatrix *A)
{
   dim_t iPtr, q, p, ib;
   const dim_t n=Paso_SystemMatrix_getNumRows(A);
   const dim_t block_size=A->block_size;
   index_t rank=A->mpi_info->rank;
   char *str1, *str2;
   str1 = TMPMEMALLOC(n*n*block_size*30+100, char);
   str2 = TMPMEMALLOC(30, char);
   
   sprintf(str1, "rank %d Main Block:\n-----------\n", rank);
   for (q=0; q< n; ++q){
      sprintf(str2, "Row %d: ",q);
      strcat(str1, str2);
      for (iPtr =A->mainBlock->pattern->ptr[q]; iPtr<A->mainBlock->pattern->ptr[q+1]; ++iPtr) {
	 sprintf(str2, "(%d ",A->mainBlock->pattern->index[iPtr]);
	 strcat(str1, str2);
	 for (ib=0; ib<block_size; ib++){
	   sprintf(str2, "%f ", A->mainBlock->val[iPtr*block_size+ib]);
	   strcat(str1, str2);
	 }
	 sprintf(str2, "),");
	 strcat(str1, str2);
      }
      sprintf(str1, "%s\n", str1);
   }
   fprintf(stderr, "%s", str1);
   if ( ( A->col_coupleBlock != NULL) && (A->mpi_info->size>1) ) {
      sprintf(str1, "rank %d Column Couple Block:\n------------------\n", rank);
      for (q=0; q< A->col_coupleBlock->pattern->numOutput; ++q){
	 sprintf(str2, "Row %d: ", q);
	 strcat(str1, str2);
	 for (iPtr =A->col_coupleBlock->pattern->ptr[q]; iPtr<A->col_coupleBlock->pattern->ptr[q+1]; ++iPtr) {
	    if (A->global_id) 
	      sprintf(str2, "(%d %f),",A->global_id[A->col_coupleBlock->pattern->index[iPtr]], A->col_coupleBlock->val[iPtr*block_size]);
	    else 
	      sprintf(str2, "(%d %f),",A->col_coupleBlock->pattern->index[iPtr], A->col_coupleBlock->val[iPtr*block_size]);
	    strcat(str1, str2);
	 }
	 sprintf(str1, "%s\n", str1);
      }
      fprintf(stderr, "%s", str1);
   }
   if ( ( A->row_coupleBlock != NULL ) && (A->mpi_info->size>1) ) {
      sprintf(str1, "rank %d Row Couple Block:\n--------------------\n", rank);
      for (p=0; p< A->row_coupleBlock->pattern->numOutput; ++p){
         sprintf(str2, "Row %d:", p);
         strcat(str1, str2);
         for (iPtr =A->row_coupleBlock->pattern->ptr[p]; iPtr<A->row_coupleBlock->pattern->ptr[p+1]; ++iPtr) {
            sprintf(str2, "(%d %f),",A->row_coupleBlock->pattern->index[iPtr], A->row_coupleBlock->val[iPtr*block_size]);
            strcat(str1, str2);
         }
         sprintf(str1, "%s\n", str1);
      }
      fprintf(stderr, "%s", str1);
   }
   if ( ( A->remote_coupleBlock != NULL ) && (A->mpi_info->size>1) ) {
      sprintf(str1, "rank %d Remote Couple Block:\n--------------------\n", rank);
      for (p=0; p< A->remote_coupleBlock->pattern->numOutput; ++p){
         sprintf(str2, "Row %d:", p);
         strcat(str1, str2);
         for (iPtr =A->remote_coupleBlock->pattern->ptr[p]; iPtr<A->remote_coupleBlock->pattern->ptr[p+1]; ++iPtr) {
            sprintf(str2, "(%d %f),",A->remote_coupleBlock->pattern->index[iPtr], A->remote_coupleBlock->val[iPtr*block_size]);
            strcat(str1, str2);
         }
         sprintf(str1, "%s\n", str1);
      }
      fprintf(stderr, "%s", str1);
   }
   TMPMEMFREE(str1);
   TMPMEMFREE(str2);
   return; 
}

