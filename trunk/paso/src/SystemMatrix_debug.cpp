
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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
   double* cols=new double[m];
   double* rows=new double[n];
   paso::Coupler_ptr col_coupler(new paso::Coupler(A->col_coupler->connector, 1));
   paso::Coupler_ptr row_coupler(new paso::Coupler(A->col_coupler->connector, 1));
   
   #pragma omp parallel for private(i)
   for (i=0; i<n; ++i) rows[i]=row_offset+i;
   col_coupler->startCollect(rows);
   
   #pragma omp parallel for private(i)
   for (i=0; i<m; ++i) cols[i]=col_offset+i;
   row_coupler->startCollect(cols);
   
   /* main block : */
   for (q=0; q< n; ++q){
      for (iPtr =A->mainBlock->pattern->ptr[q]; iPtr<A->mainBlock->pattern->ptr[q+1]; ++iPtr) {
	 p =A->mainBlock->pattern->index[iPtr];
	 for (ib=0; ib<block_size; ib++) A->mainBlock->val[iPtr*block_size+ib]=f1*rows[q]+cols[p];
      }
   }

   col_coupler->finishCollect();
   if ( A->col_coupleBlock != NULL ){
      for (q=0; q< A->col_coupleBlock->pattern->numOutput; ++q){
	 for (iPtr =A->col_coupleBlock->pattern->ptr[q]; iPtr<A->col_coupleBlock->pattern->ptr[q+1]; ++iPtr) {
	    p =A->col_coupleBlock->pattern->index[iPtr];
	    for (ib=0; ib<block_size; ib++) A->col_coupleBlock->val[iPtr*block_size+ib]=f1*rows[q]+col_coupler->recv_buffer[p];
	 }
      }
   }
   
   row_coupler->finishCollect();
   if ( A->row_coupleBlock != NULL ){
      for (p=0; p< A->row_coupleBlock->pattern->numOutput; ++p){
	 for (iPtr =A->row_coupleBlock->pattern->ptr[p]; iPtr<A->row_coupleBlock->pattern->ptr[p+1]; ++iPtr) {
	    q =A->row_coupleBlock->pattern->index[iPtr];
	    for (ib=0; ib<block_size; ib++) A->row_coupleBlock->val[iPtr*block_size+ib]=row_coupler->recv_buffer[p]*f1+cols[q];
	 }
      }      
   }
   
   delete[] cols;
   delete[] rows;
   return; 		      
}

void Paso_SystemMatrix_print(Paso_SystemMatrix *A)
{
   dim_t iPtr, q, p, ib;
   const dim_t n=Paso_SystemMatrix_getNumRows(A);
   const dim_t block_size=A->block_size;
   index_t rank=A->mpi_info->rank;
   char *str1, *str2;
   str1 = new char[n*n*block_size*30+100];
   str2 = new char[30];
   
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
   delete[] str1;
   delete[] str2;
   return; 
}

