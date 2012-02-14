
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


/**************************************************************

 Paso: SystemMatrix: copies the col_coupleBlock into 
                     row_coupleBlock. 
                     
  WARNING: function uses mpi_requests of the coupler attached to A.
  
                 No reordering on the columns received is performed.
                 In practice this means that components in
                    A->row_coupleBlock->pattern->index
                 and 
                    A->row_coupler->connector->recv->shared
                 are ordered by increasing value.
		     
                 Notice: that send and receive A->row_coupler->connectors
                 are swapping roles.

**************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "SystemMatrix.h"
#include "esysUtils/error.h"

/**************************************************************/

void Paso_SystemMatrix_copyColCoupleBlock(Paso_SystemMatrix *A)
{
   dim_t p;
   index_t z0, iPtr, rPtr;
   const dim_t block_size=A->block_size;
   const size_t block_size_size=block_size*sizeof(double);
   double * send_buffer = NULL;
   
   if ( A->mpi_info->size > 1 ) {
      if ( A->row_coupleBlock == NULL ) {
	 Esys_setError(VALUE_ERROR, "SystemMatrix_copyColCoupleBlock: creation of row_coupleBlock pattern not supported yet.");
	 return;
      }
      if ( A->row_coupler->in_use ) {
	 Esys_setError(SYSTEM_ERROR,"SystemMatrix_copyColCoupleBlock: Coupler in use.");
	 return;
      }

      /* start receiving */
      for (p=0; p<A->row_coupler->connector->recv->numNeighbors; p++) {
	    #ifdef ESYS_MPI
	    const index_t irow1= A->row_coupler->connector->recv->offsetInShared[p];
	    const index_t irow2= A->row_coupler->connector->recv->offsetInShared[p+1];
	    const index_t a = A->row_coupleBlock->pattern->ptr[irow1];
	    const index_t b = A->row_coupleBlock->pattern->ptr[irow2];
	     
	    MPI_Irecv(&(A->row_coupleBlock->val[a * block_size]), (b-a) * block_size,  MPI_DOUBLE,
	       A->row_coupler->connector->recv->neighbor[p], 
	       A->mpi_info->msg_tag_counter+A->row_coupler->connector->recv->neighbor[p],
	       A->mpi_info->comm,
	       &(A->row_coupler->mpi_requests[p]) );

	     #endif
      }
      /* start sending */
      z0=0;
      send_buffer = TMPMEMALLOC(A->col_coupleBlock->len, double);
	
      for (p=0; p<A->row_coupler->connector->send->numNeighbors; p++) {
	 /* j_min, j_max defines the range of columns to be sent to processor p*/
	 const index_t j_min = A->col_coupler->connector->recv->offsetInShared[p];
	 const index_t j_max = A->col_coupler->connector->recv->offsetInShared[p+1];
	 index_t z = z0;
	 
	 /* run over the rows to be connected to processor p */
	 for (rPtr= A->row_coupler->connector->send->offsetInShared[p];  rPtr< A->row_coupler->connector->send->offsetInShared[p+1]; ++rPtr) {
	    const index_t row=A->row_coupler->connector->send->shared[rPtr];
    
	    /* collect the entries in the col. couple block referring to columns on processor p */
	    for (iPtr =A->col_coupleBlock->pattern->ptr[row]; iPtr<A->col_coupleBlock->pattern->ptr[row+1]; ++iPtr) {
	       register index_t j =A->col_coupleBlock->pattern->index[iPtr];
	       if ( (j_min <= j) && (j < j_max) ) {
		  memcpy(&(send_buffer[z]),&(A->col_coupleBlock->val[ block_size * iPtr]), block_size_size);
		  z+=block_size;
	       }
	    }
	    
	 }
	 #ifdef ESYS_MPI
	 MPI_Issend(&(send_buffer[z0]),z-z0, MPI_DOUBLE,
		       A->row_coupler->connector->send->neighbor[p], 
		       A->mpi_info->msg_tag_counter+A->mpi_info->rank,
		       A->mpi_info->comm,
		       &(A->row_coupler->mpi_requests[p+A->row_coupler->connector->recv->numNeighbors]));
	 
	 #endif
         z0 = z;
      }
   
	 
        /* wait until everything is done */
	 #ifdef ESYS_MPI
	 MPI_Waitall(A->row_coupler->connector->send->numNeighbors+A->row_coupler->connector->recv->numNeighbors,
		     A->row_coupler->mpi_requests,
		     A->row_coupler->mpi_stati);
         #endif
         A->mpi_info->msg_tag_counter+=A->mpi_info->size;
         TMPMEMFREE(send_buffer);
   }
   return; 		      
}
