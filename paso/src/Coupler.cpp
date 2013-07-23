
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


#include "Coupler.h"
#include "esysUtils/error.h"

/***********************************************************************************
 *
 * allocates a Connector
 *
 ************************************************************************************/

Paso_Connector* Paso_Connector_alloc(Paso_SharedComponents* send,
                                     Paso_SharedComponents* recv)
{
  Paso_Connector*out=NULL;
  Esys_resetError();
  out=new Paso_Connector;
  if ( send->mpi_info != recv->mpi_info ) {
     Esys_setError(SYSTEM_ERROR,"Paso_Coupler_alloc: send and recv mpi communicators don't match.");
     return NULL;
  }
  if ( send->local_length != recv->local_length ) {
     Esys_setError(SYSTEM_ERROR,"Paso_Coupler_alloc: local length of send and recv Paso_SharedComponents must match.");
     return NULL;
  }
  
  if (!Esys_checkPtr(out)) {
      out->send=Paso_SharedComponents_getReference(send);
      out->recv= Paso_SharedComponents_getReference(recv);
      out->mpi_info = Esys_MPIInfo_getReference(send->mpi_info);
      out->reference_counter=1;

/*
{ int i;
for (i=0; i< out->recv->numNeighbors; ++i) 
   printf("Coupler: %d receive %d data at %d from %d\n",send->mpi_info->rank,out->recv->offsetInShared[i+1]- out->recv->offsetInShared[i],
out->recv->offsetInShared[i],out->recv->neighbor[i]);
for (i=0; i< out->send->numNeighbors; ++i) 
printf("Coupler: %d send %d data at %d to %d\n",send->mpi_info->rank,out->send->offsetInShared[i+1]- out->send->offsetInShared[i],
out->send->offsetInShared[i],out->send->neighbor[i]);
}
*/

  }
  if (Esys_noError()) {
     return out;
  } else {
     Paso_Connector_free(out);
     return NULL;
  }
}

/* returns a reference to Connector */

Paso_Connector* Paso_Connector_getReference(Paso_Connector* in) {
     if (in!=NULL) {
        ++(in->reference_counter);
     }
     return in;
}
  
/* deallocates a Connector: */

void Paso_Connector_free(Paso_Connector* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<=0) {
        Paso_SharedComponents_free(in->send);
        Paso_SharedComponents_free(in->recv);
        Esys_MPIInfo_free(in->mpi_info);
        delete in;
        #ifdef Paso_TRACE
        printf("Paso_Coupler_dealloc: system matrix pattern as been deallocated.\n");
        #endif
     }
   }
}

Paso_Connector* Paso_Connector_copy(Paso_Connector* in) {
  return Paso_Connector_unroll(in,1);
}

Paso_Connector* Paso_Connector_unroll(Paso_Connector* in, index_t block_size) {
     Paso_SharedComponents *new_send_shcomp=NULL, *new_recv_shcomp=NULL;
     Paso_Connector *out=NULL;
     if (Esys_noError()) {
        if (block_size>1) {
            new_send_shcomp=Paso_SharedComponents_alloc(in->send->local_length,
                                                        in->send->numNeighbors,
                                                        in->send->neighbor,
                                                        in->send->shared,
                                                        in->send->offsetInShared,
                                                        block_size,0,in->mpi_info);

            new_recv_shcomp=Paso_SharedComponents_alloc(in->recv->local_length,
                                                        in->recv->numNeighbors,
                                                        in->recv->neighbor,
                                                        in->recv->shared,
                                                        in->recv->offsetInShared,
                                                        block_size,0,in->mpi_info);
        } else {
            new_send_shcomp=Paso_SharedComponents_getReference(in->send);
            new_recv_shcomp=Paso_SharedComponents_getReference(in->recv);
        }
        if (Esys_noError()) out=Paso_Connector_alloc(new_send_shcomp,new_recv_shcomp);
     }
     Paso_SharedComponents_free(new_send_shcomp);
     Paso_SharedComponents_free(new_recv_shcomp);
     if (Esys_noError()) {
          return out;
     } else {
          Paso_Connector_free(out);
          return NULL;
     } 
}
/***********************************************************************************
 *
 * allocates a Connector
 *
 ************************************************************************************/

Paso_Coupler* Paso_Coupler_alloc(Paso_Connector* connector, dim_t block_size)
{
  Esys_MPIInfo *mpi_info = connector->mpi_info;  
  Paso_Coupler*out=NULL;
  Esys_resetError();
  out=new Paso_Coupler;
  if (!Esys_checkPtr(out)) {
      out->data=NULL;
      out->block_size=block_size;
      out->connector=Paso_Connector_getReference(connector);
      out->send_buffer=NULL;
      out->recv_buffer=NULL;
      out->mpi_requests=NULL;
      out->mpi_stati=NULL;
      out->mpi_info = Esys_MPIInfo_getReference(mpi_info);
      out->reference_counter=1;
      out->in_use = FALSE;
      
      #ifdef ESYS_MPI
         out->mpi_requests=new MPI_Request[connector->send->numNeighbors+connector->recv->numNeighbors];
         out->mpi_stati=new MPI_Status[connector->send->numNeighbors+connector->recv->numNeighbors];
         Esys_checkPtr(out->mpi_requests);
         Esys_checkPtr(out->mpi_stati);
      #endif
      if (mpi_info->size>1) {
        out->send_buffer=new double[connector->send->numSharedComponents * block_size];
        out->recv_buffer=new double[connector->recv->numSharedComponents * block_size];
        Esys_checkPtr(out->send_buffer);
        Esys_checkPtr(out->recv_buffer);
      }
  }
  if (Esys_noError()) {
     return out;
  } else {
     Paso_Coupler_free(out);
     return NULL;
  }
}

/* returns a reference to Coupler */

Paso_Coupler* Paso_Coupler_getReference(Paso_Coupler* in) {
     if (in!=NULL) {
        ++(in->reference_counter);
     }
     return in;
}
  
/* deallocates a Coupler: */

void Paso_Coupler_free(Paso_Coupler* in) {
  if (in!=NULL) {
     (in->reference_counter)--;
     if (in->reference_counter<=0) {
        Paso_Connector_free(in->connector);
        delete[] in->send_buffer;
        delete[] in->recv_buffer;
  #ifdef ESYS_MPI
        delete[] in->mpi_requests;
        delete[] in->mpi_stati;
  #endif		
        Esys_MPIInfo_free(in->mpi_info);
        delete in;
     }
   }
}


void Paso_Coupler_startCollect(Paso_Coupler* coupler,const double* in)
{
  Esys_MPIInfo *mpi_info = coupler->mpi_info;  
  dim_t block_size=coupler->block_size;
  size_t block_size_size=block_size*sizeof(double);
  dim_t i;
  coupler->data=(double*) in;
  if ( mpi_info->size>1) {
     if (coupler->in_use) {
	Esys_setError(SYSTEM_ERROR,"Paso_Coupler_startCollect: Coupler in use.");
     }
     /* start receiving input */
     {
        for (i=0; i< coupler->connector->recv->numNeighbors; ++i) {
            #ifdef ESYS_MPI
            MPI_Irecv(&(coupler->recv_buffer[coupler->connector->recv->offsetInShared[i] *  block_size]),
                      (coupler->connector->recv->offsetInShared[i+1]- coupler->connector->recv->offsetInShared[i])*block_size,
                      MPI_DOUBLE,
                      coupler->connector->recv->neighbor[i], 
                      mpi_info->msg_tag_counter+coupler->connector->recv->neighbor[i],
                      mpi_info->comm,
                      &(coupler->mpi_requests[i]));
            #endif

        }
     }
     /* collect values into buffer */
     if (block_size>1) {
        #pragma omp parallel for private(i)
        for (i=0; i < coupler->connector->send->numSharedComponents;++i) {
           memcpy(&(coupler->send_buffer[(block_size)*i]),&(in[ block_size * coupler->connector->send->shared[i]]), block_size_size);
        }
     } else {
        #pragma omp parallel for private(i)
        for (i=0; i < coupler->connector->send->numSharedComponents;++i) {
	   coupler->send_buffer[i]=in[coupler->connector->send->shared[i]];
	}
     }
     /* send buffer out */
     {
        for (i=0; i< coupler->connector->send->numNeighbors; ++i) {
             #ifdef ESYS_MPI
             MPI_Issend(&(coupler->send_buffer[coupler->connector->send->offsetInShared[i] *  block_size]),
                        (coupler->connector->send->offsetInShared[i+1]- coupler->connector->send->offsetInShared[i])*block_size,
                        MPI_DOUBLE,
                        coupler->connector->send->neighbor[i], 
                        mpi_info->msg_tag_counter+mpi_info->rank,
                        mpi_info->comm,
                        &(coupler->mpi_requests[i+ coupler->connector->recv->numNeighbors]));
             #endif 
        }
     }
     ESYS_MPI_INC_COUNTER(*mpi_info,mpi_info->size)
     coupler->in_use=TRUE;
  }
}

double* Paso_Coupler_finishCollect(Paso_Coupler* coupler)
{
  Esys_MPIInfo *mpi_info = coupler->mpi_info;  
  if ( mpi_info->size>1) {
     if (! coupler->in_use) {
	Esys_setError(SYSTEM_ERROR,"Paso_Coupler_finishCollect: Communication has not been initiated.");
	return NULL;
     }
     /* wait for receive */
     #ifdef ESYS_MPI
        MPI_Waitall(coupler->connector->recv->numNeighbors+coupler->connector->send->numNeighbors,
                    coupler->mpi_requests,
                    coupler->mpi_stati);
      #endif
      coupler->in_use=FALSE;
  }

  return coupler->recv_buffer;
}

void Paso_Coupler_copyAll(const Paso_Coupler* src, Paso_Coupler* target) 
{
   dim_t i;
   #pragma omp parallel 
   {
       #pragma omp for private(i)
       for (i =0; i< src->connector->recv->numSharedComponents * src->block_size; ++i) {
          target->recv_buffer[i]=src->recv_buffer[i];
      }
      #pragma omp for private(i)
      for (i =0; i< Paso_Coupler_getLocalLength(src) * src->block_size; ++i) {
          target->data[i]=src->data[i];
     }
  }
}

/*  */ 
void Paso_Coupler_fillOverlap(const dim_t n, double* x, Paso_Coupler *coupler)
{
   double *remote_values = NULL;
   const dim_t overlap_n = Paso_Coupler_getNumOverlapValues(coupler) ;
   const dim_t my_n= n - overlap_n;
   const dim_t block_size = coupler->block_size;
   const dim_t offset = block_size * my_n;
   dim_t i;
   
    Paso_Coupler_startCollect(coupler, x);
    Paso_Coupler_finishCollect(coupler);
    remote_values=coupler->recv_buffer;
      
    #pragma omp parallel for private(i)
    for (i=0;i<overlap_n * block_size; ++i) {
         x[offset+i]=remote_values[i];
    } 
}

/* adjusts max values across shared values x */
void Paso_Coupler_max(const dim_t n, double* x, Paso_Coupler *coupler)
{
   double *remote_values = NULL;
   const dim_t overlap_n = Paso_Coupler_getNumOverlapValues(coupler) ;
   const dim_t my_n= n - overlap_n;
   dim_t i;
   
   Paso_Coupler_startCollect(coupler, x);
   Paso_Coupler_finishCollect(coupler);
   remote_values=coupler->recv_buffer;
   #pragma omp parallel for private(i)
   for (i=0;i<overlap_n; ++i) x[my_n+i]=MAX(x[my_n+i], remote_values[i]);
}

