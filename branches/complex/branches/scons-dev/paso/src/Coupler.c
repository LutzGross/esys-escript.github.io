/* $Id: Coupler.c 1306 2007-09-18 05:51:09Z ksteube $ */

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
 **************************************************************
 *  
 * Paso: Connector and  Coupler organizes the coupling with in a pattern/matrix  
 *       across processors                                        
 *
 **************************************************************
 *
 * Author: gross@access.edu.au 
 *
 **************************************************************/

#include "Coupler.h"

/*************************************************************
 *
 * allocates a Connector
 *
 **************************************************************/

Paso_Connector* Paso_Connector_alloc(Paso_SharedComponents* send,
                                     Paso_SharedComponents* recv)
{
  Paso_Connector*out=NULL;
  Paso_resetError();
  out=MEMALLOC(1,Paso_Connector);
  if ( send->mpi_info != recv->mpi_info ) {
     Paso_setError(SYSTEM_ERROR,"Paso_Coupler_alloc: send and recv mpi communicator don't match.");
     return NULL;
  }
  if ( send->local_length != recv->local_length ) {
     Paso_setError(SYSTEM_ERROR,"Paso_Coupler_alloc: local length of send and recv Paso_SharedComponents must match.");
     return NULL;
  }
  
  if (!Paso_checkPtr(out)) {
      out->send=Paso_SharedComponents_getReference(send);
      out->recv= Paso_SharedComponents_getReference(recv);
      out->mpi_info = Paso_MPIInfo_getReference(send->mpi_info);
      out->reference_counter=1;
  }
  if (Paso_noError()) {
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
        Paso_MPIInfo_free(in->mpi_info);
        MEMFREE(in);
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
     if (Paso_noError()) {
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
        if (Paso_noError()) out=Paso_Connector_alloc(new_send_shcomp,new_recv_shcomp);
     }
     Paso_SharedComponents_free(new_send_shcomp);
     Paso_SharedComponents_free(new_recv_shcomp);
     if (Paso_noError()) {
          return out;
     } else {
          Paso_Connector_free(out);
          return NULL;
     } 
}
/*************************************************************
 *
 * allocates a Connector
 *
 **************************************************************/

Paso_Coupler* Paso_Coupler_alloc(Paso_Connector* connector, dim_t block_size)
{
  Paso_MPIInfo *mpi_info = connector->mpi_info;  
  Paso_Coupler*out=NULL;
  Paso_resetError();
  out=MEMALLOC(1,Paso_Coupler);
  if (!Paso_checkPtr(out)) {
      out->data=NULL;
      out->block_size=block_size;
      out->connector=Paso_Connector_getReference(connector);
      out->send_buffer=NULL;
      out->recv_buffer=NULL;
      out->mpi_requests=NULL;
      out->mpi_stati=NULL;
      out->mpi_info = Paso_MPIInfo_getReference(mpi_info);
      out->reference_counter=1;
      
      #ifdef PASO_MPI
         out->mpi_requests=MEMALLOC(connector->send->numNeighbors+connector->recv->numNeighbors,MPI_Request);
         out->mpi_stati=MEMALLOC(connector->send->numNeighbors+connector->recv->numNeighbors,MPI_Status);
         Paso_checkPtr(out->mpi_requests);
         Paso_checkPtr(out->mpi_stati);
      #endif
      if (mpi_info->size>1) {
        out->send_buffer=MEMALLOC(connector->send->numSharedComponents * block_size,double);
        out->recv_buffer=MEMALLOC(connector->recv->numSharedComponents * block_size,double);
        Paso_checkPtr(out->send_buffer);
        Paso_checkPtr(out->recv_buffer);
      }
  }
  if (Paso_noError()) {
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
     in->reference_counter--;
     if (in->reference_counter<=0) {
        Paso_Connector_free(in->connector);
        MEMFREE(in->send_buffer);
        MEMFREE(in->recv_buffer);
        MEMFREE(in->mpi_requests);
        MEMFREE(in->mpi_stati);
        Paso_MPIInfo_free(in->mpi_info);
        MEMFREE(in);
        #ifdef Paso_TRACE
        printf("Paso_Coupler_dealloc: system matrix pattern as been deallocated.\n");
        #endif
     }
   }
}


void Paso_Coupler_startCollect(Paso_Coupler* coupler,const double* in)
{
  Paso_MPIInfo *mpi_info = coupler->mpi_info;  
  dim_t block_size=coupler->block_size;
  size_t block_size_size=block_size*sizeof(double);
  dim_t i;
  coupler->data=(double*) in;
  if ( mpi_info->size>1) {
     /* start reveiving input */
     {
        for (i=0; i< coupler->connector->recv->numNeighbors; ++i) {
            #ifdef PASO_MPI
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
     #pragma omp parallel for private(i)
     for (i=0; i < coupler->connector->send->numSharedComponents;++i) {
        memcpy(&(coupler->send_buffer[(block_size)*i]),&(in[ block_size * coupler->connector->send->shared[i]]), block_size_size);
     }
     /* send buffer out */
     {
        for (i=0; i< coupler->connector->send->numNeighbors; ++i) {
             #ifdef PASO_MPI
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
     mpi_info->msg_tag_counter+=mpi_info->size;
  }
}

double* Paso_Coupler_finishCollect(Paso_Coupler* coupler)
{
  Paso_MPIInfo *mpi_info = coupler->mpi_info;  
  if ( mpi_info->size>1) {
     /* wait for receive */
        #ifdef PASO_MPI
        MPI_Waitall(coupler->connector->recv->numNeighbors+coupler->connector->send->numNeighbors,
                    coupler->mpi_requests,
                    coupler->mpi_stati);
        #endif
  }
  return coupler->recv_buffer;
}
dim_t Paso_Coupler_getLocalLength(const Paso_Coupler* in) {
     return in->connector->send->local_length;

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
