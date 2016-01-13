
/* $Id$ */

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
 *******************************************************/

/**************************************************************/

/*   Finley: Mesh: NodeFile                                   */

/*   gathers the NodeFile out from the NodeFile in using the entries 
/*   in index[0:out->numNodes-1] which are between min_index and max_index (exclusive) */
/*   the node index[i]

/**************************************************************/

#include "NodeFile.h"

/**************************************************************/

void Finley_NodeFile_gatherEntries(dim_t n, index_t* index, index_t min_index, index_t max_index,
                                   index_t* Id_out, index_t* Id_in, 
                                   index_t* Tag_out, index_t* Tag_in, 
                                   index_t* globalDegreesOfFreedom_out, index_t* globalDegreesOfFreedom_in, 
                                   dim_t numDim, double* Coordinates_out, double* Coordinates_in)
{
   dim_t i;
   register index_t k;
   register const index_t range=max_index-min_index;
   const  size_t numDim_size=(size_t)numDim*sizeof(double);

   #pragma omp parallel for private(i,k) schedule(static)
   for (i=0;i<n;i++) {
      k=index[i]-min_index;
      if ((k>=0) && (k <range)) {
         Id_out[i]=Id_in[k];
         Tag_out[i]=Tag_in[k];
         globalDegreesOfFreedom_out[i]=globalDegreesOfFreedom_in[k];
         memcpy(&(Coordinates_out[INDEX2(0,i,numDim)]), &(Coordinates_in[INDEX2(0,k,numDim)]), numDim_size);
      }
   }
}

void Finley_NodeFile_gather(index_t* index, Finley_NodeFile* in, Finley_NodeFile* out) 
{
   index_t min_id, max_id;
   Finley_NodeFile_setGlobalIdRange(&min_id,&max_id,in);
   Finley_NodeFile_gatherEntries(out->numNodes, index, min_id, max_id,
                                 out->Id, in->Id, 
                                 out->Tag, in->Tag, 
                                 out->globalDegreesOfFreedom, in->globalDegreesOfFreedom, 
                                 out->numDim, out->Coordinates, in->Coordinates);
}

void Finley_NodeFile_gather_global(index_t* index, Finley_NodeFile* in, Finley_NodeFile* out)
{
  index_t min_id, max_id, undefined_node;
  Paso_MPI_rank buffer_rank, dest, source, *distribution=NULL;
  index_t  *Id_buffer=NULL, *Tag_buffer=NULL, *globalDegreesOfFreedom_buffer=NULL;
  double* Coordinates_buffer=NULL;
  dim_t p, buffer_len,n;
  char error_msg[100];
  #ifdef PASO_MPI
  MPI_Status status;
  #endif

  /* get the global range of node ids */
  Finley_NodeFile_setGlobalIdRange(&min_id,&max_id,in);
  undefined_node=min_id-1;

  distribution=TMPMEMALLOC(in->MPIInfo->size+1, index_t);

  if ( !Finley_checkPtr(distribution) ) {
      /* distribute the range of node ids */
      buffer_len=Paso_MPIInfo_setDistribution(in->MPIInfo,min_id,max_id,distribution);
      /* allocate buffers */
      Id_buffer=TMPMEMALLOC(buffer_len,index_t);
      Tag_buffer=TMPMEMALLOC(buffer_len,index_t);
      globalDegreesOfFreedom_buffer=TMPMEMALLOC(buffer_len,index_t);
      Coordinates_buffer=TMPMEMALLOC(buffer_len*out->numDim,double);
      if (! (Finley_checkPtr(Id_buffer) || Finley_checkPtr(Tag_buffer) || 
                     Finley_checkPtr(globalDegreesOfFreedom_buffer) || Finley_checkPtr(Coordinates_buffer) ) ) {
            /* fill Id_buffer by the undefined_node marker to check if nodes are defined */
            #pragma omp parallel for private(n) schedule(static)
            for (n=0;n<buffer_len;n++) Id_buffer[n]=undefined_node;
            
            /* fill the buffer by sending portions around in a circle */
            dest=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
            source=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
            buffer_rank=in->MPIInfo->rank;
            for (p=0; p< in->MPIInfo->size; ++p) {
                 if (p>0) {  /* the initial send can be skipped */
                     #ifdef PASO_MPI
                     MPI_Sendrecv_replace(Id_buffer, buffer_len, MPI_INT,
                                          dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
                                          in->MPIInfo->comm,&status);
                     MPI_Sendrecv_replace(Tag_buffer, buffer_len, MPI_INT,
                                          dest, in->MPIInfo->msg_tag_counter+1, source, in->MPIInfo->msg_tag_counter+1,
                                          in->MPIInfo->comm,&status);
                     MPI_Sendrecv_replace(globalDegreesOfFreedom_buffer, buffer_len, MPI_INT,
                                          dest, in->MPIInfo->msg_tag_counter+2, source, in->MPIInfo->msg_tag_counter+2,
                                          in->MPIInfo->comm,&status);
                     MPI_Sendrecv_replace(Coordinates_buffer, buffer_len*out->numDim, MPI_DOUBLE,
                                          dest, in->MPIInfo->msg_tag_counter+3, source, in->MPIInfo->msg_tag_counter+3,
                                          in->MPIInfo->comm,&status);
                     #endif
                     in->MPIInfo->msg_tag_counter+=4;
                 }
                 buffer_rank=Paso_MPIInfo_mod(in->MPIInfo->size, buffer_rank-1);
                 Finley_NodeFile_scatterEntries(in->numNodes, in->Id, 
                                                distribution[buffer_rank], distribution[buffer_rank+1],
                                                Id_buffer, in->Id,
                                                Tag_buffer, in->Tag, 
                                                globalDegreesOfFreedom_buffer, in->globalDegreesOfFreedom,
                                                out->numDim, Coordinates_buffer, in->Coordinates);
            }
            /* now entries are collected from the buffer again by sending the entries around in a circle */
            dest=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
            source=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
            buffer_rank=in->MPIInfo->rank;
            for (p=0; p< in->MPIInfo->size; ++p) {
                 Finley_NodeFile_gatherEntries(out->numNodes, index, 
                                               distribution[buffer_rank], distribution[buffer_rank+1],
                                               out->Id, Id_buffer, 
                                               out->Tag, Tag_buffer, 
                                               out->globalDegreesOfFreedom, globalDegreesOfFreedom_buffer, 
                                               out->numDim, out->Coordinates, Coordinates_buffer);
                 if (p<in->MPIInfo->size-1) {  /* the last send can be skipped */
                     #ifdef PASO_MPI
                     MPI_Sendrecv_replace(Id_buffer, buffer_len, MPI_INT,
                                          dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
                                          in->MPIInfo->comm,&status);
                     MPI_Sendrecv_replace(Tag_buffer, buffer_len, MPI_INT,
                                          dest, in->MPIInfo->msg_tag_counter+1, source, in->MPIInfo->msg_tag_counter+1,
                                          in->MPIInfo->comm,&status);
                     MPI_Sendrecv_replace(globalDegreesOfFreedom_buffer, buffer_len, MPI_INT,
                                          dest, in->MPIInfo->msg_tag_counter+2, source, in->MPIInfo->msg_tag_counter+2,
                                          in->MPIInfo->comm,&status);
                     MPI_Sendrecv_replace(Coordinates_buffer, buffer_len*out->numDim, MPI_DOUBLE,
                                          dest, in->MPIInfo->msg_tag_counter+3, source, in->MPIInfo->msg_tag_counter+3,
                                          in->MPIInfo->comm,&status);
                     #endif
                     in->MPIInfo->msg_tag_counter+=4;
                 }
                 buffer_rank=Paso_MPIInfo_mod(in->MPIInfo->size, buffer_rank-1);
            }
            /* check if all nodes are set: */
            #pragma omp parallel for private(n) schedule(static)
            for (n=0; n< out->numNodes; ++n) {
                if (out->Id[n] == undefined_node ) {
                 sprintf(error_msg,"Finley_NodeFile_gather_global: Node id %d is referenced but is not defined.",out->Id[n]);
                 Finley_setError(VALUE_ERROR,error_msg);
               }
             }

      }
      TMPMEMFREE(Id_buffer);
      TMPMEMFREE(Tag_buffer);
      TMPMEMFREE(globalDegreesOfFreedom_buffer);
      TMPMEMFREE(Coordinates_buffer);
  }
  TMPMEMFREE(distribution);
  /* make sure that the error is global */
  Paso_MPIInfo_noError(in->MPIInfo);
}
