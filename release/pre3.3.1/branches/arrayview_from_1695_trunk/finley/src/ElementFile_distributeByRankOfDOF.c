
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

/*   Finley: ElementFile: this will redistribute the Elements including overlap by */

/**************************************************************/

#include "ElementFile.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/

void Finley_ElementFile_distributeByRankOfDOF(Finley_ElementFile* self, Paso_MPI_rank* mpiRankOfDOF, index_t* Id) {
     size_t size_size;
     Paso_MPI_rank myRank, p, *Owner_buffer=NULL, loc_proc_mask_max;
     dim_t e, j, i, size, *send_count=NULL, *recv_count=NULL, *newOwner=NULL, *loc_proc_mask=NULL, *loc_send_count=NULL,
           newNumElements, numElementsInBuffer, numNodes, numRequests, NN;
     index_t *send_offset=NULL, *recv_offset=NULL, *Id_buffer=NULL, *Tag_buffer=NULL, *Nodes_buffer=NULL, k;
     bool_t *proc_mask=NULL;
     #ifdef PASO_MPI
        MPI_Request* mpi_requests=NULL;
        MPI_Status* mpi_stati=NULL;
     #endif
     if (self==NULL) return;
     myRank=self->MPIInfo->rank;
     size=self->MPIInfo->size;
     size_size=size*sizeof(dim_t);
     numNodes=self->numNodes;
     NN=self->numNodes;
     if (size>1) {
         #ifdef PASO_MPI
            mpi_requests=TMPMEMALLOC(8*size, MPI_Request);
            mpi_stati=TMPMEMALLOC(8*size, MPI_Status);
            Finley_checkPtr(mpi_requests);
            Finley_checkPtr(mpi_stati);
         #endif

        /* count the number elements that have to be send to each processor (send_count) 
           and define a new element owner as the processor with the largest number of DOFs and the smallest id */
        send_count=TMPMEMALLOC(size,dim_t);
        recv_count=TMPMEMALLOC(size,dim_t);
        newOwner=TMPMEMALLOC(self->numElements,Paso_MPI_rank);
        if ( !( Finley_checkPtr(send_count) || Finley_checkPtr(recv_count) || Finley_checkPtr(newOwner) ) ) {
           memset(send_count, 0, size_size);
           #pragma omp parallel private(p,loc_proc_mask,loc_send_count)
           {
               loc_proc_mask=THREAD_MEMALLOC(size,dim_t);
               loc_send_count=THREAD_MEMALLOC(size,dim_t);
               memset(loc_send_count, 0, size_size); 
               #pragma omp for private(e,j,loc_proc_mask_max) schedule(static)
               for (e=0;e<self->numElements;e++) {
                  if (self->Owner[e] == myRank) {
                    newOwner[e]=myRank;
                    memset(loc_proc_mask, 0, size_size);
                    for(j=0;j<numNodes;j++) {
                        p=mpiRankOfDOF[self->Nodes[INDEX2(j,e,NN)]];
                        loc_proc_mask[p]++;
                    }
                    loc_proc_mask_max=0;
                    for (p=0;p<size;++p) {
                       if (loc_proc_mask[p]>0) loc_send_count[p]++;
                       if (loc_proc_mask[p]>loc_proc_mask_max) {
                          newOwner[e]=p;
                          loc_proc_mask_max=loc_proc_mask[p];
                       }
                    }
                  } else {
                     newOwner[e]=-1;
                  }
               }
               #pragma omp critical
               {
                 for (p=0;p<size;++p) send_count[p]+=loc_send_count[p];
               }
               THREAD_MEMFREE(loc_proc_mask);
               THREAD_MEMFREE(loc_send_count);
           }
           #ifdef PASO_MPI
              MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,self->MPIInfo->comm);
           #else
              for (p=0;p<size;++p) recv_count[p]=send_count[p];
           #endif
           /* get the new number of elements for this processor */
           newNumElements=0;
           for (p=0;p<size;++p) newNumElements+=recv_count[p];

           /* get the new number of elements for this processor */
           numElementsInBuffer=0;
           for (p=0;p<size;++p) numElementsInBuffer+=send_count[p];
           /* allocate buffers */
           Id_buffer=TMPMEMALLOC(numElementsInBuffer,index_t);
           Tag_buffer=TMPMEMALLOC(numElementsInBuffer,index_t);
           Owner_buffer=TMPMEMALLOC(numElementsInBuffer,Paso_MPI_rank);
           Nodes_buffer=TMPMEMALLOC(numElementsInBuffer*NN,index_t);
           send_offset=TMPMEMALLOC(size,index_t);
           recv_offset=TMPMEMALLOC(size,index_t);
           proc_mask=TMPMEMALLOC(size,bool_t);
           if ( !( Finley_checkPtr(Id_buffer) || Finley_checkPtr(Tag_buffer) || Finley_checkPtr(Owner_buffer) ||
                   Finley_checkPtr(Nodes_buffer) || Finley_checkPtr(send_offset) || Finley_checkPtr(recv_offset) || 
                   Finley_checkPtr(proc_mask) )) {

              /* callculate the offsets for the processor buffers */
              recv_offset[0]=0;
              for (p=0;p<size-1;++p) recv_offset[p+1]=recv_offset[p]+recv_count[p];
              send_offset[0]=0;
              for (p=0;p<size-1;++p) send_offset[p+1]=send_offset[p]+send_count[p];

              memset(send_count, 0, size_size);
              /* copy element into buffers. proc_mask makes sure that an element is copied once only for each processor */
              for (e=0;e<self->numElements;e++) {
                 if (self->Owner[e] == myRank) {
                    memset(proc_mask, TRUE, size_size);
                    for(j=0;j<numNodes;j++) {
                         p=mpiRankOfDOF[self->Nodes[INDEX2(j,e,NN)]];
                         if (proc_mask[p]) {
                            k=send_offset[p]+send_count[p];
                            Id_buffer[k]=self->Id[e];
                            Tag_buffer[k]=self->Tag[e];
                            Owner_buffer[k]=newOwner[e];
                            for (i=0;i<numNodes;i++) Nodes_buffer[INDEX2(i,k,NN)]=Id[self->Nodes[INDEX2(i,e,NN)]];
                            send_count[p]++;
                            proc_mask[p]=FALSE;
                         }
                    }
                 }
              }
              /* allocate new tables */
              Finley_ElementFile_allocTable(self,newNumElements);

              /* start to receive new elements */
              numRequests=0;
              for (p=0;p<size;++p) {
                 if (recv_count[p]>0) {
                    #ifdef PASO_MPI
                    MPI_Irecv(&(self->Id[recv_offset[p]]), recv_count[p], 
                              MPI_INT, p, self->MPIInfo->msg_tag_counter+myRank,
                              self->MPIInfo->comm, &mpi_requests[numRequests]);
                    numRequests++;
                    MPI_Irecv(&(self->Tag[recv_offset[p]]), recv_count[p], 
                              MPI_INT, p, self->MPIInfo->msg_tag_counter+size+myRank,
                              self->MPIInfo->comm, &mpi_requests[numRequests]);
                    numRequests++;
                    MPI_Irecv(&(self->Owner[recv_offset[p]]), recv_count[p], 
                              MPI_INT, p, self->MPIInfo->msg_tag_counter+2*size+myRank,
                              self->MPIInfo->comm, &mpi_requests[numRequests]);
                    numRequests++;
                    MPI_Irecv(&(self->Nodes[recv_offset[p]*NN]), recv_count[p]*NN, 
                              MPI_INT, p, self->MPIInfo->msg_tag_counter+3*size+myRank,
                              self->MPIInfo->comm, &mpi_requests[numRequests]);
                    numRequests++;
                    #endif
                 }
              }
              /* now the buffers can be send away */
              for (p=0;p<size;++p) {
                 if (send_count[p]>0) {
                   #ifdef PASO_MPI
                   MPI_Issend(&(Id_buffer[send_offset[p]]), send_count[p], 
                              MPI_INT, p, self->MPIInfo->msg_tag_counter+p,
                              self->MPIInfo->comm, &mpi_requests[numRequests]);
                   numRequests++;
                   MPI_Issend(&(Tag_buffer[send_offset[p]]), send_count[p], 
                              MPI_INT, p, self->MPIInfo->msg_tag_counter+size+p,
                              self->MPIInfo->comm, &mpi_requests[numRequests]);
                   numRequests++;
                   MPI_Issend(&(Owner_buffer[send_offset[p]]), send_count[p], 
                              MPI_INT, p, self->MPIInfo->msg_tag_counter+2*size+p,
                              self->MPIInfo->comm, &mpi_requests[numRequests]);
                   numRequests++;
                   MPI_Issend(&(Nodes_buffer[send_offset[p]*NN]), send_count[p]*NN, 
                              MPI_INT, p, self->MPIInfo->msg_tag_counter+3*size+p,
                              self->MPIInfo->comm, &mpi_requests[numRequests]);
                   numRequests++;
                   #endif

                 }
              }
              self->MPIInfo->msg_tag_counter+=4*size;
              /* wait for the requests to be finalized */
              #ifdef PASO_MPI
              MPI_Waitall(numRequests,mpi_requests,mpi_stati);
              #endif
           }
           /* clear buffer */
           TMPMEMFREE(Id_buffer);
           TMPMEMFREE(Tag_buffer);
           TMPMEMFREE(Owner_buffer);
           TMPMEMFREE(Nodes_buffer);
           TMPMEMFREE(send_offset);
           TMPMEMFREE(recv_offset);
           TMPMEMFREE(proc_mask);
        }
        #ifdef PASO_MPI
            TMPMEMFREE(mpi_requests);
            TMPMEMFREE(mpi_stati);
        #endif
        TMPMEMFREE(send_count);
        TMPMEMFREE(recv_count);
        TMPMEMFREE(newOwner);
     } else {
        #pragma omp for private(e,i) schedule(static)
        for (e=0;e<self->numElements;e++) {
            self->Owner[e]=myRank;
            for (i=0;i<numNodes;i++) self->Nodes[INDEX2(i,e,NN)]=Id[self->Nodes[INDEX2(i,e,NN)]];
        }
     }
     return;
}

