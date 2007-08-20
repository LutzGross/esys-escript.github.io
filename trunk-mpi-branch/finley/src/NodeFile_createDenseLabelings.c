/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

/**************************************************************/

/*   Finley: Mesh: NodeFile                                   */

/*   creates a dense labeling of the global degrees of freedom  */
/*   and returns the new number of  global degrees of freedom  */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id:$ */

/**************************************************************/

#include "NodeFile.h"

/**************************************************************/

dim_t Finley_NodeFile_createDenseDOFLabeling(Finley_NodeFile* in) 
{
  index_t min_dof, max_dof, unset_dof=-1,set_dof=1, dof_0, dof_1, *DOF_buffer=NULL, k;
  Paso_MPI_rank buffer_rank, dest, source, *distribution=NULL;
  dim_t p, buffer_len,n, myDOFs, *offsets=NULL, *loc_offsets=NULL, new_numGlobalDOFs=0, myNewDOFs;
  bool_t *set_new_DOF=NULL;
  #ifdef PASO_MPI
  MPI_Status status;
  #endif

  /* get the global range of node ids */
  Finley_NodeFile_setGlobalDOFRange(&min_dof,&max_dof,in);

  distribution=TMPMEMALLOC(in->MPIInfo->size+1, index_t);
  offsets=TMPMEMALLOC(in->MPIInfo->size, dim_t);
  loc_offsets=TMPMEMALLOC(in->MPIInfo->size, dim_t);
  set_new_DOF=TMPMEMALLOC(in->numNodes, bool_t);

  if ( ! (Finley_checkPtr(distribution) || Finley_checkPtr(offsets) || Finley_checkPtr(loc_offsets) || Finley_checkPtr(set_new_DOF)) ) {
      /* distribute the range of node ids */
      buffer_len=Paso_MPIInfo_setDistribution(in->MPIInfo,min_dof,max_dof,distribution);
      myDOFs=distribution[in->MPIInfo->rank+1]-distribution[in->MPIInfo->rank];
      /* allocate buffers */
      DOF_buffer=TMPMEMALLOC(buffer_len,index_t);
      if (! Finley_checkPtr(DOF_buffer)) {
            /* fill DOF_buffer by the unset_dof marker to check if nodes are defined */
            #pragma omp parallel for private(n) schedule(static)
            for (n=0;n<buffer_len;n++) DOF_buffer[n]=unset_dof;
            
            /* fill the buffer by sending portions around in a circle */
            dest=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
            source=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
            buffer_rank=in->MPIInfo->rank;
            for (p=0; p< in->MPIInfo->size; ++p) {
                 if (p>0) {  /* the initial send can be skipped */
                     #ifdef PASO_MPI
                     MPI_Sendrecv_replace(DOF_buffer, buffer_len, MPI_INT,
                                          dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
                                          in->MPIInfo->comm,&status);
                     #endif
                     in->MPIInfo->msg_tag_counter++;
                 }
                 buffer_rank=Paso_MPIInfo_mod(in->MPIInfo->size, buffer_rank-1);
                 dof_0=distribution[buffer_rank];
                 dof_1=distribution[buffer_rank+1];
                 #pragma omp parallel for private(n,k) schedule(static)
                 for (n=0;n<in->numNodes;n++) {
                     k=in->globalDegreesOfFreedom[n];
                     if ((dof_0<=k) && (k<dof_1)) {
                         DOF_buffer[k-dof_0] = set_dof;
                     }
                 }
            }
            /* count the entries in the DOF_buffer */ 
            /* TODO: OMP parallel */
            myNewDOFs=0;
            for (n=0; n<myDOFs; ++n) {
                if ( DOF_buffer[n] == set_dof) {
                      DOF_buffer[n]=myNewDOFs;
                      myNewDOFs++;
                }
            }
            memset(loc_offsets,0,in->MPIInfo->size*sizeof(dim_t));
            loc_offsets[in->MPIInfo->rank]=myNewDOFs;
            #ifdef PASO_MPI
               MPI_Allreduce(loc_offsets,offsets,in->MPIInfo->size, MPI_INT, MPI_SUM, in->MPIInfo->comm );
               new_numGlobalDOFs=0;
               for (n=0; n< in->MPIInfo->size; ++n) {
                      loc_offsets[n]=new_numGlobalDOFs;
                      new_numGlobalDOFs+=offsets[n];
               }
            #else
               new_numGlobalDOFs=loc_offsets[0];
               loc_offsets[0]=0;
            #endif
            #pragma omp parallel for private(n) schedule(static)
            for (n=0; n<myDOFs; ++n) DOF_buffer[n]+=loc_offsets[in->MPIInfo->rank];
            /* now entries are collected from the buffer again by sending the entries around in a circle */
            #pragma omp parallel for private(n) schedule(static)
            for (n=0; n<in->numNodes; ++n) set_new_DOF[n]=TRUE;
            dest=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
            source=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
            buffer_rank=in->MPIInfo->rank;
            for (p=0; p< in->MPIInfo->size; ++p) {
                 dof_0=distribution[buffer_rank];
                 dof_1=distribution[buffer_rank+1];
                 #pragma omp parallel for private(n,k) schedule(static)
                 for (n=0;n<in->numNodes;n++) {
                      k=in->globalDegreesOfFreedom[n];
                      if ( set_new_DOF[n] && (dof_0<=k) && (k<dof_1)) {
                           in->globalDegreesOfFreedom[n]=DOF_buffer[k-dof_0];
                           set_new_DOF[n]=FALSE;
                      }
                 }
                 if (p<in->MPIInfo->size-1) {  /* the last send can be skipped */
                     #ifdef PASO_MPI
                     MPI_Sendrecv_replace(DOF_buffer, buffer_len, MPI_INT,
                                          dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
                                          in->MPIInfo->comm,&status);
                     #endif
                     in->MPIInfo->msg_tag_counter+=1;
                 }
                 buffer_rank=Paso_MPIInfo_mod(in->MPIInfo->size, buffer_rank-1);
            }
      }
      TMPMEMFREE(DOF_buffer);
  }
  TMPMEMFREE(distribution);
  TMPMEMFREE(loc_offsets);
  TMPMEMFREE(offsets);
  TMPMEMFREE(set_new_DOF);
  in->isPrepared=FINLEY_UNPREPARED;
  return new_numGlobalDOFs;
}

void Finley_NodeFile_assignMPIRankToDOFs(Finley_NodeFile* in,Paso_MPI_rank* mpiRankOfDOF, index_t *distribution){
  index_t min_DOF,max_DOF, k;
  dim_t n;
  Paso_MPI_rank p, p_min=in->MPIInfo->size, p_max=-1;
  /* first we calculate the min and max dof on this processor to reduce costs for seraching */
  Finley_NodeFile_setDOFRange(&min_DOF,&max_DOF,in);

  for (p=0; p<in->MPIInfo->size; ++p) {
       if (distribution[p]<=min_DOF) p_min=p;
       if (distribution[p]<=max_DOF) p_max=p;
  }
  #pragma omp parallel for private(n,k,p) schedule(static)
  for (n=0; n<in->numNodes; ++n) {
     k=in->globalDegreesOfFreedom[n];
     for (p=p_min; p<=p_max; ++p) {
        if (k<distribution[p+1]) {
              mpiRankOfDOF[n]=p;
              break;
        }
     }
  }
} 
dim_t Finley_NodeFile_createDenseReducedDOFLabeling(Finley_NodeFile* in,index_t* reducedNodeMask) 
{
  index_t min_dof, max_dof, unset_dof=-1,set_dof=1, dof_0, dof_1, *DOF_buffer=NULL, k;
  Paso_MPI_rank buffer_rank, dest, source, *distribution=NULL;
  dim_t p, buffer_len,n, myDOFs, *offsets=NULL, *loc_offsets=NULL, globalNumReducedDOFs=0, myNewDOFs;
  #ifdef PASO_MPI
  MPI_Status status;
  #endif

  /* get the global range of node ids */
  Finley_NodeFile_setGlobalDOFRange(&min_dof,&max_dof,in);

  distribution=TMPMEMALLOC(in->MPIInfo->size+1, index_t);
  offsets=TMPMEMALLOC(in->MPIInfo->size, dim_t);
  loc_offsets=TMPMEMALLOC(in->MPIInfo->size, dim_t);

  if ( ! (Finley_checkPtr(distribution) || Finley_checkPtr(offsets) || Finley_checkPtr(loc_offsets) ) ) {
      /* distribute the range of node ids */
      buffer_len=Paso_MPIInfo_setDistribution(in->MPIInfo,min_dof,max_dof,distribution);
      myDOFs=distribution[in->MPIInfo->rank+1]-distribution[in->MPIInfo->rank];
      /* allocate buffers */
      DOF_buffer=TMPMEMALLOC(buffer_len,index_t);
      if (! Finley_checkPtr(DOF_buffer)) {
            /* fill DOF_buffer by the unset_dof marker to check if nodes are defined */
            #pragma omp parallel for private(n) schedule(static)
            for (n=0;n<buffer_len;n++) DOF_buffer[n]=unset_dof;
            
            /* fill the buffer by sending portions around in a circle */
            dest=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
            source=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
            buffer_rank=in->MPIInfo->rank;
            for (p=0; p< in->MPIInfo->size; ++p) {
                 if (p>0) {  /* the initial send can be skipped */
                     #ifdef PASO_MPI
                     MPI_Sendrecv_replace(DOF_buffer, buffer_len, MPI_INT,
                                          dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
                                          in->MPIInfo->comm,&status);
                     #endif
                     in->MPIInfo->msg_tag_counter++;
                 }
                 buffer_rank=Paso_MPIInfo_mod(in->MPIInfo->size, buffer_rank-1);
                 dof_0=distribution[buffer_rank];
                 dof_1=distribution[buffer_rank+1];
                 #pragma omp parallel for private(n,k) schedule(static)
                 for (n=0;n<in->numNodes;n++) {
                     if (reducedNodeMask[n] >-1) {
                        k=in->globalDegreesOfFreedom[n];
                        if ((dof_0<=k) && (k<dof_1)) {
                            DOF_buffer[k-dof_0] = set_dof;
                        }
                     }
                 }
            }
            /* count the entries in the DOF_buffer */ 
            /* TODO: OMP parallel */
            myNewDOFs=0;
            for (n=0; n<myDOFs; ++n) {
                if ( DOF_buffer[n] == set_dof) {
                      DOF_buffer[n]=myNewDOFs;
                      myNewDOFs++;
                }
            }
            memset(loc_offsets,0,in->MPIInfo->size*sizeof(dim_t));
            loc_offsets[in->MPIInfo->rank]=myNewDOFs;
            #ifdef PASO_MPI
               MPI_Allreduce(loc_offsets,offsets,in->MPIInfo->size, MPI_INT, MPI_SUM, in->MPIInfo->comm );
               globalNumReducedDOFs=0;
               for (n=0; n< in->MPIInfo->size; ++n) {
                      loc_offsets[n]=globalNumReducedDOFs;
                      globalNumReducedDOFs+=offsets[n];
               }
            #else
               globalNumReducedDOFs=loc_offsets[0];
               loc_offsets[n]=0;
            #endif
            #pragma omp parallel for private(n) schedule(static)
            for (n=0; n<myDOFs; ++n) DOF_buffer[n]+=loc_offsets[in->MPIInfo->rank];
            /* now entries are collected from the buffer again by sending the entries around in a circle */
            #pragma omp parallel for private(n) schedule(static)
            for (n=0; n<in->numNodes; ++n) in->globalReducedDOFIndex[n]=-1;
            dest=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
            source=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
            buffer_rank=in->MPIInfo->rank;
            for (p=0; p< in->MPIInfo->size; ++p) {
                 dof_0=distribution[buffer_rank];
                 dof_1=distribution[buffer_rank+1];
                 #pragma omp parallel for private(n,k) schedule(static)
                 for (n=0;n<in->numNodes;n++) {
                      if (reducedNodeMask[n] >-1) {
                         k=in->globalDegreesOfFreedom[n];
                         if ( (dof_0<=k) && (k<dof_1)) in->globalReducedDOFIndex[n]=DOF_buffer[k-dof_0];
                      }
                 }
                 if (p<in->MPIInfo->size-1) {  /* the last send can be skipped */
                     #ifdef PASO_MPI
                     MPI_Sendrecv_replace(DOF_buffer, buffer_len, MPI_INT,
                                          dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
                                          in->MPIInfo->comm,&status);
                     #endif
                     in->MPIInfo->msg_tag_counter+=1;
                 }
                 buffer_rank=Paso_MPIInfo_mod(in->MPIInfo->size, buffer_rank-1);
            }
      }
      TMPMEMFREE(DOF_buffer);
  }
  TMPMEMFREE(distribution);
  TMPMEMFREE(loc_offsets);
  TMPMEMFREE(offsets);
  in->isPrepared=FINLEY_UNPREPARED;
  return globalNumReducedDOFs;
}
dim_t Finley_NodeFile_createDenseNodeLabeling(Finley_NodeFile* in) 
{
  index_t min_nodeID, max_nodeID, unset_nodeID=-1,set_nodeID=1, nodeID_0, nodeID_1, *Node_buffer=NULL, k;
  Paso_MPI_rank buffer_rank, dest, source, *distribution=NULL;
  dim_t p, buffer_len,n, myNodes, *offsets=NULL, *loc_offsets=NULL, globalNumNodes=0, myNewNodes;
  #ifdef PASO_MPI
  MPI_Status status;
  #endif

  /* get the global range of node ids */
  Finley_NodeFile_setGlobalIdRange(&min_nodeID,&max_nodeID,in);

  distribution=TMPMEMALLOC(in->MPIInfo->size+1, index_t);
  offsets=TMPMEMALLOC(in->MPIInfo->size, dim_t);
  loc_offsets=TMPMEMALLOC(in->MPIInfo->size, dim_t);

  if ( ! (Finley_checkPtr(distribution) || Finley_checkPtr(offsets) || Finley_checkPtr(loc_offsets)) ) {
      /* distribute the range of node ids */
      buffer_len=Paso_MPIInfo_setDistribution(in->MPIInfo,min_nodeID,max_nodeID,distribution);
      myNodes=distribution[in->MPIInfo->rank+1]-distribution[in->MPIInfo->rank];
      /* allocate buffers */
      Node_buffer=TMPMEMALLOC(buffer_len,index_t);
      if (! Finley_checkPtr(Node_buffer)) {
            /* fill Node_buffer by the unset_nodeID marker to check if nodes are defined */
            #pragma omp parallel for private(n) schedule(static)
            for (n=0;n<buffer_len;n++) Node_buffer[n]=unset_nodeID;
            
            /* fill the buffer by sending portions around in a circle */
            dest=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
            source=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
            buffer_rank=in->MPIInfo->rank;
            for (p=0; p< in->MPIInfo->size; ++p) {
                 if (p>0) {  /* the initial send can be skipped */
                     #ifdef PASO_MPI
                     MPI_Sendrecv_replace(Node_buffer, buffer_len, MPI_INT,
                                          dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
                                          in->MPIInfo->comm,&status);
                     #endif
                     in->MPIInfo->msg_tag_counter++;
                 }
                 buffer_rank=Paso_MPIInfo_mod(in->MPIInfo->size, buffer_rank-1);
                 nodeID_0=distribution[buffer_rank];
                 nodeID_1=distribution[buffer_rank+1];
                 #pragma omp parallel for private(n,k) schedule(static)
                 for (n=0;n<in->numNodes;n++) {
                     k=in->Id[n];
                     if ((nodeID_0<=k) && (k<nodeID_1)) {
                         Node_buffer[k-nodeID_0] = set_nodeID;
                     }
                 }
            }
            /* count the entries in the Node_buffer */ 
            /* TODO: OMP parallel */
            myNewNodes=0;
            for (n=0; n<myNodes; ++n) {
                if ( Node_buffer[n] == set_nodeID) {
                      Node_buffer[n]=myNewNodes;
                      myNewNodes++;
                }
            }
            memset(loc_offsets,0,in->MPIInfo->size*sizeof(dim_t));
            loc_offsets[in->MPIInfo->rank]=myNewNodes;
            #ifdef PASO_MPI
               MPI_Allreduce(loc_offsets,offsets,in->MPIInfo->size, MPI_INT, MPI_SUM, in->MPIInfo->comm );
               globalNumNodes=0;
               for (n=0; n< in->MPIInfo->size; ++n) {
                      loc_offsets[n]=globalNumNodes;
                      globalNumNodes+=offsets[n];
               }
            #else
               globalNumNodes=loc_offsets[0];
               loc_offsets[n]=0;
            #endif
            #pragma omp parallel for private(n) schedule(static)
            for (n=0; n<myNodes; ++n) Node_buffer[n]+=loc_offsets[in->MPIInfo->rank];
            /* now entries are collected from the buffer again by sending the entries around in a circle */
            #pragma omp parallel for private(n) schedule(static)
            for (n=0; n<in->numNodes; ++n) in->globalNodesIndex[n]=-1;
            dest=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
            source=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
            buffer_rank=in->MPIInfo->rank;
            for (p=0; p< in->MPIInfo->size; ++p) {
                 nodeID_0=distribution[buffer_rank];
                 nodeID_1=distribution[buffer_rank+1];
                 #pragma omp parallel for private(n,k) schedule(static)
                 for (n=0;n<in->numNodes;n++) {
                      k=in->Id[n];
                      if ( ( nodeID_0<=k) && (k<nodeID_1) ) in->globalNodesIndex[n]=Node_buffer[k-nodeID_0];
                 }
                 if (p<in->MPIInfo->size-1) {  /* the last send can be skipped */
                     #ifdef PASO_MPI
                     MPI_Sendrecv_replace(Node_buffer, buffer_len, MPI_INT,
                                          dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
                                          in->MPIInfo->comm,&status);
                     #endif
                     in->MPIInfo->msg_tag_counter+=1;
                 }
                 buffer_rank=Paso_MPIInfo_mod(in->MPIInfo->size, buffer_rank-1);
            }
      }
      TMPMEMFREE(Node_buffer);
  }
  TMPMEMFREE(distribution);
  TMPMEMFREE(loc_offsets);
  TMPMEMFREE(offsets);
  in->isPrepared=FINLEY_UNPREPARED;
  return globalNumNodes;
}

dim_t Finley_NodeFile_createDenseReducedNodeLabeling(Finley_NodeFile* in, index_t* reducedNodeMask)
{
  index_t min_nodeID, max_nodeID, unset_nodeID=-1,set_nodeID=1, nodeID_0, nodeID_1, *Node_buffer=NULL, k;
  Paso_MPI_rank buffer_rank, dest, source, *distribution=NULL;
  dim_t p, buffer_len,n, myNodes, *offsets=NULL, *loc_offsets=NULL, globalNumReducedNodes=0, myNewNodes;
  #ifdef PASO_MPI
  MPI_Status status;
  #endif

  /* get the global range of node ids */
  Finley_NodeFile_setGlobalIdRange(&min_nodeID,&max_nodeID,in);

  distribution=TMPMEMALLOC(in->MPIInfo->size+1, index_t);
  offsets=TMPMEMALLOC(in->MPIInfo->size, dim_t);
  loc_offsets=TMPMEMALLOC(in->MPIInfo->size, dim_t);

  if ( ! (Finley_checkPtr(distribution) || Finley_checkPtr(offsets) || Finley_checkPtr(loc_offsets) )) {
      /* distribute the range of node ids */
      buffer_len=Paso_MPIInfo_setDistribution(in->MPIInfo,min_nodeID,max_nodeID,distribution);
      myNodes=distribution[in->MPIInfo->rank+1]-distribution[in->MPIInfo->rank];
      /* allocate buffers */
      Node_buffer=TMPMEMALLOC(buffer_len,index_t);
      if (! Finley_checkPtr(Node_buffer)) {
            /* fill Node_buffer by the unset_nodeID marker to check if nodes are defined */
            #pragma omp parallel for private(n) schedule(static)
            for (n=0;n<buffer_len;n++) Node_buffer[n]=unset_nodeID;
            
            /* fill the buffer by sending portions around in a circle */
            dest=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
            source=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
            buffer_rank=in->MPIInfo->rank;
            for (p=0; p< in->MPIInfo->size; ++p) {
                 if (p>0) {  /* the initial send can be skipped */
                     #ifdef PASO_MPI
                     MPI_Sendrecv_replace(Node_buffer, buffer_len, MPI_INT,
                                          dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
                                          in->MPIInfo->comm,&status);
                     #endif
                     in->MPIInfo->msg_tag_counter++;
                 }
                 buffer_rank=Paso_MPIInfo_mod(in->MPIInfo->size, buffer_rank-1);
                 nodeID_0=distribution[buffer_rank];
                 nodeID_1=distribution[buffer_rank+1];
                 #pragma omp parallel for private(n,k) schedule(static)
                 for (n=0;n<in->numNodes;n++) {
                     if (reducedNodeMask[n] >-1) {
                        k=in->Id[n];
                        if ((nodeID_0<=k) && (k<nodeID_1)) {
                            Node_buffer[k-nodeID_0] = set_nodeID;
                        }
                     }
                 }
            }
            /* count the entries in the Node_buffer */ 
            /* TODO: OMP parallel */
            myNewNodes=0;
            for (n=0; n<myNodes; ++n) {
                if ( Node_buffer[n] == set_nodeID) {
                      Node_buffer[n]=myNewNodes;
                      myNewNodes++;
                }
            }
            memset(loc_offsets,0,in->MPIInfo->size*sizeof(dim_t));
            loc_offsets[in->MPIInfo->rank]=myNewNodes;
            #ifdef PASO_MPI
               MPI_Allreduce(loc_offsets,offsets,in->MPIInfo->size, MPI_INT, MPI_SUM, in->MPIInfo->comm );
               globalNumReducedNodes=0;
               for (n=0; n< in->MPIInfo->size; ++n) {
                      loc_offsets[n]=globalNumReducedNodes;
                      globalNumReducedNodes+=offsets[n];
               }
            #else
               globalNumReducedNodes=loc_offsets[0];
               loc_offsets[n]=0;
            #endif
            #pragma omp parallel for private(n) schedule(static)
            for (n=0; n<myNodes; ++n) Node_buffer[n]+=loc_offsets[in->MPIInfo->rank];
            /* now entries are collected from the buffer again by sending the entries around in a circle */
            #pragma omp parallel for private(n) schedule(static)
            for (n=0; n<in->numNodes; ++n) in->globalReducedNodesIndex[n]=TRUE;
            dest=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank + 1);
            source=Paso_MPIInfo_mod(in->MPIInfo->size, in->MPIInfo->rank - 1);
            buffer_rank=in->MPIInfo->rank;
            for (p=0; p< in->MPIInfo->size; ++p) {
                 nodeID_0=distribution[buffer_rank];
                 nodeID_1=distribution[buffer_rank+1];
                 #pragma omp parallel for private(n,k) schedule(static)
                 for (n=0;n<in->numNodes;n++) {
                     if (reducedNodeMask[n] >-1) {
                        k=in->Id[n];
                        if ( (nodeID_0<=k) && (k<nodeID_1)) in->globalReducedNodesIndex[n]=Node_buffer[k-nodeID_0];
                     }
                 }
                 if (p<in->MPIInfo->size-1) {  /* the last send can be skipped */
                     #ifdef PASO_MPI
                     MPI_Sendrecv_replace(Node_buffer, buffer_len, MPI_INT,
                                          dest, in->MPIInfo->msg_tag_counter, source, in->MPIInfo->msg_tag_counter,
                                          in->MPIInfo->comm,&status);
                     #endif
                     in->MPIInfo->msg_tag_counter+=1;
                 }
                 buffer_rank=Paso_MPIInfo_mod(in->MPIInfo->size, buffer_rank-1);
            }
      }
      TMPMEMFREE(Node_buffer);
  }
  TMPMEMFREE(distribution);
  TMPMEMFREE(loc_offsets);
  TMPMEMFREE(offsets);
  in->isPrepared=FINLEY_UNPREPARED;
  return globalNumReducedNodes;
}
