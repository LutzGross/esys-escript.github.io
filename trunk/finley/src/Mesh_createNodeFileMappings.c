
/* $Id:$ */

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

/*   Finley: NodeFile : creates the mappings using the indexReducedNodes */
/*                 no distribution is happening                          */

/**************************************************************/

#include "Mesh.h"
#define UNUSED -1

#define BOUNDS_CHECK 1

/**************************************************************/

void Mesh_createDOFMappingAndCoupling(Finley_Mesh* in, bool_t use_reduced_elements) 
{
  index_t min_DOF, max_DOF, *shared=NULL, *offsetInShared=NULL, *locDOFMask=NULL, i, k, myFirstDOF, myLastDOF, *nodeMask=NULL, firstDOF, lastDOF, *globalDOFIndex, *wanted_DOFs=NULL;
  dim_t mpiSize, len_loc_dof, numNeighbors, n, lastn, numNodes,*rcv_len=NULL, *snd_len=NULL, count;
  Paso_MPI_rank myRank,p,p_min,p_max, *neighbor=NULL;
  Paso_SharedComponents *rcv_shcomp=NULL, *snd_shcomp=NULL;
  Finley_NodeMapping *this_mapping=NULL;
  Paso_Connector* this_connector=NULL;
  Paso_Distribution* dof_distribution;
  Paso_MPIInfo *mpi_info = in->MPIInfo;
  #ifdef PASO_MPI
      MPI_Request* mpi_requests=NULL;
      MPI_Status* mpi_stati=NULL;
  #else
      int *mpi_requests=NULL, *mpi_stati=NULL;
  #endif

  numNodes=in->Nodes->numNodes;
  if (use_reduced_elements) {
    dof_distribution=in->Nodes->reducedDegreesOfFreedomDistribution;
    globalDOFIndex=in->Nodes->globalReducedDOFIndex;
  } else {
    dof_distribution=in->Nodes->degreesOfFreedomDistribution;
    globalDOFIndex=in->Nodes->globalDegreesOfFreedom;
  }
  myFirstDOF=Paso_Distribution_getFirstComponent(dof_distribution);
  myLastDOF=Paso_Distribution_getLastComponent(dof_distribution);


  mpiSize=mpi_info->size;
  myRank=mpi_info->rank;

  min_DOF=Finley_Util_getFlaggedMinInt(1,numNodes,globalDOFIndex,-1);
  max_DOF=Finley_Util_getFlaggedMaxInt(1,numNodes,globalDOFIndex,-1);

  if (max_DOF < min_DOF) {
      min_DOF=myFirstDOF;
      max_DOF=myLastDOF-1;
  }

  p_min=mpiSize;
  p_max=-1;
  if (max_DOF >= min_DOF) {
      for (p=0; p<mpiSize; ++p) {
         if (dof_distribution->first_component[p]<=min_DOF) p_min=p;
         if (dof_distribution->first_component[p]<=max_DOF) p_max=p;
     }
   }

  len_loc_dof=max_DOF-min_DOF+1;
  if (! ((min_DOF<=myFirstDOF) && (myLastDOF-1<=max_DOF)) ) {
      Finley_setError(SYSTEM_ERROR,"Local elements do not span local degrees of freedom.");
      return;
  }
  rcv_len=TMPMEMALLOC(mpiSize,dim_t);
  snd_len=TMPMEMALLOC(mpiSize,dim_t);
  #ifdef PASO_MPI
    mpi_requests=MEMALLOC(mpiSize*2,MPI_Request);
    mpi_stati=MEMALLOC(mpiSize*2,MPI_Status);
  #else
    mpi_requests=MEMALLOC(mpiSize*2,int);
    mpi_stati=MEMALLOC(mpiSize*2,int);
  #endif
  wanted_DOFs=TMPMEMALLOC(numNodes,index_t);
  nodeMask=TMPMEMALLOC(numNodes,index_t);
  neighbor=TMPMEMALLOC(mpiSize,Paso_MPI_rank);
  shared=TMPMEMALLOC(numNodes*(p_max-p_min+1),index_t);
  offsetInShared=TMPMEMALLOC(mpiSize+1,index_t);
  locDOFMask=TMPMEMALLOC(len_loc_dof, index_t);
  if (! ( Finley_checkPtr(neighbor) || Finley_checkPtr(shared) || Finley_checkPtr(offsetInShared) || Finley_checkPtr(locDOFMask) || 
     Finley_checkPtr(nodeMask) || Finley_checkPtr(rcv_len) || Finley_checkPtr(snd_len) || Finley_checkPtr(mpi_requests) || Finley_checkPtr(mpi_stati) || 
     Finley_checkPtr(mpi_stati) )) {

    memset(rcv_len,0,sizeof(dim_t)*mpiSize);
    #pragma omp parallel 
    {
        #pragma omp for private(i) schedule(static)
        for (i=0;i<len_loc_dof;++i) locDOFMask[i]=UNUSED;
        #pragma omp for private(i) schedule(static)
        for (i=0;i<numNodes;++i) nodeMask[i]=UNUSED;
        #pragma omp for private(i,k) schedule(static)
        for (i=0;i<numNodes;++i) {
           k=globalDOFIndex[i];
           if (k>-1) {
              locDOFMask[k-min_DOF]=UNUSED-1;
              #ifdef BOUNDS_CHECK
              if ((k-min_DOF) >= len_loc_dof) { printf("BOUNDS_CHECK %s %d i=%d k=%d min_DOF=%d\n", __FILE__, __LINE__, i, k, min_DOF); exit(1); }
              #endif
           }
        }

        #pragma omp for private(i) schedule(static)
        for (i=myFirstDOF-min_DOF;i<myLastDOF-min_DOF;++i) {
          locDOFMask[i]=i-myFirstDOF+min_DOF;
          #ifdef BOUNDS_CHECK
          if (i < 0 || i >= len_loc_dof) { printf("BOUNDS_CHECK %s %d i=%d\n", __FILE__, __LINE__, i); exit(1); }
          #endif
        }
    }

    numNeighbors=0;
    n=0;
    lastn=n;
    for (p=p_min;p<=p_max;++p) {
       firstDOF=MAX(min_DOF,dof_distribution->first_component[p]);
       lastDOF=MIN(max_DOF+1,dof_distribution->first_component[p+1]);
       if (p != myRank) {
           for (i=firstDOF-min_DOF;i<lastDOF-min_DOF;++i) {
                #ifdef BOUNDS_CHECK
                if (i < 0 || i >= len_loc_dof) { printf("BOUNDS_CHECK %s %d p=%d i=%d\n", __FILE__, __LINE__, p, i); exit(1); }
                #endif
                if (locDOFMask[i] == UNUSED-1) {
                   locDOFMask[i]=myLastDOF-myFirstDOF+n;
                   wanted_DOFs[n]=i+min_DOF;
                   ++n;
                }
           }
           if (n>lastn) {
              rcv_len[p]=n-lastn;
              neighbor[numNeighbors]=p;
              #ifdef BOUNDS_CHECK
              if (numNeighbors < 0 || numNeighbors >= mpiSize+1) { printf("BOUNDS_CHECK %s %d p=%d numNeighbors=%d n=%d\n", __FILE__, __LINE__, p, numNeighbors, n); exit(1); }
              #endif
              offsetInShared[numNeighbors]=lastn;
              numNeighbors++;
              lastn=n;
           }
        }
    }
#ifdef BOUNDS_CHECK
    if (numNeighbors < 0 || numNeighbors >= mpiSize+1) { printf("BOUNDS_CHECK %s %d numNeighbors=%d\n", __FILE__, __LINE__, numNeighbors); exit(1); }
#endif
    offsetInShared[numNeighbors]=lastn;

    /* assign new DOF labels to nodes */
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<numNodes;++i) {
       k=globalDOFIndex[i];
       if (k>-1) nodeMask[i]=locDOFMask[k-min_DOF];
    }

    /* now we can set the mapping from nodes to local DOFs */
    this_mapping=Finley_NodeMapping_alloc(numNodes,nodeMask,UNUSED);
    /* define how to get DOF values for controlled bu other processors */
    #ifdef BOUNDS_CHECK
    for (i=0;i<offsetInShared[numNeighbors];++i) {
      if (i < 0 || i >= numNodes*(p_max-p_min+1)) { printf("BOUNDS_CHECK %s %d i=%d\n", __FILE__, __LINE__, i); exit(1); }
    }
    #endif
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<offsetInShared[numNeighbors];++i) shared[i]=myLastDOF-myFirstDOF+i;

    rcv_shcomp=Paso_SharedComponents_alloc(myLastDOF-myFirstDOF,numNeighbors,neighbor,shared,offsetInShared,1,0,mpi_info);

    /*
     *    now we build the sender
     */
    #ifdef PASO_MPI
         MPI_Alltoall(rcv_len,1,MPI_INT,snd_len,1,MPI_INT,mpi_info->comm);
    #else
        for (p=0;p<mpiSize;++p) snd_len[p]=rcv_len[p];
    #endif
    count=0;
    for (p=0;p<rcv_shcomp->numNeighbors;p++) {
       #ifdef PASO_MPI
       MPI_Isend(&(wanted_DOFs[rcv_shcomp->offsetInShared[p]]), rcv_shcomp->offsetInShared[p+1]-rcv_shcomp->offsetInShared[p],
                  MPI_INT,rcv_shcomp->neighbor[p],mpi_info->msg_tag_counter+myRank,mpi_info->comm,&mpi_requests[count]);
        #endif
        count++;
    }
    n=0;
    numNeighbors=0;
    for (p=0;p<mpiSize;p++) {
         if (snd_len[p] > 0) {
            #ifdef PASO_MPI
            MPI_Irecv(&(shared[n]),snd_len[p],
                      MPI_INT, p, mpi_info->msg_tag_counter+p, mpi_info->comm, &mpi_requests[count]);
            #endif
            count++;
            neighbor[numNeighbors]=p;
            offsetInShared[numNeighbors]=n;
            numNeighbors++;
            n+=snd_len[p];
         }
    }
    mpi_info->msg_tag_counter+=mpi_info->size;
    offsetInShared[numNeighbors]=n;
    #ifdef PASO_MPI
    MPI_Waitall(count,mpi_requests,mpi_stati);
    #endif
    /* map global ids to local id's */
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<offsetInShared[numNeighbors];++i) {
        shared[i]=locDOFMask[shared[i]-min_DOF];
    }

    snd_shcomp=Paso_SharedComponents_alloc(myLastDOF-myFirstDOF,numNeighbors,neighbor,shared,offsetInShared,1,0,dof_distribution->mpi_info);

    if (Finley_noError()) this_connector=Paso_Connector_alloc(snd_shcomp,rcv_shcomp);
    /* assign new DOF labels to nodes */
    Paso_SharedComponents_free(rcv_shcomp);
    Paso_SharedComponents_free(snd_shcomp);
  }
  TMPMEMFREE(rcv_len);
  TMPMEMFREE(snd_len);
  TMPMEMFREE(mpi_requests);
  TMPMEMFREE(mpi_stati);
  TMPMEMFREE(wanted_DOFs);
  TMPMEMFREE(nodeMask);
  TMPMEMFREE(neighbor);
  TMPMEMFREE(shared);
  TMPMEMFREE(offsetInShared);
  TMPMEMFREE(locDOFMask);
  if (Finley_noError()) {
     if (use_reduced_elements) {
        in->Nodes->reducedDegreesOfFreedomMapping=this_mapping;
        in->Nodes->reducedDegreesOfFreedomConnector=this_connector;
     } else {
        in->Nodes->degreesOfFreedomMapping=this_mapping;
        in->Nodes->degreesOfFreedomConnector=this_connector;
    }
  } else {
     Finley_NodeMapping_free(this_mapping);
     Paso_Connector_free(this_connector);

  }
}
 
void Finley_Mesh_createMappings(Finley_Mesh* mesh, index_t* distribution) {
  int i;
  index_t *maskReducedNodes=NULL, *indexReducedNodes=NULL;
  dim_t numReducedNodes;

  maskReducedNodes=TMPMEMALLOC(mesh->Nodes->numNodes,index_t);
  indexReducedNodes=TMPMEMALLOC(mesh->Nodes->numNodes,index_t);

  if (! ( Finley_checkPtr(maskReducedNodes) || Finley_checkPtr(indexReducedNodes) ) ) {
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<mesh->Nodes->numNodes;++i) maskReducedNodes[i]=-1;
    Finley_Mesh_markNodes(maskReducedNodes,0,mesh,TRUE);
    numReducedNodes=Finley_Util_packMask(mesh->Nodes->numNodes,maskReducedNodes,indexReducedNodes);
    if (Finley_noError()) Finley_Mesh_createNodeFileMappings(mesh,numReducedNodes,indexReducedNodes,distribution);
  }

  TMPMEMFREE(maskReducedNodes);
  TMPMEMFREE(indexReducedNodes);
}

void Finley_Mesh_createNodeFileMappings(Finley_Mesh* in, dim_t numReducedNodes, index_t* indexReducedNodes, index_t* dof_first_component) {


  index_t myFirstDOF, myLastDOF, myFirstNode, myLastNode, *reduced_dof_first_component=NULL, *nodeMask=NULL,
         *reduced_nodes_first_component=NULL, *nodes_first_component=NULL,k,
         *maskMyReducedDOF=NULL, *indexMyReducedDOF=NULL, *maskMyReducedNodes=NULL, *indexMyReducedNodes=NULL;
  dim_t myNumDOF, myNumNodes, myNumReducedNodes, myNumReducedDOF, globalNumReducedNodes, globalNumReducedDOF,i,mpiSize, minGlobalNodeIndex,maxGlobalNodeIndex;
  Paso_MPI_rank myRank;

  mpiSize=in->Nodes->MPIInfo->size;
  myRank=in->Nodes->MPIInfo->rank;
  /* mark the nodes used by the reduced mesh */

  reduced_dof_first_component=TMPMEMALLOC(mpiSize+1,index_t);
  reduced_nodes_first_component=TMPMEMALLOC(mpiSize+1,index_t);
  nodes_first_component=TMPMEMALLOC(mpiSize+1,index_t);

  if (! ( Finley_checkPtr(reduced_dof_first_component) || Finley_checkPtr(reduced_nodes_first_component) || Finley_checkPtr(nodes_first_component)  ) ) {

     Finley_NodeFile_setGlobalNodeIDIndexRange(&minGlobalNodeIndex,&maxGlobalNodeIndex,in->Nodes);
     Paso_MPIInfo_setDistribution(in->Nodes->MPIInfo,minGlobalNodeIndex,maxGlobalNodeIndex,nodes_first_component);

     myFirstDOF=dof_first_component[myRank];
     myLastDOF=dof_first_component[myRank+1];
     myNumDOF=myLastDOF-myFirstDOF;
     myFirstNode=nodes_first_component[myRank];
     myLastNode=nodes_first_component[myRank+1];
     myNumNodes=myLastNode-myFirstNode;

     maskMyReducedDOF=TMPMEMALLOC(myNumDOF,index_t);
     indexMyReducedDOF=TMPMEMALLOC(myNumDOF,index_t);
     maskMyReducedNodes=TMPMEMALLOC(myNumNodes,index_t);
     indexMyReducedNodes=TMPMEMALLOC(myNumNodes,index_t);

     if (! ( Finley_checkPtr(maskMyReducedDOF) || Finley_checkPtr(indexMyReducedDOF) || Finley_checkPtr(maskMyReducedNodes) || Finley_checkPtr(indexMyReducedNodes)  ) ) {

        #pragma omp parallel private(i)
        {
            #pragma omp for schedule(static)
            for (i=0;i<myNumNodes;++i) maskMyReducedNodes[i]=-1;
            #pragma omp for schedule(static)
            for (i=0;i<myNumDOF;++i) maskMyReducedDOF[i]=-1;
            #pragma omp for private(k) schedule(static)
            for (i=0;i<numReducedNodes;++i) {
               k=in->Nodes->globalNodesIndex[indexReducedNodes[i]];
               if ( (k>=myFirstNode) && (myLastNode>k) ) maskMyReducedNodes[k-myFirstNode]=i;
               k=in->Nodes->globalDegreesOfFreedom[indexReducedNodes[i]];
               if ( (k>=myFirstDOF) && (myLastDOF>k) ) {
                  maskMyReducedDOF[k-myFirstDOF]=i;
               }
            }
        }
        myNumReducedNodes=Finley_Util_packMask(myNumNodes,maskMyReducedNodes,indexMyReducedNodes);
        myNumReducedDOF=Finley_Util_packMask(myNumDOF,maskMyReducedDOF,indexMyReducedDOF);
        
        #ifdef PASO_MPI
           MPI_Allgather(&myNumReducedNodes,1,MPI_INT,reduced_nodes_first_component,1,MPI_INT,in->Nodes->MPIInfo->comm);
           MPI_Allgather(&myNumReducedDOF,1,MPI_INT,reduced_dof_first_component,1,MPI_INT,in->Nodes->MPIInfo->comm);
        #else
           reduced_nodes_first_component[0]=myNumReducedNodes;
           reduced_dof_first_component[0]=myNumReducedDOF;
        #endif
        globalNumReducedNodes=0;
        globalNumReducedDOF=0;
        for (i=0;i<mpiSize;++i) {
            k=reduced_nodes_first_component[i];
            reduced_nodes_first_component[i]=globalNumReducedNodes;
            globalNumReducedNodes+=k;
    
            k=reduced_dof_first_component[i];
            reduced_dof_first_component[i]=globalNumReducedDOF;
            globalNumReducedDOF+=k;
        }
        reduced_nodes_first_component[mpiSize]=globalNumReducedNodes;
        reduced_dof_first_component[mpiSize]=globalNumReducedDOF;
        /* ==== distribution of Nodes ===============================*/
        Paso_MPIInfo_setDistribution(in->Nodes->MPIInfo,minGlobalNodeIndex,maxGlobalNodeIndex,nodes_first_component);
        in->Nodes->nodesDistribution=Paso_Distribution_alloc(in->Nodes->MPIInfo,nodes_first_component,1,0);
    
        /* ==== distribution of Nodes ===============================*/
        in->Nodes->degreesOfFreedomDistribution=Paso_Distribution_alloc(in->Nodes->MPIInfo,dof_first_component,1,0);
    
        /* ==== distribution of reduced Nodes ===============================*/
        reduced_nodes_first_component[mpiSize]=globalNumReducedNodes;
        in->Nodes->reducedNodesDistribution=Paso_Distribution_alloc(in->Nodes->MPIInfo,reduced_nodes_first_component,1,0);
    
        /* ==== distribution of reduced DOF ===============================*/
        in->Nodes->reducedDegreesOfFreedomDistribution=Paso_Distribution_alloc(in->Nodes->MPIInfo,reduced_dof_first_component,1,0);
     }
     TMPMEMFREE(maskMyReducedDOF);
     TMPMEMFREE(indexMyReducedDOF);
     TMPMEMFREE(maskMyReducedNodes);
     TMPMEMFREE(indexMyReducedNodes);
  }
  TMPMEMFREE(reduced_dof_first_component);
  TMPMEMFREE(reduced_nodes_first_component);
  TMPMEMFREE(nodes_first_component);

  nodeMask=TMPMEMALLOC(in->Nodes->numNodes,index_t);
  if (! Finley_checkPtr(nodeMask) && Finley_noError()) {

    /* ==== nodes mapping which is a dummy structure ======== */
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<in->Nodes->numNodes;++i) nodeMask[i]=i;
    in->Nodes->nodesMapping=Finley_NodeMapping_alloc(in->Nodes->numNodes,nodeMask,UNUSED);

    /* ==== mapping between nodes and reduced nodes ========== */
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<in->Nodes->numNodes;++i) nodeMask[i]=UNUSED;
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<numReducedNodes;++i) nodeMask[indexReducedNodes[i]]=i;
    in->Nodes->reducedNodesMapping=Finley_NodeMapping_alloc(in->Nodes->numNodes,nodeMask,UNUSED);

  }
  TMPMEMFREE(nodeMask);
  /* ==== mapping between nodes and DOFs + DOF connector ========== */
  if ( Finley_noError()) Mesh_createDOFMappingAndCoupling(in,FALSE);
  /* ==== mapping between nodes and reduced DOFs + reduced DOF connector ========== */
  if ( Finley_noError()) Mesh_createDOFMappingAndCoupling(in,TRUE);

  /* get the Ids for DOFs and reduced nodes */
  if (Finley_noError()) {
     #pragma omp parallel private(i)
     {
         #pragma omp for
         for (i=0;i<in->Nodes->reducedNodesMapping->numTargets;++i) in->Nodes->reducedNodesId[i]=in->Nodes->Id[in->Nodes->reducedNodesMapping->map[i]];
         #pragma omp for
         for (i=0;i<in->Nodes->degreesOfFreedomMapping->numTargets;++i) in->Nodes->degreesOfFreedomId[i]=in->Nodes->Id[in->Nodes->degreesOfFreedomMapping->map[i]];
         #pragma omp for
         for (i=0;i<in->Nodes->reducedDegreesOfFreedomMapping->numTargets;++i) in->Nodes->reducedDegreesOfFreedomId[i]=in->Nodes->Id[in->Nodes->reducedDegreesOfFreedomMapping->map[i]];
     }
  } else {
    Finley_NodeMapping_free(in->Nodes->nodesMapping);
    Finley_NodeMapping_free(in->Nodes->reducedNodesMapping);
    Finley_NodeMapping_free(in->Nodes->degreesOfFreedomMapping);
    Finley_NodeMapping_free(in->Nodes->reducedDegreesOfFreedomMapping);
    Paso_Distribution_free(in->Nodes->nodesDistribution);
    Paso_Distribution_free(in->Nodes->reducedNodesDistribution);
    Paso_Distribution_free(in->Nodes->degreesOfFreedomDistribution);
    Paso_Distribution_free(in->Nodes->reducedDegreesOfFreedomDistribution);
    Paso_Connector_free(in->Nodes->degreesOfFreedomConnector);
    Paso_Connector_free(in->Nodes->reducedDegreesOfFreedomConnector);
    in->Nodes->nodesMapping=NULL;
    in->Nodes->reducedNodesMapping=NULL;
    in->Nodes->degreesOfFreedomMapping=NULL;
    in->Nodes->reducedDegreesOfFreedomMapping=NULL;
    in->Nodes->nodesDistribution=NULL;
    in->Nodes->reducedNodesDistribution=NULL;
    in->Nodes->degreesOfFreedomDistribution=NULL;
    in->Nodes->reducedDegreesOfFreedomDistribution=NULL;
    in->Nodes->degreesOfFreedomConnector=NULL;
    in->Nodes->reducedDegreesOfFreedomConnector=NULL;
  }
}

