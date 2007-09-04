/*
 ************************************************************
 *          Copyright 2007 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

/**************************************************************/

/*   Finley: NodeFile : creates the mappings using the indexReducedNodes */
/*                 no distribution is happening                          */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id:$ */

/**************************************************************/

#include "Mesh.h"
#define UNUSED -1
#define BOUNDS_CHECK 1 

/**************************************************************/

void Mesh_createDOFMappingAndCoupling(Finley_Mesh* in, bool_t use_reduced_elements) 
{
  index_t min_DOF, max_DOF, *shared=NULL, *offsetInShared=NULL, *locDOFMask=NULL, i, k, myFirstDOF, myLastDOF, *nodeMask=NULL, firstDOF, lastDOF, *globalDOFIndex;
  dim_t mpiSize, len_loc_dof, numNeighbors, n, lastn, numNodes;
  Paso_MPI_rank myRank,p,p_min,p_max, *neighbor=NULL;
  Paso_SharedComponents *rcv_shcomp=NULL, *snd_shcomp=NULL;
  Finley_NodeMapping *this_mapping=NULL;
  Paso_Coupler* this_coupler=NULL;
  Paso_Distribution* dof_distribution;
  
  numNodes=in->Nodes->numNodes;
  if (use_reduced_elements) {
    dof_distribution=in->Nodes->reducedDegreesOfFreedomDistribution;
    globalDOFIndex=in->Nodes->globalReducedDOFIndex;
    Finley_NodeFile_setReducedDOFRange(&min_DOF, &max_DOF,in->Nodes);
  } else {
    dof_distribution=in->Nodes->degreesOfFreedomDistribution;
    globalDOFIndex=in->Nodes->globalDegreesOfFreedom;
    Finley_NodeFile_setDOFRange(&min_DOF, &max_DOF,in->Nodes);
  }

  mpiSize=dof_distribution->mpi_info->size;
  myRank=dof_distribution->mpi_info->rank;

  min_DOF=Finley_Util_getFlaggedMinInt(1,numNodes,globalDOFIndex,-1);
  max_DOF=Finley_Util_getFlaggedMaxInt(1,numNodes,globalDOFIndex,-1);

  p_min=mpiSize;
  p_max=-1;

  for (p=0; p<mpiSize; ++p) {
      if (dof_distribution->first_component[p]<=min_DOF) p_min=p;
      if (dof_distribution->first_component[p]<=max_DOF) p_max=p;
  }

  len_loc_dof=max_DOF-min_DOF+1;
  myFirstDOF=Paso_Distribution_getFirstComponent(dof_distribution);
  myLastDOF=Paso_Distribution_getLastComponent(dof_distribution);

 
  nodeMask=TMPMEMALLOC(numNodes,index_t);
  neighbor=TMPMEMALLOC(mpiSize,Paso_MPI_rank);
  shared=TMPMEMALLOC(numNodes*(p_max-p_min+1),index_t);
  offsetInShared=TMPMEMALLOC(mpiSize+1,index_t);
  locDOFMask=TMPMEMALLOC(len_loc_dof, index_t);
  if (! ( Finley_checkPtr(neighbor) || Finley_checkPtr(shared) || Finley_checkPtr(offsetInShared) || Finley_checkPtr(locDOFMask) || Finley_checkPtr(nodeMask) ))  {


    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<len_loc_dof;++i) locDOFMask[i]=UNUSED;
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<numNodes;++i) nodeMask[i]=UNUSED;

    for (i=0;i<numNodes;++i) {
       k=globalDOFIndex[i];
#ifdef BOUNDS_CHECK
       if ((k-min_DOF) >= len_loc_dof) { printf("BOUNDS_CHECK %s:%i i=%i k=%i min_DOF=%i\n", __FILE__, __LINE__, i, k, min_DOF); exit(1); }
#endif
       if (k>-1) locDOFMask[k-min_DOF]=UNUSED-1;
    }

    for (i=myFirstDOF-min_DOF;i<myLastDOF-min_DOF;++i) {
      locDOFMask[i]=i-myFirstDOF+min_DOF;
#ifdef BOUNDS_CHECK
      if (i < 0 || i >= len_loc_dof) { printf("BOUNDS_CHECK %s:%i i=%i\n", __FILE__, __LINE__, i); exit(1); }
#endif
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
                   if (i < 0 || i >= len_loc_dof) { printf("BOUNDS_CHECK %s:%i p=%i i=%i\n", __FILE__, __LINE__, p, i); exit(1); }
#endif
                if (locDOFMask[i] == UNUSED-1) {
                   locDOFMask[i]=myLastDOF-myFirstDOF+n;
                   ++n;
                }
           }
           if (n>lastn) {
              neighbor[numNeighbors]=p;
#ifdef BOUNDS_CHECK
              if (numNeighbors < 0 || numNeighbors >= mpiSize+1) { printf("BOUNDS_CHECK %s:%i p=%i numNeighbors=%i n=%i\n", __FILE__, __LINE__, p, numNeighbors, n); exit(1); }
#endif
              offsetInShared[numNeighbors]=lastn;
              numNeighbors++;
              lastn=n;
           }
        }
    }
#ifdef BOUNDS_CHECK
    if (numNeighbors < 0 || numNeighbors >= mpiSize+1) { printf("BOUNDS_CHECK %s:%i numNeighbors=%i\n", __FILE__, __LINE__, numNeighbors); exit(1); }
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
      if (i < 0 || i >= numNodes*(p_max-p_min+1)) { printf("BOUNDS_CHECK %s:%i i=%i\n", __FILE__, __LINE__, i); exit(1); }
    }
#endif
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<offsetInShared[numNeighbors];++i) shared[i]=myLastDOF-myFirstDOF+i;

    rcv_shcomp=Paso_SharedComponents_alloc(numNeighbors,neighbor,shared,offsetInShared,1,0,dof_distribution->mpi_info);

    /* now it is determined which DOFs needs to be send off:*/
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<len_loc_dof;++i) locDOFMask[i]=UNUSED;
    n=0;
    numNeighbors=0;
    lastn=n;
    for (p=p_min;p<=p_max;++p) {
       firstDOF=dof_distribution->first_component[p];
       lastDOF=dof_distribution->first_component[p+1];
       if (p != myRank) {
           /* mark a DOF by p if it will be requested by processor p */
           Finley_Mesh_markDOFsConnectedToRange(locDOFMask,min_DOF,p,firstDOF,lastDOF,in,use_reduced_elements);

           for (i=myFirstDOF-min_DOF;i<myLastDOF-min_DOF;++i) {
                if (locDOFMask[i] == p) {
#ifdef BOUNDS_CHECK
		   if (n < 0 || n >= numNodes*(p_max-p_min+1)) { printf("BOUNDS_CHECK %s:%i p=%i i=%i n=%i\n", __FILE__, __LINE__, p, i, n); exit(1); }
#endif
                   shared[n]=i-myFirstDOF+min_DOF;
                   ++n;
                }
           }
           if (n>lastn) {
              neighbor[numNeighbors]=p;
#ifdef BOUNDS_CHECK
              if (numNeighbors < 0 || numNeighbors >= mpiSize+1) { printf("BOUNDS_CHECK %s:%i p=%i n=%i numNeighbors=%i\n", __FILE__, __LINE__, p, n, numNeighbors); exit(1); }
#endif
              offsetInShared[numNeighbors]=lastn;
              numNeighbors++;
              lastn=n;
           }
        }
    }
#ifdef BOUNDS_CHECK
    if (numNeighbors < 0 || numNeighbors >= mpiSize+1) { printf("BOUNDS_CHECK %s:%i numNeighbors=%i\n", __FILE__, __LINE__, numNeighbors); exit(1); }
#endif
    offsetInShared[numNeighbors]=lastn;
    snd_shcomp=Paso_SharedComponents_alloc(numNeighbors,neighbor,shared,offsetInShared,1,0,dof_distribution->mpi_info);

    if (Finley_noError()) this_coupler=Paso_Coupler_alloc(snd_shcomp,rcv_shcomp);
    /* assign new DOF labels to nodes */
    Paso_SharedComponents_free(rcv_shcomp);
    Paso_SharedComponents_free(snd_shcomp);
  }
  TMPMEMFREE(nodeMask);
  TMPMEMFREE(neighbor);
  TMPMEMFREE(shared);
  TMPMEMFREE(offsetInShared);
  TMPMEMFREE(locDOFMask);
  if (Finley_noError()) {
     if (use_reduced_elements) {
        in->Nodes->reducedDegreesOfFreedomMapping=this_mapping;
        in->Nodes->reducedDegreesOfFreedomCoupler=this_coupler;
     } else {
        in->Nodes->degreesOfFreedomMapping=this_mapping;
        in->Nodes->degreesOfFreedomCoupler=this_coupler;
    }
  } else {
     Finley_NodeMapping_free(this_mapping);
     Paso_Coupler_free(this_coupler);

  }
}
void Finley_Mesh_createNodeFileMappings(Finley_Mesh* in, dim_t numReducedNodes, index_t* indexReducedNodes, index_t* dof_first_component) {


  index_t myFirstDOF, myLastDOF, myFirstNode, myLastNode, *reduced_dof_first_component=NULL, *nodeMask=NULL,
         *reduced_nodes_first_component=NULL, *nodes_first_component=NULL,k,
         *maskMyReducedDOF=NULL, *indexMyReducedDOF=NULL, *maskMyReducedNodes=NULL, *indexMyReducedNodes=NULL;
  dim_t myNumDOF, myNumNodes, myNumReducedNodes, myNumReducedDOF, globalNumReducedNodes, globalNumReducedDOF,i,mpiSize, globalNumNodes, n, lastn,n0, numNeighbors;
  Paso_MPI_rank myRank;

  mpiSize=in->Nodes->MPIInfo->size;
  myRank=in->Nodes->MPIInfo->rank;
  /* mark the nodes used by the reduced mesh */

  reduced_dof_first_component=TMPMEMALLOC(mpiSize+1,index_t);
  reduced_nodes_first_component=TMPMEMALLOC(mpiSize+1,index_t);
  nodes_first_component=TMPMEMALLOC(mpiSize+1,index_t);

  if (! ( Finley_checkPtr(reduced_dof_first_component) || Finley_checkPtr(reduced_nodes_first_component) || Finley_checkPtr(nodes_first_component)  ) ) {

     globalNumNodes=Finley_NodeFile_maxGlobalNodeIDIndex(in->Nodes)+1;
     Paso_MPIInfo_setDistribution(in->Nodes->MPIInfo,0,globalNumNodes-1,nodes_first_component);

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
               if ( (k>=myFirstDOF) && (myLastDOF>k) ) maskMyReducedDOF[k-myFirstDOF]=i;
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
        Paso_MPIInfo_setDistribution(in->Nodes->MPIInfo,0,globalNumNodes-1,nodes_first_component);
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
  /* ==== mapping between nodes and DOFs + DOF coupler ========== */
  if ( Finley_noError()) Mesh_createDOFMappingAndCoupling(in,FALSE);
  /* ==== mapping between nodes and reduced DOFs + reduced DOF coupler ========== */
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
    Paso_Coupler_free(in->Nodes->degreesOfFreedomCoupler);
    Paso_Coupler_free(in->Nodes->reducedDegreesOfFreedomCoupler);
  }
}

