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

/**************************************************************/

void Finley_NodeFile_createMappings(Finley_NodeFile* in, dim_t numReducedNodes, index_t* indexReducedNodes, index_t* dof_first_component) {


  index_t myFirstDOF, myLastDOF, myFirstNode, myLastNode, *reduced_dof_first_component=NULL,
             *reduced_nodes_first_component=NULL, *nodes_first_component=NULL,k,
             *maskMyReducedDOF=NULL, *indexMyReducedDOF=NULL, *maskMyReducedNodes=NULL, *indexMyReducedNodes=NULL,
             firstDOF, lastDOF, myFirstReducedDOF, myLastReducedDOF, *nodeMask=NULL, min_DOF, max_DOF;
  dim_t myNumDOF, myNumNodes, myNumReducedNodes, myNumReducedDOF, globalNumReducedNodes, globalNumReducedDOF,i,mpiSize, globalNumNodes, n;
  Paso_MPI_rank myRank,p,p_min,p_max;

  mpiSize=in->MPIInfo->size;
  myRank=in->MPIInfo->rank;
  /* mark the nodes used by the reduced mesh */

  reduced_dof_first_component=TMPMEMALLOC(mpiSize+1,index_t);
  reduced_nodes_first_component=TMPMEMALLOC(mpiSize+1,index_t);
  nodes_first_component=TMPMEMALLOC(mpiSize+1,index_t);

  if (! ( Finley_checkPtr(reduced_dof_first_component) || Finley_checkPtr(reduced_nodes_first_component) || Finley_checkPtr(nodes_first_component)  ) ) {

     globalNumNodes=Finley_NodeFile_maxGlobalNodeIDIndex(in);
     Paso_MPIInfo_setDistribution(in->MPIInfo,0,globalNumNodes-1,nodes_first_component);

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
               k=in->globalNodesIndex[indexReducedNodes[i]];
               if ( (k>=myFirstNode) && (myLastNode>k) ) maskMyReducedNodes[k-myFirstNode]=i;
               k=in->globalDegreesOfFreedom[indexReducedNodes[i]];
               if ( (k>=myFirstDOF) && (myLastDOF>k) ) maskMyReducedDOF[k-myFirstDOF]=i;
            }
        }
        myNumReducedNodes=Finley_Util_packMask(myNumNodes,maskMyReducedNodes,indexMyReducedNodes);
        myNumReducedDOF=Finley_Util_packMask(myNumDOF,maskMyReducedDOF,indexMyReducedDOF);
        
        #ifdef PASO_MPI
           MPI_Allgather(&myNumReducedNodes,1,MPI_INT,reduced_nodes_first_component,1,MPI_INT,in->MPIInfo->comm);
           MPI_Allgather(&myNumReducedDOF,1,MPI_INT,reduced_dof_first_component,1,MPI_INT,in->MPIInfo->comm);
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
        Paso_MPIInfo_setDistribution(in->MPIInfo,0,globalNumNodes-1,nodes_first_component);
        in->nodesDistribution=Paso_Distribution_alloc(in->MPIInfo,nodes_first_component,1,0);
    
        /* ==== distribution of reduced Nodes ===============================*/
        reduced_nodes_first_component[mpiSize]=globalNumReducedNodes;
        in->reducedDegreesOfFreedomDistribution=Paso_Distribution_alloc(in->MPIInfo,reduced_nodes_first_component,1,0);
    
        /* ==== distribution of Nodes ===============================*/
        in->degreesOfFreedomDistribution=Paso_Distribution_alloc(in->MPIInfo,dof_first_component,1,0);
    
        /* ==== distribution of reduced DOF ===============================*/
        in->reducedDegreesOfFreedomDistribution=Paso_Distribution_alloc(in->MPIInfo,reduced_dof_first_component,1,0);
     }
     TMPMEMFREE(maskMyReducedDOF);
     TMPMEMFREE(indexMyReducedDOF);
     TMPMEMFREE(maskMyReducedNodes);
     TMPMEMFREE(indexMyReducedNodes);
  }
  TMPMEMFREE(reduced_dof_first_component);
  TMPMEMFREE(reduced_nodes_first_component);
  TMPMEMFREE(nodes_first_component);

  nodeMask=TMPMEMALLOC(in->numNodes,index_t);
  if (! Finley_checkPtr(nodeMask) && Finley_noError()) {

    /* ==== nodes mapping which is a dummy structure ======== */
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<in->numNodes;++i) nodeMask[i]=i;
    in->nodesMapping=Finley_NodeMapping_alloc(in->numNodes,nodeMask,UNUSED);

    /* ==== mapping between nodes and reduced nodes ========== */
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<in->numNodes;++i) nodeMask[i]=UNUSED;
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<numReducedNodes;++i) nodeMask[indexReducedNodes[i]]=i;
    in->reducedNodesMapping=Finley_NodeMapping_alloc(in->numNodes,nodeMask,UNUSED);

    /* ==== mapping between nodes and DOFs ========== */
    /* the trick is to count the local DOFs first and then the romote DOFs grouped by processor */
    p_min=mpiSize;
    p_max=-1;
    Finley_NodeFile_setDOFRange(&min_DOF,&max_DOF,in);

    for (p=0; p<in->MPIInfo->size; ++p) {
         if (in->degreesOfFreedomDistribution->first_component[p]<=min_DOF) p_min=p;
         if (in->degreesOfFreedomDistribution->first_component[p]<=max_DOF) p_max=p;
    }
    n=0;
    for (i=0;i<in->numNodes;++i) {
       k=in->globalDegreesOfFreedom[i];
       if ( (k>=myFirstDOF) && (myLastDOF>k) ) {
          nodeMask[i]=n;
          n++;
       }
    }
    for (p=p_min;p<=p_max;++p) {
       firstDOF=in->degreesOfFreedomDistribution->first_component[p];
       lastDOF=in->degreesOfFreedomDistribution->first_component[p+1];
       if (p != myRank) {
           for (i=0;i<in->numNodes;++i) {
               k=in->globalDegreesOfFreedom[i];
               if ( (k>=firstDOF) && (lastDOF>k) ) {
                 nodeMask[i]=n;
                 n++;
               }
           }
       }
    }
    in->degreesOfFreedomMapping=Finley_NodeMapping_alloc(in->numNodes,nodeMask,UNUSED);

    /* ==== mapping between nodes and DOFs ========== */
    /* the trick is to count the local DOFs first and then the romote DOFs grouped by processor */
    p_min=mpiSize;
    p_max=-1;
    Finley_NodeFile_setReducedDOFRange(&min_DOF,&max_DOF,in); 


    for (p=0; p<in->MPIInfo->size; ++p) {
         if (in->reducedDegreesOfFreedomDistribution->first_component[p]<=min_DOF) p_min=p;
         if (in->reducedDegreesOfFreedomDistribution->first_component[p]<=max_DOF) p_max=p;
    }
    myFirstReducedDOF=in->reducedDegreesOfFreedomDistribution->first_component[myRank];
    myLastReducedDOF=in->reducedDegreesOfFreedomDistribution->first_component[myRank+1];
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<in->numNodes;++i) nodeMask[i]=UNUSED;

    n=0;
    for (i=0;i<in->numNodes;++i) {
       k=in->globalReducedDOFIndex[i];
       if ( (k>=myFirstReducedDOF) && (myLastReducedDOF>k) ) {
          nodeMask[i]=n;
          n++;
       }
    }
    for (p=p_min;p<=p_max;++p) {
       firstDOF=in->reducedDegreesOfFreedomDistribution->first_component[p];
       lastDOF=in->reducedDegreesOfFreedomDistribution->first_component[p+1];
       if (p != myRank) {
           for (i=0;i<in->numNodes;++i) {
               k=in->globalReducedDOFIndex[i];
               if ( (k>=firstDOF) && (lastDOF>k) ) {
                 nodeMask[i]=n;
                 n++;
               }
           }
       }
    }
    in->reducedDegreesOfFreedomMapping=Finley_NodeMapping_alloc(in->numNodes,nodeMask,UNUSED);
  }
  TMPMEMFREE(nodeMask);
   
/*
 Paso_Coupler* degreesOfFreedomCoupler;
 Paso_Coupler *reducedDegreesOfFreedomCoupler;
*/
  if (Finley_noError()) {
     #pragma omp parallel private(i)
     {
         #pragma omp for
         for (i=0;i<in->reducedNodesMapping->numTargets;++i) in->reducedNodesId[i]=in->Id[in->reducedNodesMapping->map[i]];
         #pragma omp for
         for (i=0;i<in->degreesOfFreedomMapping->numTargets;++i) in->degreesOfFreedomId[i]=in->Id[in->degreesOfFreedomMapping->map[i]];
         #pragma omp for
         for (i=0;i<in->reducedDegreesOfFreedomMapping->numTargets;++i) in->reducedDegreesOfFreedomId[i]=in->Id[in->reducedDegreesOfFreedomMapping->map[i]];
     }
  } else {
    Finley_NodeMapping_free(in->nodesMapping);
    Finley_NodeMapping_free(in->reducedNodesMapping);
    Finley_NodeMapping_free(in->degreesOfFreedomMapping);
    Finley_NodeMapping_free(in->reducedDegreesOfFreedomMapping);
    Paso_Distribution_free(in->nodesDistribution);
    Paso_Distribution_free(in->reducedNodesDistribution);
    Paso_Distribution_free(in->degreesOfFreedomDistribution);
    Paso_Distribution_free(in->reducedDegreesOfFreedomDistribution);
    Paso_Coupler_free(in->degreesOfFreedomCoupler);
    Paso_Coupler_free(in->reducedDegreesOfFreedomCoupler);
  }
}
