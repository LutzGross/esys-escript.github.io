
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Finley: Mesh */

/*   at input the element nodes refers to the numbering defined the global Id assigned to the nodes in the */
/*   NodeFile. It is also not ensured that all nodes refered by an element is actually available */
/*   on the process.  At the output, a local node labeling is used and all nodes are available */
/*   In particular the numbering of the element nodes is between 0 and in->NodeFile->numNodes */
/*   The function does not create a distribution of the degrees of freedom. */

/**************************************************************/

#include "Mesh.h"
#include "Util.h"

/**************************************************************/

void  Finley_Mesh_resolveNodeIds(Finley_Mesh* in) {

  index_t min_id, max_id,  min_id2, max_id2, global_min_id, global_max_id, 
          *globalToNewLocalNodeLabels=NULL, *newLocalToGlobalNodeLabels=NULL;
  dim_t len, n, newNumNodes, numDim;
  Finley_NodeFile *newNodeFile=NULL;
  #ifdef PASO_MPI
  index_t id_range[2], global_id_range[2];
  #endif 
  numDim=Finley_Mesh_getDim(in);
  /*  find the minimum and maximum id used by elements: */
  min_id=INDEX_T_MAX;
  max_id=-INDEX_T_MAX;
  Finley_ElementFile_setNodeRange(&min_id2,&max_id2,in->Elements);
  max_id=MAX(max_id,max_id2);
  min_id=MIN(min_id,min_id2);
  Finley_ElementFile_setNodeRange(&min_id2,&max_id2,in->FaceElements);
  max_id=MAX(max_id,max_id2);
  min_id=MIN(min_id,min_id2);
  Finley_ElementFile_setNodeRange(&min_id2,&max_id2,in->ContactElements);
  max_id=MAX(max_id,max_id2);
  min_id=MIN(min_id,min_id2);
  Finley_ElementFile_setNodeRange(&min_id2,&max_id2,in->Points);
  max_id=MAX(max_id,max_id2);
  min_id=MIN(min_id,min_id2);
  #ifdef PASO_MPI
     id_range[0]=-min_id;
     id_range[1]=max_id;
     MPI_Allreduce( id_range, global_id_range, 2, MPI_INT, MPI_MAX, in->MPIInfo->comm );
     global_min_id=-global_id_range[0];
     global_max_id=global_id_range[1];
  #else
     global_min_id=min_id;
     global_max_id=max_id;
  #endif
  #ifdef Finley_TRACE
  printf("Node id range used by elements is %d:%d\n",global_min_id,global_max_id);
  #endif
  if (min_id>max_id) {
     max_id=-1;
     min_id=0;
  }
  
  /* allocate mappings for new local node labeling to global node labeling (newLocalToGlobalNodeLabels)
     and global node labeling to the new local node labeling (globalToNewLocalNodeLabels[i-min_id] is the 
     new local id of global node i) */
  len=(max_id>=min_id) ? max_id-min_id+1 : 0 ;
  globalToNewLocalNodeLabels=TMPMEMALLOC(len,index_t); /* local mask for used nodes */
  newLocalToGlobalNodeLabels=TMPMEMALLOC(len,index_t);
  if (! ( (Finley_checkPtr(globalToNewLocalNodeLabels) && Finley_checkPtr(newLocalToGlobalNodeLabels) ) ) ) {

       #pragma omp parallel
       {
           #pragma omp for private(n) schedule(static)
           for (n=0;n<len;n++) newLocalToGlobalNodeLabels[n]=-1;
           #pragma omp for private(n) schedule(static)
           for (n=0;n<len;n++) globalToNewLocalNodeLabels[n]=-1;
       }

       /*  mark the nodes referred by elements in globalToNewLocalNodeLabels which is currently used as a mask: */
       Finley_Mesh_markNodes(globalToNewLocalNodeLabels,min_id,in,FALSE);

       /* create a local labeling newLocalToGlobalNodeLabels of the local nodes by packing the mask globalToNewLocalNodeLabels*/

       newNumNodes=Finley_Util_packMask(len,globalToNewLocalNodeLabels,newLocalToGlobalNodeLabels);

       /* invert the new labeling and shift the index newLocalToGlobalNodeLabels to global node ids */
       #pragma omp parallel for private(n) schedule(static)
       for (n=0;n<newNumNodes;n++) {
              #ifdef BOUNDS_CHECK
                     if (n >= len || n < 0) { printf("BOUNDS_CHECK %s %d n=%d\n", __FILE__, __LINE__, n); exit(1); }
                     if (newLocalToGlobalNodeLabels[n] >= len || newLocalToGlobalNodeLabels[n] < 0) { printf("BOUNDS_CHECK %s %d n=%d\n", __FILE__, __LINE__, n); exit(1); }
              #endif
              globalToNewLocalNodeLabels[newLocalToGlobalNodeLabels[n]]=n;
              newLocalToGlobalNodeLabels[n]+=min_id;
        }
        /* create a new table */
        newNodeFile=Finley_NodeFile_alloc(numDim,in->MPIInfo);
        if (Finley_noError()) {
           Finley_NodeFile_allocTable(newNodeFile,newNumNodes);
        }
        if (Finley_noError()) {
            Finley_NodeFile_gather_global(newLocalToGlobalNodeLabels,in->Nodes, newNodeFile);
        }
        if (Finley_noError()) {
           Finley_NodeFile_free(in->Nodes);
           in->Nodes=newNodeFile;
           /*  relable nodes of the elements: */
           Finley_Mesh_relableElementNodes(globalToNewLocalNodeLabels,min_id,in);
        }
  }
  TMPMEMFREE(globalToNewLocalNodeLabels);
  TMPMEMFREE(newLocalToGlobalNodeLabels);
  if (! Finley_noError()) {
       Finley_NodeFile_free(newNodeFile);
  }
}
