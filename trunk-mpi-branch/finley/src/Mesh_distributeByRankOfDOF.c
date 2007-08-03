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

/*   Finley: Mesh: this will redistribute the Nodes and Elements including overlap and */
/*   will create an element coloring but will not create any mappings. The distribution is done */
/*   according to the rank assigned to each node/DOF via mpiRankOfDOF (length in->Nodes->numNodes) */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id:$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Finley_Mesh_distributeByRankOfDOF(Finley_Mesh* self, Paso_MPI_rank* mpiRankOfDOF) {

     if (self==NULL) return;

     /* first the elements are redistributed according to mpiRankOfDOF */
     /* at the input the Node tables refering to a the local labeling of the nodes */
     /* while at the output they refer to the global labeling which is rectified in the next step */
     Finley_ElementFile_distributeByRankOfDOF(self->Elements,mpiRankOfDOF, self->Nodes->Id);
     if (!Finley_noError()) return;

     Finley_ElementFile_distributeByRankOfDOF(self->FaceElements,mpiRankOfDOF, self->Nodes->Id);
     if (!Finley_noError()) return;

     Finley_ElementFile_distributeByRankOfDOF(self->ContactElements,mpiRankOfDOF, self->Nodes->Id);
     if (!Finley_noError()) return;

     Finley_ElementFile_distributeByRankOfDOF(self->Points,mpiRankOfDOF, self->Nodes->Id);
     if (!Finley_noError()) return;
     /* resolve the node ids */
     Finley_Mesh_resolveNodeIds(self);
     if (!Finley_noError()) return;

     /* create a local labeling of the DOFs */
  /* allocate mappings for new local node labeling to global node labeling (newLocalToGlobalNodeLabels)
     and global node labeling to the new local node labeling (globalToNewLocalNodeLabels[i-min_id] is the 
     new local id of global node i) */
  
  len=max_id-min_id+1;
  globalToNewLocalNodeLabels=TMPMEMALLOC(len,index_t); /* local mask for used nodes */
  newLocalToGlobalNodeLabels=TMPMEMALLOC(len,index_t);
  if (! ( (Finley_checkPtr(globalToNewLocalNodeLabels) && Finley_checkPtr(newLocalToGlobalNodeLabels) ) ) ) {

       #pragma omp parallel
       {
           #pragma omp for private(n) schedule(static)
           for (n=0;n<in->Nodes->numNodes;n++) newLocalToGlobalNodeLabels[n]=-1;
           #pragma omp for private(n) schedule(static)
           for (n=0;n<len;n++) globalToNewLocalNodeLabels[n]=-1;
       }

       /*  mark the nodes referred by elements in globalToNewLocalNodeLabels which is currently used as a mask: */

       Finley_Mesh_markNodes(globalToNewLocalNodeLabels,min_id,in,FALSE);


     
     /* create a element pre-coloring */

     /* refine element coloring */

     return;
}

