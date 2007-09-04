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

/*   Finley: Mesh: prepares the mesh for further calculations  */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Finley_Mesh_prepare(Finley_Mesh* in, bool_t optimize) {
     dim_t newGlobalNumDOFs=0, numReducedNodes=0,i;
     index_t* distribution=NULL, *maskReducedNodes=NULL, *indexReducedNodes=NULL;
     if (in==NULL) return;
     if (in->Nodes == NULL) return;

     /* first step is to distribute the elements according to a global distribution of DOF */

     distribution=TMPMEMALLOC(in->MPIInfo->size+1,index_t);
     if (! Finley_checkPtr(distribution)) {
        /* first we create dense labeling for the DOFs */
        newGlobalNumDOFs=Finley_NodeFile_createDenseDOFLabeling(in->Nodes);

        /* create a distribution of the global DOFs and determine
           the MPI_rank controling the DOFs on this processor      */
        Paso_MPIInfo_setDistribution(in->MPIInfo,0,newGlobalNumDOFs-1,distribution);

        /* now the mesh is re-distributed according to the mpiRankOfDOF vector */
        /* this will redistribute the Nodes and Elements including overlap and will create an element coloring 
           but will not create any mappings (see later in this function)                                   */
        if (Finley_noError()) Finley_Mesh_distributeByRankOfDOF(in,distribution);
     }

     /* at this stage we are able to start an optimization of the DOF distribution using ParaMetis */
     /* on return distribution is altered and new DOF ids have been assigned */
     if (Finley_noError() && optimize && in->MPIInfo->size>1) {

         Finley_Mesh_optimizeDOFDistribution(in,distribution);

         if (Finley_noError()) Finley_Mesh_distributeByRankOfDOF(in,distribution); 

     }
     /* now a local labeling of the DOF is introduced */
     if (Finley_noError() && optimize) {
       Finley_Mesh_optimizeDOFLabeling(in,distribution); 
     }
     /* rearrange elements with the attempt to bring elements closer to memory locations of the nodes (distributed shared memory!): */
     if (Finley_noError()) Finley_Mesh_optimizeElementOrdering(in);


     /* create the global indices */
     if (Finley_noError()) {
        maskReducedNodes=TMPMEMALLOC(in->Nodes->numNodes,index_t);
        indexReducedNodes=TMPMEMALLOC(in->Nodes->numNodes,index_t);
        if (! ( Finley_checkPtr(maskReducedNodes) ||  Finley_checkPtr(indexReducedNodes) ) ) {

          #pragma omp parallel for private(i) schedule(static)
          for (i=0;i<in->Nodes->numNodes;++i) maskReducedNodes[i]=-1;

          Finley_Mesh_markNodes(maskReducedNodes,0,in,1);
   
          numReducedNodes=Finley_Util_packMask(in->Nodes->numNodes,maskReducedNodes,indexReducedNodes);

          Finley_NodeFile_createDenseNodeLabeling(in->Nodes); 
          Finley_NodeFile_createDenseReducedNodeLabeling(in->Nodes,maskReducedNodes); 
          Finley_NodeFile_createDenseReducedDOFLabeling(in->Nodes,maskReducedNodes); 
          /* create the missing mappings */
          if (Finley_noError()) Finley_Mesh_createNodeFileMappings(in,numReducedNodes,indexReducedNodes,distribution);
        }

        TMPMEMFREE(maskReducedNodes);
        TMPMEMFREE(indexReducedNodes);
     }

     TMPMEMFREE(distribution);

     return;
}

/*                                                      */
/*  tries to reduce the coloring for all element files: */
/*                                                      */
void Finley_Mesh_createColoring(Finley_Mesh* in, index_t *node_localDOF_map) {
  if (Finley_noError()) Finley_ElementFile_createColoring(in->Elements,in->Nodes->numNodes,node_localDOF_map);
  if (Finley_noError()) Finley_ElementFile_createColoring(in->FaceElements,in->Nodes->numNodes,node_localDOF_map);
  if (Finley_noError()) Finley_ElementFile_createColoring(in->Points,in->Nodes->numNodes,node_localDOF_map);
  if (Finley_noError()) Finley_ElementFile_createColoring(in->ContactElements,in->Nodes->numNodes,node_localDOF_map);
}
/*                                                                    */
/*  redistribute elements to minimize communication during assemblage */
/*                                                                    */
void Finley_Mesh_optimizeElementOrdering(Finley_Mesh* in) {
  if (Finley_noError()) Finley_ElementFile_optimizeOrdering(&(in->Elements));
  if (Finley_noError()) Finley_ElementFile_optimizeOrdering(&(in->FaceElements));
  if (Finley_noError()) Finley_ElementFile_optimizeOrdering(&(in->Points));
  if (Finley_noError()) Finley_ElementFile_optimizeOrdering(&(in->ContactElements));
}
