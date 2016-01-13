
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/************************************************************************************/

/*   Finley: Mesh: prepares the mesh for further calculations  */

/************************************************************************************/

#include "Mesh.h"

using namespace finley;

void Finley_Mesh_prepare(Finley_Mesh* in, bool_t optimize) {
     dim_t newGlobalNumDOFs=0, numReducedNodes=0,i;
     index_t* distribution=NULL, *maskReducedNodes=NULL, *indexReducedNodes=NULL, *node_distribution=NULL;
     if (in==NULL) return;
     if (in->Nodes == NULL) return;

     Finley_Mesh_setOrders(in);

     /* first step is to distribute the elements according to a global distribution of DOF */

     distribution=new index_t[in->MPIInfo->size+1];
     node_distribution=new index_t[in->MPIInfo->size+1];
     if (! (Finley_checkPtr(distribution) || Finley_checkPtr(node_distribution))) {
        /* first we create dense labeling for the DOFs */

        newGlobalNumDOFs=in->Nodes->createDenseDOFLabeling();

        /* create a distribution of the global DOFs and determine
           the MPI_rank controlling the DOFs on this processor      */
        Esys_MPIInfo_setDistribution(in->MPIInfo,0,newGlobalNumDOFs-1,distribution);

        /* now the mesh is re-distributed according to the mpiRankOfDOF vector */
        /* this will redistribute the Nodes and Elements including overlap and will create an element coloring 
           but will not create any mappings (see later in this function) */
        if (Finley_noError()) Finley_Mesh_distributeByRankOfDOF(in,distribution);
     }

     /* at this stage we are able to start an optimization of the DOF distribution using ParMetis */
     /* on return distribution is altered and new DOF ids have been assigned */
     if (Finley_noError() && optimize && in->MPIInfo->size>1) {
         Finley_Mesh_optimizeDOFDistribution(in,distribution); 
         if (Finley_noError()) Finley_Mesh_distributeByRankOfDOF(in,distribution); 
     }
     /* the local labelling of the degrees of freedom is optimized */
     if (Finley_noError() && optimize) {
       Finley_Mesh_optimizeDOFLabeling(in,distribution); 
     }
     /* rearrange elements with the attempt to bring elements closer to memory locations of the nodes (distributed shared memory!): */
     if (Finley_noError()) Finley_Mesh_optimizeElementOrdering(in);


     /* create the global indices */
     if (Finley_noError()) {


        maskReducedNodes=new index_t[in->Nodes->numNodes];
        indexReducedNodes=new index_t[in->Nodes->numNodes];
        if (! ( Finley_checkPtr(maskReducedNodes) ||  Finley_checkPtr(indexReducedNodes) ) ) {

/* useful DEBUG:
{index_t MIN_id,MAX_id;
printf("Mesh_prepare: global DOF : %d\n",newGlobalNumDOFs);
in->Nodes->setGlobalIdRange(&MIN_id,&MAX_id);
printf("Mesh_prepare: global node id range = %d :%d\n", MIN_id,MAX_id);
in->Nodes->setIdRange(&MIN_id,&MAX_id);
printf("Mesh_prepare: local node id range = %d :%d\n", MIN_id,MAX_id);
}
*/
          #pragma omp parallel for private(i) schedule(static)
          for (i=0;i<in->Nodes->numNodes;++i) maskReducedNodes[i]=-1;

          Finley_Mesh_markNodes(maskReducedNodes,0,in,TRUE);
   
          numReducedNodes=util::packMask(in->Nodes->numNodes,maskReducedNodes,indexReducedNodes);

          in->Nodes->createDenseNodeLabeling(node_distribution, distribution); 
          // created reduced DOF labeling
          in->Nodes->createDenseReducedLabeling(maskReducedNodes, false); 
          // created reduced node labeling
          in->Nodes->createDenseReducedLabeling(maskReducedNodes, true);

          /* create the missing mappings */
          if (Finley_noError()) Finley_Mesh_createNodeFileMappings(in,numReducedNodes,indexReducedNodes,distribution, node_distribution);
        }

        delete[] maskReducedNodes;
        delete[] indexReducedNodes;
     }

     delete[] distribution;
     delete[] node_distribution;

     Finley_Mesh_setTagsInUse(in);
}

/*                                                      */
/*  tries to reduce the coloring for all element files: */
/*                                                      */
void Finley_Mesh_createColoring(Finley_Mesh* in, index_t *node_localDOF_map) {
  if (Finley_noError()) in->Elements->createColoring(in->Nodes->numNodes,node_localDOF_map);
  if (Finley_noError()) in->FaceElements->createColoring(in->Nodes->numNodes,node_localDOF_map);
  if (Finley_noError()) in->Points->createColoring(in->Nodes->numNodes,node_localDOF_map);
  if (Finley_noError()) in->ContactElements->createColoring(in->Nodes->numNodes,node_localDOF_map);
}
/*                                                                    */
/*  redistribute elements to minimize communication during assemblage */
/*                                                                    */
void Finley_Mesh_optimizeElementOrdering(Finley_Mesh* in) {
  if (Finley_noError()) in->Elements->optimizeOrdering();
  if (Finley_noError()) in->FaceElements->optimizeOrdering();
  if (Finley_noError()) in->Points->optimizeOrdering();
  if (Finley_noError()) in->ContactElements->optimizeOrdering();
}

void Finley_Mesh_setTagsInUse(Finley_Mesh* in)
{
    if (Finley_noError()) in->Nodes->updateTagList();
    if (Finley_noError()) in->Elements->updateTagList();
    if (Finley_noError()) in->FaceElements->updateTagList();
    if (Finley_noError()) in->Points->updateTagList();
    if (Finley_noError()) in->ContactElements->updateTagList();
}

