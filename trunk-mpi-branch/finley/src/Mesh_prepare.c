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

void Finley_Mesh_prepare(Finley_Mesh* in) {

     dim_t newGlobalNumDOFs=0;
     Paso_MPI_rank* mpiRankOfDOF=NULL;
     index_t* distribution=NULL;
     if (in==NULL) return;
     if (in->Nodes == NULL) return;

     /* first step is to distribute the elements according to a global distribution of DOF */

     mpiRankOfDOF=TMPMEMALLOC(in->Nodes->numNodes,Paso_MPI_rank);
     distribution=TMPMEMALLOC(in->MPIInfo->size+1,index_t);
     if (!(Finley_checkPtr(mpiRankOfDOF) || Finley_checkPtr(distribution))) {
        /* first we create dense labeling for the DOFs */
        newGlobalNumDOFs=Finley_NodeFile_createDenseDOFLabeling(in->Nodes);
        /* create a distribution of the global DOFs and determine
           the MPI_rank controling the DOFs on this processor      */
        Paso_MPIInfo_setDistribution(in->MPIInfo,0,newGlobalNumDOFs-1,distribution);
        Finley_NodeFile_assignMPIRankToDOFs(in->Nodes,mpiRankOfDOF,distribution);

        /* now the mesh is re-distributed according to the mpiRankOfDOF vector */
        /* this will redistribute the Nodes and Elements including overlap and will create an element coloring 
           but will not create any mappings (see later in this function)                                   */
        if (Finley_noError()) Finley_Mesh_distributeByRankOfDOF(in,mpiRankOfDOF);
     }
     TMPMEMFREE(mpiRankOfDOF);
     TMPMEMFREE(distribution);
     if (!Finley_noError()) return;

     /* at this stage we are able to start an optimization of the DOF distribution using ParaMetis */
     
     /* now a local labeling of the DOF is introduced */
return;

      /* set the labeling vectors in node files: */
       Finley_Mesh_prepareNodes(in);

       /* rearrange elements: */
       Finley_Mesh_optimizeElementOrdering(in);
}

bool_t Finley_Mesh_isPrepared(Finley_Mesh* in) {
     /* returns true if nodes and elements have been prepared for calculation */
     index_t out = TRUE;
     if ( in->Nodes != NULL) out = MIN(out , in->Nodes->isPrepared);
     if ( in->Elements != NULL) out = MIN(out ,  in->Elements->isPrepared);
     if ( in->FaceElements != NULL) out = MIN(out ,  in->FaceElements->isPrepared);
     if ( in->ContactElements != NULL) out = MIN(out ,  in->ContactElements->isPrepared);
     if ( in->Points != NULL) out = MIN(out ,  in->Points->isPrepared);

     return out == FINLEY_PREPARED;
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
