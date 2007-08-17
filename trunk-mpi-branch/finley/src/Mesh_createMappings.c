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

/*   Finley: Mesh: build up the mappings and distributions    */
/*                 no distribution is happening               */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id:$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Finley_Mesh_createMappings(Finley_Mesh* in, index_t *dof_first_component) {
  dim_t i, numReducedNodes;
  index_t *maskReducedNodes=NULL, *indexReducedNodes=NULL;

  maskReducedNodes=TMPMEMALLOC(in->Nodes->numNodes,index_t);
  indexReducedNodes=TMPMEMALLOC(in->Nodes->numNodes,index_t);
  if (! ( Finley_checkPtr(maskReducedNodes) ||  Finley_checkPtr(indexReducedNodes) ) ) {


    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<in->Nodes->numNodes;++i) maskReducedNodes[i]=-1;

    Finley_Mesh_markNodes(maskReducedNodes,0,in,TRUE);

    numReducedNodes=Finley_Util_packMask(in->Nodes->numNodes,maskReducedNodes,indexReducedNodes);

    Finley_Mesh_createNodeFileMappings(in,numReducedNodes,indexReducedNodes,dof_first_component);

  }

  TMPMEMFREE(maskReducedNodes);
  TMPMEMFREE(indexReducedNodes);
}
