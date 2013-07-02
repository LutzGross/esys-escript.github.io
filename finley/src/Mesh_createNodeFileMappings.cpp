
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

/*   Finley: NodeFile : creates the mappings using by indexReducedNodes */
/*                 no distribution is happening                          */

/************************************************************************************/

#include "Mesh.h"

/************************************************************************************/

using namespace finley;

void Finley_Mesh_createMappings(Finley_Mesh* mesh, const int* dof_distribution, const int* node_distribution)
{
    std::vector<int> maskReducedNodes(mesh->Nodes->numNodes, -1);
    std::vector<int> indexReducedNodes(mesh->Nodes->numNodes);

    Finley_Mesh_markNodes(&maskReducedNodes[0], 0, mesh, true);
    int numReducedNodes=util::packMask(mesh->Nodes->numNodes, &maskReducedNodes[0], &indexReducedNodes[0]);
    if (Finley_noError())
        mesh->Nodes->createNodeMappings(numReducedNodes, indexReducedNodes, dof_distribution, node_distribution);
}


