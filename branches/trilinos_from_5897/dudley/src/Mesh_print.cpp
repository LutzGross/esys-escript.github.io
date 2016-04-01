
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include "Mesh.h"

namespace dudley {

/// prints the mesh details to standard output
void Mesh::print()
{
    // write header
    std::cout << "Mesh name: " << m_name << std::endl;

    // write nodes
    if (Nodes) {
        const index_t* dofs = Nodes->borrowTargetDegreesOfFreedom();
        const index_t* nodes = Nodes->borrowTargetNodes();
        const int numDim = Nodes->numDim;
        printf("=== %1dD-Nodes:\nnumber of nodes=%d\n", numDim, Nodes->getNumNodes());
        std::cout << "Id,Tag,globalDegreesOfFreedom,degreesOfFreedom,reducedDegreesOfFeedom,node,reducedNode,Coordinates\n";
        for (index_t i = 0; i < Nodes->getNumNodes(); i++) {
            printf("%d,%d,%d,%d,%d,%d,%d ",
                   Nodes->Id[i], Nodes->Tag[i], Nodes->globalDegreesOfFreedom[i],
                   dofs[i], dofs[i], nodes[i], nodes[i]);
            for (int j = 0; j < numDim; j++)
                printf(" %20.15e", Nodes->Coordinates[INDEX2(j, i, numDim)]);
            std::cout << std::endl;
        }
    }

    // write elements
    if (Elements) {
        std::cout << "=== "
                 << Elements->ename << ":\nnumber of elements="
                 << Elements->numElements << "\ncolor range=["
                 << Elements->minColor << "," << Elements->maxColor << "]\n";

        if (Elements->numElements > 0) {
            const int NN = Elements->numNodes;
            std::cout << "Id,Tag,Owner,Color,Nodes" << std::endl;
            for (index_t i = 0; i < Elements->numElements; i++) {
                std::cout << Elements->Id[i] << "," << Elements->Tag[i] << ","
                    << Elements->Owner[i] << "," << Elements->Color[i] << ",";
                for (int j = 0; j < NN; j++)
                    std::cout << " " << Nodes->Id[Elements->Nodes[INDEX2(j, i, NN)]];
                std::cout << std::endl;
            }
        }
    }

    // write face elements
    if (FaceElements) {
        std::cout << "=== "
                 << FaceElements->ename << ":\nnumber of elements="
                 << FaceElements->numElements << "\ncolor range=["
                 << FaceElements->minColor << "," << FaceElements->maxColor
                 << "]\n";

        if (FaceElements->numElements > 0) {
            const int NN = FaceElements->numNodes;
            std::cout << "Id,Tag,Owner,Color,Nodes" << std::endl;
            for (index_t i = 0; i < FaceElements->numElements; i++) {
                std::cout << FaceElements->Id[i] << "," << FaceElements->Tag[i]
                    << "," << FaceElements->Owner[i] << ","
                    << FaceElements->Color[i] << ",";
                for (int j=0; j<NN; j++)
                    std::cout << " " << Nodes->Id[FaceElements->Nodes[INDEX2(j,i,NN)]];
                std::cout << std::endl;
            }
        }
    }

    // write points
    if (Points) {
        std::cout << "=== "
                 << Points->ename << ":\nnumber of elements="
                 << Points->numElements << "\ncolor range=["
                 << Points->minColor << "," << Points->maxColor << "]\n";

        if (Points->numElements > 0) {
            const int NN = Points->numNodes;
            std::cout << "Id,Tag,Owner,Color,Nodes" << std::endl;
            for (index_t i=0; i<Points->numElements; i++) {
                std::cout << Points->Id[i] << "," << Points->Tag[i] << ","
                    << Points->Owner[i] << "," << Points->Color[i] << ",";
                for (int j=0; j<NN; j++)
                    std::cout << " " << Nodes->Id[Points->Nodes[INDEX2(j,i,NN)]];
                std::cout << std::endl;
            }
        }
    }
}

} // namespace dudley

