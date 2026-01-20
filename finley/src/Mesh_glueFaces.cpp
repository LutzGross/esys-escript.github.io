
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


/****************************************************************************

  Finley: Mesh

  removes matching face elements from the mesh

*****************************************************************************/

#include "FinleyDomain.h"

#include <escript/index.h>

namespace finley {

void FinleyDomain::glueFaces(double safetyFactor, double tolerance, bool optimize)
{
    if (m_mpiInfo->size > 1) {
        throw escript::NotImplementedError("glueFaces: MPI is not supported yet.");
    }
    if (!m_faceElements)
        return;

    const_ReferenceElement_ptr faceRefElement(m_faceElements->
                        referenceElementSet->borrowReferenceElement(false));
    const int NNFace = faceRefElement->Type->numNodesOnFace;
    const int NN = m_faceElements->numNodes;
    const int numDim = m_nodes->numDim;
    const int* faceNodes = faceRefElement->Type->faceNodes;
   
    if (NNFace <= 0) {
        std::stringstream ss;
        ss << "Mesh::glueFaces: glueing faces cannot be applied to face "
            "elements of type " << faceRefElement->Type->Name;
        throw escript::ValueError(ss.str());
    }

    // allocate work arrays
    int* elem1 = new int[m_faceElements->numElements];
    int* elem0 = new int[m_faceElements->numElements];
    IndexVector elem_mask(m_faceElements->numElements, 0);
    int* matching_nodes_in_elem1 = new int[m_faceElements->numElements*NN];
    IndexVector new_node_label(m_nodes->getNumNodes());
    // find the matching face elements
    int numPairs;
    findMatchingFaces(safetyFactor, tolerance, &numPairs, elem0, elem1,
                      matching_nodes_in_elem1);
    for (index_t n = 0; n < m_nodes->getNumNodes(); n++)
        new_node_label[n] = n;
    // mark matching face elements to be removed
    for (int e = 0; e < numPairs; e++) {
        elem_mask[elem0[e]] = 1;
        elem_mask[elem1[e]] = 1;
        for (int i = 0; i < NNFace; i++) {
            const int face_node = faceNodes[i];
            new_node_label[matching_nodes_in_elem1[INDEX2(face_node,e,NN)]] =
                    m_faceElements->Nodes[INDEX2(face_node,elem0[e],NN)];
        }
    }
    // create an index of face elements
    dim_t new_numFaceElements = 0;
    for (index_t e = 0; e < m_faceElements->numElements; e++) {
        if (elem_mask[e] < 1) {
            elem_mask[new_numFaceElements] = e;
            new_numFaceElements++;
        }
    }
    // get the new number of nodes
    IndexVector new_node_mask(m_nodes->getNumNodes(), -1);
    IndexVector new_node_list;
    dim_t newNumNodes = 0;
    for (index_t n = 0; n < m_nodes->getNumNodes(); n++)
        new_node_mask[new_node_label[n]] = 1;
    for (index_t n = 0; n < m_nodes->getNumNodes(); n++) {
        if (new_node_mask[n] > 0) {
            new_node_mask[n] = newNumNodes;
            new_node_list.push_back(n);
            newNumNodes++;
        }
    }

    for (index_t n = 0; n < m_nodes->getNumNodes(); n++)
        new_node_label[n] = new_node_mask[new_node_label[n]];


    // allocate new node and element files
    NodeFile* newNodeFile = new NodeFile(numDim, m_mpiInfo); 
    newNodeFile->allocTable(newNumNodes);
    
    ElementFile* newFaceElementsFile = new ElementFile(
            m_faceElements->referenceElementSet, m_mpiInfo);
    newFaceElementsFile->allocTable(new_numFaceElements);
    // get the new nodes
    newNodeFile->gather(&new_node_list[0], m_nodes);

    // they are the new nodes
    delete m_nodes;
    m_nodes = newNodeFile;
    // get the face elements which are still in use
    newFaceElementsFile->gather(&elem_mask[0], m_faceElements);
    
    
    // they are the new face elements
    delete m_faceElements;
    m_faceElements = newFaceElementsFile;

    // assign new node ids to elements
    


    relabelElementNodes(new_node_label, 0);
    prepare(optimize);
    delete[] elem1;
    delete[] elem0;
    delete[] matching_nodes_in_elem1;
}

} // namespace finley

