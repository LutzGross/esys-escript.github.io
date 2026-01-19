
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
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

  detects faces in the mesh that match and replaces it by step elements

*****************************************************************************/

#include "FinleyDomain.h"

#include <escript/index.h>

namespace finley {

void FinleyDomain::joinFaces(double safety_factor, double tolerance, bool optimize)
{
    if (m_mpiInfo->size > 1) {
        throw escript::NotImplementedError("Mesh::joinFaces: MPI is not supported yet.");
    }
    if (!m_contactElements) {
        throw escript::ValueError("Mesh::joinFaces: No contact elements present.");
    }
    if (!m_faceElements)
        return;

    const_ReferenceElement_ptr faceRefElement(m_faceElements->referenceElementSet->borrowReferenceElement(false));
    const_ReferenceElement_ptr contactRefElement(m_contactElements->referenceElementSet->borrowReferenceElement(false));

    if (faceRefElement->Type->numNodesOnFace <= 0) {
        std::stringstream ss;
        ss << "Mesh::joinFaces: joining faces cannot be applied to face "
            "elements of type " << faceRefElement->Type->Name;
        throw escript::ValueError(ss.str());
    }

    if (contactRefElement->Type->numNodes != 2*faceRefElement->Type->numNodes) {
        std::stringstream ss;
        ss << "Mesh::joinFaces: contact element file for "
            << contactRefElement->Type->Name << " needs to hold elements "
            "created from face elements " << faceRefElement->Type->Name;
        throw escript::ValueError(ss.str());
    }

    const int NN = m_faceElements->numNodes;
    const int NN_Contact = m_contactElements->numNodes;

    // allocate work arrays
    int* elem1 = new int[m_faceElements->numElements];
    int* elem0 = new int[m_faceElements->numElements];
    index_t* elem_mask = new index_t[m_faceElements->numElements];
    int* matching_nodes_in_elem1 = new int[m_faceElements->numElements*NN];

    // find the matching face elements
    int numPairs;
    findMatchingFaces(safety_factor, tolerance, &numPairs, elem0, elem1, matching_nodes_in_elem1);
    // get a list of the face elements to be kept
#pragma omp parallel for
    for (index_t e = 0; e < m_faceElements->numElements; e++)
        elem_mask[e] = 1;

    for (int e = 0; e < numPairs; e++) {
        elem_mask[elem0[e]] = 0;
        elem_mask[elem1[e]] = 0;
    }
    dim_t new_numFaceElements = 0;
    // OMP
    for (index_t e = 0; e < m_faceElements->numElements; e++) {
        if (elem_mask[e] > 0) {
            elem_mask[new_numFaceElements] = e;
            new_numFaceElements++;
        }
    }
    // allocate new face element and Contact element files
    ElementFile* newContactElementsFile = new ElementFile(m_contactElements->referenceElementSet, m_mpiInfo);
    ElementFile* newFaceElementsFile = new ElementFile(m_faceElements->referenceElementSet, m_mpiInfo);
    newContactElementsFile->allocTable(numPairs+m_contactElements->numElements);
    newFaceElementsFile->allocTable(new_numFaceElements);
    // copy the old elements over
    // get the face elements which are still in use
    newFaceElementsFile->gather(elem_mask, m_faceElements);
    // get the contact elements which are still in use
    newContactElementsFile->copyTable(0, 0, 0, m_contactElements);
    dim_t c = m_contactElements->numElements;
    // OMP
    for (int e = 0; e < numPairs; e++) {
        const int e0 = elem0[e];
        const int e1 = elem1[e];
        newContactElementsFile->Id[c] = std::min(m_faceElements->Id[e0], m_faceElements->Id[e1]);
        newContactElementsFile->Tag[c] = std::min(m_faceElements->Tag[e0], m_faceElements->Tag[e1]);
        newContactElementsFile->Color[c] = e;
        for (int i = 0; i < NN; i++)
            newContactElementsFile->Nodes[INDEX2(i,c,NN_Contact)] =
                                        m_faceElements->Nodes[INDEX2(i,e0,NN)];
        for (int i = 0; i < NN; i++)
            newContactElementsFile->Nodes[INDEX2(i+NN,c,NN_Contact)] =
                                       matching_nodes_in_elem1[INDEX2(i,e,NN)];
        c++;
    }
    newContactElementsFile->minColor = 0;
    newContactElementsFile->maxColor = numPairs-1;
    // set new face and Contact elements
    
    delete m_faceElements;
    m_faceElements = newFaceElementsFile;

    delete m_contactElements;
    m_contactElements = newContactElementsFile;

    prepare(optimize);

    delete[] elem1;
    delete[] elem0;
    delete[] matching_nodes_in_elem1;
    delete[] elem_mask;
}

} // namespace finley

