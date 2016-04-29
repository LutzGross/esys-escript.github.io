
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


/****************************************************************************

  Finley: Mesh

  detects faces in the mesh that match and replaces it by step elements

*****************************************************************************/

#include "Mesh.h"

#include <escript/index.h>

namespace finley {

void Mesh::joinFaces(double safety_factor, double tolerance, bool optimize)
{
    if (MPIInfo->size > 1) {
        throw escript::NotImplementedError("Mesh::joinFaces: MPI is not supported yet.");
    }
    if (!ContactElements) {
        throw escript::ValueError("Mesh::joinFaces: No contact elements present.");
    }
    if (!FaceElements)
        return;

    const_ReferenceElement_ptr faceRefElement(FaceElements->referenceElementSet->borrowReferenceElement(false));
    const_ReferenceElement_ptr contactRefElement(ContactElements->referenceElementSet->borrowReferenceElement(false));

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

    const int NN = FaceElements->numNodes;
    const int NN_Contact = ContactElements->numNodes;

    // allocate work arrays
    int* elem1 = new int[FaceElements->numElements];
    int* elem0 = new int[FaceElements->numElements];
    index_t* elem_mask = new index_t[FaceElements->numElements];
    int* matching_nodes_in_elem1 = new int[FaceElements->numElements*NN];

    // find the matching face elements
    int numPairs;
    findMatchingFaces(safety_factor, tolerance, &numPairs, elem0, elem1, matching_nodes_in_elem1);
    // get a list of the face elements to be kept
#pragma omp parallel for
    for (index_t e = 0; e < FaceElements->numElements; e++)
        elem_mask[e] = 1;

    for (int e = 0; e < numPairs; e++) {
        elem_mask[elem0[e]] = 0;
        elem_mask[elem1[e]] = 0;
    }
    dim_t new_numFaceElements = 0;
    // OMP
    for (index_t e = 0; e < FaceElements->numElements; e++) {
        if (elem_mask[e] > 0) {
            elem_mask[new_numFaceElements] = e;
            new_numFaceElements++;
        }
    }
    // allocate new face element and Contact element files
    ElementFile *newFaceElementsFile, *newContactElementsFile;
    newContactElementsFile = new ElementFile(ContactElements->referenceElementSet, MPIInfo);
    newFaceElementsFile = new ElementFile(FaceElements->referenceElementSet, MPIInfo);
    newContactElementsFile->allocTable(numPairs+ContactElements->numElements);
    newFaceElementsFile->allocTable(new_numFaceElements);
    // copy the old elements over
    // get the face elements which are still in use
    newFaceElementsFile->gather(elem_mask, FaceElements);
    // get the contact elements which are still in use
    newContactElementsFile->copyTable(0, 0, 0, ContactElements);
    dim_t c = ContactElements->numElements;
    // OMP
    for (int e = 0; e < numPairs; e++) {
        const int e0 = elem0[e];
        const int e1 = elem1[e];
        newContactElementsFile->Id[c] = std::min(FaceElements->Id[e0],FaceElements->Id[e1]);
        newContactElementsFile->Tag[c] = std::min(FaceElements->Tag[e0],FaceElements->Tag[e1]);
        newContactElementsFile->Color[c] = e;
        for (int i = 0; i < NN; i++)
            newContactElementsFile->Nodes[INDEX2(i,c,NN_Contact)]=FaceElements->Nodes[INDEX2(i,e0,NN)];
        for (int i = 0; i < NN; i++)
            newContactElementsFile->Nodes[INDEX2(i+NN,c,NN_Contact)]=matching_nodes_in_elem1[INDEX2(i,e,NN)];
        c++;
    }
    newContactElementsFile->minColor = 0;
    newContactElementsFile->maxColor = numPairs-1;
    // set new face and Contact elements
    delete FaceElements;
    FaceElements=newFaceElementsFile;
    delete ContactElements;
    ContactElements=newContactElementsFile;
    prepare(optimize);
    delete[] elem1;
    delete[] elem0;
    delete[] matching_nodes_in_elem1;
    delete[] elem_mask;
}

} // namespace finley

