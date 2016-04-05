
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

  takes nodes, elements, etc. of all in put meshes and copies them into
  a new mesh. Ids of output are shifted by the maximum Id of input.

*****************************************************************************/

#include "Mesh.h"
#include "Util.h"

namespace finley {

Mesh* Mesh_merge(const std::vector<Mesh*>& msh)
{
    if (msh.size()==0) {
        throw escript::ValueError("Mesh_merge: Empty mesh list");
    }
    for (int i=0; i<msh.size(); i++) {
        if (msh[i]->MPIInfo->size > 1) {
            throw escript::NotImplementedError("Mesh_merge: more than 1 processor is not supported yet.");
        }
    }

    const_ReferenceElementSet_ptr refPoints, refContactElements, refFaceElements, refElements;
    dim_t numNodes = 0;
    dim_t numElements = 0;
    dim_t numFaceElements = 0;
    dim_t numContactElements = 0;
    dim_t numPoints = 0;
    index_t maxNodeID = 0;
    index_t maxDOF = 0;
    index_t maxElementID = 0;
    index_t maxElementID2 = 0;
    ElementTypeId elementTypeId = NoRef;
    ElementTypeId faceElementTypeId = NoRef;
    ElementTypeId pointTypeId = NoRef;
    ElementTypeId contactTypeId = NoRef;

    int order = msh[0]->integrationOrder;
    int reduced_order = msh[0]->reducedIntegrationOrder;
    const int numDim = msh[0]->Nodes->numDim;
    escript::JMPI& mpiInfo = msh[0]->MPIInfo;
    std::stringstream newName;

    for (int i=0; i < msh.size(); i++) {
        // check if all meshes have the same type and dimensions
        order=std::max(order, msh[i]->integrationOrder);
        reduced_order=std::min(reduced_order, msh[i]->reducedIntegrationOrder);
        numNodes+=msh[i]->Nodes->numNodes;
        if (mpiInfo->comm != msh[i]->MPIInfo->comm) {
            throw escript::ValueError("Mesh_merge: MPI communicators of meshes don't match.");
        }
        if (numDim != msh[i]->Nodes->numDim) {
            throw escript::ValueError("Mesh_merge: Spatial dimensions of meshes don't match.");
        }

        if (msh[i]->Elements) {
            numElements+=msh[i]->Elements->numElements;
            if (elementTypeId==NoRef) {
                elementTypeId=msh[i]->Elements->referenceElementSet->referenceElement->Type->TypeId;
            } else {
                if (elementTypeId != msh[i]->Elements->referenceElementSet->referenceElement->Type->TypeId) {
                    throw escript::ValueError("Mesh_merge: element types of meshes don't match.");
                }
            }
        }

        if (msh[i]->FaceElements) {
            numFaceElements+=msh[i]->FaceElements->numElements;
            if (faceElementTypeId==NoRef) {
                faceElementTypeId=msh[i]->FaceElements->referenceElementSet->referenceElement->Type->TypeId;
            } else {
                if (faceElementTypeId != msh[i]->FaceElements->referenceElementSet->referenceElement->Type->TypeId) {
                    throw escript::ValueError("Mesh_merge: face element types of meshes don't match.");
                }
            }
        }

        if (msh[i]->ContactElements) {
            numContactElements+=msh[i]->ContactElements->numElements;
            if (contactTypeId==NoRef) {
                contactTypeId=msh[i]->ContactElements->referenceElementSet->referenceElement->Type->TypeId;
            } else {
                if (contactTypeId != msh[i]->ContactElements->referenceElementSet->referenceElement->Type->TypeId) {
                    throw escript::ValueError("Mesh_merge: contact element types of meshes don't match.");
                }
            }
        }

        if (msh[i]->Points) {
            numPoints+=msh[i]->Points->numElements;
            if (pointTypeId==NoRef) {
                pointTypeId=msh[i]->Points->referenceElementSet->referenceElement->Type->TypeId;
            } else {
                if (pointTypeId != msh[i]->Points->referenceElementSet->referenceElement->Type->TypeId ) {
                    throw escript::ValueError("Mesh_merge: point element types of meshes don't match.");
                }
            }
        }

        if (i>0)
            newName << "+";
        newName << msh[i]->m_name;
    }

    // allocate
    Mesh* out = new Mesh(newName.str(), numDim, mpiInfo);
    refElements.reset(new ReferenceElementSet(elementTypeId, order, reduced_order));
    refFaceElements.reset(new ReferenceElementSet(faceElementTypeId, order, reduced_order));
    refContactElements.reset(new ReferenceElementSet(contactTypeId, order, reduced_order));
    refPoints.reset(new ReferenceElementSet(pointTypeId, order, reduced_order));
    out->Elements = new ElementFile(refElements, mpiInfo);
    out->FaceElements = new ElementFile(refFaceElements, mpiInfo);
    out->Points = new ElementFile(refPoints, mpiInfo);
    out->ContactElements = new ElementFile(refContactElements, mpiInfo);

    // allocate new tables
    out->Nodes->allocTable(numNodes);
    out->Elements->allocTable(numElements);
    out->FaceElements->allocTable(numFaceElements);
    out->ContactElements->allocTable(numContactElements);
    out->Points->allocTable(numPoints);

    // copy tables
    numNodes = 0;
    numElements = 0;
    numFaceElements = 0;
    numContactElements = 0;
    numPoints = 0;

    for (int i=0; i<msh.size(); i++) {
        out->Nodes->copyTable(numNodes, maxNodeID, maxDOF, msh[i]->Nodes);
        out->Elements->copyTable(numElements, numNodes, maxElementID, msh[i]->Elements);
        out->FaceElements->copyTable(numFaceElements, numNodes, maxElementID, msh[i]->FaceElements);
        out->ContactElements->copyTable(numContactElements, numNodes, maxElementID,msh[i]->ContactElements);
        out->Points->copyTable(numPoints, numNodes, maxElementID, msh[i]->Points);

        numNodes += msh[i]->Nodes->numNodes;
        numElements += msh[i]->Elements->numElements;
        numFaceElements += msh[i]->FaceElements->numElements;
        numContactElements += msh[i]->ContactElements->numElements;
        numPoints += msh[i]->Points->numElements;

        if (msh[i]->Nodes->numNodes > 0)
            maxNodeID += util::getMaxInt(1, msh[i]->Nodes->numNodes, msh[i]->Nodes->Id)+1;
        maxDOF += util::getMaxInt(1, msh[i]->Nodes->numNodes, msh[i]->Nodes->globalDegreesOfFreedom)+1;
        maxElementID2 = 0;
        if (msh[i]->Elements->numElements > 0)
            maxElementID2 = std::max(maxElementID2, util::getMaxInt(1, msh[i]->Elements->numElements, msh[i]->Elements->Id));
        if (msh[i]->FaceElements->numElements>0)
            maxElementID2 = std::max(maxElementID2, util::getMaxInt(1, msh[i]->FaceElements->numElements, msh[i]->FaceElements->Id));
        if (msh[i]->ContactElements->numElements>0)
            maxElementID2 = std::max(maxElementID2, util::getMaxInt(1, msh[i]->ContactElements->numElements, msh[i]->ContactElements->Id));
        if (msh[i]->Points->numElements)
            maxElementID2 = std::max(maxElementID2, util::getMaxInt(1, msh[i]->Points->numElements, msh[i]->Points->Id));
        maxElementID += maxElementID2+1;
    }

    // all done
    out->prepare(false);
    return out;
}

} // namespace finley

