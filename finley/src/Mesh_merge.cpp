
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "FinleyDomain.h"
#include "Util.h"

using escript::ValueError;

namespace finley {

FinleyDomain* FinleyDomain::merge(const std::vector<const FinleyDomain*>& msh)
{
    if (msh.empty()) {
        throw ValueError("merge: Empty mesh list");
    }
    for (int i = 0; i < msh.size(); i++) {
        if (msh[i]->getMPISize() > 1) {
            throw escript::NotImplementedError("merge: more than 1 processor is not supported yet.");
        }
    }

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
    int reducedOrder = msh[0]->reducedIntegrationOrder;
    const int numDim = msh[0]->getDim();
    escript::JMPI mpiInfo = msh[0]->getMPI();
    std::stringstream newName;

    for (int i=0; i < msh.size(); i++) {
        // check if all meshes have the same type and dimensions
        order = std::max(order, msh[i]->integrationOrder);
        reducedOrder = std::min(reducedOrder, msh[i]->reducedIntegrationOrder);
        numNodes += msh[i]->getNodes()->getNumNodes();
        if (mpiInfo->comm != msh[i]->getMPIComm()) {
            throw ValueError("merge: MPI communicators of meshes don't match.");
        }
        if (numDim != msh[i]->getDim()) {
            throw ValueError("merge: Spatial dimensions of meshes don't match.");
        }

        if (msh[i]->getElements()) {
            numElements += msh[i]->getElements()->numElements;
            if (elementTypeId == NoRef) {
                elementTypeId = msh[i]->getElements()->referenceElementSet->referenceElement->Type->TypeId;
            } else {
                if (elementTypeId != msh[i]->getElements()->referenceElementSet->referenceElement->Type->TypeId) {
                    throw ValueError("merge: element types of meshes don't match.");
                }
            }
        }

        if (msh[i]->getFaceElements()) {
            numFaceElements += msh[i]->getFaceElements()->numElements;
            if (faceElementTypeId == NoRef) {
                faceElementTypeId = msh[i]->getFaceElements()->referenceElementSet->referenceElement->Type->TypeId;
            } else {
                if (faceElementTypeId != msh[i]->getFaceElements()->referenceElementSet->referenceElement->Type->TypeId) {
                    throw ValueError("merge: face element types of meshes don't match.");
                }
            }
        }

        if (msh[i]->getContactElements()) {
            numContactElements += msh[i]->getContactElements()->numElements;
            if (contactTypeId == NoRef) {
                contactTypeId = msh[i]->getContactElements()->referenceElementSet->referenceElement->Type->TypeId;
            } else {
                if (contactTypeId != msh[i]->getContactElements()->referenceElementSet->referenceElement->Type->TypeId) {
                    throw ValueError("merge: contact element types of meshes don't match.");
                }
            }
        }

        if (msh[i]->getPoints()) {
            numPoints += msh[i]->getPoints()->numElements;
            if (pointTypeId == NoRef) {
                pointTypeId = msh[i]->getPoints()->referenceElementSet->referenceElement->Type->TypeId;
            } else {
                if (pointTypeId != msh[i]->getPoints()->referenceElementSet->referenceElement->Type->TypeId ) {
                    throw ValueError("merge: point element types of meshes don't match.");
                }
            }
        }

        if (i > 0)
            newName << "+";
        newName << msh[i]->m_name;
    }

    // allocate
    FinleyDomain* out = new FinleyDomain(newName.str(), numDim, mpiInfo);
    const_ReferenceElementSet_ptr refElements(new ReferenceElementSet(elementTypeId, order, reducedOrder));
    const_ReferenceElementSet_ptr refFaceElements(new ReferenceElementSet(faceElementTypeId, order, reducedOrder));
    const_ReferenceElementSet_ptr refContactElements(new ReferenceElementSet(contactTypeId, order, reducedOrder));
    const_ReferenceElementSet_ptr refPoints(new ReferenceElementSet(pointTypeId, order, reducedOrder));

    NodeFile* nodes = out->getNodes();
    out->setElements(new ElementFile(refElements, mpiInfo));
    out->setFaceElements(new ElementFile(refFaceElements, mpiInfo));
    out->setPoints(new ElementFile(refPoints, mpiInfo));
    out->setContactElements(new ElementFile(refContactElements, mpiInfo));

    // allocate new tables
    nodes->allocTable(numNodes);
    out->getElements()->allocTable(numElements);
    out->getFaceElements()->allocTable(numFaceElements);
    out->getContactElements()->allocTable(numContactElements);
    out->getPoints()->allocTable(numPoints);

    // copy tables
    numNodes = 0;
    numElements = 0;
    numFaceElements = 0;
    numContactElements = 0;
    numPoints = 0;

    for (int i = 0; i < msh.size(); i++) {
        nodes->copyTable(numNodes, maxNodeID, maxDOF, msh[i]->getNodes());
        out->getElements()->copyTable(numElements, numNodes, maxElementID, msh[i]->getElements());
        out->getFaceElements()->copyTable(numFaceElements, numNodes, maxElementID, msh[i]->getFaceElements());
        out->getContactElements()->copyTable(numContactElements, numNodes, maxElementID, msh[i]->getContactElements());
        out->getPoints()->copyTable(numPoints, numNodes, maxElementID, msh[i]->getPoints());

        numNodes += msh[i]->getNodes()->getNumNodes();
        numElements += msh[i]->getElements()->numElements;
        numFaceElements += msh[i]->getFaceElements()->numElements;
        numContactElements += msh[i]->getContactElements()->numElements;
        numPoints += msh[i]->getPoints()->numElements;

        if (msh[i]->getNodes()->getNumNodes() > 0)
            maxNodeID += util::getMaxInt(1, msh[i]->getNodes()->getNumNodes(), msh[i]->getNodes()->Id) + 1;
        maxDOF += util::getMaxInt(1, msh[i]->getNodes()->getNumNodes(), msh[i]->getNodes()->globalDegreesOfFreedom) + 1;
        maxElementID2 = 0;
        if (msh[i]->getElements()->numElements > 0)
            maxElementID2 = std::max(maxElementID2, util::getMaxInt(1, msh[i]->getElements()->numElements, msh[i]->getElements()->Id));
        if (msh[i]->getFaceElements()->numElements > 0)
            maxElementID2 = std::max(maxElementID2, util::getMaxInt(1, msh[i]->getFaceElements()->numElements, msh[i]->getFaceElements()->Id));
        if (msh[i]->getContactElements()->numElements > 0)
            maxElementID2 = std::max(maxElementID2, util::getMaxInt(1, msh[i]->getContactElements()->numElements, msh[i]->getContactElements()->Id));
        if (msh[i]->getPoints()->numElements > 0)
            maxElementID2 = std::max(maxElementID2, util::getMaxInt(1, msh[i]->getPoints()->numElements, msh[i]->getPoints()->Id));
        maxElementID += maxElementID2 + 1;
    }

    // all done
    out->prepare(false);
    return out;
}

} // namespace finley

