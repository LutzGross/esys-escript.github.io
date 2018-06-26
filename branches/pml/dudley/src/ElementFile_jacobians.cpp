
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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

#include "ElementFile.h"
#include "Assemble.h"
#include "ShapeTable.h"

namespace dudley {

ElementFile_Jacobians::ElementFile_Jacobians() :
    status(DUDLEY_INITIAL_STATUS - 1),
    numDim(0),
    numQuad(0),
    numElements(0),
    absD(NULL),
    quadweight(0),
    DSDX(NULL)
{
}

ElementFile_Jacobians::~ElementFile_Jacobians()
{
    delete[] DSDX;
    delete[] absD;
}

ElementFile_Jacobians* ElementFile::borrowJacobians(const NodeFile* nodes,
                                           bool reducedIntegrationOrder) const
{
    ElementFile_Jacobians* out =
                (reducedIntegrationOrder ? jacobians_reducedQ : jacobians);

    if (out->status < nodes->status) {
        out->numDim = nodes->numDim;
        out->numQuad = QuadNums[numDim][!reducedIntegrationOrder];
        out->numShapes = numDim + 1;
        out->numElements = numElements;
        if (out->DSDX == NULL)
            out->DSDX = new double[out->numElements * out->numShapes * out->numDim * out->numQuad];
        if (out->absD == NULL)
            out->absD = new double[out->numElements];

        /*========================== dim = 2 ============================= */
        if (out->numDim == 2) {
            if (numLocalDim == 1) {
                Assemble_jacobians_2D_M1D_E1D(nodes->Coordinates, out->numQuad,
                        numElements, numNodes, Nodes, out->DSDX, out->absD,
                        &out->quadweight, Id);
            } else if (numLocalDim == 2) {
                Assemble_jacobians_2D(nodes->Coordinates, out->numQuad,
                        numElements, numNodes, Nodes, out->DSDX, out->absD,
                        &out->quadweight, Id);
            } else {
                throw DudleyException("ElementFile::borrowJacobians: local "
                                "dimension in a 2D domain has to be 1 or 2.");
            }
        /*========================== dim = 3 ============================= */
        } else if (out->numDim == 3) {
            if (numLocalDim == 2) {
                Assemble_jacobians_3D_M2D_E2D(nodes->Coordinates, out->numQuad,
                        numElements, numNodes, Nodes, out->DSDX, out->absD,
                        &out->quadweight, Id);
            } else if (numLocalDim == 3) {
                Assemble_jacobians_3D(nodes->Coordinates, out->numQuad,
                        numElements, numNodes, Nodes, out->DSDX, out->absD,
                        &out->quadweight, Id);
            } else {
                throw DudleyException("ElementFile::borrowJacobians: local "
                                "dimension in a 3D domain has to be 2 or 3.");
            }
        } else {
            throw DudleyException("ElementFile::borrowJacobians: number of "
                                  "spatial dimensions has to be 2 or 3.");
        }
        out->status = nodes->status;
    }
    return out;
}

} // namespace dudley

