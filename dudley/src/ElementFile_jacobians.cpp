
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

#include "ElementFile.h"
#include "Assemble.h"
#include "ShapeTable.h"

namespace dudley {

Dudley_ElementFile_Jacobians *Dudley_ElementFile_Jacobians_alloc(void)
{
    Dudley_ElementFile_Jacobians *out = new Dudley_ElementFile_Jacobians;
    out->status = DUDLEY_INITIAL_STATUS - 1;
    out->numDim = 0;
    out->numQuad = 0;
    out->numElements = 0;
    out->absD = NULL;
    out->quadweight = 0;
    out->DSDX = NULL;
    return out;
}

void Dudley_ElementFile_Jacobians_dealloc(Dudley_ElementFile_Jacobians * in)
{
    if (in != NULL)
    {
        delete[] in->DSDX;
        delete[] in->absD;
        delete in;
    }
}

Dudley_ElementFile_Jacobians *Dudley_ElementFile_borrowJacobians(
                 const Dudley_ElementFile* self, const Dudley_NodeFile* nodes,
                 bool reducedIntegrationOrder)
{
    Dudley_ElementFile_Jacobians *out = NULL;

    dim_t numNodes = self->numNodes;

    if (reducedIntegrationOrder)
    {
        out = self->jacobians_reducedQ;
    }
    else
    {
        out = self->jacobians;
    }
    if (out->status < nodes->status)
    {
        out->numDim = nodes->numDim;
        out->numQuad = QuadNums[self->numDim][!reducedIntegrationOrder];
        out->numShapes = self->numDim + 1;
        out->numElements = self->numElements;
        if (out->DSDX == NULL)
            out->DSDX = new  double[(out->numElements) * (out->numShapes) * (out->numDim) * (out->numQuad)];
        if (out->absD == NULL)
            out->absD = new  double[out->numElements];

        /*========================== dim = 1 ============================= */
        if (out->numDim == 1) {
            throw DudleyException("Dudley does not support 1D domains.");
        /*========================== dim = 2 ============================= */
        } else if (out->numDim == 2) {
            if (self->numLocalDim == 1) {
                dudley::Assemble_jacobians_2D_M1D_E1D(nodes->Coordinates, out->numQuad, self->numElements, numNodes,
                                              self->Nodes, out->DSDX, out->absD, &(out->quadweight), self->Id);
            } else if (self->numLocalDim == 2) {
                dudley::Assemble_jacobians_2D(nodes->Coordinates, out->numQuad, self->numElements, numNodes, self->Nodes,
                                      out->DSDX, out->absD, &(out->quadweight), self->Id);
            } else {
                throw DudleyException("Dudley_ElementFile_borrowJacobians: local dimension in a 2D domain has to be 1 or 2.");
            }
        /*========================== dim = 3 ============================= */
        } else if (out->numDim == 3) {
            if (self->numLocalDim == 2) {
                dudley::Assemble_jacobians_3D_M2D_E2D(nodes->Coordinates, out->numQuad, self->numElements, numNodes,
                                              self->Nodes, out->DSDX, out->absD, &(out->quadweight), self->Id);
            } else if (self->numLocalDim == 3) {
                dudley::Assemble_jacobians_3D(nodes->Coordinates, out->numQuad, self->numElements, numNodes, self->Nodes,
                                      out->DSDX, out->absD, &(out->quadweight), self->Id);
            } else {
                throw DudleyException("Dudley_ElementFile_borrowJacobians: local dimension in a 3D domain has to be 2 or 3.");
            }
        } else {
            throw DudleyException("Dudley_ElementFile_borrowJacobians: spatial dimension has to be 1, 2 or 3.");
        }
        out->status = nodes->status;
    }
    return out;
}

} // namespace dudley

