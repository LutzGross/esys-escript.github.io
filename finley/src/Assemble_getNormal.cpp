
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


/****************************************************************************

  Assemblage routines: calculates the normal vector at quadrature points on
  face elements

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <escript/index.h>

namespace finley {

void Assemble_getNormal(const NodeFile* nodes, const ElementFile* elements,
                        escript::Data& normal)
{
    if (!nodes || !elements)
        return;

    const_ReferenceElement_ptr refElement(elements->referenceElementSet->
            borrowReferenceElement(util::hasReducedIntegrationOrder(normal)));
    const int NN = elements->numNodes;
    const int numDim = nodes->numDim;
    const int numQuad = refElement->Parametrization->numQuadNodes;
    const int numDim_local = refElement->Parametrization->Type->numDim;
    const int NS = refElement->Parametrization->Type->numShapes;
  
    int sign, node_offset;
    if (normal.getFunctionSpace().getTypeCode() == FINLEY_CONTACT_ELEMENTS_2) {
        node_offset = refElement->Type->offsets[1];
        sign = -1;
    } else {
        node_offset = refElement->Type->offsets[0];
        sign = 1;
    }

    // check the dimensions of normal
    if (!(numDim == numDim_local || numDim - 1 == numDim_local)) {
        throw escript::ValueError("Assemble_setNormal: Cannot calculate normal vector");
    } else if (!normal.numSamplesEqual(numQuad, elements->numElements)) {
        throw escript::ValueError("Assemble_setNormal: illegal number of samples of normal Data object");
    } else if (!normal.isDataPointShapeEqual(1, &numDim)) {
        throw escript::ValueError("Assemble_setNormal: illegal point data shape of normal Data object");
    } else if (!normal.actsExpanded()) {
        throw escript::ValueError("Assemble_setNormal: expanded Data object is expected for normal.");
    }

    normal.requireWrite();
#pragma omp parallel
    {
        std::vector<double> local_X(NS * numDim); 
        std::vector<double> dVdv(numQuad * numDim * numDim_local);
#pragma omp for
        for (index_t e = 0; e < elements->numElements; e++) {
            // gather local coordinates of nodes into local_X
            util::gather(NS, &elements->Nodes[INDEX2(node_offset, e, NN)],
                         numDim, nodes->Coordinates, &local_X[0]);
            // calculate dVdv(i,j,q)=local_X(i,n)*DSDv(n,j,q)
            util::smallMatMult(numDim, numDim_local * numQuad, &dVdv[0], NS,
                    local_X, refElement->Parametrization->dSdv);
            double* normal_array = normal.getSampleDataRW(e);
            util::normalVector(numQuad, numDim, numDim_local, &dVdv[0],
                               normal_array);
            for (int q = 0; q < numQuad * numDim; q++)
                normal_array[q] *= sign;
        }
    } // parallel section
}

} // namespace finley
