
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

  Assemblage routines: calculates the minimum distance between two vertices
  of elements and assigns the value to each quadrature point in out.

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <escript/index.h>

namespace finley {

void Assemble_getSize(const NodeFile* nodes, const ElementFile* elements,
                      escript::Data& out)
{
    if (!nodes || !elements)
        return;

    const_ReferenceElement_ptr refElement(elements->referenceElementSet->
            borrowReferenceElement(
                util::hasReducedIntegrationOrder(out)));

    const int numDim = nodes->numDim;
    const int numQuad = refElement->Parametrization->numQuadNodes;
    const int NN = elements->numNodes;
    const int NS = refElement->Parametrization->Type->numShapes;
    const int NVertices = refElement->Parametrization->Type->numVertices;

    // check the dimensions of out
    if (!out.numSamplesEqual(numQuad, elements->numElements)) {
        throw escript::ValueError("Assemble_getSize: illegal number of samples of out Data object");
    } else if (!out.isDataPointShapeEqual(0, &numDim)) {
        throw escript::ValueError("Assemble_getSize: illegal data point shape of out Data object");
    } else if (!out.actsExpanded()) {
        throw escript::ValueError("Assemble_getSize: expanded Data object is expected for element size.");
    }

    // now we can start
    int node_offset;
    if (out.getFunctionSpace().getTypeCode() == FINLEY_CONTACT_ELEMENTS_2) {
        node_offset = refElement->Type->offsets[1];
    } else {
        node_offset = refElement->Type->offsets[0];
    }
    const double f = pow(0.5, pow((double)(refElement->Type->numSubElements),
                1./(double)(numDim))-1);

    out.requireWrite();
#pragma omp parallel
    {
        std::vector<double> local_X(NN * numDim);
#pragma omp for
        for (index_t e = 0; e < elements->numElements; e++) {
            // gather local coordinates of nodes into local_X(numDim,NN)
            util::gather(NS, &elements->Nodes[INDEX2(node_offset, e, NN)],
                         numDim, nodes->Coordinates, &local_X[0]);
            // calculate minimal differences:
            double max_diff = 0.;
            for (int n0 = 0; n0 < NVertices; n0++) {
                for (int n1 = n0 + 1; n1 < NVertices; n1++) {
                    double diff = 0;
                    for (int i = 0; i < numDim; i++) {
                        const double d = local_X[INDEX2(i,n0,numDim)] - local_X[INDEX2(i,n1,numDim)];
                        diff += d * d;
                    }
                    max_diff = std::max(max_diff, diff);
                }
            }
            max_diff = sqrt(max_diff) * f;
            // set all values to max_diff
            double* out_array = out.getSampleDataRW(e);
            for (int q = 0; q < numQuad; q++)
                out_array[q] = max_diff;
        }
    } // end of parallel region
}

} // namespace finley

