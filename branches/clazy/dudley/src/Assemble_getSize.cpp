
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

#include "Assemble.h"
#include "Util.h"

#include <escript/index.h>

namespace dudley {

void Assemble_getSize(const NodeFile* nodes, const ElementFile* elements,
                      escript::Data& out)
{
    if (!nodes || !elements)
        return;

    if (out.isComplex())
    {
        throw DudleyException("Assemble_getSize: complex arguments are not supported.");      
    }
    const int numDim = nodes->numDim;

    // now we look up what type of elements we need based on the function space
    // of out. If it is DUDLEY_REDUCED_ELEMENTS or
    // DUDLEY_REDUCED_FACE_ELEMENTS then we have single quad point
    int numQuad = (hasReducedIntegrationOrder(out) ? 1 : elements->numNodes);
    const int NN = elements->numNodes;
    const int NS = elements->numDim + 1;
    const int NVertices = elements->numDim + 1;

    // check the dimensions of out
    if (!out.numSamplesEqual(numQuad, elements->numElements)) {
        throw DudleyException("Assemble_getSize: illegal number of samples of element size Data object");
    } else if (!out.isDataPointShapeEqual(0, &numDim)) {
        throw DudleyException("Assemble_getSize: illegal data point shape of element size Data object");
    } else if (!out.actsExpanded()) {
        throw DudleyException("Assemble_getSize: expanded Data object is expected for element size.");
    }

    // now we can start
    out.requireWrite();
#pragma omp parallel
    {
        std::vector<double> local_X(NN * numDim);
#pragma omp for
        for (index_t e = 0; e < elements->numElements; e++) {
            // gather local coordinates of nodes into local_X(numDim,NN)
            util::gather(NS, &elements->Nodes[INDEX2(0, e, NN)], numDim,
                         nodes->Coordinates, &local_X[0]);
            // calculate minimal differences
            double max_diff = 0;
            for (int n0 = 0; n0 < NVertices; n0++) {
                for (int n1 = n0 + 1; n1 < NVertices; n1++) {
                    double diff = 0;
                    for (int i = 0; i < numDim; i++) {
                        const double d = local_X[INDEX2(i, n0, numDim)] - local_X[INDEX2(i, n1, numDim)];
                        diff += d * d;
                    }

                    max_diff = std::max(max_diff, diff);
                }
            }
            max_diff = sqrt(max_diff);
            // set all values to max_diff
            double* out_array = out.getSampleDataRW(e, static_cast<escript::DataTypes::real_t>(0));
            for (int q = 0; q < numQuad; q++)
                out_array[q] = max_diff;
        }
    } // end of parallel region
}

} // namespace dudley

