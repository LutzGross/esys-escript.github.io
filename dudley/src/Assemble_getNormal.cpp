
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
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
#include "ShapeTable.h"
#include "Util.h"

#include <escript/index.h>

namespace dudley {

void Assemble_getNormal(const NodeFile* nodes, const ElementFile* elements,
                        escript::Data& normal)
{
    if (!nodes || !elements)
        return;

    if (normal.isComplex())
    {
        throw DudleyException("Assemble_setNormal: complex arguments not supported.");
    }
    const int NN = elements->numNodes;
    const int numDim = nodes->numDim;
    const int numQuad = (hasReducedIntegrationOrder(normal) ? 1 : NN);
    const int numDim_local = elements->numLocalDim;
    const int NS = elements->numDim + 1;

    const double *dSdv = NULL;
    switch (elements->numDim) {
        case 2:
            dSdv = &DTDV_2D[0][0];
        break;
        case 3:
            dSdv = &DTDV_3D[0][0];
        break;
        default:
            dSdv = &DTDV_1D[0][0];
        break;
    }

    // check the dimensions of normal
    if (!(numDim == numDim_local || numDim - 1 == numDim_local)) {
        throw DudleyException("Assemble_setNormal: Cannot calculate normal vector");
    } else if (!normal.isDataPointShapeEqual(1, &numDim)) {
        throw DudleyException("Assemble_setNormal: illegal point data shape of normal Data object");
    } else if (!normal.numSamplesEqual(numQuad, elements->numElements)) {
        throw DudleyException("Assemble_setNormal: illegal number of samples of normal Data object");
    } else if (!normal.actsExpanded()) {
        throw DudleyException("Assemble_setNormal: expanded Data object is expected for normal.");
    }

    normal.requireWrite();
#pragma omp parallel
    {
        std::vector<double> local_X(NS * numDim);
        std::vector<double> dVdv(numQuad * numDim * numDim_local);
#pragma omp for
        for (index_t e = 0; e < elements->numElements; e++) {
            // gather local coordinates of nodes into local_X
            util::gather(NS, &elements->Nodes[INDEX2(0, e, NN)], numDim,
                         nodes->Coordinates, &local_X[0]);

            // calculate dVdv(i,j,q)=local_X(i,n)*DSDv(n,j,q)
            util::smallMatMult(numDim, numDim_local * numQuad,
                                     &dVdv[0], NS, &local_X[0], dSdv);
            // get normalized vector
            double* normal_array = normal.getSampleDataRW(e, static_cast<escript::DataTypes::real_t>(0));
            util::normalVector(numQuad, numDim, numDim_local, &dVdv[0], normal_array);
        }
    }
}

} // namespace dudley

