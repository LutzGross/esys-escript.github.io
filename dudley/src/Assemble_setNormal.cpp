
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

/****************************************************************************

  Assemblage routines: calculates the normal vector at quadrature points on
  face elements

*****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Assemble.h"
#include "ShapeTable.h"
#include "Util.h"

namespace dudley {

void Assemble_setNormal(Dudley_NodeFile* nodes, Dudley_ElementFile* elements, escript::Data* normal)
{
    Dudley_resetError();
    if (!nodes || !elements)
        return;

    const int NN = elements->numNodes;
    const int numDim = nodes->numDim;
    bool reduced_integration = Assemble_reducedIntegrationOrder(normal);
    const int numQuad = (!reduced_integration) ? (elements->numDim + 1) : 1;
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
        Dudley_setError(TYPE_ERROR, "Assemble_setNormal: Cannot calculate normal vector");
    } else if (!normal->isDataPointShapeEqual(1, &numDim)) {
        Dudley_setError(TYPE_ERROR, "Assemble_setNormal: illegal point data shape of normal Data object");
    } else if (!normal->numSamplesEqual(numQuad, elements->numElements)) {
        Dudley_setError(TYPE_ERROR, "Assemble_setNormal: illegal number of samples of normal Data object");
    } else if (!normal->actsExpanded()) {
        Dudley_setError(TYPE_ERROR, "Assemble_setNormal: expanded Data object is expected for normal.");
    }

    if (Dudley_noError()) {
        normal->requireWrite();
#pragma omp parallel
        {
            std::vector<double> local_X(NS * numDim);
            std::vector<double> dVdv(numQuad * numDim * numDim_local);
#pragma omp for
            for (index_t e = 0; e < elements->numElements; e++) {
                // gather local coordinates of nodes into local_X
                Dudley_Util_Gather_double(NS,
                        &elements->Nodes[INDEX2(0, e, NN)], numDim,
                        nodes->Coordinates, &local_X[0]);

                // calculate dVdv(i,j,q)=local_X(i,n)*DSDv(n,j,q)
                Dudley_Util_SmallMatMult(numDim, numDim_local * numQuad,
                                         &dVdv[0], NS, &local_X[0], dSdv);
                // get normalized vector
                double* normal_array = normal->getSampleDataRW(e);
                Dudley_NormalVector(numQuad, numDim, numDim_local, &dVdv[0], normal_array);
            }
        }
    }
}

} // namespace dudley

