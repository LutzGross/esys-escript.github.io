
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

void Assemble_interpolate(const NodeFile* nodes, const ElementFile* elements,
                          const escript::Data& data,
                          escript::Data& interpolated_data)
{
    if (!nodes || !elements)
        return;

    if (data.isComplex() || interpolated_data.isComplex())
    {
        throw DudleyException("Assemble_interpolate: complex arguments are not supported.");
    }
    const int data_type = data.getFunctionSpace().getTypeCode();
    const bool reduced_integration = hasReducedIntegrationOrder(interpolated_data);

    dim_t numNodes = 0;
    const index_t* map = NULL;

    if (data_type == DUDLEY_NODES) {
        numNodes = nodes->getNumNodes();
        map = nodes->borrowTargetNodes();
    } else if (data_type == DUDLEY_DEGREES_OF_FREEDOM) {
        if (elements->MPIInfo->size > 1) {
            throw DudleyException("Assemble_interpolate: for more than one processor DEGREES_OF_FREEDOM data are not accepted as input.");
        }
        numNodes = nodes->getNumDegreesOfFreedom();
        map = nodes->borrowTargetDegreesOfFreedom();
    } else {
        throw DudleyException("Assemble_interpolate: Cannot interpolate data");
    }

    const int numComps = data.getDataPointSize();
    const int NN = elements->numNodes;
    const int numQuad = reduced_integration ? 1 : elements->numNodes;
    const int NS_DOF = elements->numDim + 1;
    const double *shapeFns = NULL;

    // check the dimensions of interpolated_data and data
    if (!interpolated_data.numSamplesEqual(numQuad, elements->numElements)) {
        throw DudleyException("Assemble_interpolate: illegal number of samples of output Data object");
    } else if (!data.numSamplesEqual(1, numNodes)) {
        throw DudleyException("Assemble_interpolate: illegal number of samples of input Data object");
    } else if (numComps != interpolated_data.getDataPointSize()) {
        throw DudleyException("Assemble_interpolate: number of components of input and interpolated Data do not match.");
    } else if (!interpolated_data.actsExpanded()) {
        throw DudleyException("Assemble_interpolate: expanded Data object is expected for output data.");
    }

    if (!getQuadShape(elements->numDim, reduced_integration, &shapeFns)) {
        throw DudleyException("Assemble_interpolate: unable to locate shape function.");
    }

    interpolated_data.requireWrite();
#pragma omp parallel
    {
        std::vector<double> local_data(NS_DOF * numComps);
        const size_t numComps_size = numComps * sizeof(double);
        // open the element loop
#pragma omp for
        for (index_t e = 0; e < elements->numElements; e++) {
            for (int q = 0; q < NS_DOF; q++) {
                const index_t i = elements->Nodes[INDEX2(q, e, NN)];
                const double* data_array = data.getSampleDataRO(map[i], static_cast<escript::DataTypes::real_t>(0));
                memcpy(&local_data[INDEX3(0, q, 0, numComps, NS_DOF)],
                       data_array, numComps_size);
            }
            // calculate interpolated_data=local_data*S
            util::smallMatSetMult1(1, numComps, numQuad,
                            interpolated_data.getSampleDataRW(e, static_cast<escript::DataTypes::real_t>(0)), NS_DOF,
                            &local_data[0], shapeFns);
        } // end of element loop
    } // end of parallel region
}

} // namespace dudley

