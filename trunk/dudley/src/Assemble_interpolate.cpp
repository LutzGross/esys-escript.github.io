
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

  Assemblage routines: interpolates nodal data in a data array onto elements
  (=integration points)

*****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Assemble.h"
#include "ShapeTable.h"
#include "Util.h"

namespace dudley {

void Assemble_interpolate(Dudley_NodeFile* nodes,
                          Dudley_ElementFile* elements,
                          const escript::Data* data,
                          escript::Data* interpolated_data)
{
    if (!nodes || !elements)
        return;

    const int data_type = data->getFunctionSpace().getTypeCode();
    const bool reduced_integration = Assemble_reducedIntegrationOrder(interpolated_data);

    dim_t numNodes = 0;
    index_t *map = NULL;

    if (data_type == DUDLEY_NODES) {
        numNodes = Dudley_NodeFile_getNumNodes(nodes);
        map = Dudley_NodeFile_borrowTargetNodes(nodes);
    } else if (data_type == DUDLEY_REDUCED_NODES) {
        numNodes = Dudley_NodeFile_getNumReducedNodes(nodes);
        map = Dudley_NodeFile_borrowTargetReducedNodes(nodes);
    } else if (data_type == DUDLEY_DEGREES_OF_FREEDOM) {
        if (elements->MPIInfo->size > 1) {
            throw DudleyException("Assemble_interpolate: for more than one processor DEGREES_OF_FREEDOM data are not accepted as input.");
        }
        numNodes = Dudley_NodeFile_getNumDegreesOfFreedom(nodes);
        map = Dudley_NodeFile_borrowTargetDegreesOfFreedom(nodes);
    } else if (data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (elements->MPIInfo->size > 1) {
            throw DudleyException("Assemble_interpolate: for more than one processor REDUCED_DEGREES_OF_FREEDOM data are not accepted as input.");
        }
        numNodes = Dudley_NodeFile_getNumReducedDegreesOfFreedom(nodes);
        map = Dudley_NodeFile_borrowTargetReducedDegreesOfFreedom(nodes);
    } else {
        throw DudleyException("Assemble_interpolate: Cannot interpolate data");
    }

    const dim_t numComps = data->getDataPointSize();
    const int NN = elements->numNodes;
    const int numQuad = reduced_integration ? 1 : (elements->numDim + 1);
    const int NS_DOF = elements->numDim + 1;
    const double *shapeFns = NULL;

    // check the dimensions of interpolated_data and data
    if (!interpolated_data->numSamplesEqual(numQuad, elements->numElements)) {
        throw DudleyException("Assemble_interpolate: illegal number of samples of output Data object");
    } else if (!data->numSamplesEqual(1, numNodes)) {
        throw DudleyException("Assemble_interpolate: illegal number of samples of input Data object");
    } else if (numComps != interpolated_data->getDataPointSize()) {
        throw DudleyException("Assemble_interpolate: number of components of input and interpolated Data do not match.");
    } else if (!interpolated_data->actsExpanded()) {
        throw DudleyException("Assemble_interpolate: expanded Data object is expected for output data.");
    }

    if (!getQuadShape(elements->numDim, reduced_integration, &shapeFns)) {
        throw DudleyException("Assemble_interpolate: unable to locate shape function.");
    }

    interpolated_data->requireWrite();
#pragma omp parallel
    {
        std::vector<double> local_data(NS_DOF * numComps);
        const size_t numComps_size = numComps *sizeof(double);
        /* open the element loop */
#pragma omp for
        for (index_t e = 0; e < elements->numElements; e++) {
            for (int q = 0; q < NS_DOF; q++) {
                const index_t i = elements->Nodes[INDEX2(q, e, NN)];
                const double* data_array = data->getSampleDataRO(map[i]);
                memcpy(&local_data[INDEX3(0, q, 0, numComps, NS_DOF)],
                       data_array, numComps_size);
            }
            // calculate interpolated_data=local_data*S
            Dudley_Util_SmallMatSetMult1(1, numComps, numQuad,
                            interpolated_data->getSampleDataRW(e), NS_DOF,
                            &local_data[0], shapeFns);
        } // end of element loop
    } // end of parallel region
}

} // namespace dudley

