
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
#include "Util.h"

#include <escript/index.h>

namespace dudley {

void Assemble_integrate(const NodeFile* nodes, const ElementFile* elements,
                        const escript::Data& data, std::vector<double>& out)
{
    if (!nodes || !elements)
        return;

    if (data.isComplex())
    {
        throw DudleyException("Assemble_integrate: complex arguments not supported.");
    }
    const int my_mpi_rank = nodes->MPIInfo->rank;
    const ElementFile_Jacobians* jac = elements->borrowJacobians(nodes,
                                         hasReducedIntegrationOrder(data));

    const int numQuadTotal = jac->numQuad;
    // check the shape of the data
    if (!data.numSamplesEqual(numQuadTotal, elements->numElements)) {
        throw DudleyException("Assemble_integrate: illegal number of samples of integrant kernel Data object");
    }

    const int numComps = data.getDataPointSize();

    for (int q = 0; q < numComps; q++)
        out[q] = 0;

#pragma omp parallel
    {
        std::vector<double> out_local(numComps);

        if (data.actsExpanded()) {
#pragma omp for
            for (index_t e = 0; e < elements->numElements; e++) {
                if (elements->Owner[e] == my_mpi_rank) {
                    const double vol = jac->absD[e] * jac->quadweight;
                    const double* data_array = data.getSampleDataRO(e, static_cast<escript::DataTypes::real_t>(0));
                    for (int q = 0; q < numQuadTotal; q++) {
                        for (int i = 0; i < numComps; i++)
                            out_local[i] += data_array[INDEX2(i, q, numComps)] * vol;
                    }
                }
            }
        } else {
#pragma omp for
            for (index_t e = 0; e < elements->numElements; e++) {
                if (elements->Owner[e] == my_mpi_rank) {
                    const double vol = jac->absD[e] * jac->quadweight;
                    const double* data_array = data.getSampleDataRO(e, static_cast<escript::DataTypes::real_t>(0));
                    double rtmp = 0.;
                    for (int q = 0; q < numQuadTotal; q++)
                        rtmp += vol;
                    for (int i = 0; i < numComps; i++)
                        out_local[i] += data_array[i] * rtmp;
                }
            }
        }
        // add local results to global result
#pragma omp critical
        for (int i = 0; i < numComps; i++)
            out[i] += out_local[i];
    }
}

} // namespace dudley

