
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

/****************************************************************************

  Assemblage routines: integrates data on quadrature points

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

namespace dudley {

void Assemble_integrate(Dudley_NodeFile* nodes, Dudley_ElementFile* elements, const escript::Data* data, double* out)
{
    if (!nodes || !elements)
        return;

    const int my_mpi_rank = nodes->MPIInfo->rank;
    Dudley_ElementFile_Jacobians* jac = Dudley_ElementFile_borrowJacobians(
            elements, nodes, Assemble_reducedIntegrationOrder(data));

    const dim_t numQuadTotal = jac->numQuad;
    // check the shape of the data
    if (!data->numSamplesEqual(numQuadTotal, elements->numElements)) {
        throw DudleyException("Assemble_integrate: illegal number of samples of integrant kernel Data object");
    }

    const int numComps = data->getDataPointSize();

    for (int q = 0; q < numComps; q++)
        out[q] = 0;

#pragma omp parallel
    {
        std::vector<double> out_local(numComps);

        if (data->actsExpanded()) {
#pragma omp for
            for (index_t e = 0; e < elements->numElements; e++) {
                if (elements->Owner[e] == my_mpi_rank) {
                    const double vol = jac->absD[e] * jac->quadweight;
                    const double* data_array = data->getSampleDataRO(e);
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
                    const double* data_array = data->getSampleDataRO(e);
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

