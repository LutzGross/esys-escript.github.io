
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


/****************************************************************************

  Assemblage routines: integrates data on quadrature points

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <escript/index.h>
#include <escript/Utils.h>

namespace finley {

template<typename Scalar>
void Assemble_integrate(const NodeFile* nodes, const ElementFile* elements,
                        const escript::Data& data, Scalar* out)
{
    if (!nodes || !elements)
        return;

    const int my_mpi_rank = nodes->MPIInfo->rank;
    ElementFile_Jacobians* jac = elements->borrowJacobians(nodes, false,
                                    util::hasReducedIntegrationOrder(data));

    const int numQuadTotal = jac->numQuadTotal;
    // check the shape of the data
    if (!data.numSamplesEqual(numQuadTotal, elements->numElements)) {
        throw escript::ValueError("Assemble_integrate: illegal number of samples of integrant kernel Data object");
    }

    const int numComps = data.getDataPointSize();
    const Scalar zero = static_cast<Scalar>(0);

    for (int q = 0; q < numComps; q++)
        out[q] = zero;

#pragma omp parallel
    {
        std::vector<Scalar> out_local(numComps);
        {
            if (data.actsExpanded()) {
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    if (elements->Owner[e] == my_mpi_rank) {
                        const Scalar* data_array = data.getSampleDataRO(e, zero);
                        for (int q = 0; q < numQuadTotal; q++) {
                            for (int i = 0; i < numComps; i++)
                                out_local[i] += data_array[INDEX2(i,q,numComps)]*jac->volume[INDEX2(q,e,numQuadTotal)];
                        }
                    }
                }
            } else {
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    if (elements->Owner[e] == my_mpi_rank) {
                        const Scalar* data_array = data.getSampleDataRO(e, zero);
                        double rtmp = 0.;
                        for (int q = 0; q < numQuadTotal; q++)
                            rtmp += jac->volume[INDEX2(q, e, numQuadTotal)];
                        for (int i = 0; i < numComps; i++)
                            out_local[i] += data_array[i] * rtmp;
                    }
                }
            }
        }
        // add local results to global result
#pragma omp critical
        for (int i = 0; i < numComps; i++)
            out[i] += out_local[i];
} // parallel section
}
// instantiate our two supported versions
template void Assemble_integrate<escript::DataTypes::real_t>(
                    const NodeFile* nodes, const ElementFile* elements,
                    const escript::Data& data, escript::DataTypes::real_t* out);
template void Assemble_integrate<escript::DataTypes::cplx_t>(
                    const NodeFile* nodes, const ElementFile* elements,
                    const escript::Data& data, escript::DataTypes::cplx_t* out);

} // namespace finley
