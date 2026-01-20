
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

  Assemblage routines: integrates data on quadrature points

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <escript/index.h>
#include <escript/Utils.h>

namespace finley {

template<typename Scalar>
void Assemble_integrate_points(const ElementFile* points,
                        const escript::Data& data, Scalar* out)
{
    if (!points)
        return;

    const int my_mpi_rank = points->MPIInfo->rank;
    // check the shape of the data
    if (!data.numSamplesEqual(1, points->numElements)) {
        throw escript::ValueError("Assemble_integrate_points: illegal number of samples of integrant kernel Data object");
    }

    const int numComps = data.getDataPointSize();
    const Scalar zero = static_cast<Scalar>(0);

    for (int q = 0; q < numComps; q++)
        out[q] = zero;

#pragma omp parallel
    {
        std::vector<Scalar> out_local(numComps);
        {
#pragma omp for
                for (index_t e = 0; e < points->numElements; e++) {
                    if (points->Owner[e] == my_mpi_rank) {
                        const Scalar* data_array = data.getSampleDataRO(e, zero);
                        for (int i = 0; i < numComps; i++)
                                out_local[i] += data_array[INDEX2(i,0,numComps)];
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
template void Assemble_integrate_points<escript::DataTypes::real_t>(
                    const ElementFile* points,
                    const escript::Data& data, escript::DataTypes::real_t* out);
template void Assemble_integrate_points<escript::DataTypes::cplx_t>(
                    const ElementFile* points,
                    const escript::Data& data, escript::DataTypes::cplx_t* out);

} // namespace finley
