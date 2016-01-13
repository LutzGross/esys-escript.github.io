
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

  Assemblage routines: integrates data on quadrature points

*****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"


#include "Assemble.h"
#include "Util.h"

namespace finley {

void Assemble_integrate(const NodeFile* nodes, const ElementFile* elements,
                        const escript::Data& data, double* out)
{
    resetError();
    if (!nodes || !elements)
        return;

    Esys_MPI_rank my_mpi_rank = nodes->MPIInfo->rank;
    ElementFile_Jacobians *jac = elements->borrowJacobians(nodes, FALSE,
                                    util::hasReducedIntegrationOrder(data));
    if (noError()) {
        const int numQuadTotal = jac->numQuadTotal;
        // check the shape of the data
        if (!data.numSamplesEqual(numQuadTotal, elements->numElements)) {
            setError(TYPE_ERROR, "Assemble_integrate: illegal number of samples of integrant kernel Data object");
            return;
        }

        const int numComps = data.getDataPointSize();

        for (int q=0; q<numComps; q++)
            out[q]=0;

#pragma omp parallel
        {
            std::vector<double> out_local(numComps);

            if (data.actsExpanded()) {
#pragma omp for
                for (int e=0; e<elements->numElements; e++) {
                    if (elements->Owner[e] == my_mpi_rank) {
                        const double *data_array=data.getSampleDataRO(e);
                        for (int q=0; q<numQuadTotal; q++) {
                            for (int i=0; i<numComps; i++)
                                out_local[i]+=data_array[INDEX2(i,q,numComps)]*jac->volume[INDEX2(q,e,numQuadTotal)];
                        }
                    }
                }
            } else {
#pragma omp for
                for (int e=0; e<elements->numElements; e++) {
                    if (elements->Owner[e] == my_mpi_rank) {
                        const double *data_array=data.getSampleDataRO(e);
                        double rtmp=0.;
                        for (int q=0; q<numQuadTotal; q++)
                            rtmp+=jac->volume[INDEX2(q,e,numQuadTotal)];
                        for (int i=0; i<numComps; i++)
                            out_local[i]+=data_array[i]*rtmp;
                    }
                }
            }
            // add local results to global result
#pragma omp critical
            for (int i=0; i<numComps; i++)
                out[i]+=out_local[i];
        }
    }
}

} // namespace finley

