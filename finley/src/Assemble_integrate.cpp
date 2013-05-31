
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/****************************************************************************

  Assemblage routines: integrates data on quadrature points

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <vector>

void Finley_Assemble_integrate(NodeFile* nodes,
                               ElementFile* elements,
                               escriptDataC* data, double* out)
{
    Finley_resetError();
    if (!nodes || !elements)
        return;

    Esys_MPI_rank my_mpi_rank = nodes->MPIInfo->rank;
    ElementFile_Jacobians *jac = elements->borrowJacobians(nodes, FALSE,
            Finley_Assemble_reducedIntegrationOrder(data));
    if (Finley_noError()) {
        const int numQuadTotal = jac->numQuadTotal;
        // check the shape of the data
        if (!numSamplesEqual(data,numQuadTotal,elements->numElements)) {
            Finley_setError(TYPE_ERROR, "Finley_Assemble_integrate: illegal number of samples of integrant kernel Data object");
            return;
        }

        const int numComps = getDataPointSize(data);

        // now we can start
        for (int q=0; q<numComps; q++)
            out[q]=0;

#pragma omp parallel
        {
            // initialize local result
            std::vector<double> out_local(numComps);

            // open the element loop
            if (isExpanded(data)) {
#pragma omp for
                for (int e=0; e<elements->numElements; e++) {
                    if (elements->Owner[e] == my_mpi_rank) {
                        const double *data_array=getSampleDataRO(data, e);
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
                        const double *data_array=getSampleDataRO(data,e);
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

