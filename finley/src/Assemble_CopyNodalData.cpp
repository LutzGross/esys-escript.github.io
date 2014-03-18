
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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

  Assemblage routines: copies data between different types of nodal
  representations

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

namespace finley {

void Assemble_CopyNodalData(const NodeFile* nodes, escript::Data& out,
                            const escript::Data& in)
{
    resetError();
    if (!nodes)
        return;

    const int mpiSize = nodes->MPIInfo->size;
    const int numComps = out.getDataPointSize();
    const int in_data_type=in.getFunctionSpace().getTypeCode();
    const int out_data_type=out.getFunctionSpace().getTypeCode();

    // check out and in
    if (numComps != in.getDataPointSize()) {
        setError(TYPE_ERROR,"Assemble_CopyNodalData: number of components of input and output Data do not match.");
    } else if (!out.actsExpanded()) {
        setError(TYPE_ERROR,"Assemble_CopyNodalData: expanded Data object is expected for output data.");
    }

    // more sophisticated test needed for overlapping node/DOF counts
    if (in_data_type == FINLEY_NODES) {
        if (!in.numSamplesEqual(1, nodes->getNumNodes())) {
            setError(TYPE_ERROR,"Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
    } else if (in_data_type == FINLEY_REDUCED_NODES) {
        if (!in.numSamplesEqual(1, nodes->getNumReducedNodes())) {
            setError(TYPE_ERROR,"Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
    } else if (in_data_type == FINLEY_DEGREES_OF_FREEDOM) {
        if (!in.numSamplesEqual(1, nodes->getNumDegreesOfFreedom())) {
            setError(TYPE_ERROR,"Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
        if (((out_data_type == FINLEY_NODES) || (out_data_type == FINLEY_DEGREES_OF_FREEDOM)) && !in.actsExpanded() && (mpiSize>1)) {
            setError(TYPE_ERROR,"Assemble_CopyNodalData: FINLEY_DEGREES_OF_FREEDOM to FINLEY_NODES or FINLEY_DEGREES_OF_FREEDOM requires expanded input data on more than one processor.");
        }
    } else if (in_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (!in.numSamplesEqual(1, nodes->getNumReducedDegreesOfFreedom())) {
            setError(TYPE_ERROR,"Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
        if ((out_data_type == FINLEY_DEGREES_OF_FREEDOM) && !in.actsExpanded() && (mpiSize>1)) {
            setError(TYPE_ERROR,"Assemble_CopyNodalData: FINLEY_REDUCED_DEGREES_OF_FREEDOM to FINLEY_DEGREES_OF_FREEDOM requires expanded input data on more than one processor.");
        }
    } else {
        setError(TYPE_ERROR, "Assemble_CopyNodalData: illegal function space type for target object");
    }

    int numOut=0;
    switch (out_data_type) {
        case FINLEY_NODES:
            numOut=nodes->getNumNodes();
            break;

        case FINLEY_REDUCED_NODES:
            numOut=nodes->getNumReducedNodes();
            break;

        case FINLEY_DEGREES_OF_FREEDOM:
            numOut=nodes->getNumDegreesOfFreedom();
            break;

        case FINLEY_REDUCED_DEGREES_OF_FREEDOM:
            numOut=nodes->getNumReducedDegreesOfFreedom();
            break;

        default:
            setError(TYPE_ERROR,"Assemble_CopyNodalData: illegal function space type for source object");
    }

    if (!out.numSamplesEqual(1, numOut)) {
        setError(TYPE_ERROR,"Assemble_CopyNodalData: illegal number of samples of output Data object");
    }

    if (!noError())
        return;

    const size_t numComps_size = numComps*sizeof(double);

    /*********************** FINLEY_NODES ********************************/
    if (in_data_type == FINLEY_NODES) {
        out.requireWrite();
        if (out_data_type == FINLEY_NODES) {
#pragma omp parallel for
            for (int n=0; n<numOut; n++) {
                memcpy(out.getSampleDataRW(n), in.getSampleDataRO(n), numComps_size);
            }
        } else if (out_data_type == FINLEY_REDUCED_NODES) {
            const std::vector<int>& map = nodes->borrowReducedNodesTarget();
#pragma omp parallel for
            for (int n=0; n<map.size(); n++) {
                memcpy(out.getSampleDataRW(n), in.getSampleDataRO(map[n]),
                       numComps_size);
            }
        } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
            const std::vector<int>& map = nodes->borrowDegreesOfFreedomTarget();
#pragma omp parallel for
            for (int n=0; n<numOut; n++) {
                memcpy(out.getSampleDataRW(n), in.getSampleDataRO(map[n]),
                       numComps_size);
            }
        } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            const std::vector<int>& map = nodes->borrowReducedDegreesOfFreedomTarget();
#pragma omp parallel for
            for (int n=0; n<numOut; n++) {
                memcpy(out.getSampleDataRW(n), in.getSampleDataRO(map[n]),
                       numComps_size);
            }
        }

    /*********************** FINLEY_REDUCED_NODES ***************************/
    } else if (in_data_type == FINLEY_REDUCED_NODES) {
        if (out_data_type == FINLEY_NODES) {
            setError(TYPE_ERROR,"Assemble_CopyNodalData: cannot copy from reduced nodes to nodes.");
        } else if (out_data_type == FINLEY_REDUCED_NODES) {
            out.requireWrite();
#pragma omp parallel for
            for (int n=0; n<nodes->getNumNodes(); n++) {
                memcpy(out.getSampleDataRW(n), in.getSampleDataRO(n), numComps_size);
            }
       } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
            setError(TYPE_ERROR,"Assemble_CopyNodalData: cannot copy from reduced nodes to degrees of freedom.");
       } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            out.requireWrite();
            const int* target = nodes->borrowTargetReducedNodes();
            const std::vector<int>& map = nodes->borrowReducedDegreesOfFreedomTarget();
#pragma omp parallel for
            for (int n=0; n<numOut; n++) {
               memcpy(out.getSampleDataRW(n),
                      in.getSampleDataRO(target[map[n]]), numComps_size);
            }
        }

    /******************** FINLEY_DEGREES_OF_FREEDOM *********************/
    } else if (in_data_type == FINLEY_DEGREES_OF_FREEDOM) {
        out.requireWrite();
        if (out_data_type == FINLEY_NODES) {
            Paso_Coupler *coupler = Paso_Coupler_alloc(nodes->degreesOfFreedomConnector, numComps);
            if (Esys_noError()) {
                // Coupler holds the pointer but it doesn't appear to get
                // used so RO should work.
                Paso_Coupler_startCollect(coupler, in.getDataRO());
                const double *recv_buffer=Paso_Coupler_finishCollect(coupler);
                const int upperBound=nodes->getNumDegreesOfFreedom();
                const int* target = nodes->borrowTargetDegreesOfFreedom();
#pragma omp parallel for
                for (int n=0; n<nodes->numNodes; n++) {
                    const int k=target[n];
                    if (k < upperBound) {
                        memcpy(out.getSampleDataRW(n), in.getSampleDataRO(k),
                               numComps_size);
                    } else {
                        memcpy(out.getSampleDataRW(n),
                               &recv_buffer[(k-upperBound)*numComps],
                               numComps_size);
                    }
                }
            }
            Paso_Coupler_free(coupler);
        } else if  (out_data_type == FINLEY_REDUCED_NODES) {
            Paso_Coupler *coupler = Paso_Coupler_alloc(nodes->degreesOfFreedomConnector, numComps);
            if (Esys_noError()) {
                Paso_Coupler_startCollect(coupler, in.getDataRO());
                const double *recv_buffer=Paso_Coupler_finishCollect(coupler);
                const int upperBound=nodes->getNumDegreesOfFreedom();
                const std::vector<int>& map = nodes->borrowReducedNodesTarget();
                const int* target = nodes->borrowTargetDegreesOfFreedom();

#pragma omp parallel for
                for (int n=0; n<map.size(); n++) {
                    const int k=target[map[n]];
                    if (k < upperBound) {
                        memcpy(out.getSampleDataRW(n), in.getSampleDataRO(k),
                               numComps_size);
                    } else {
                        memcpy(out.getSampleDataRW(n),
                               &recv_buffer[(k-upperBound)*numComps],
                               numComps_size);
                    }
                }
            }
            Paso_Coupler_free(coupler);
        } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
#pragma omp parallel for
            for (int n=0; n<numOut; n++) {
                memcpy(out.getSampleDataRW(n), in.getSampleDataRO(n),
                       numComps_size);
            }
        } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            const std::vector<int>& map = nodes->borrowReducedDegreesOfFreedomTarget();
            const int* target = nodes->borrowTargetDegreesOfFreedom();
#pragma omp parallel for
            for (int n=0; n<numOut; n++) {
                memcpy(out.getSampleDataRW(n),
                       in.getSampleDataRO(target[map[n]]), numComps_size);
            }
        }

    /**************** FINLEY_REDUCED_DEGREES_OF_FREEDOM *****************/
    } else if (in_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (out_data_type == FINLEY_NODES) {
            setError(TYPE_ERROR, "Assemble_CopyNodalData: cannot copy from reduced degrees of freedom to nodes.");
        } else if (out_data_type == FINLEY_REDUCED_NODES) {
            Paso_Coupler *coupler=Paso_Coupler_alloc(nodes->reducedDegreesOfFreedomConnector,numComps);
            if (Esys_noError()) {
                Paso_Coupler_startCollect(coupler, in.getDataRO());
                out.requireWrite();
                const int upperBound=nodes->getNumReducedDegreesOfFreedom();
                const std::vector<int>& map=nodes->borrowReducedNodesTarget();
                const int* target=nodes->borrowTargetReducedDegreesOfFreedom();
                const double *recv_buffer=Paso_Coupler_finishCollect(coupler);
#pragma omp parallel for
                for (int n=0; n<map.size(); n++) {
                    const int k=target[map[n]];
                    if (k < upperBound) {
                        memcpy(out.getSampleDataRW(n), in.getSampleDataRO(k),
                               numComps_size);
                    } else {
                        memcpy(out.getSampleDataRW(n),
                               &recv_buffer[(k-upperBound)*numComps],
                               numComps_size);
                    }
                }
            }
            Paso_Coupler_free(coupler);
        } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            out.requireWrite();
#pragma omp parallel for
            for (int n=0; n<numOut; n++) {
                memcpy(out.getSampleDataRW(n), in.getSampleDataRO(n), numComps_size);
            }
        } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM ) {
            setError(TYPE_ERROR, "Assemble_CopyNodalData: cannot copy from reduced degrees of freedom to degrees of freedom.");
        }
    } // in_data_type
}

} // namespace finley

