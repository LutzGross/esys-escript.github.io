
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

  Assemblage routines: copies data between different types nodal
  representations

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

namespace dudley {

void Assemble_CopyNodalData(Dudley_NodeFile* nodes, escript::Data* out, const escript::Data* in)
{
    if (!nodes)
        return;

    const int mpiSize = nodes->MPIInfo->size;
    const int numComps = out->getDataPointSize();
    const int in_data_type = in->getFunctionSpace().getTypeCode();
    const int out_data_type = out->getFunctionSpace().getTypeCode();

    // check out and in
    if (numComps != in->getDataPointSize()) {
        throw DudleyException("Assemble_CopyNodalData: number of components of input and output Data do not match.");
    } else if (!out->actsExpanded()) {
        throw DudleyException("Assemble_CopyNodalData: expanded Data object is expected for output data.");
    }

    // more sophisticated test needed for overlapping node/DOF counts
    if (in_data_type == DUDLEY_NODES) {
        if (!in->numSamplesEqual(1, Dudley_NodeFile_getNumNodes(nodes))) {
            throw DudleyException("Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
    } else if (in_data_type == DUDLEY_REDUCED_NODES) {
        if (!in->numSamplesEqual(1, Dudley_NodeFile_getNumReducedNodes(nodes))) {
            throw DudleyException("Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
    } else if (in_data_type == DUDLEY_DEGREES_OF_FREEDOM) {
        if (!in->numSamplesEqual(1, Dudley_NodeFile_getNumDegreesOfFreedom(nodes))) {
            throw DudleyException("Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
        if ((((out_data_type == DUDLEY_NODES) || (out_data_type == DUDLEY_DEGREES_OF_FREEDOM)) && !in->actsExpanded() && (mpiSize > 1))) {

            throw DudleyException("Assemble_CopyNodalData: DUDLEY_DEGREES_OF_FREEDOM to DUDLEY_NODES or DUDLEY_DEGREES_OF_FREEDOM requires expanded input data on more than one processor.");
        }
    } else if (in_data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (!in->numSamplesEqual(1, Dudley_NodeFile_getNumReducedDegreesOfFreedom(nodes))) {
            throw DudleyException("Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
        if ((out_data_type == DUDLEY_DEGREES_OF_FREEDOM) && !in->actsExpanded() && (mpiSize > 1)) {
            throw DudleyException("Assemble_CopyNodalData: DUDLEY_REDUCED_DEGREES_OF_FREEDOM to DUDLEY_DEGREES_OF_FREEDOM requires expanded input data on more than one processor.");
        }
    } else {
        throw DudleyException("Assemble_CopyNodalData: illegal function space type for target object");
    }

    if (out_data_type == DUDLEY_NODES) {
        if (!out->numSamplesEqual(1, Dudley_NodeFile_getNumNodes(nodes))) {
            throw DudleyException("Assemble_CopyNodalData: illegal number of samples of output Data object");
        }
    } else if (out_data_type == DUDLEY_REDUCED_NODES) {
        if (!out->numSamplesEqual(1, Dudley_NodeFile_getNumReducedNodes(nodes))) {
            throw DudleyException("Assemble_CopyNodalData: illegal number of samples of output Data object");
        }
    } else if (out_data_type == DUDLEY_DEGREES_OF_FREEDOM) {
        if (!out->numSamplesEqual(1, Dudley_NodeFile_getNumDegreesOfFreedom(nodes))) {
            throw DudleyException("Assemble_CopyNodalData: illegal number of samples of output Data object");
        }
    } else if (out_data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (!out->numSamplesEqual(1, Dudley_NodeFile_getNumReducedDegreesOfFreedom(nodes))) {
            throw DudleyException("Assemble_CopyNodalData: illegal number of samples of output Data object");
        }
    } else {
        throw DudleyException("Assemble_CopyNodalData: illegal function space type for source object");
    }

    size_t numComps_size = numComps * sizeof(double);

    /**************************** DUDLEY_NODES ******************************/
    if (in_data_type == DUDLEY_NODES) {
        out->requireWrite();
        if (out_data_type == DUDLEY_NODES) {
#pragma omp parallel for
            for (index_t n = 0; n < nodes->nodesMapping->numNodes; n++) {
                memcpy(out->getSampleDataRW(n), in->getSampleDataRO(n), numComps_size);
            }
        } else if (out_data_type == DUDLEY_REDUCED_NODES) {
#pragma omp parallel for
            for (index_t n = 0; n < nodes->reducedNodesMapping->numTargets; n++) {
                memcpy(out->getSampleDataRW(n),
                       in->getSampleDataRO(nodes->reducedNodesMapping->map[n]),
                       numComps_size);
            }
        } else if (out_data_type == DUDLEY_DEGREES_OF_FREEDOM) {
            const dim_t nComps = nodes->degreesOfFreedomDistribution->getMyNumComponents();
#pragma omp parallel for
            for (index_t n = 0; n < nComps; n++) {
                memcpy(out->getSampleDataRW(n),
                       in->getSampleDataRO(nodes->degreesOfFreedomMapping->map[n]), numComps_size);
            }
        } else if (out_data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM) {
            const dim_t nComps = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
#pragma omp parallel for
            for (index_t n = 0; n < nComps; n++) {
                memcpy(out->getSampleDataRW(n),
                       in->getSampleDataRO(nodes->reducedDegreesOfFreedomMapping->map[n]), numComps_size);
            }
        }
    /************************ DUDLEY_REDUCED_NODES **************************/
    } else if (in_data_type == DUDLEY_REDUCED_NODES) {
        out->requireWrite();
        if (out_data_type == DUDLEY_NODES) {
            throw DudleyException("Assemble_CopyNodalData: cannot copy from reduced nodes to nodes.");
        } else if (out_data_type == DUDLEY_REDUCED_NODES) {
#pragma omp parallel for
            for (index_t n = 0; n < nodes->reducedNodesMapping->numNodes; n++) {
                memcpy(out->getSampleDataRW(n), in->getSampleDataRO(n), numComps_size);
            }
        } else if (out_data_type == DUDLEY_DEGREES_OF_FREEDOM) {
            throw DudleyException("Assemble_CopyNodalData: cannot copy from reduced nodes to degrees of freedom.");
        } else if (out_data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM) {
            const dim_t nComps = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
#pragma omp parallel for
            for (index_t n = 0; n < nComps; n++) {
                const dim_t k = nodes->reducedDegreesOfFreedomMapping->map[n];
                memcpy(out->getSampleDataRW(n),
                       in->getSampleDataRO(nodes->reducedNodesMapping->target[k]), numComps_size);
            }
        }
    /********************** DUDLEY_DEGREES_OF_FREEDOM ***********************/
    } else if (in_data_type == DUDLEY_DEGREES_OF_FREEDOM) {
        out->requireWrite();
        if (out_data_type == DUDLEY_NODES) {
            paso::Coupler_ptr coupler(new paso::Coupler(nodes->degreesOfFreedomConnector, numComps));
            // safe provided coupler->copyAll is called before the pointer
            // in "in" is invalidated
            const_cast<escript::Data*>(in)->resolve();
            coupler->startCollect(in->getDataRO());  
            const double* recv_buffer = coupler->finishCollect();
            const index_t upperBound = nodes->degreesOfFreedomDistribution->getMyNumComponents();
#pragma omp parallel for
            for (index_t n = 0; n < nodes->numNodes; n++) {
                const index_t k = nodes->degreesOfFreedomMapping->target[n];
                if (k < upperBound) {
                    memcpy(out->getSampleDataRW(n), in->getSampleDataRO(k), numComps_size);
                } else {
                    memcpy(out->getSampleDataRW(n),
                           &recv_buffer[(k - upperBound) * numComps], numComps_size);
                }
            }
        } else if (out_data_type == DUDLEY_REDUCED_NODES) {
            paso::Coupler_ptr coupler(new paso::Coupler(nodes->degreesOfFreedomConnector, numComps));
            // safe provided coupler->copyAll is called before the pointer
            // in "in" is invalidated
            const_cast<escript::Data*>(in)->resolve();
            coupler->startCollect(in->getDataRO());  
            const double* recv_buffer = coupler->finishCollect();
            const index_t upperBound = nodes->degreesOfFreedomDistribution->getMyNumComponents();

#pragma omp parallel for
            for (index_t n = 0; n < nodes->reducedNodesMapping->numTargets; n++) {
                const index_t l = nodes->reducedNodesMapping->map[n];
                const index_t k = nodes->degreesOfFreedomMapping->target[l];
                if (k < upperBound) {
                    memcpy(out->getSampleDataRW(n), in->getSampleDataRO(k), numComps_size);
                } else {
                    memcpy(out->getSampleDataRW(n),
                           &recv_buffer[(k - upperBound)*numComps], numComps_size);
                }
            }
        } else if (out_data_type == DUDLEY_DEGREES_OF_FREEDOM) {
            const dim_t nComps = nodes->degreesOfFreedomDistribution->getMyNumComponents();
#pragma omp parallel for
            for (index_t n = 0; n < nComps; n++) {
                memcpy(out->getSampleDataRW(n), in->getSampleDataRO(n), numComps_size);
            }
        } else if (out_data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM) {
            const dim_t nComps = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
#pragma omp parallel for
            for (index_t n = 0; n < nComps; n++) {
                const dim_t k = nodes->reducedDegreesOfFreedomMapping->map[n];
                memcpy(out->getSampleDataRW(n),
                       in->getSampleDataRO(nodes->degreesOfFreedomMapping->target[k]), numComps_size);
            }
        }

    /****************** DUDLEY_REDUCED_DEGREES_OF_FREEDOM *******************/
    } else if (in_data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (out_data_type == DUDLEY_NODES) {
            throw DudleyException("Assemble_CopyNodalData: cannot copy from reduced degrees of freedom to nodes.");
        } else if (out_data_type == DUDLEY_REDUCED_NODES) {
            paso::Coupler_ptr coupler(new paso::Coupler(nodes->reducedDegreesOfFreedomConnector, numComps));
            const dim_t upperBound = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
            // safe provided coupler->copyAll is called before the pointer
            // in "in" is invalidated
            const_cast<escript::Data*>(in)->resolve();
            coupler->startCollect(in->getDataRO());  
            out->requireWrite();
            const double* recv_buffer = coupler->finishCollect();
#pragma omp parallel for
            for (index_t n = 0; n < nodes->reducedNodesMapping->numTargets; n++) {
                const index_t l = nodes->reducedNodesMapping->map[n];
                const index_t k = nodes->reducedDegreesOfFreedomMapping->target[l];
                if (k < upperBound) {
                    memcpy(out->getSampleDataRW(n), in->getSampleDataRO(k), numComps_size);
                } else {
                    memcpy(out->getSampleDataRW(n),
                           &recv_buffer[(k - upperBound)*numComps], numComps_size);
                }
            }
        } else if (out_data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM) {
            const dim_t nComps = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
            out->requireWrite();
#pragma omp parallel for
            for (index_t n = 0; n < nComps; n++) {
                memcpy(out->getSampleDataRW(n), in->getSampleDataRO(n), numComps_size);
            }
        } else if (out_data_type == DUDLEY_DEGREES_OF_FREEDOM) {
            throw DudleyException("Assemble_CopyNodalData: cannot copy from reduced degrees of freedom to degrees of freedom.");
        }
    }
}

} // namespace dudley

