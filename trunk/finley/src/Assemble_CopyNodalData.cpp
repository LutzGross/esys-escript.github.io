
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

  Assemblage routines: copies data between different types of nodal
  representations

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

void Finley_Assemble_CopyNodalData(Finley_NodeFile* nodes,
                                   escriptDataC* out, escriptDataC* in)
{
    Finley_resetError();
    if (!nodes)
        return;

    const int mpiSize = nodes->MPIInfo->size;
    const int numComps = getDataPointSize(out);
    const type_t in_data_type=getFunctionSpaceType(in);
    const type_t out_data_type=getFunctionSpaceType(out);

    // check out and in
    if (numComps != getDataPointSize(in)) {
        Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: number of components of input and output Data do not match.");
    } else if (!isExpanded(out)) {
        Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: expanded Data object is expected for output data.");
    }

    // more sophisticated test needed for overlapping node/DOF counts
    if (in_data_type == FINLEY_NODES) {
        if (!numSamplesEqual(in, 1, Finley_NodeFile_getNumNodes(nodes))) {
            Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
    } else if (in_data_type == FINLEY_REDUCED_NODES) {
        if (! numSamplesEqual(in,1,Finley_NodeFile_getNumReducedNodes(nodes))) {
            Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
    } else if (in_data_type == FINLEY_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(in,1,Finley_NodeFile_getNumDegreesOfFreedom(nodes))) {
            Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
        if ( (((out_data_type == FINLEY_NODES) || (out_data_type == FINLEY_DEGREES_OF_FREEDOM)) && !isExpanded(in) && (mpiSize>1))) {
            Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: FINLEY_DEGREES_OF_FREEDOM to FINLEY_NODES or FINLEY_DEGREES_OF_FREEDOM requires expanded input data on more than one processor.");
        }
    } else if (in_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(in,1,Finley_NodeFile_getNumReducedDegreesOfFreedom(nodes))) {
            Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
        if ( (out_data_type == FINLEY_DEGREES_OF_FREEDOM) && !isExpanded(in) && (mpiSize>1)) {

            Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: FINLEY_REDUCED_DEGREES_OF_FREEDOM to FINLEY_DEGREES_OF_FREEDOM requires expanded input data on more than one processor.");
        }
    } else {
        Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal function space type for target object");
    }

    if (out_data_type == FINLEY_NODES) {
        if (! numSamplesEqual(out,1,Finley_NodeFile_getNumNodes(nodes))) {
            Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of output Data object");
        }
    } else if (out_data_type == FINLEY_REDUCED_NODES) {
        if (! numSamplesEqual(out,1,Finley_NodeFile_getNumReducedNodes(nodes))) {
            Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of output Data object");
        }
    } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(out,1,Finley_NodeFile_getNumDegreesOfFreedom(nodes))) {
            Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of output Data object");
        }
    } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(out,1,Finley_NodeFile_getNumReducedDegreesOfFreedom(nodes))) {
            Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of output Data object");
        }
    } else {
        Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal function space type for source object");
    }

    if (!Finley_noError())
        return;

    const size_t numComps_size = numComps*sizeof(double);

    /*********************** FINLEY_NODES ********************************/
    if (in_data_type == FINLEY_NODES) {
        requireWrite(out);
        if (out_data_type == FINLEY_NODES) {
#pragma omp parallel for
            for (int n=0; n<nodes->nodesMapping->numNodes; n++) {
                memcpy(getSampleDataRWFast(out,n), getSampleDataROFast(in,n), numComps_size);
            }
        } else if (out_data_type == FINLEY_REDUCED_NODES) {
#pragma omp parallel for
            for (int n=0; n<nodes->reducedNodesMapping->numTargets; n++) {
                memcpy(getSampleDataRWFast(out,n),
                       getSampleDataROFast(in,nodes->reducedNodesMapping->map[n]),
                       numComps_size);
            }
        } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
            int nComps = Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution);
#pragma omp parallel for
            for (int n=0; n<nComps; n++) {
                memcpy(getSampleDataRWFast(out,n),
                       getSampleDataROFast(in,nodes->degreesOfFreedomMapping->map[n]),
                       numComps_size);
            }
        } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            int nComps = Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
#pragma omp parallel for
            for (int n=0; n<nComps; n++) {
                memcpy(getSampleDataRWFast(out,n),
                       getSampleDataROFast(in,nodes->reducedDegreesOfFreedomMapping->map[n]),
                       numComps_size);
            }
        }

    /*********************** FINLEY_REDUCED_NODES ***************************/
    } else if (in_data_type == FINLEY_REDUCED_NODES) {
        requireWrite(out);
        if (out_data_type == FINLEY_NODES) {
            Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: cannot copy from reduced nodes to nodes.");
        } else if (out_data_type == FINLEY_REDUCED_NODES) {
#pragma omp parallel for
            for (int n=0; n<nodes->reducedNodesMapping->numNodes; n++) {
                memcpy(getSampleDataRWFast(out,n),getSampleDataROFast(in,n),numComps_size);
            }
       } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
            Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: cannot copy from reduced nodes to degrees of freedom.");
       } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            int nComps = Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
#pragma omp parallel for
            for (int n=0; n<nComps; n++) {
               const int k = nodes->reducedDegreesOfFreedomMapping->map[n];
               memcpy(getSampleDataRWFast(out,n),
                      getSampleDataROFast(in,nodes->reducedNodesMapping->target[k]),
                      numComps_size);
            }
        }

    /******************** FINLEY_DEGREES_OF_FREEDOM *********************/
    } else if (in_data_type == FINLEY_DEGREES_OF_FREEDOM) {
        requireWrite(out);
        if (out_data_type == FINLEY_NODES) {
            Paso_Coupler *coupler = Paso_Coupler_alloc(nodes->degreesOfFreedomConnector, numComps);
            if (Esys_noError()) {
                // It is not immediately clear whether coupler can be
                // trusted with constant data so I'll assume RW.
                // Also, it holds pointers so it might not be safe to use
                // on lazy data anyway?
                requireWrite(in);
                Paso_Coupler_startCollect(coupler, getDataRW(in));
                const double *recv_buffer=Paso_Coupler_finishCollect(coupler);
                const index_t upperBound=Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution);
#pragma omp parallel for
                for (int n=0; n<nodes->numNodes; n++) {
                    const int k=nodes->degreesOfFreedomMapping->target[n];
                    if (k < upperBound) {
                        memcpy(getSampleDataRWFast(out,n),
                               getSampleDataROFast(in,k),
                               numComps_size);
                    } else {
                        memcpy(getSampleDataRWFast(out,n),
                               &recv_buffer[(k-upperBound)*numComps],
                               numComps_size);
                    }
                }
            }
            Paso_Coupler_free(coupler);
        } else if  (out_data_type == FINLEY_REDUCED_NODES) {
            Paso_Coupler *coupler = Paso_Coupler_alloc(nodes->degreesOfFreedomConnector, numComps);
            if (Esys_noError()) {
                requireWrite(in); // See comment above about coupler and const
                Paso_Coupler_startCollect(coupler, getDataRW(in));
                const double *recv_buffer=Paso_Coupler_finishCollect(coupler);
                const index_t upperBound=Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution);
                requireWrite(out);

#pragma omp parallel for
                for (int n=0; n<nodes->reducedNodesMapping->numTargets; n++) {
                    const int l=nodes->reducedNodesMapping->map[n];
                    const int k=nodes->degreesOfFreedomMapping->target[l];
                    if (k < upperBound) {
                        memcpy(getSampleDataRWFast(out,n),
                               getSampleDataROFast(in,k),
                               numComps_size);
                    } else {
                        memcpy(getSampleDataRWFast(out,n),
                               &recv_buffer[(k-upperBound)*numComps],
                               numComps_size);
                    }
                }
            }
            Paso_Coupler_free(coupler);
        } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
            const int nComps = Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution);
            requireWrite(out);
#pragma omp parallel for
            for (int n=0; n<nComps; n++) {
                memcpy(getSampleDataRWFast(out,n), getSampleDataROFast(in,n),
                       numComps_size);
            }
        } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            const int nComps = Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
            requireWrite(out);
#pragma omp parallel for
            for (int n=0; n<nComps; n++) {
                const int k=nodes->reducedDegreesOfFreedomMapping->map[n];
                memcpy(getSampleDataRWFast(out,n),
                       getSampleDataROFast(in,nodes->degreesOfFreedomMapping->target[k]),
                       numComps_size);
            }
        }

    /**************** FINLEY_REDUCED_DEGREES_OF_FREEDOM *****************/
    } else if (in_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (out_data_type == FINLEY_NODES) {
            Finley_setError(TYPE_ERROR, "Finley_Assemble_CopyNodalData: cannot copy from reduced degrees of freedom to nodes.");
        } else if (out_data_type == FINLEY_REDUCED_NODES) {
            Paso_Coupler *coupler=Paso_Coupler_alloc(nodes->reducedDegreesOfFreedomConnector,numComps);
            if (Esys_noError()) {
                const index_t upperBound=Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
                requireWrite(in); // See comment about coupler and const
                Paso_Coupler_startCollect(coupler, getDataRW(in));
                const double *recv_buffer=Paso_Coupler_finishCollect(coupler);
                requireWrite(out);
#pragma omp parallel for
                for (int n=0; n<nodes->reducedNodesMapping->numTargets; n++) {
                    const int l=nodes->reducedNodesMapping->map[n];
                    const int k=nodes->reducedDegreesOfFreedomMapping->target[l];
                    if (k < upperBound) {
                        memcpy(getSampleDataRWFast(out,n),
                               getSampleDataROFast(in,k),
                               numComps_size);
                    } else {
                        memcpy(getSampleDataRWFast(out,n),
                               &recv_buffer[(k-upperBound)*numComps],
                               numComps_size);
                    }
                }
            }
            Paso_Coupler_free(coupler);
        } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            const int nComps = Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
            requireWrite(out);
#pragma omp parallel for
            for (int n=0; n<nComps; n++) {
                memcpy(getSampleDataRWFast(out,n), getSampleDataROFast(in,n), numComps_size);
            }
        } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM ) {
            Finley_setError(TYPE_ERROR, "Finley_Assemble_CopyNodalData: cannot copy from reduced degrees of freedom to degrees of freedom.");
        }
    } // in_data_type
}

