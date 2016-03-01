
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

  Assemblage of Jacobians: calculate the gradient of nodal data at quadrature
  points

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

// Unless the loops in here get complicated again this file should be compiled
// with loop unrolling

namespace dudley {

void Assemble_gradient(Dudley_NodeFile* nodes, Dudley_ElementFile* elements,
                              escript::Data* grad_data, const escript::Data* data)
{
    if (!nodes || !elements)
        return;

    const int numComps = data->getDataPointSize();
    const int NN = elements->numNodes;
    const bool reducedIntegrationOrder = Assemble_reducedIntegrationOrder(grad_data);
    const int data_type = data->getFunctionSpace().getTypeCode();

    dim_t numNodes = 0;
    if (data_type == DUDLEY_NODES) {
        numNodes = nodes->nodesMapping->numTargets;
    } else if (data_type == DUDLEY_REDUCED_NODES) {
        numNodes = nodes->reducedNodesMapping->numTargets;
    } else if (data_type == DUDLEY_DEGREES_OF_FREEDOM) {
        if (elements->MPIInfo->size > 1) {
            throw DudleyException("Assemble_gradient: for more than one "
                "processor DEGREES_OF_FREEDOM data are not accepted as input.");
        }
        numNodes = nodes->degreesOfFreedomMapping->numTargets;
    } else if (data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (elements->MPIInfo->size > 1) {
            throw DudleyException("Assemble_gradient: for more than one "
                "processor REDUCED_DEGREES_OF_FREEDOM data are not accepted as input.");
        }
        numNodes = nodes->reducedDegreesOfFreedomMapping->numTargets;
    } else {
        throw DudleyException("Assemble_gradient: Cannot calculate gradient "
               "of data because of unsuitable input data representation.");
    }

    Dudley_ElementFile_Jacobians* jac = Dudley_ElementFile_borrowJacobians(
            elements, nodes, reducedIntegrationOrder);
    const int numDim = jac->numDim;
    const int numShapesTotal = jac->numShapes;
    const int numQuad = jac->numQuad;
    const size_t localGradSize = sizeof(double) * numDim * numQuad * numComps;

    // check the dimensions of data
    if (!grad_data->numSamplesEqual(numQuad, elements->numElements)) {
        throw DudleyException("Assemble_gradient: illegal number of samples in gradient Data object");
    } else if (!data->numSamplesEqual(1, numNodes)) {
        throw DudleyException("Assemble_gradient: illegal number of samples of input Data object");
    } else if (numDim * numComps != grad_data->getDataPointSize()) {
        throw DudleyException("Assemble_gradient: illegal number of components in gradient data object.");
    } else if (!grad_data->actsExpanded()) {
        throw DudleyException("Assemble_gradient: expanded Data object is expected for output data.");
    }

    grad_data->requireWrite();
#pragma omp parallel
    {
        if (data_type == DUDLEY_NODES) {
            if (numDim == 1) {
                const dim_t numShapes = 2;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    double* grad_data_e = grad_data->getSampleDataRW(e);
                    memset(grad_data_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const double* data_array = data->getSampleDataRO(n);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                grad_data_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            } else if (numDim == 2) {
                const dim_t numShapes = 3;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    double* grad_data_e = grad_data->getSampleDataRW(e);
                    memset(grad_data_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const double* data_array = data->getSampleDataRO(n);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                grad_data_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                grad_data_e[INDEX4(l, 1, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 1, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            } else if (numDim == 3) {
                const dim_t numShapes = 4;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    double* grad_data_e = grad_data->getSampleDataRW(e);
                    memset(grad_data_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const double* data_array = data->getSampleDataRO(n);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                grad_data_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                grad_data_e[INDEX4(l, 1, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 1, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                grad_data_e[INDEX4(l, 2, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 2, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            }
        } else if (data_type == DUDLEY_REDUCED_NODES) {
            if (numDim == 1) {
                const dim_t numShapes = 2;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    double* grad_data_e = grad_data->getSampleDataRW(e);
                    memset(grad_data_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const double* data_array = data->getSampleDataRO(nodes->reducedNodesMapping->target[n]);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                grad_data_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            } else if (numDim == 2) {
                const dim_t numShapes = 3;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    double* grad_data_e = grad_data->getSampleDataRW(e);
                    memset(grad_data_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const double* data_array = data->getSampleDataRO(nodes->reducedNodesMapping->target[n]);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                grad_data_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                grad_data_e[INDEX4(l, 1, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 1, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            } else if (numDim == 3) {
                const dim_t numShapes = 4;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    double* grad_data_e = grad_data->getSampleDataRW(e);
                    memset(grad_data_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const double* data_array = data->getSampleDataRO(nodes->reducedNodesMapping->target[n]);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                grad_data_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                grad_data_e[INDEX4(l, 1, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 1, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                grad_data_e[INDEX4(l, 2, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 2, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            }
        } else if (data_type == DUDLEY_DEGREES_OF_FREEDOM) {
            if (numDim == 1)
            {
                const dim_t numShapes = 2;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    double* grad_data_e = grad_data->getSampleDataRW(e);
                    memset(grad_data_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const double* data_array = data->getSampleDataRO(nodes->degreesOfFreedomMapping->target[n]);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                grad_data_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            } else if (numDim == 2) {
                const dim_t numShapes = 3;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    double* grad_data_e = grad_data->getSampleDataRW(e);
                    memset(grad_data_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const double* data_array = data->getSampleDataRO(nodes->degreesOfFreedomMapping->target[n]);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                grad_data_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                grad_data_e[INDEX4(l, 1, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 1, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            } else if (numDim == 3) {
                const dim_t numShapes = 4;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    double* grad_data_e = grad_data->getSampleDataRW(e);
                    memset(grad_data_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const double* data_array = data->getSampleDataRO(nodes->degreesOfFreedomMapping->target[n]);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                grad_data_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                grad_data_e[INDEX4(l, 1, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 1, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                grad_data_e[INDEX4(l, 2, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 2, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            }
        } else if (data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM) {
            if (numDim == 1) {
                const dim_t numShapes = 2;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    double* grad_data_e = grad_data->getSampleDataRW(e);
                    memset(grad_data_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const double* data_array = data->getSampleDataRO(nodes->reducedDegreesOfFreedomMapping->target[n]);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                grad_data_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            } else if (numDim == 2) {
                const dim_t numShapes = 3;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    double* grad_data_e = grad_data->getSampleDataRW(e);
                    memset(grad_data_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const double* data_array = data->getSampleDataRO(nodes->reducedDegreesOfFreedomMapping->target[n]);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                grad_data_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                grad_data_e[INDEX4(l, 1, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 1, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            } else if (numDim == 3) {
                const dim_t numShapes = 4;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    double* grad_data_e = grad_data->getSampleDataRW(e);
                    memset(grad_data_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const double* data_array = data->getSampleDataRO(nodes->reducedDegreesOfFreedomMapping->target[n]);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                grad_data_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                grad_data_e[INDEX4(l, 1, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 1, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                grad_data_e[INDEX4(l, 2, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 2, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            }
        }
    } // end parallel region
}

} // namespace dudley

