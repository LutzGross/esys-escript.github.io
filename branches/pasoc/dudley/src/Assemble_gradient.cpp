
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

// Unless the loops in here get complicated again this file should be compiled
// with loop unrolling

namespace dudley {

template<typename Scalar>
void Assemble_gradient(const NodeFile* nodes, const ElementFile* elements,
                       escript::Data& out, const escript::Data& data)
{
    if (!nodes || !elements)
        return;

    const int numComps = data.getDataPointSize();
    const int NN = elements->numNodes;
    const bool reducedIntegrationOrder = hasReducedIntegrationOrder(out);
    const int data_type = data.getFunctionSpace().getTypeCode();

    dim_t numNodes = 0;
    if (data_type == DUDLEY_NODES) {
        numNodes = nodes->getNumNodes();
    } else if (data_type == DUDLEY_DEGREES_OF_FREEDOM) {
        if (elements->MPIInfo->size > 1) {
            throw DudleyException("Assemble_gradient: for more than one "
                "processor DEGREES_OF_FREEDOM data are not accepted as input.");
        }
        numNodes = nodes->getNumDegreesOfFreedom();
    } else {
        throw DudleyException("Assemble_gradient: Cannot calculate gradient "
               "of data because of unsuitable input data representation.");
    }

    ElementFile_Jacobians* jac = elements->borrowJacobians(nodes,
                                                     reducedIntegrationOrder);
    const int numDim = jac->numDim;
    const int numShapesTotal = jac->numShapes;
    const int numQuad = jac->numQuad;

    // check the dimensions of data
    if (!out.numSamplesEqual(numQuad, elements->numElements)) {
        throw DudleyException("Assemble_gradient: illegal number of samples in gradient Data object");
    } else if (!data.numSamplesEqual(1, numNodes)) {
        throw DudleyException("Assemble_gradient: illegal number of samples of input Data object");
    } else if (numDim * numComps != out.getDataPointSize()) {
        throw DudleyException("Assemble_gradient: illegal number of components in gradient data object.");
    } else if (!out.actsExpanded()) {
        throw DudleyException("Assemble_gradient: expanded Data object is expected for output data.");
    }

    const Scalar zero = static_cast<Scalar>(0);
    const size_t localGradSize = sizeof(Scalar) * numDim * numQuad * numComps;
    out.requireWrite();
#pragma omp parallel
    {
        if (data_type == DUDLEY_NODES) {
            if (numDim == 1) {
                const int numShapes = 2;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    Scalar* gradData_e = out.getSampleDataRW(e, zero);
                    memset(gradData_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const Scalar* data_array = data.getSampleDataRO(n, zero);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                gradData_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            } else if (numDim == 2) {
                const int numShapes = 3;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    Scalar* gradData_e = out.getSampleDataRW(e, zero);
                    memset(gradData_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const Scalar* data_array = data.getSampleDataRO(n, zero);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                gradData_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                gradData_e[INDEX4(l, 1, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 1, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            } else if (numDim == 3) {
                const int numShapes = 4;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    Scalar* gradData_e = out.getSampleDataRW(e, zero);
                    memset(gradData_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const Scalar* data_array = data.getSampleDataRO(n, zero);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                gradData_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                gradData_e[INDEX4(l, 1, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 1, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                gradData_e[INDEX4(l, 2, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 2, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            }
        } else if (data_type == DUDLEY_DEGREES_OF_FREEDOM) {
            const index_t* target = nodes->borrowTargetDegreesOfFreedom();
            if (numDim == 1) {
                const int numShapes = 2;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    Scalar* gradData_e = out.getSampleDataRW(e, zero);
                    memset(gradData_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const Scalar* data_array = data.getSampleDataRO(target[n], zero);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                gradData_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            } else if (numDim == 2) {
                const int numShapes = 3;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    Scalar* gradData_e = out.getSampleDataRW(e, zero);
                    memset(gradData_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const Scalar* data_array = data.getSampleDataRO(target[n], zero);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                gradData_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                gradData_e[INDEX4(l, 1, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 1, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                            }
                        }
                    }
                }
            } else if (numDim == 3) {
                const int numShapes = 4;
#pragma omp for
                for (index_t e = 0; e < elements->numElements; e++) {
                    Scalar* gradData_e = out.getSampleDataRW(e, zero);
                    memset(gradData_e, 0, localGradSize);
                    for (int s = 0; s < numShapes; s++) {
                        const index_t n = elements->Nodes[INDEX2(s, e, NN)];
                        const Scalar* data_array = data.getSampleDataRO(target[n], zero);
                        for (int q = 0; q < numQuad; q++) {
#pragma ivdep
                            for (int l = 0; l < numComps; l++) {
                                gradData_e[INDEX4(l, 0, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 0, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                gradData_e[INDEX4(l, 1, q, 0, numComps, numDim, numQuad)] +=
                                    data_array[l] *
                                    jac->DSDX[INDEX5(s, 1, q, 0, e, numShapesTotal, numDim, numQuad, 1)];
                                gradData_e[INDEX4(l, 2, q, 0, numComps, numDim, numQuad)] +=
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

// instantiate our two supported versions
template void Assemble_gradient<escript::DataTypes::real_t>(
                       const NodeFile* nodes, const ElementFile* elements,
                       escript::Data& out, const escript::Data& data);
template void Assemble_gradient<escript::DataTypes::cplx_t>(
                       const NodeFile* nodes, const ElementFile* elements,
                       escript::Data& out, const escript::Data& data);


} // namespace dudley

