
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

  Assemblage of jacobians: calculates the gradient of nodal data at
  quadrature points

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

void Finley_Assemble_gradient(Finley_NodeFile* nodes,
                              Finley_ElementFile* elements,
                              escriptDataC* grad_data,
                              escriptDataC* data)
{
    Finley_resetError();
    if (!nodes || !elements)
        return;

    const dim_t numComps=getDataPointSize(data);
    const dim_t NN=elements->numNodes;
    const bool_t reducedIntegrationOrder=Finley_Assemble_reducedIntegrationOrder(grad_data);

    const type_t data_type=getFunctionSpaceType(data);
    bool_t reducedShapefunction = FALSE;
    dim_t numNodes = 0;
    if (data_type == FINLEY_NODES) {
        numNodes = nodes->nodesMapping->numTargets;
    } else if (data_type==FINLEY_REDUCED_NODES) { 
        reducedShapefunction = TRUE;
        numNodes = nodes->reducedNodesMapping->numTargets;
    } else if (data_type==FINLEY_DEGREES_OF_FREEDOM) {
        if (elements->MPIInfo->size > 1) {
            Finley_setError(TYPE_ERROR, "Finley_Assemble_gradient: for more than one processor DEGREES_OF_FREEDOM data are not accepted as input.");
            return;
        }
        numNodes = nodes->degreesOfFreedomMapping->numTargets;
    } else if (data_type==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (elements->MPIInfo->size > 1) {
            Finley_setError(TYPE_ERROR, "Finley_Assemble_gradient: for more than one processor REDUCED_DEGREES_OF_FREEDOM data are not accepted as input.");
            return;
        }
        reducedShapefunction = TRUE;
        numNodes = nodes->reducedDegreesOfFreedomMapping->numTargets;
    } else {
        Finley_setError(TYPE_ERROR, "Finley_Assemble_gradient: Cannot calculate gradient of data because of unsuitable input data representation.");
        return;
    }

    Finley_ElementFile_Jacobeans *jac = Finley_ElementFile_borrowJacobeans(
            elements, nodes, reducedShapefunction, reducedIntegrationOrder);
    Finley_ReferenceElement *refElement =
        Finley_ReferenceElementSet_borrowReferenceElement(
                elements->referenceElementSet, reducedIntegrationOrder);
    const dim_t numDim=jac->numDim;
    const dim_t numShapes=jac->BasisFunctions->Type->numShapes;
    const dim_t numShapesTotal=jac->numShapesTotal;
    const dim_t numSub=jac->numSub;
    const dim_t numQuad=jac->numQuadTotal/numSub;
    dim_t numShapesTotal2=0;
    index_t s_offset=0, *nodes_selector=NULL;
  
    if (Finley_noError()) {
        const type_t grad_data_type=getFunctionSpaceType(grad_data);
        if (grad_data_type==FINLEY_CONTACT_ELEMENTS_2 || grad_data_type==FINLEY_REDUCED_CONTACT_ELEMENTS_2)  {
            s_offset=jac->offsets[1];
            s_offset=jac->offsets[1];
        } else {
            s_offset=jac->offsets[0];
            s_offset=jac->offsets[0];
        }
        if (data_type==FINLEY_REDUCED_NODES || data_type==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            nodes_selector=refElement->Type->linearNodes;
            numShapesTotal2=refElement->LinearBasisFunctions->Type->numShapes * refElement->Type->numSides;
        } else { 
            nodes_selector=refElement->Type->subElementNodes;
            numShapesTotal2=refElement->BasisFunctions->Type->numShapes * refElement->Type->numSides;
        }

        // check the dimensions of data
        if (!numSamplesEqual(grad_data, numQuad*numSub, elements->numElements)) {
            Finley_setError(TYPE_ERROR, "Finley_Assemble_gradient: illegal number of samples in gradient Data object");
        } else if (!numSamplesEqual(data, 1, numNodes)) {
            Finley_setError(TYPE_ERROR, "Finley_Assemble_gradient: illegal number of samples of input Data object");
        } else if (numDim*numComps != getDataPointSize(grad_data)) {
            Finley_setError(TYPE_ERROR, "Finley_Assemble_gradient: illegal number of components in gradient data object.");
        } else if (!isExpanded(grad_data)) {
            Finley_setError(TYPE_ERROR, "Finley_Assemble_gradient: expanded Data object is expected for output data.");
        } else if (!(s_offset+numShapes <= numShapesTotal)) {
            Finley_setError(SYSTEM_ERROR, "Finley_Assemble_gradient: nodes per element is inconsistent with number of jacobians.");
        }
    }

    // now we can start
    if (!Finley_noError())
        return;

    const size_t localGradSize=sizeof(double)*numDim*numQuad*numSub*numComps;
    requireWrite(grad_data);
#pragma omp parallel
    {
        if (data_type==FINLEY_NODES) {
            if (numDim==1) {
#pragma omp for schedule(static)
                for (dim_t e=0; e<elements->numElements; e++) {
                    double *grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e, 0, localGradSize);
                    for (dim_t isub=0; isub<numSub; isub++) {
                        for (dim_t s=0; s<numShapes; s++) {
                            const dim_t n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
                            const double *data_array=getSampleDataRO(data,n);
                            for (dim_t q=0; q<numQuad; q++) {
#pragma ivdep
                                for (dim_t l=0; l<numComps; l++) {
                                    grad_data_e[INDEX4(l,0,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                }
                            }
                        }
                    }
                }
            } else if (numDim==2) {
#pragma omp for schedule(static)
                for (dim_t e=0; e<elements->numElements; e++) {
                    double *grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e, 0, localGradSize);
                    for (dim_t isub=0; isub<numSub; isub++) {
                        for (dim_t s=0; s<numShapes; s++) {
                            const dim_t n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
                            const double *data_array=getSampleDataRO(data,n);
                            for (dim_t q=0; q<numQuad; q++) {
#pragma ivdep
                                for (dim_t l=0; l<numComps; l++) {
                                    grad_data_e[INDEX4(l,0,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                    grad_data_e[INDEX4(l,1,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                }
                            }
                        }
                    }
                }
            } else if (numDim==3) {
#pragma omp for schedule(static)
                for (dim_t e=0; e<elements->numElements; e++) {
                    double *grad_data_e=getSampleDataRW(grad_data,e); 
                    memset(grad_data_e,0, localGradSize);
                    for (dim_t isub=0; isub<numSub; isub++) {
                        for (dim_t s=0; s<numShapes; s++) {
                            const dim_t n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
                            const double *data_array=getSampleDataRO(data,n);
                            for (dim_t q=0; q<numQuad; q++) {
#pragma ivdep
                                for (dim_t l=0; l<numComps; l++) {
                                    grad_data_e[INDEX4(l,0,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                    grad_data_e[INDEX4(l,1,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                    grad_data_e[INDEX4(l,2,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,2,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                }
                            }
                        }
                    }
                }
            }
        } else if (data_type==FINLEY_REDUCED_NODES) {
            if (numDim==1) {
#pragma omp for schedule(static)
                for (dim_t e=0; e<elements->numElements; e++) {
                    double *grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e, 0, localGradSize);
                    for (dim_t isub=0; isub<numSub; isub++) {
                        for (dim_t s=0;s<numShapes;s++) {
                            const dim_t n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
                            const double *data_array=getSampleDataRO(data,nodes->reducedNodesMapping->target[n]);            
                            for (dim_t q=0; q<numQuad; q++) {
#pragma ivdep
                                for (dim_t l=0; l<numComps; l++) {                              
                                    grad_data_e[INDEX4(l,0,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                }
                            }
                        }
                    }
                }
            } else if (numDim==2) {
#pragma omp for schedule(static)
                for (dim_t e=0; e<elements->numElements; e++) {
                    double *grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e, 0, localGradSize);
                    for (dim_t isub=0; isub<numSub; isub++) {
                        for (dim_t s=0; s<numShapes; s++) {
                            const dim_t n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
                            const double *data_array=getSampleDataRO(data,nodes->reducedNodesMapping->target[n]);
                            for (dim_t q=0; q<numQuad; q++) {
#pragma ivdep
                                for (dim_t l=0; l<numComps; l++) {
                                    grad_data_e[INDEX4(l,0,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                    grad_data_e[INDEX4(l,1,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                }
                            }
                        }
                    }
                }
            } else if (numDim==3) {
#pragma omp for schedule(static)
                for (dim_t e=0;e<elements->numElements;e++) {
                    double *grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e, 0, localGradSize);
                    for (dim_t isub=0; isub<numSub; isub++) {
                        for (dim_t s=0; s<numShapes; s++) {
                            const dim_t n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
                            const double *data_array=getSampleDataRO(data,nodes->reducedNodesMapping->target[n]);
                            for (dim_t q=0; q<numQuad; q++) {   
#pragma ivdep
                                for (dim_t l=0;l<numComps;l++) {
                                    grad_data_e[INDEX4(l,0,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                    grad_data_e[INDEX4(l,1,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                    grad_data_e[INDEX4(l,2,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,2,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                }
                            }
                        }
                    }
                }
            }
        } else if (data_type==FINLEY_DEGREES_OF_FREEDOM) {

            if (numDim==1) {
#pragma omp for schedule(static)
                for (dim_t e=0; e<elements->numElements; e++) {
                    double *grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e, 0, localGradSize);
                    for (dim_t isub=0; isub<numSub; isub++) {
                        for (dim_t s=0; s<numShapes; s++) {
                            const dim_t n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
                            const double *data_array=getSampleDataRO(data,nodes->degreesOfFreedomMapping->target[n]);
                            for (dim_t q=0; q<numQuad; q++) {
#pragma ivdep
                                for (dim_t l=0; l<numComps; l++) {
                                    grad_data_e[INDEX4(l,0,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                }
                            }
                        }
                    }
                }
            } else if (numDim==2) {
#pragma omp for schedule(static)
                for (dim_t e=0; e<elements->numElements; e++) {
                    double *grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e, 0, localGradSize);
                    for (dim_t isub=0; isub<numSub; isub++) {
                        for (dim_t s=0; s<numShapes; s++) {
                            const dim_t n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
                            const double *data_array=getSampleDataRO(data,nodes->degreesOfFreedomMapping->target[n]);
                            for (dim_t q=0; q<numQuad; q++) {
#pragma ivdep
                                for (dim_t l=0; l<numComps; l++) {
                                    grad_data_e[INDEX4(l,0,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                    grad_data_e[INDEX4(l,1,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                }
                            }
                        }
                    }
                }
            } else if (numDim==3) {
#pragma omp for schedule(static)
                for (dim_t e=0; e<elements->numElements; e++) {
                    double *grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e, 0, localGradSize);
                    for (dim_t isub=0; isub<numSub; isub++) {
                        for (dim_t s=0; s<numShapes; s++) {
                            const dim_t n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
                            const double *data_array=getSampleDataRO(data,nodes->degreesOfFreedomMapping->target[n]);
                            for (dim_t q=0; q<numQuad; q++) {
#pragma ivdep
                                for (dim_t l=0; l<numComps; l++) {
                                    grad_data_e[INDEX4(l,0,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                    grad_data_e[INDEX4(l,1,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                    grad_data_e[INDEX4(l,2,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,2,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                }
                            }
                        }
                    }
                }
            }
        } else if (data_type==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            if (numDim==1) {
#pragma omp for schedule(static)
                for (dim_t e=0; e<elements->numElements; e++) {
                    double *grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
                    for (dim_t isub=0; isub<numSub; isub++) {
                        for (dim_t s=0; s<numShapes; s++) {
                            const dim_t n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
                            const double *data_array=getSampleDataRO(data,nodes->reducedDegreesOfFreedomMapping->target[n]);
                            for (dim_t q=0; q<numQuad; q++) {
#pragma ivdep
                                for (dim_t l=0; l<numComps; l++) {
                                    grad_data_e[INDEX4(l,0,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                }
                            }
                        }
                    }
                }
            } else if (numDim==2) {
#pragma omp for schedule(static)
                for (dim_t e=0;e<elements->numElements;e++) {
                    double *grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e, 0, localGradSize);
                    for (dim_t isub=0; isub<numSub; isub++) {
                        for (dim_t s=0; s<numShapes; s++) {
                            const dim_t n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
                            const double *data_array=getSampleDataRO(data,nodes->reducedDegreesOfFreedomMapping->target[n]);
                            for (dim_t q=0; q<numQuad; q++) {
#pragma ivdep
                                for (dim_t l=0; l<numComps; l++) {
                                    grad_data_e[INDEX4(l,0,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                    grad_data_e[INDEX4(l,1,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                }
                            }
                        }
                    }
                }

            } else if (numDim==3) {
#pragma omp for schedule(static)
                for (dim_t e=0; e<elements->numElements; e++) {
                    double *grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
                    for (dim_t isub=0; isub<numSub; isub++) {
                        for (dim_t s=0; s<numShapes; s++) {
                            const dim_t n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
                            const double *data_array=getSampleDataRO(data,nodes->reducedDegreesOfFreedomMapping->target[n]);
                            for (dim_t q=0; q<numQuad; q++) {
#pragma ivdep
                                for (dim_t l=0; l<numComps; l++) {
                                    grad_data_e[INDEX4(l,0,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                    grad_data_e[INDEX4(l,1,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                    grad_data_e[INDEX4(l,2,q,isub,numComps,numDim,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,2,q,isub,e,numShapesTotal,numDim,numQuad,numSub)];
                                }
                            }
                        }
                    }
                }
            } // numDim
        } // data_type
    } // end parallel region
}
