
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

#include "ElementFile.h"
#include "Assemble.h"

namespace finley {

ElementFile_Jacobians::ElementFile_Jacobians(const_ShapeFunction_ptr basis) :
    status(FINLEY_INITIAL_STATUS-1),
    numDim(0),
    BasisFunctions(basis),
    numQuadTotal(0),
    numElements(0),
    volume(NULL),
    DSDX(NULL)
{
}

ElementFile_Jacobians::~ElementFile_Jacobians()
{
    delete[] volume;
    delete[] DSDX;
}


ElementFile_Jacobians* ElementFile::borrowJacobians(const NodeFile* nodefile, 
        bool reducedShapefunction, bool reducedIntegrationOrder) const
{
    ElementFile_Jacobians *out = NULL;
  
    if (reducedShapefunction) {
       if (reducedIntegrationOrder) {
           out=jacobians_reducedS_reducedQ;
       } else {
           out=jacobians_reducedS;
       }
    } else {
       if (reducedIntegrationOrder) {
           out=jacobians_reducedQ;
       } else {
           out=jacobians;
       }
    }

    if (out->status < nodefile->status) {
        const_ShapeFunction_ptr basis(out->BasisFunctions);
        const_ShapeFunction_ptr shape(referenceElementSet->
                            borrowParametrization(reducedIntegrationOrder));
        const_ReferenceElement_ptr refElement(referenceElementSet->
                            borrowReferenceElement(reducedIntegrationOrder));

        out->numDim=nodefile->numDim;
        out->numQuadTotal=shape->numQuadNodes; 
        out->numSides=refElement->Type->numSides;
        out->numShapesTotal=basis->Type->numShapes * out->numSides; 
        out->numElements=numElements;
        const double *dBdv;

        if (reducedShapefunction) {
            out->numSub=1;
            out->node_selection=refElement->Type->linearNodes;
            out->offsets=refElement->LinearType->offsets;
            dBdv=&basis->dSdv[0];
        } else {
            out->numSub=refElement->Type->numSubElements;
            out->node_selection=refElement->Type->subElementNodes; 
            out->offsets=refElement->Type->offsets;
            dBdv=&refElement->DBasisFunctionDv[0];
        }
     
        if (out->numQuadTotal != out->numSub*basis->numQuadNodes) {
            throw FinleyException("ElementFile::borrowJacobians: Incorrect total number of quadrature points.");
        }
        if (refElement->numNodes > numNodes) {
            throw FinleyException("ElementFile::borrowJacobians: Too many nodes expected.");
        }

        if (out->volume==NULL)
            out->volume=new double[out->numElements*out->numQuadTotal];
        if (out->DSDX==NULL)
            out->DSDX=new double[out->numElements
                                  *out->numShapesTotal
                                  *out->numDim
                                  *out->numQuadTotal];

        /*========================== dim = 1 ============================== */
        if (out->numDim==1) {
            if (refElement->numLocalDim==0) {
                // nothing to be done
            } else if (refElement->numLocalDim==1) {
                if (out->numSides==1) {
                    Assemble_jacobians_1D(nodefile->Coordinates,
                            out->numQuadTotal, &shape->QuadWeights[0],
                            shape->Type->numShapes, numElements, numNodes,
                            Nodes, &shape->dSdv[0], basis->Type->numShapes,
                            dBdv, out->DSDX, out->volume, Id);
                } else {
                    throw FinleyException("ElementFile::borrowJacobians: only one-sided elements supported in 1D.");
                }
            } else {
                throw escript::ValueError("ElementFile::borrowJacobians: local dimension in a 1D domain has to be 0 or 1.");
            }
        /*========================== dim = 2 ============================== */
        } else if (out->numDim==2) {
            if (refElement->numLocalDim==0) {
                // nothing to be done
            } else if (refElement->numLocalDim==1) {
                if (out->BasisFunctions->Type->numDim==2) {
                    if (out->numSides==1) {
                        Assemble_jacobians_2D_M1D_E2D(
                                nodefile->Coordinates, out->numQuadTotal,
                                &shape->QuadWeights[0], shape->Type->numShapes,
                                numElements, numNodes, Nodes, &shape->dSdv[0],
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else if (out->numSides==2) {
                        Assemble_jacobians_2D_M1D_E2D_C(
                                nodefile->Coordinates, out->numQuadTotal,
                                &shape->QuadWeights[0], shape->Type->numShapes,
                                numElements, numNodes, Nodes, &shape->dSdv[0],
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else {
                        throw escript::ValueError("ElementFile::borrowJacobians: elements must be one- or two-sided.");
                    }
                } else if (out->BasisFunctions->Type->numDim==1) {
                    if (out->numSides==1) {
                        Assemble_jacobians_2D_M1D_E1D(
                                nodefile->Coordinates, out->numQuadTotal,
                                &shape->QuadWeights[0], shape->Type->numShapes,
                                numElements, numNodes, Nodes, &shape->dSdv[0],
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else if (out->numSides==2) {
                        Assemble_jacobians_2D_M1D_E1D_C(
                                nodefile->Coordinates, out->numQuadTotal,
                                &shape->QuadWeights[0], shape->Type->numShapes,
                                numElements, numNodes, Nodes, &shape->dSdv[0],
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else {
                        throw escript::ValueError("ElementFile::borrowJacobians: elements must be one- or two-sided.");
                    }
                } else {
                    throw escript::ValueError("ElementFile::borrowJacobians: element dimension for local dimension 1 in a 2D domain has to be 1 or 2.");
                }
            } else if (refElement->numLocalDim==2) {
                if (out->numSides==1) {
                    Assemble_jacobians_2D(nodefile->Coordinates,
                            out->numQuadTotal, &shape->QuadWeights[0],
                            shape->Type->numShapes, numElements, numNodes,
                            Nodes, &shape->dSdv[0], basis->Type->numShapes,
                            dBdv, out->DSDX, out->volume, Id);
                } else {
                    throw escript::ValueError("ElementFile::borrowJacobians: 2D volume supports one-sided elements only.");
                }
            } else {
                throw escript::ValueError("ElementFile::borrowJacobians: local dimension in a 2D domain has to be 1 or 2.");
            }
        /*========================== dim = 3 ============================== */
        } else if (out->numDim==3) {
            if (refElement->numLocalDim==0) {
                // nothing to be done
            } else if (refElement->numLocalDim==2) {
                if (out->BasisFunctions->Type->numDim==3) {
                    if (out->numSides==1) {
                        Assemble_jacobians_3D_M2D_E3D(
                                nodefile->Coordinates, out->numQuadTotal,
                                &shape->QuadWeights[0], shape->Type->numShapes,
                                numElements, numNodes, Nodes, &shape->dSdv[0],
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else if (out->numSides==2) {
                        Assemble_jacobians_3D_M2D_E3D_C(
                                nodefile->Coordinates, out->numQuadTotal,
                                &shape->QuadWeights[0], shape->Type->numShapes,
                                numElements, numNodes, Nodes, &shape->dSdv[0],
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else {
                        throw escript::ValueError("ElementFile::borrowJacobians: elements must be one- or two-sided.");
                    }
                } else if (out->BasisFunctions->Type->numDim==2) {
                    if (out->numSides==1) {
                        Assemble_jacobians_3D_M2D_E2D(
                                nodefile->Coordinates, out->numQuadTotal,
                                &shape->QuadWeights[0], shape->Type->numShapes,
                                numElements, numNodes, Nodes, &shape->dSdv[0],
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else if (out->numSides==2) {
                        Assemble_jacobians_3D_M2D_E2D_C(
                                nodefile->Coordinates, out->numQuadTotal,
                                &shape->QuadWeights[0], shape->Type->numShapes,
                                numElements, numNodes, Nodes, &shape->dSdv[0],
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else {
                        throw escript::ValueError("ElementFile::borrowJacobians: elements must be one- or two-sided.");
                    }
                } else {
                    throw escript::ValueError("ElementFile::borrowJacobians: element dimension for local dimension 2 in a 3D domain has to be 2 or 3.");
                }
            } else if (refElement->numLocalDim==3) {
                if (out->numSides==1) {
                    Assemble_jacobians_3D(nodefile->Coordinates,
                            out->numQuadTotal, &shape->QuadWeights[0],
                            shape->Type->numShapes, numElements, numNodes,
                            Nodes, &shape->dSdv[0], basis->Type->numShapes,
                            dBdv, out->DSDX, out->volume, Id);
                } else {
                    throw escript::ValueError("ElementFile::borrowJacobians: 3D volume supports one sided elements only..");
                }
            } else {
                throw escript::ValueError("ElementFile::borrowJacobians: local dimension in a 3D domain has to be 2 or 3.");
            }
        } else {
            throw escript::ValueError("ElementFile::borrowJacobians: number of spatial dimensions has to be 1, 2 or 3.");
        }

        out->status = nodefile->status;
    }
    return out;
}

} // namespace

