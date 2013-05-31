
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

#include "ElementFile.h"
#include "Assemble.h"

namespace finley {

ElementFile_Jacobians::ElementFile_Jacobians(Finley_ShapeFunction* basis)
{
    status=FINLEY_INITIAL_STATUS-1;
    BasisFunctions=Finley_ShapeFunction_reference(basis);
    numDim=0;
    numQuadTotal=0;
    numElements=0;
    volume=NULL;
    DSDX=NULL;
}

ElementFile_Jacobians::~ElementFile_Jacobians()
{
    Finley_ShapeFunction_dealloc(BasisFunctions);  
    delete[] volume;
    delete[] DSDX;
}


ElementFile_Jacobians* ElementFile::borrowJacobians(NodeFile* nodefile, 
        bool_t reducedShapefunction, bool_t reducedIntegrationOrder)
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
        Finley_ShapeFunction *basis = out->BasisFunctions;
        Finley_ShapeFunction *shape =
            Finley_ReferenceElementSet_borrowParametrization(referenceElementSet, reducedIntegrationOrder);
        Finley_ReferenceElement *refElement =
            Finley_ReferenceElementSet_borrowReferenceElement(referenceElementSet, reducedIntegrationOrder);

        out->numDim=nodefile->numDim;
        out->numQuadTotal=shape->numQuadNodes; 
        out->numSides=refElement->Type->numSides;
        out->numShapesTotal=basis->Type->numShapes * out->numSides; 
        out->numElements=numElements;
        double *dBdv;

        if (reducedShapefunction) {
            out->numSub=1;
            out->node_selection=refElement->Type->linearNodes;
            out->offsets=refElement->LinearType->offsets;
            dBdv=basis->dSdv;
        } else {
            out->numSub=refElement->Type->numSubElements;
            out->node_selection=refElement->Type->subElementNodes; 
            out->offsets=refElement->Type->offsets;
            dBdv=refElement->DBasisFunctionDv;
        }
     
        if (out->numQuadTotal != out->numSub*basis->numQuadNodes) {
            Finley_setError(SYSTEM_ERROR, "ElementFile::borrowJacobians: Incorrect total number of quadrature points.");
            return NULL;
        }
        if (refElement->numNodes > numNodes) {
            Finley_setError(SYSTEM_ERROR, "ElementFile::borrowJacobians: Too many nodes expected.");
            return NULL;
        }

        if (out->volume==NULL)
            out->volume=new double[(out->numElements)*(out->numQuadTotal)];
        if (out->DSDX==NULL)
            out->DSDX=new double[(out->numElements)
                                  *(out->numShapesTotal)
                                  *(out->numDim)
                                  *(out->numQuadTotal)];

        /*========================== dim = 1 ============================== */
        if (out->numDim==1) {
            if (refElement->numLocalDim==0) {
                // nothing to be done
            } else if (refElement->numLocalDim==1) {
                if (out->numSides==1) {
                    Finley_Assemble_jacobians_1D(nodefile->Coordinates,
                            out->numQuadTotal, shape->QuadWeights,
                            shape->Type->numShapes, numElements, numNodes,
                            Nodes, shape->dSdv, basis->Type->numShapes, dBdv,
                            out->DSDX, out->volume, Id);
                } else {
                    Finley_setError(SYSTEM_ERROR, "ElementFile::borrowJacobians: only one-sided elements supported in 1D.");
                }
            } else {
                Finley_setError(SYSTEM_ERROR, "ElementFile::borrowJacobians: local dimension in a 1D domain has to be 0 or 1.");
            }
        /*========================== dim = 2 ============================== */
        } else if (out->numDim==2) {
            if (refElement->numLocalDim==0) {
                // nothing to be done
            } else if (refElement->numLocalDim==1) {
                if (out->BasisFunctions->Type->numDim==2) {
                    if (out->numSides==1) {
                        Finley_Assemble_jacobians_2D_M1D_E2D(
                                nodefile->Coordinates, out->numQuadTotal,
                                shape->QuadWeights, shape->Type->numShapes,
                                numElements, numNodes, Nodes, shape->dSdv,
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else if (out->numSides==2) {
                        Finley_Assemble_jacobians_2D_M1D_E2D_C(
                                nodefile->Coordinates, out->numQuadTotal,
                                shape->QuadWeights, shape->Type->numShapes,
                                numElements, numNodes, Nodes, shape->dSdv,
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else {
                        Finley_setError(SYSTEM_ERROR, "ElementFile::borrowJacobians: elements must be one- or two-sided.");
                    }
                } else if (out->BasisFunctions->Type->numDim==1) {
                    if (out->numSides==1) {
                        Finley_Assemble_jacobians_2D_M1D_E1D(
                                nodefile->Coordinates, out->numQuadTotal,
                                shape->QuadWeights, shape->Type->numShapes,
                                numElements, numNodes, Nodes, shape->dSdv,
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else if (out->numSides==2) {
                        Finley_Assemble_jacobians_2D_M1D_E1D_C(
                                nodefile->Coordinates, out->numQuadTotal,
                                shape->QuadWeights, shape->Type->numShapes,
                                numElements, numNodes, Nodes, shape->dSdv,
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else {
                        Finley_setError(SYSTEM_ERROR, "ElementFile::borrowJacobians: elements must be one- or two-sided.");
                    }
                } else {
                    Finley_setError(SYSTEM_ERROR, "ElementFile::borrowJacobians: element dimension for local dimension 1 in a 2D domain has to be 1 or 2.");
                }
            } else if (refElement->numLocalDim==2) {
                if (out->numSides==1) {
                    Finley_Assemble_jacobians_2D(nodefile->Coordinates,
                            out->numQuadTotal, shape->QuadWeights,
                            shape->Type->numShapes, numElements, numNodes,
                            Nodes, shape->dSdv, basis->Type->numShapes, dBdv,
                            out->DSDX, out->volume, Id);
                } else {
                    Finley_setError(SYSTEM_ERROR, "ElementFile::borrowJacobians: 2D volume supports one-sided elements only.");
                }
            } else {
                Finley_setError(SYSTEM_ERROR, "ElementFile::borrowJacobians: local dimension in a 2D domain has to be 1 or 2.");
            }
        /*========================== dim = 3 ============================== */
        } else if (out->numDim==3) {
            if (refElement->numLocalDim==0) {
                // nothing to be done
            } else if (refElement->numLocalDim==2) {
                if (out->BasisFunctions->Type->numDim==3) {
                    if (out->numSides==1) {
                        Finley_Assemble_jacobians_3D_M2D_E3D(
                                nodefile->Coordinates, out->numQuadTotal,
                                shape->QuadWeights, shape->Type->numShapes,
                                numElements, numNodes, Nodes, shape->dSdv,
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else if (out->numSides==2) {
                        Finley_Assemble_jacobians_3D_M2D_E3D_C(
                                nodefile->Coordinates, out->numQuadTotal,
                                shape->QuadWeights, shape->Type->numShapes,
                                numElements, numNodes, Nodes, shape->dSdv,
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else {
                        Finley_setError(SYSTEM_ERROR, "ElementFile::borrowJacobians: elements must be one- or two-sided.");
                    }
                } else if (out->BasisFunctions->Type->numDim==2) {
                    if (out->numSides==1) {
                        Finley_Assemble_jacobians_3D_M2D_E2D(
                                nodefile->Coordinates, out->numQuadTotal,
                                shape->QuadWeights, shape->Type->numShapes,
                                numElements, numNodes, Nodes, shape->dSdv,
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else if (out->numSides==2) {
                        Finley_Assemble_jacobians_3D_M2D_E2D_C(
                                nodefile->Coordinates, out->numQuadTotal,
                                shape->QuadWeights, shape->Type->numShapes,
                                numElements, numNodes, Nodes, shape->dSdv,
                                basis->Type->numShapes, dBdv, out->DSDX,
                                out->volume, Id);
                    } else {
                        Finley_setError(SYSTEM_ERROR,"ElementFile::borrowJacobians: elements must be one- or two-sided.");
                    }
                } else {
                    Finley_setError(SYSTEM_ERROR, "ElementFile::borrowJacobians: element dimension for local dimension 2 in a 3D domain has to be 2 or 3.");
                }
            } else if (refElement->numLocalDim==3) {
                if (out->numSides==1) {
                    Finley_Assemble_jacobians_3D(nodefile->Coordinates,
                            out->numQuadTotal, shape->QuadWeights,
                            shape->Type->numShapes, numElements, numNodes,
                            Nodes, shape->dSdv, basis->Type->numShapes, dBdv,
                            out->DSDX, out->volume, Id);
                } else {
                    Finley_setError(SYSTEM_ERROR, "ElementFile::borrowJacobians: 3D volume supports one sided elements only..");
                }
            } else {
                Finley_setError(SYSTEM_ERROR, "ElementFile::borrowJacobians: local dimension in a 3D domain has to be 2 or 3.");
            }
        } else {
            Finley_setError(SYSTEM_ERROR, "ElementFile::borrowJacobians: number of spatial dimensions has to be 1, 2 or 3.");
        }

        if (Finley_noError()) {
            out->status = nodefile->status;
        } else {
            out=NULL;
        }
    }
    return out;
}

} // namespace

