
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "ElementFile.h"
#include "Assemble.h"
#ifdef _OPENMP
#include <omp.h>
#endif


/**************************************************************/

Finley_ElementFile_Jacobeans* Finley_ElementFile_Jacobeans_alloc(Finley_ShapeFunction* BasisFunctions)
{
  Finley_ElementFile_Jacobeans* out=MEMALLOC(1,Finley_ElementFile_Jacobeans);
  if (Finley_checkPtr(out)) {
     return NULL;
  } else {
     out->status=FINLEY_INITIAL_STATUS-1;
     out->BasisFunctions=Finley_ShapeFunction_reference(BasisFunctions);
	 out->numDim=0;
	 out->numQuadTotal=0;
	 out->numElements=0;
     out->volume=NULL;
     out->DSDX=NULL;
     return out;
  }
}

/**************************************************************/

void Finley_ElementFile_Jacobeans_dealloc(Finley_ElementFile_Jacobeans* in)
{
  if (in!=NULL) {
	Finley_ShapeFunction_dealloc(in->BasisFunctions);  
    MEMFREE(in->volume);
    MEMFREE(in->DSDX);
    MEMFREE(in);
  }
}

/**************************************************************/


Finley_ElementFile_Jacobeans* Finley_ElementFile_borrowJacobeans(Finley_ElementFile* self, Finley_NodeFile* nodes, 
                                                                 bool_t reducedShapefunction, bool_t reducedIntegrationOrder) {
  Finley_ElementFile_Jacobeans *out = NULL;
  Finley_ShapeFunction *shape=NULL, *basis;
  Finley_ReferenceElement*  refElement=NULL;
  double *dBdv;
  
  dim_t numNodes=self->numNodes;
  
  if (reducedShapefunction) {
       if (reducedIntegrationOrder) {
           out=self->jacobeans_reducedS_reducedQ;
       } else {
           out=self->jacobeans_reducedS;
       }
  } else {
       if (reducedIntegrationOrder) {
           out=self->jacobeans_reducedQ;
       } else {
           out=self->jacobeans;
       }
  }
  if (out->status < nodes->status) {
    
     basis=out->BasisFunctions;
     shape=Finley_ReferenceElementSet_borrowParametrization(self->referenceElementSet, reducedIntegrationOrder);
     refElement= Finley_ReferenceElementSet_borrowReferenceElement(self->referenceElementSet, reducedIntegrationOrder);

     out->numDim=nodes->numDim;
     out->numQuadTotal=shape->numQuadNodes; 
     out->numSides=refElement->Type->numSides;
     out->numShapesTotal=basis->Type->numShapes * out->numSides; 
     out->numElements=self->numElements;
     
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
     
     if (out->numQuadTotal != out->numSub * basis->numQuadNodes) {
        Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: Incorrect total number of quadrature points.");
        return NULL;
     }
     if (refElement->numNodes> numNodes) {
           Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: Too many nodes expected.");
           return NULL;
     }

     if (out->volume==NULL) out->volume=MEMALLOC((out->numElements)*(out->numQuadTotal),double);
     if (out->DSDX==NULL) out->DSDX=MEMALLOC((out->numElements)
                                            *(out->numShapesTotal)
                                            *(out->numDim)
                                            *(out->numQuadTotal),double);
     if (! (Finley_checkPtr(out->volume) || Finley_checkPtr(out->DSDX)) ) {
          /*========================== dim = 1 ============================================== */
          if (out->numDim==1) {
             if (refElement->numLocalDim==0) {
		  	          /* all done */
             } else if (refElement->numLocalDim==1) {
                  if (out->numSides==2) {
				       Finley_Assemble_jacobeans_1D(nodes->Coordinates,out->numQuadTotal,shape->QuadWeights,
                                            shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                            shape->dSdv,basis->Type->numShapes,dBdv, out->DSDX,out->volume,self->Id);

                  } else {
                      Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: 1D supports one sided elements only.");
                  }
             } else {
                  Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: local dimenion in a 1D domain has to be 0 or 1.");
             }
          /*========================== dim = 2 ============================================== */
          } else if (out->numDim==2) {
             if (refElement->numLocalDim==0) {
	          /* all done */
             } else if (refElement->numLocalDim==1) {
                  if (out->BasisFunctions->Type->numDim==2) {
                     if (out->numSides==1) {
                        Finley_Assemble_jacobeans_2D_M1D_E2D(nodes->Coordinates,out->numQuadTotal,shape->QuadWeights,
                                                      shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                                      shape->dSdv,basis->Type->numShapes,dBdv,
                                                      out->DSDX,out->volume,self->Id);
                     } else if (out->numSides==2) {
                        Finley_Assemble_jacobeans_2D_M1D_E2D_C(nodes->Coordinates,out->numQuadTotal,shape->QuadWeights,
                                                        shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                                        shape->dSdv,basis->Type->numShapes,dBdv,
                                                        out->DSDX,out->volume,self->Id);
                     } else {
                          Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: 2D supports one or two sided elements only.");
                      }
                  }  else if (out->BasisFunctions->Type->numDim==1) {
                     if (out->numSides==1) {
                        Finley_Assemble_jacobeans_2D_M1D_E1D(nodes->Coordinates,out->numQuadTotal,shape->QuadWeights,
                                                      shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                                      shape->dSdv,basis->Type->numShapes,dBdv,
                                                      out->DSDX,out->volume,self->Id);
                     } else if (out->numSides==2) {
                        Finley_Assemble_jacobeans_2D_M1D_E1D_C(nodes->Coordinates,out->numQuadTotal,shape->QuadWeights,
                                                        shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                                        shape->dSdv,basis->Type->numShapes,dBdv,
                                                        out->DSDX,out->volume,self->Id);
                     } else {
                          Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: 2D supports one or two sided elements only.");
                      }
                  } else {
                    Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: element dimension for local dimenion 1 in a 2D domain has to be 1 or 2.");
                  }
             } else if (refElement->numLocalDim==2) {
                  if (out->numSides==1) {
                     Finley_Assemble_jacobeans_2D(nodes->Coordinates,out->numQuadTotal,shape->QuadWeights,
                                           shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                           shape->dSdv,basis->Type->numShapes,dBdv,
                                           out->DSDX,out->volume,self->Id);
                  } else {
                      Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: 2D volume supports one sided elements only.");
                  }

             } else {
               Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: local dimenion in a 2D domain has to be  1 or 2.");
             }
          /*========================== dim = 3 ============================================== */
          } else if (out->numDim==3) {
             if (refElement->numLocalDim==0) {
		  	          /* all done */
             } else if (refElement->numLocalDim==2) {
                  if (out->BasisFunctions->Type->numDim==3) {
                     if (out->numSides==1) {
                        Finley_Assemble_jacobeans_3D_M2D_E3D(nodes->Coordinates,out->numQuadTotal,shape->QuadWeights,
                                                      shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                                      shape->dSdv,basis->Type->numShapes,dBdv,
                                                      out->DSDX,out->volume,self->Id);
                     } else if (out->numSides==2) {
                        Finley_Assemble_jacobeans_3D_M2D_E3D_C(nodes->Coordinates,out->numQuadTotal,shape->QuadWeights,
                                                        shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                                        shape->dSdv,basis->Type->numShapes,dBdv,
                                                        out->DSDX,out->volume,self->Id);
                     } else {
                          Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: 3D supports one or two sided elements only.");
                      }
                  }  else if (out->BasisFunctions->Type->numDim==2) {
                     if (out->numSides==1) {
                        Finley_Assemble_jacobeans_3D_M2D_E2D(nodes->Coordinates,out->numQuadTotal,shape->QuadWeights,
                                                      shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                                      shape->dSdv,basis->Type->numShapes,dBdv,
                                                      out->DSDX,out->volume,self->Id);
                     } else if (out->numSides==2) {
                        Finley_Assemble_jacobeans_3D_M2D_E2D_C(nodes->Coordinates,out->numQuadTotal,shape->QuadWeights,
                                                        shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                                        shape->dSdv,basis->Type->numShapes,dBdv,
                                                        out->DSDX,out->volume,self->Id);
                     } else {
                          Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: 3D supports one or two sided elements only.");
                      }
                  } else {
                    Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: element dimension for local dimenion 2 in a 3D domain has to be 3 or 2.");
                  }
             } else if (refElement->numLocalDim==3) {
                  if (out->numSides==1) {
                     Finley_Assemble_jacobeans_3D(nodes->Coordinates,out->numQuadTotal,shape->QuadWeights,
                                           shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                           shape->dSdv,basis->Type->numShapes,dBdv,
                                           out->DSDX,out->volume,self->Id);
                  } else {
                      Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: 3D volume supports one sided elements only..");
                  }
             } else {
               Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: local dimenion in a 3D domain has to be 2 or 3.");
             }
          } else {
            Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: spatial dimension has to be 1, 2 or 3.");
          }
     }
     if (Finley_noError()) {
         out->status = nodes->status;
     } else {
         out=NULL;
     }

   }

   return out;
}
