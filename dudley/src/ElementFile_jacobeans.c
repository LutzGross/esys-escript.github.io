
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

Dudley_ElementFile_Jacobeans* Dudley_ElementFile_Jacobeans_alloc(Dudley_ShapeFunction* BasisFunctions)
{
  Dudley_ElementFile_Jacobeans* out=MEMALLOC(1,Dudley_ElementFile_Jacobeans);
  if (Dudley_checkPtr(out)) {
     return NULL;
  } else {
     out->status=DUDLEY_INITIAL_STATUS-1;
     out->BasisFunctions=Dudley_ShapeFunction_reference(BasisFunctions);
	 out->numDim=0;
	 out->numQuadTotal=0;
	 out->numElements=0;
     out->volume=NULL;
     out->DSDX=NULL;
     return out;
  }
}

/**************************************************************/

void Dudley_ElementFile_Jacobeans_dealloc(Dudley_ElementFile_Jacobeans* in)
{
  if (in!=NULL) {
	Dudley_ShapeFunction_dealloc(in->BasisFunctions);  
    MEMFREE(in->volume);
    MEMFREE(in->DSDX);
    MEMFREE(in);
  }
}

/**************************************************************/


Dudley_ElementFile_Jacobeans* Dudley_ElementFile_borrowJacobeans(Dudley_ElementFile* self, Dudley_NodeFile* nodes, 
                                                                 bool_t reducedShapefunction, bool_t reducedIntegrationOrder) {
  Dudley_ElementFile_Jacobeans *out = NULL;
  Dudley_ShapeFunction *shape=NULL, *basis;
  Dudley_ReferenceElement*  refElement=NULL;
  double *dBdv;
  
  dim_t numNodes=self->numNodes;
  
  if (reducedIntegrationOrder)
  {
	out=self->jacobeans_reducedQ;
  }
  else
  {
	out=self->jacobeans;
  }
  if (out->status < nodes->status)
  {
     basis=out->BasisFunctions;
     shape=Dudley_ReferenceElementSet_borrowParametrization(self->referenceElementSet, reducedIntegrationOrder);
     refElement= Dudley_ReferenceElementSet_borrowReferenceElement(self->referenceElementSet, reducedIntegrationOrder);

     out->numDim=nodes->numDim;
     out->numQuadTotal=shape->numQuadNodes; 
     out->numShapesTotal=basis->Type->numShapes; 
     out->numElements=self->numElements;
     
     if (reducedShapefunction) {
        dBdv=basis->dSdv;
     } else {
        dBdv=refElement->DBasisFunctionDv;
     }
     
     if (out->numQuadTotal != basis->numQuadNodes) {
        Dudley_setError(SYSTEM_ERROR,"Dudley_ElementFile_borrowJacobeans: Incorrect total number of quadrature points.");
        return NULL;
     }
     if (refElement->numNodes> numNodes) {
           Dudley_setError(SYSTEM_ERROR,"Dudley_ElementFile_borrowJacobeans: Too many nodes expected.");
           return NULL;
     }

     if (out->volume==NULL) out->volume=MEMALLOC((out->numElements)*(out->numQuadTotal),double);
     if (out->DSDX==NULL) out->DSDX=MEMALLOC((out->numElements)
                                            *(out->numShapesTotal)
                                            *(out->numDim)
                                            *(out->numQuadTotal),double);
     if (! (Dudley_checkPtr(out->volume) || Dudley_checkPtr(out->DSDX)) ) {
          /*========================== dim = 1 ============================================== */
          if (out->numDim==1) {
	     Dudley_setError(SYSTEM_ERROR, "Dudley does not support 1D domains.");
          /*========================== dim = 2 ============================================== */
          } else if (out->numDim==2) {
             if (refElement->numLocalDim==0) {
				 Dudley_setError(SYSTEM_ERROR,"Dudley_ElementFile_borrowJacobeans: 2D does not support local dimension 0.");
             } else if (refElement->numLocalDim==1) {
                  if (out->BasisFunctions->Type->numDim==2) {

                        Assemble_jacobeans_2D_M1D_E2D(nodes->Coordinates,out->numQuadTotal,shape->QuadWeights,
                                                      shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                                      shape->dSdv,basis->Type->numShapes,dBdv,
                                                      out->DSDX,out->volume,self->Id);
                  }  else if (out->BasisFunctions->Type->numDim==1) {

                        Assemble_jacobeans_2D_M1D_E1D(nodes->Coordinates, out->numQuadTotal, self->numElements, numNodes,self->Nodes,
                                                      out->DSDX,out->volume,self->Id);
                  } else {
                    Dudley_setError(SYSTEM_ERROR,"Dudley_ElementFile_borrowJacobeans: element dimension for local dimenion 1 in a 2D domain has to be 1 or 2.");
                  }
             } else if (refElement->numLocalDim==2) {
                     Assemble_jacobeans_2D(nodes->Coordinates,out->numQuadTotal, self->numElements,numNodes,self->Nodes,
                                           out->DSDX, out->volume, self->Id);
             } else {
               Dudley_setError(SYSTEM_ERROR,"Dudley_ElementFile_borrowJacobeans: local dimenion in a 2D domain has to be  1 or 2.");
             }
          /*========================== dim = 3 ============================================== */
          } else if (out->numDim==3) {
             if (refElement->numLocalDim==0) {
		  Dudley_setError(SYSTEM_ERROR,"Dudley_ElementFile_borrowJacobeans: 3D does not support local dimension 0.");
             } else if (refElement->numLocalDim==2) {
                  if (out->BasisFunctions->Type->numDim==3) {
                        Assemble_jacobeans_3D_M2D_E3D(nodes->Coordinates,out->numQuadTotal,shape->QuadWeights,
                                                      shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                                      shape->dSdv,basis->Type->numShapes,dBdv,
                                                      out->DSDX,out->volume,self->Id);
                  }  else if (out->BasisFunctions->Type->numDim==2) {
                        Assemble_jacobeans_3D_M2D_E2D(nodes->Coordinates,out->numQuadTotal,shape->QuadWeights,
                                                      shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                                      shape->dSdv,basis->Type->numShapes,dBdv,
                                                      out->DSDX,out->volume,self->Id);
                  } else {
                    Dudley_setError(SYSTEM_ERROR,"Dudley_ElementFile_borrowJacobeans: element dimension for local dimenion 2 in a 3D domain has to be 3 or 2.");
                  }
             } else if (refElement->numLocalDim==3) {
                     Assemble_jacobeans_3D(nodes->Coordinates,out->numQuadTotal,shape->QuadWeights,
                                           shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                           shape->dSdv,basis->Type->numShapes,dBdv,
                                           out->DSDX,out->volume,self->Id);
             } else {
               Dudley_setError(SYSTEM_ERROR,"Dudley_ElementFile_borrowJacobeans: local dimenion in a 3D domain has to be 2 or 3.");
             }
          } else {
            Dudley_setError(SYSTEM_ERROR,"Dudley_ElementFile_borrowJacobeans: spatial dimension has to be 1, 2 or 3.");
          }
     }
     if (Dudley_noError()) {
         out->status = nodes->status;
     } else {
         out=NULL;
     }

   }

   return out;
}
