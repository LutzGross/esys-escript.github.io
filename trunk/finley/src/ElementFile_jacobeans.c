/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/


/**************************************************************/

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "ElementFile.h"
#include "Assemble.h"
#ifdef _OPENMP
#include <omp.h>
#endif


/**************************************************************/

Finley_ElementFile_Jacobeans* Finley_ElementFile_Jacobeans_alloc(Finley_RefElement* ReferenceElement)
{
  Finley_ElementFile_Jacobeans* out=MEMALLOC(1,Finley_ElementFile_Jacobeans);
  if (Finley_checkPtr(out)) {
     return NULL;
  } else {
     out->status=FINLEY_INITIAL_STATUS-1;
     out->ReferenceElement=ReferenceElement;
     out->volume=NULL;
     out->DSDX=NULL;
     return out;
  }
}

/**************************************************************/

void Finley_ElementFile_Jacobeans_dealloc(Finley_ElementFile_Jacobeans* in)
{
  if (in!=NULL) {
    if (in->volume!=NULL) MEMFREE(in->volume);
    if (in->DSDX!=NULL) MEMFREE(in->DSDX);  
    MEMFREE(in);
  }
}

/**************************************************************/


/**************************************************************/

Finley_ElementFile_Jacobeans* Finley_ElementFile_borrowJacobeans(Finley_ElementFile* self, Finley_NodeFile* nodes, 
                                                                 bool_t reducedShapefunction, bool_t reducedIntegrationOrder) {
  Finley_ElementFile_Jacobeans *out = NULL;
  
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
     dim_t numNodes=self->ReferenceElement->Type->numShapes;
     Finley_RefElement *shape;
     if (reducedIntegrationOrder) {
           shape=self->ReferenceElement;
     } else {
           shape=self->ReferenceElement;
     }
     if (out->volume==NULL) out->volume=MEMALLOC((self->numElements)*(out->ReferenceElement->numQuadNodes),double);
     if (out->DSDX==NULL) out->DSDX=MEMALLOC((self->numElements)
                                            *(out->ReferenceElement->Type->numShapes)
                                            *(nodes->numDim)
                                            *(out->ReferenceElement->numQuadNodes),double);
     if (! (Finley_checkPtr(out->volume) || Finley_checkPtr(out->DSDX)) ) {
          if (nodes->numDim==1) {
             if (out->ReferenceElement->Type->numLocalDim==0) {

             } else if (out->ReferenceElement->Type->numLocalDim==1) {
                  Assemble_jacobeans_1D(nodes->Coordinates,out->ReferenceElement->numQuadNodes,out->ReferenceElement->QuadWeights,
                                        shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                        shape->dSdv,out->ReferenceElement->Type->numShapes,out->ReferenceElement->dSdv,
                                        out->DSDX,out->volume,self->Id);
             } else {
                  Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: local dimenion in a 1D domain has to be 0 or 1.");
             }
          } else if (nodes->numDim==2) {
             if (out->ReferenceElement->Type->numLocalDim==0) {

             } else if (out->ReferenceElement->Type->numLocalDim==1) {
                  if (out->ReferenceElement->Type->numDim==2) {
                      Assemble_jacobeans_2D_M1D_E2D(nodes->Coordinates,out->ReferenceElement->numQuadNodes,out->ReferenceElement->QuadWeights,
                                                    shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                                    shape->dSdv,out->ReferenceElement->Type->numShapes,out->ReferenceElement->dSdv,
                                                    out->DSDX,out->volume,self->Id);

                  }  else if (out->ReferenceElement->Type->numDim==1) {
                      Assemble_jacobeans_2D_M1D_E1D(nodes->Coordinates,out->ReferenceElement->numQuadNodes,out->ReferenceElement->QuadWeights,
                                                    shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                                    shape->dSdv,out->ReferenceElement->Type->numShapes,out->ReferenceElement->dSdv,
                                                    out->DSDX,out->volume,self->Id);
                  } else {
                    Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: element dimension for local dimenion 1 in a 2D domain has to be 1 or 2.");
                  }
             } else if (out->ReferenceElement->Type->numLocalDim==2) {
                  Assemble_jacobeans_2D(nodes->Coordinates,out->ReferenceElement->numQuadNodes,out->ReferenceElement->QuadWeights,
                                        shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                        shape->dSdv,out->ReferenceElement->Type->numShapes,out->ReferenceElement->dSdv,
                                        out->DSDX,out->volume,self->Id);
             } else {
               Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: local dimenion in a 2D domain has to be  1 or 2.");
             }
          } else if (nodes->numDim==3) {
             if (out->ReferenceElement->Type->numLocalDim==0) {

             } else if (out->ReferenceElement->Type->numLocalDim==2) {
                  if (out->ReferenceElement->Type->numDim==3) {
                      Assemble_jacobeans_3D_M2D_E3D(nodes->Coordinates,out->ReferenceElement->numQuadNodes,out->ReferenceElement->QuadWeights,
                                                    shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                                    shape->dSdv,out->ReferenceElement->Type->numShapes,out->ReferenceElement->dSdv,
                                                    out->DSDX,out->volume,self->Id);

                  }  else if (out->ReferenceElement->Type->numDim==2) {
                      Assemble_jacobeans_3D_M2D_E2D(nodes->Coordinates,out->ReferenceElement->numQuadNodes,out->ReferenceElement->QuadWeights,
                                                    shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                                    shape->dSdv,out->ReferenceElement->Type->numShapes,out->ReferenceElement->dSdv,
                                                    out->DSDX,out->volume,self->Id);
                  } else {
                    Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: element dimension for local dimenion 2 in a 3D domain has to be 3 or 2.");
                  }
             } else if (out->ReferenceElement->Type->numLocalDim==3) {
                  Assemble_jacobeans_3D(nodes->Coordinates,out->ReferenceElement->numQuadNodes,out->ReferenceElement->QuadWeights,
                                        shape->Type->numShapes,self->numElements,numNodes,self->Nodes,
                                        shape->dSdv,out->ReferenceElement->Type->numShapes,out->ReferenceElement->dSdv,
                                        out->DSDX,out->volume,self->Id);
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
/*
 * $Log$
 *
 */
