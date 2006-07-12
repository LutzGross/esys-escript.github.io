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

#include "Assemble.h"
#ifdef _OPENMP
#include <omp.h>
#endif


/**************************************************************/

Finley_ElementFile_Jacobeans* Finley_ElementFile_Jacobeans_alloc(void)
{
  Finley_ElementFile_Jacobeans* out=MEMALLOC(1,Finley_ElementFile_Jacobeans);
  if (Finley_checkPtr(out)) {
     return NULL;
  } else {
     out->status=FINLEY_INITIAL_STATUS-1;
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
     Finley_RefElement *shape, *test;
     if (reducedShapefunction) {
       if (reducedIntegrationOrder) {
           shape=self->LinearReferenceElement;
           test=self->LinearReferenceElement;
       } else {
           shape=self->LinearReferenceElement;
           test=self->LinearReferenceElement;
       }
     } else {
       if (reducedIntegrationOrder) {
           shape=self->ReferenceElement;
           test=self->ReferenceElement;
       } else {
           shape=self->ReferenceElement;
           test=self->ReferenceElement;
       }
     }

     if (out->volume==NULL) out->volume=MEMALLOC((self->numElements)*(shape->numQuadNodes),double);
     if (out->DSDX==NULL) out->DSDX=MEMALLOC((self->numElements)*(test->Type->numShapes)*(shape->Type->numDim)*(shape->numQuadNodes),double);
     if (! (Finley_checkPtr(out->volume) || Finley_checkPtr(out->DSDX)) ) {
          if (nodes->numDim==1) {
             if (shape->Type->numDim==0) {

             } else if (shape->Type->numDim==1) {
                  Assemble_jacobeans_1D(nodes->Coordinates,shape->numQuadNodes,shape->QuadWeights,shape->Type->numShapes,
                                        self->numElements,self->Nodes,
                                        shape->dSdv,test->Type->numShapes,test->dSdv,
                                        out->DSDX,out->volume,self->Id);
             } else {
                  Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: local dimenion in a 1D domain has to be less or equal 1.");
             }
          } else if (nodes->numDim==2) {
             if (shape->Type->numDim==0) {

             } else if (shape->Type->numDim==1) {
                  Assemble_jacobeans_2D_M1D(nodes->Coordinates,shape->numQuadNodes,shape->QuadWeights,shape->Type->numShapes,
                                            self->numElements,self->Nodes,
                                            shape->dSdv,test->Type->numShapes,test->dSdv,
                                            out->DSDX,out->volume,self->Id);
             } else if (shape->Type->numDim==2) {
                  Assemble_jacobeans_2D(nodes->Coordinates,shape->numQuadNodes,shape->QuadWeights,shape->Type->numShapes,
                                        self->numElements,self->Nodes,
                                        shape->dSdv,test->Type->numShapes,test->dSdv,
                                        out->DSDX,out->volume,self->Id);
             } else {
               Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: local dimenion in a 2D domain has to be less or equal 2.");
             }
          } else if (nodes->numDim==3) {
             if (shape->Type->numDim==0) {

             } else if (shape->Type->numDim==1) {
                  Assemble_jacobeans_3D_M1D(nodes->Coordinates,shape->numQuadNodes,shape->QuadWeights,shape->Type->numShapes,
                                            self->numElements,self->Nodes,
                                            shape->dSdv,test->Type->numShapes,test->dSdv,
                                            out->DSDX,out->volume,self->Id);
             } else if (shape->Type->numDim==2) {
                  Assemble_jacobeans_3D_M2D(nodes->Coordinates,shape->numQuadNodes,shape->QuadWeights,shape->Type->numShapes,
                                            self->numElements,self->Nodes,
                                            shape->dSdv,test->Type->numShapes,test->dSdv,
                                            out->DSDX,out->volume,self->Id);
             } else if (shape->Type->numDim==3) {
                  Assemble_jacobeans_3D(nodes->Coordinates,shape->numQuadNodes,shape->QuadWeights,shape->Type->numShapes,
                                        self->numElements,self->Nodes,
                                        shape->dSdv,test->Type->numShapes,test->dSdv,
                                        out->DSDX,out->volume,self->Id);
             } else {
               Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: local dimenion in a 3D domain has to be less or equal 3.");
             }
          } else {
            Finley_setError(SYSTEM_ERROR,"Finley_ElementFile_borrowJacobeans: spatial dimension has to be less or equal 3.");
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
