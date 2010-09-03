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


/**************************************************************/

/*   Dudley: Mesh */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

/*   allocates a Mesh with name name for elements of type id using an integration order. If order is negative, */
/*   the most appropriate order is selected indepently. */

Dudley_Mesh* Dudley_Mesh_alloc(char* name,dim_t numDim, Paso_MPIInfo *mpi_info) 
{
  Dudley_Mesh *out;
  
  /*  allocate the return value */
  
  out=MEMALLOC(1,Dudley_Mesh);
  if (Dudley_checkPtr(out)) return NULL;
  out->Name=NULL;  
  out->Nodes=NULL; 
  out->Elements=NULL;   
  out->FaceElements=NULL; 
  out->Points=NULL;      
  out->TagMap=NULL;      
  out->reference_counter=0;

  out->FullFullPattern=NULL;
  out->FullReducedPattern=NULL;
  out->ReducedFullPattern=NULL;
  out->ReducedReducedPattern=NULL;
  out->MPIInfo = Paso_MPIInfo_getReference( mpi_info );
  if (! Dudley_noError()) {
      Dudley_Mesh_free(out);
      return NULL;
  }
  /*   copy name: */
  
  out->Name=MEMALLOC(strlen(name)+1,char);
  if (Dudley_checkPtr(out->Name)) {
      Dudley_Mesh_free(out);
      return NULL;
  }
  strcpy(out->Name,name);
  
  /*   allocate node table: */
  out->Nodes=Dudley_NodeFile_alloc( numDim, mpi_info );
  if (! Dudley_noError()) {
      Dudley_Mesh_free(out);
      return NULL;
  }
  out->approximationOrder=-1;
  out->reducedApproximationOrder=-1;
  out->integrationOrder=-1;
  out->reducedIntegrationOrder=-1;

  out->Elements=NULL;
  out->FaceElements=NULL;
  out->Points=NULL;
  out->reference_counter++;
  return out;
}

/* returns a reference to Dudley_Mesh in */

Dudley_Mesh* Dudley_Mesh_reference(Dudley_Mesh* in) {
     if (in!=NULL) ++(in->reference_counter);
     return in;
}

/*   freeates a mesh: */

void Dudley_Mesh_free(Dudley_Mesh* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<1) {
       MEMFREE(in->Name);
       Dudley_NodeFile_free(in->Nodes);
       Dudley_ElementFile_free(in->FaceElements);
       Dudley_ElementFile_free(in->Elements);   
       Dudley_ElementFile_free(in->Points);
       Dudley_TagMap_free(in->TagMap);
       Paso_SystemMatrixPattern_free(in->FullFullPattern);
       Paso_SystemMatrixPattern_free(in->FullReducedPattern);
       Paso_SystemMatrixPattern_free(in->ReducedFullPattern);
       Paso_SystemMatrixPattern_free(in->ReducedReducedPattern);
       Paso_MPIInfo_free( in->MPIInfo );
       MEMFREE(in);      
     }
  }
}

/**************************************************************/

/*  returns the spatial dimension of the mesh: */

dim_t Dudley_Mesh_getDim(Dudley_Mesh *in) {
  return in->Nodes->numDim;
}

void Dudley_Mesh_setElements(Dudley_Mesh* self,Dudley_ElementFile *elements) {
    Dudley_ElementFile_free(self->Elements);
    self->Elements=elements;
}
void Dudley_Mesh_setFaceElements(Dudley_Mesh* self,Dudley_ElementFile *elements) {
    Dudley_ElementFile_free(self->FaceElements);
    self->FaceElements=elements;
}
void Dudley_Mesh_setPoints(Dudley_Mesh* self,Dudley_ElementFile *elements) {
    Dudley_ElementFile_free(self->Points);
    self->Points=elements;
}
int  Dudley_Mesh_getStatus(Dudley_Mesh* in) {
   if  (in == NULL) {
        return -1;
   } else if (in->Nodes == NULL) {
        return -1;
   } else {
        return in->Nodes->status;
   }
}

void Mesh_setOrders(Dudley_Mesh *in) 
{
   const dim_t order_max=9999999;
   dim_t locals[3];
   #ifdef PASO_MPI
       dim_t globals[4];
   #endif
   locals[0]=order_max; locals[1]=order_max; locals[2]=order_max;

  if ( in->Elements!=NULL) {
     if (in->Elements->numElements > 0) {
         locals[0]=MIN(locals[0], in->Elements->referenceElementSet->referenceElement->BasisFunctions->Type->numOrder);
         locals[1]=MIN(locals[1], in->Elements->referenceElementSet->referenceElement->integrationOrder);
         locals[2]=MIN(locals[2], in->Elements->referenceElementSet->referenceElementReducedQuadrature->integrationOrder);
     }
  }
  if ( in->FaceElements!=NULL) {
     if (in->FaceElements->numElements > 0) {
         locals[0]=MIN(locals[0], in->FaceElements->referenceElementSet->referenceElement->BasisFunctions->Type->numOrder);
         locals[1]=MIN(locals[1], in->FaceElements->referenceElementSet->referenceElement->integrationOrder);
         locals[2]=MIN(locals[2], in->FaceElements->referenceElementSet->referenceElementReducedQuadrature->integrationOrder);
     }


  }

   #ifdef PASO_MPI
       MPI_Allreduce( locals, globals, 4, MPI_INT, MPI_MIN, in->MPIInfo->comm );
       in->approximationOrder=(globals[0] < order_max ? globals[0] : -1 );
       in->reducedApproximationOrder=(globals[1] < order_max ? globals[1] : -1 );
       in->integrationOrder=(globals[2] < order_max ? globals[2] : -1 );
       in->reducedIntegrationOrder=(globals[3] < order_max ? globals[3] : -1 );
   #else
       in->approximationOrder=(locals[0] < order_max ? locals[0] : -1 );
       in->reducedApproximationOrder=(locals[0] < order_max ? locals[0] : -1 );
       in->integrationOrder=(locals[1] < order_max ? locals[1] : -1 );
       in->reducedIntegrationOrder=(locals[2] < order_max ? locals[2] : -1 );
   #endif


}
  
