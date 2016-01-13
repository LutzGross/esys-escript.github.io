
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/*   Finley: Mesh */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

/*   allocates a Mesh with name name for elements of type id using an integration order. If order is negative, */
/*   the most appropriate order is selected indepently. */

extern Finley_RefElementInfo Finley_RefElement_InfoList[];

Finley_Mesh* Finley_Mesh_alloc(char* name,dim_t numDim, index_t order, index_t reduced_order, Paso_MPIInfo *mpi_info) 
{
  Finley_Mesh *out;
  
  /*  allocate the return value */
  
  out=MEMALLOC(1,Finley_Mesh);
  if (Finley_checkPtr(out)) return NULL;
  out->Name=NULL;  
  out->Nodes=NULL; 
  out->Elements=NULL;   
  out->FaceElements=NULL; 
  out->Points=NULL;      
  out->ContactElements=NULL;      
  out->TagMap=NULL;      
  out->reference_counter=0;

  out->FullFullPattern=NULL;
  out->FullReducedPattern=NULL;
  out->ReducedFullPattern=NULL;
  out->ReducedReducedPattern=NULL;
  out->MPIInfo = Paso_MPIInfo_getReference( mpi_info );
  if (! Finley_noError()) {
      Finley_Mesh_free(out);
      return NULL;
  }
  /*   copy name: */
  
  out->Name=MEMALLOC(strlen(name)+1,char);
  if (Finley_checkPtr(out->Name)) {
      Finley_Mesh_free(out);
      return NULL;
  }
  strcpy(out->Name,name);
  
  /*   allocate node table: */
  out->Nodes=Finley_NodeFile_alloc( numDim, mpi_info );
  if (! Finley_noError()) {
      Finley_Mesh_free(out);
      return NULL;
  }
  out->order=order;
  out->reduced_order=reduced_order;
  out->Elements=NULL;
  out->FaceElements=NULL;
  out->Points=NULL;
  out->ContactElements=NULL;
  out->reference_counter++;
  return out;
}

/* returns a reference to Finley_Mesh in */

Finley_Mesh* Finley_Mesh_reference(Finley_Mesh* in) {
     if (in!=NULL) ++(in->reference_counter);
     return in;
}

/*   freeates a mesh: */

void Finley_Mesh_free(Finley_Mesh* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<1) {
       #ifdef Finley_TRACE
       if (in->Name!=NULL) {
           printf("Finley_Mesh_free: mesh %s is freed.\n",in->Name);
       } else {
           printf("Finley_Mesh_free\n");
       }
       #endif
       MEMFREE(in->Name);
       Finley_NodeFile_free(in->Nodes);
       Finley_ElementFile_free(in->Elements);   
       Finley_ElementFile_free(in->FaceElements);
       Finley_ElementFile_free(in->ContactElements);
       Finley_ElementFile_free(in->Points);
       Finley_TagMap_free(in->TagMap);
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

dim_t Finley_Mesh_getDim(Finley_Mesh *in) {
  return in->Nodes->numDim;
}

void Finley_Mesh_setElements(Finley_Mesh* self,Finley_ElementFile *elements) {
    Finley_ElementFile_free(self->Elements);
    self->Elements=elements;
}
void Finley_Mesh_setFaceElements(Finley_Mesh* self,Finley_ElementFile *elements) {
    Finley_ElementFile_free(self->FaceElements);
    self->FaceElements=elements;
}
void Finley_Mesh_setContactElements(Finley_Mesh* self,Finley_ElementFile *elements) {
    Finley_ElementFile_free(self->ContactElements);
    self->ContactElements=elements;
}
void Finley_Mesh_setPoints(Finley_Mesh* self,Finley_ElementFile *elements) {
    Finley_ElementFile_free(self->Points);
    self->Points=elements;
}

