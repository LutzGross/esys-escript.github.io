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

/*   Finley: Mesh */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

/*   allocates a Mesh with name name for elements of type id using an integration order. If order is negative, */
/*   the most appropriate order is selected indepently. */

extern Finley_RefElementInfo Finley_RefElement_InfoList[];

#ifndef PASO_MPI
Finley_Mesh* Finley_Mesh_alloc(char* name,dim_t numDim, index_t order) 
#else
Finley_Mesh* Finley_Mesh_alloc(char* name,dim_t numDim, index_t order, Paso_MPIInfo *mpi_info) 
#endif
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
  out->reference_counter=0;

  out->FullFullPattern=NULL;
  out->FullReducedPattern=NULL;
  out->ReducedFullPattern=NULL;
  out->ReducedReducedPattern=NULL;

#ifdef PASO_MPI 
  out->MPIInfo = NULL;
 
  /* get MPI info */
  out->MPIInfo = Paso_MPIInfo_getReference( mpi_info );
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
#endif
  /*   copy name: */
  
  out->Name=MEMALLOC(strlen(name)+1,char);
  if (Finley_checkPtr(out->Name)) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  strcpy(out->Name,name);
  
  /*   allocate node table: */
#ifdef PASO_MPI
  out->Nodes=Finley_NodeFile_alloc( numDim, mpi_info );
#else
  out->Nodes=Finley_NodeFile_alloc(numDim);
#endif
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  out->order=order;
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

/*   deallocates a mesh: */

void Finley_Mesh_dealloc(Finley_Mesh* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<1) {
       #ifdef Finley_TRACE
       if (in->Name!=NULL) {
           printf("Finley_Mesh_dealloc: mesh %s is deallocated.\n",in->Name);
       } else {
           printf("Finley_Mesh_dealloc\n");
       }
       #endif
       MEMFREE(in->Name);
       Finley_NodeFile_dealloc(in->Nodes);
       Finley_ElementFile_dealloc(in->Elements);   
       Finley_ElementFile_dealloc(in->FaceElements);
       Finley_ElementFile_dealloc(in->ContactElements);
       Finley_ElementFile_dealloc(in->Points);
       Paso_SystemMatrixPattern_dealloc(in->FullFullPattern);
       Paso_SystemMatrixPattern_dealloc(in->FullReducedPattern);
       Paso_SystemMatrixPattern_dealloc(in->ReducedFullPattern);
       Paso_SystemMatrixPattern_dealloc(in->ReducedReducedPattern);
#ifdef PASO_MPI
       Paso_MPIInfo_dealloc( in->MPIInfo );
#endif
       MEMFREE(in);      
     }
  }
}

/**************************************************************/

/*  returns the spatial dimension of the mesh: */

dim_t Finley_Mesh_getDim(Finley_Mesh *in) {
  return in->Nodes->numDim;
}

/**************************************************************/

/*  returns the number of nodes in the mesh: */

dim_t Finley_Mesh_getNumNodes(Finley_Mesh *in) {
  return in->Nodes->numNodes;
}
/**************************************************************/

/*  returns the number of degrees of freedom in the mesh: */

dim_t Finley_Mesh_getNumDegreesOfFreedom(Finley_Mesh *in) {
  return in->Nodes->numDegreesOfFreedom;
}
/**************************************************************/

/*  returns the number of degrees of freedom in the mesh: */

dim_t Finley_Mesh_getReducedNumDegreesOfFreedom(Finley_Mesh *in) {
  return in->Nodes->reducedNumDegreesOfFreedom;
}
/* 
* $Log$
* Revision 1.6  2005/09/15 03:44:22  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.5.2.1  2005/09/07 06:26:19  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.5  2005/07/08 04:07:51  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.4  2004/12/15 07:08:32  jgs
* *** empty log message ***
* Revision 1.1.1.1.2.3  2005/06/29 02:34:51  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1.2.2  2004/11/24 01:37:13  gross
* some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
*
*
*
*/
