/**************************************************************/

/*   Finley: Mesh */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003,04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Finley.h"
#include "Mesh.h"

/**************************************************************/

/*   allocates a Mesh with name name for elements of type id using an integration order. If order is negative, */
/*   the most appropriate order is selected indepently. */

extern Finley_RefElementInfo Finley_RefElement_InfoList[];

Finley_Mesh* Finley_Mesh_alloc(char* name,int numDim, int order) {
  Finley_Mesh *out;
  
  /*  allocate the return value */
  
  out=(Finley_Mesh*)MEMALLOC(sizeof(Finley_Mesh));
  if (Finley_checkPtr(out)) return NULL;
  out->Name=NULL;  
  out->Nodes=NULL; 
  out->Elements=NULL;   
  out->FaceElements=NULL; 
  out->Points=NULL;      
  out->ContactElements=NULL;      
  out->reference_counter=0;
  
  /*   copy name: */
  
  out->Name=(char*)MEMALLOC((strlen(name)+1)*sizeof(char));
  if (Finley_checkPtr(out->Name)) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  strcpy(out->Name,name);
  
  /*   allocate node table: */
  
  out->Nodes=Finley_NodeFile_alloc(numDim);
  if (Finley_ErrorCode!=NO_ERROR) {
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
       MEMFREE(in);      
     }
  }
}

/**************************************************************/

/*  returns the spatial dimension of the mesh: */

int Finley_Mesh_getDim(Finley_Mesh *in) {
  return in->Nodes->numDim;
}

/**************************************************************/

/*  returns the number of nodes in the mesh: */

int Finley_Mesh_getNumNodes(Finley_Mesh *in) {
  return in->Nodes->numNodes;
}
/**************************************************************/

/*  returns the number of degrees of freedom in the mesh: */

int Finley_Mesh_getNumDegreesOfFreedom(Finley_Mesh *in) {
  return in->Nodes->numDegreesOfFreedom;
}
/**************************************************************/

/*  returns the number of degrees of freedom in the mesh: */

int Finley_Mesh_getReducedNumDegreesOfFreedom(Finley_Mesh *in) {
  return in->Nodes->reducedNumDegreesOfFreedom;
}
/* 
* $Log$
* Revision 1.1  2004/10/26 06:53:57  jgs
* Initial revision
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/
