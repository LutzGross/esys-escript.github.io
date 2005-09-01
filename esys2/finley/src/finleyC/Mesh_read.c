/**************************************************************/

/*   Finley: read mesh */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003/04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Finley.h"
#include "Mesh.h"

/**************************************************************/

/*  reads a mesh from a Finley file of name fname */

Finley_Mesh* Finley_Mesh_read(char* fname,index_t order) {

  dim_t numNodes, numDim, numEle, i0, i1;
  char name[LenString_MAX],element_type[LenString_MAX];
  Finley_NodeFile *nodes_p=NULL;
  double time0=Finley_timer();

  /* get file handle */
  FILE * fileHandle_p = fopen(fname, "r");
  if (fileHandle_p==NULL) {
    Finley_ErrorCode=IO_ERROR;
    sprintf(Finley_ErrorMsg,"Opening file %s for reading failed.",fname);
    return NULL;
  }

  /* read header */
  fscanf(fileHandle_p, "%63[^\n]", name);

  /* get the nodes */

  fscanf(fileHandle_p, "%1d%*s %d\n", &numDim,&numNodes);
  nodes_p=Finley_NodeFile_alloc(numDim);
  if (Finley_ErrorCode!=NO_ERROR) return NULL;
  Finley_NodeFile_allocTable(nodes_p, numNodes);
  if (Finley_ErrorCode!=NO_ERROR) return NULL;

  if (1 == numDim) {
      for (i0 = 0; i0 < numNodes; i0++)
	fscanf(fileHandle_p, "%d %d %d %le\n", &nodes_p->Id[i0],
	       &nodes_p->degreeOfFreedom[i0], &nodes_p->Tag[i0],
	       &nodes_p->Coordinates[INDEX2(0,i0,numDim)]);
  } else if (2 == numDim) {
      for (i0 = 0; i0 < numNodes; i0++)
	fscanf(fileHandle_p, "%d %d %d %le %le\n", &nodes_p->Id[i0],
	       &nodes_p->degreeOfFreedom[i0], &nodes_p->Tag[i0],
	       &nodes_p->Coordinates[INDEX2(0,i0,numDim)],
	       &nodes_p->Coordinates[INDEX2(1,i0,numDim)]);
  } else if (3 == numDim) {
      for (i0 = 0; i0 < numNodes; i0++)
	fscanf(fileHandle_p, "%d %d %d %le %le %le\n", &nodes_p->Id[i0],
	       &nodes_p->degreeOfFreedom[i0], &nodes_p->Tag[i0],
	       &nodes_p->Coordinates[INDEX2(0,i0,numDim)],
	       &nodes_p->Coordinates[INDEX2(1,i0,numDim)],
	       &nodes_p->Coordinates[INDEX2(2,i0,numDim)]);
  } /* if else else */

  /* get the element type */

  fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
  ElementTypeId typeID=Finley_RefElement_getTypeId(element_type);
  if (typeID==NoType) {
    Finley_ErrorCode=VALUE_ERROR;
    sprintf(Finley_ErrorMsg,"Unidentified element type %s",element_type);
    return NULL;
  }

  /* allocate mesh */

  Finley_Mesh * mesh_p =Finley_Mesh_alloc(name,numDim,order);
  if (Finley_ErrorCode!=NO_ERROR) return NULL;
  mesh_p->Nodes=nodes_p;

  /* read the elements */
  mesh_p->Elements=Finley_ElementFile_alloc(typeID,mesh_p->order);
  Finley_ElementFile_allocTable(mesh_p->Elements, numEle);
  mesh_p->Elements->minColor=0;
  mesh_p->Elements->maxColor=numEle-1;
  for (i0 = 0; i0 < numEle; i0++) {
    fscanf(fileHandle_p, "%d %d", &mesh_p->Elements->Id[i0], &mesh_p->Elements->Tag[i0]);
    mesh_p->Elements->Color[i0]=i0;
    for (i1 = 0; i1 < mesh_p->Elements->ReferenceElement->Type->numNodes; i1++) {
         fscanf(fileHandle_p, " %d",
            &mesh_p->Elements->Nodes[INDEX2(i1, i0, mesh_p->Elements->ReferenceElement->Type->numNodes)]);
    }	/* for i1 */
    fscanf(fileHandle_p, "\n");
  } /* for i0 */

  /* get the face elements */
  fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
  ElementTypeId faceTypeID=Finley_RefElement_getTypeId(element_type);
  if (faceTypeID==NoType) {
    Finley_ErrorCode=VALUE_ERROR;
    sprintf(Finley_ErrorMsg,"Unidentified element type %s for face elements",element_type);
    return NULL;
  }
  mesh_p->FaceElements=Finley_ElementFile_alloc(faceTypeID,mesh_p->order);
  Finley_ElementFile_allocTable(mesh_p->FaceElements, numEle);
  mesh_p->FaceElements->minColor=0;
  mesh_p->FaceElements->maxColor=numEle-1;
  for (i0 = 0; i0 < numEle; i0++) {
    fscanf(fileHandle_p, "%d %d", &mesh_p->FaceElements->Id[i0], &mesh_p->FaceElements->Tag[i0]);
    mesh_p->FaceElements->Color[i0]=i0;
    for (i1 = 0; i1 < mesh_p->FaceElements->ReferenceElement->Type->numNodes; i1++) {
         fscanf(fileHandle_p, " %d",
            &mesh_p->FaceElements->Nodes[INDEX2(i1, i0, mesh_p->FaceElements->ReferenceElement->Type->numNodes)]);
    }	/* for i1 */
    fscanf(fileHandle_p, "\n");
  } /* for i0 */

  /* get the Contact face element */
  fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
  ElementTypeId contactTypeID=Finley_RefElement_getTypeId(element_type);
  if (contactTypeID==NoType) {
    Finley_ErrorCode=VALUE_ERROR;
    sprintf(Finley_ErrorMsg,"Unidentified element type %s for contact elements",element_type);
    return NULL;
  }
  mesh_p->ContactElements=Finley_ElementFile_alloc(contactTypeID,mesh_p->order);
  Finley_ElementFile_allocTable(mesh_p->ContactElements, numEle);
  mesh_p->ContactElements->minColor=0;
  mesh_p->ContactElements->maxColor=numEle-1;
  for (i0 = 0; i0 < numEle; i0++) {
    fscanf(fileHandle_p, "%d %d", &mesh_p->ContactElements->Id[i0], &mesh_p->ContactElements->Tag[i0]);
    mesh_p->ContactElements->Color[i0]=i0;
    for (i1 = 0; i1 < mesh_p->ContactElements->ReferenceElement->Type->numNodes; i1++) {
        fscanf(fileHandle_p, " %d",
           &mesh_p->ContactElements->Nodes[INDEX2(i1, i0, mesh_p->ContactElements->ReferenceElement->Type->numNodes)]);
    }	/* for i1 */
    fscanf(fileHandle_p, "\n");
  } /* for i0 */

  /* get the nodal element */
  fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
  ElementTypeId pointTypeID=Finley_RefElement_getTypeId(element_type);
  if (pointTypeID==NoType) {
    Finley_ErrorCode=VALUE_ERROR;
    sprintf(Finley_ErrorMsg,"Unidentified element type %s for points",element_type);
    return NULL;
  }
  mesh_p->Points=Finley_ElementFile_alloc(pointTypeID,mesh_p->order);

  Finley_ElementFile_allocTable(mesh_p->Points, numEle);
  mesh_p->Points->minColor=0;
  mesh_p->Points->maxColor=numEle-1;
  for (i0 = 0; i0 < numEle; i0++) {
    fscanf(fileHandle_p, "%d %d", &mesh_p->Points->Id[i0], &mesh_p->Points->Tag[i0]);
    mesh_p->Points->Color[i0]=i0;
    for (i1 = 0; i1 < mesh_p->Points->ReferenceElement->Type->numNodes; i1++) {
        fscanf(fileHandle_p, " %d",
           &mesh_p->Points->Nodes[INDEX2(i1, i0, mesh_p->Points->ReferenceElement->Type->numNodes)]);
    }	/* for i1 */
    fscanf(fileHandle_p, "\n");
  } /* for i0 */


  /* close file */

  fclose(fileHandle_p);

  /*   resolve id's : */

  Finley_Mesh_resolveNodeIds(mesh_p);

  /* rearrange elements: */

  Finley_Mesh_prepare(mesh_p);
  printf ("nodes read!\n");

  /* that's it */
  #ifdef Finley_TRACE
  printf("timing: reading mesh: %.4e sec\n",Finley_timer()-time0);
  #endif
  if (Finley_ErrorCode!=NO_ERROR) Finley_Mesh_dealloc(mesh_p);
  return mesh_p;
}
/*
* $Log$
* Revision 1.4  2005/09/01 03:31:36  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-01
*
* Revision 1.3.2.1  2005/08/24 02:02:18  gross
* timing output switched off. solver output can be swiched through getSolution(verbose=True) now.
*
* Revision 1.3  2005/07/22 03:53:08  jgs
* Merge of development branch back to main trunk on 2005-07-22
*
* Revision 1.2  2005/07/08 04:07:54  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.1.1.1.2.2  2005/07/18 10:34:54  gross
* some informance improvements when reading meshes
*
* Revision 1.1.1.1.2.1  2005/06/29 02:34:53  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/

