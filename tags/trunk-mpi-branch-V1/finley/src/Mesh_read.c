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

/*   Finley: read mesh */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

/*  reads a mesh from a Finley file of name fname */

Finley_Mesh* Finley_Mesh_read(char* fname,index_t order, index_t reduced_order,  bool_t optimize_labeling) {

  dim_t numNodes, numDim, numEle, i0, i1;
  index_t tag_key;
  Finley_Mesh *mesh_p=NULL;
  char name[LenString_MAX],element_type[LenString_MAX],frm[20];
  char error_msg[LenErrorMsg_MAX];
  double time0=Finley_timer();
  FILE *fileHandle_p = NULL;
  ElementTypeId typeID, faceTypeID, contactTypeID, pointTypeID;

  Finley_resetError();
#ifdef PASO_MPI
  /* TODO */
  Finley_setError(SYSTEM_ERROR,"Finley_Mesh_read: MPI is not suporrted yet.");
#endif

  /* get file handle */
  fileHandle_p = fopen(fname, "r");
  if (fileHandle_p==NULL) {
    sprintf(error_msg,"Finley_Mesh_read: Opening file %s for reading failed.",fname);
    Finley_setError(IO_ERROR,error_msg);
    return NULL;
  }

  /* read header */
  sprintf(frm,"%%%d[^\n]",LenString_MAX-1);
  fscanf(fileHandle_p, frm, name);

  /* get the nodes */

  fscanf(fileHandle_p, "%1d%*s %d\n", &numDim,&numNodes);
  /* allocate mesh */
#ifdef PASO_MPI
  mesh_p = Finley_Mesh_alloc(name,numDim,order,reduced_order,NULL); /* NULL should to be replaced by mpi_info before use */
#else
  mesh_p = Finley_Mesh_alloc(name,numDim,order,reduced_order);
#endif
  if (! Finley_noError()) return NULL;

  Finley_NodeFile_allocTable(mesh_p->Nodes, numNodes);
  if (! Finley_noError()) return NULL;

  if (1 == numDim) {
      for (i0 = 0; i0 < numNodes; i0++)
	fscanf(fileHandle_p, "%d %d %d %le\n", &mesh_p->Nodes->Id[i0],
	       &mesh_p->Nodes->degreeOfFreedom[i0], &mesh_p->Nodes->Tag[i0],
	       &mesh_p->Nodes->Coordinates[INDEX2(0,i0,numDim)]);
  } else if (2 == numDim) {
      for (i0 = 0; i0 < numNodes; i0++)
	fscanf(fileHandle_p, "%d %d %d %le %le\n", &mesh_p->Nodes->Id[i0],
	       &mesh_p->Nodes->degreeOfFreedom[i0], &mesh_p->Nodes->Tag[i0],
	       &mesh_p->Nodes->Coordinates[INDEX2(0,i0,numDim)],
	       &mesh_p->Nodes->Coordinates[INDEX2(1,i0,numDim)]);
  } else if (3 == numDim) {
      for (i0 = 0; i0 < numNodes; i0++)
	fscanf(fileHandle_p, "%d %d %d %le %le %le\n", &mesh_p->Nodes->Id[i0],
	       &mesh_p->Nodes->degreeOfFreedom[i0], &mesh_p->Nodes->Tag[i0],
	       &mesh_p->Nodes->Coordinates[INDEX2(0,i0,numDim)],
	       &mesh_p->Nodes->Coordinates[INDEX2(1,i0,numDim)],
	       &mesh_p->Nodes->Coordinates[INDEX2(2,i0,numDim)]);
  } /* if else else */

  /* get the element type */

  fscanf(fileHandle_p, "%s %d\n", element_type, &numEle);
  typeID=Finley_RefElement_getTypeId(element_type);
  if (typeID==NoType) {
    sprintf(error_msg,"Finley_Mesh_read :Unidentified element type %s",element_type);
    Finley_setError(VALUE_ERROR,error_msg);
    return NULL;
  }
  /* read the elements */
#ifdef PASO_MPI
  mesh_p->Elements=Finley_ElementFile_alloc(typeID,mesh_p->order, mesh_p->reduced_order, NULL); /* NULL should be replaced by mpi_info before use */
#else
  mesh_p->Elements=Finley_ElementFile_alloc(typeID,mesh_p->order, mesh_p->reduced_order);
#endif
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
  faceTypeID=Finley_RefElement_getTypeId(element_type);
  faceTypeID=Finley_RefElement_getTypeId(element_type);
  if (faceTypeID==NoType) {
    sprintf(error_msg,"Finley_Mesh_read :Unidentified element type %s for face elements",element_type);
    Finley_setError(VALUE_ERROR,error_msg);
    return NULL;
  }
#ifdef PASO_MPI
  mesh_p->FaceElements=Finley_ElementFile_alloc(faceTypeID,mesh_p->order, mesh_p->reduced_order, NULL); /* NULL should be replaced by mpi_info before use */
#else
  mesh_p->FaceElements=Finley_ElementFile_alloc(faceTypeID,mesh_p->order, mesh_p->reduced_order);
#endif
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
  contactTypeID=Finley_RefElement_getTypeId(element_type);
  if (contactTypeID==NoType) {
    sprintf(error_msg,"Finley_Mesh_read: Unidentified element type %s for contact elements",element_type);
    Finley_setError(VALUE_ERROR,error_msg);
    return NULL;
  }
#ifdef PASO_MPI
  mesh_p->ContactElements=Finley_ElementFile_alloc(contactTypeID,mesh_p->order, mesh_p->reduced_order, NULL); /* NULL should be replaced by mpi_info before use */
#else
  mesh_p->ContactElements=Finley_ElementFile_alloc(contactTypeID,mesh_p->order, mesh_p->reduced_order);
#endif
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
  pointTypeID=Finley_RefElement_getTypeId(element_type);
  if (pointTypeID==NoType) {
    sprintf(error_msg,"Finley_Mesh_read: Unidentified element type %s for points",element_type);
    Finley_setError(VALUE_ERROR,error_msg);
    return NULL;
  }
#ifdef PASO_MPI
  mesh_p->Points=Finley_ElementFile_alloc(pointTypeID,mesh_p->order, mesh_p->reduced_order, NULL); /* NULL should be replaced by mpi_info before use */
#else
  mesh_p->Points=Finley_ElementFile_alloc(pointTypeID,mesh_p->order, mesh_p->reduced_order);
#endif
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
  /* get the name tags */
  if (feof(fileHandle_p) == 0) {
     fscanf(fileHandle_p, "%s\n", name);
     while (feof(fileHandle_p) == 0) {
       fscanf(fileHandle_p, "%s %d\n", name, &tag_key);
       Finley_Mesh_addTagMap(mesh_p,name,tag_key);
     }
  }
  /* close file */

  fclose(fileHandle_p);

  /*   resolve id's : */

  Finley_Mesh_resolveNodeIds(mesh_p);

  /* rearrange elements: */

  Finley_Mesh_prepare(mesh_p);

  /* that's it */
  #ifdef Finley_TRACE
  printf("timing: reading mesh: %.4e sec\n",Finley_timer()-time0);
  #endif
  if (Finley_noError()) {
       if ( ! Finley_Mesh_isPrepared(mesh_p)) {
          Finley_setError(SYSTEM_ERROR,"Mesh is not prepared for calculation. Contact the programmers.");
       }
  } else {
       Finley_Mesh_dealloc(mesh_p);
  }
  return mesh_p;
}