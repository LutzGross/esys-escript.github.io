/* $Id$ */

/**************************************************************/

/*   writes data and mesh in an opendx file */

/**************************************************************/

/*   Copyrights by ACcESS, Australia 2004 */
/*   Author: Lutz Gross, gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "Common.h"
#include "Mesh.h"
#include "escript/Data/DataC.h"

void Finley_Mesh_saveDX(const char * filename_p, Finley_Mesh *mesh_p, escriptDataC* data_p) {
  /* if there is no mesh we just return */
  if (mesh_p==NULL) return;
  /* some tables needed for reordering */
  int resort[6][9]={
                    {0,1},   /* line */
                    {0,1,2},  /* triagle */
                    {0,1,2,3}, /* tetrahedron */
                    {0,3,1,2}, /* quadrilateral */
                    {3,0,7,4,2,1,6,5}, /* hexahedron */
                   };
  Finley_ElementFile* elements=NULL;
  char elemTypeStr[32];
  int i,j,k,numDXNodesPerElement,*resortIndex,isCellCentered=TRUE,nodetype=FINLEY_DEGREES_OF_FREEDOM;
  double* values,rtmp;
  int nDim = mesh_p->Nodes->numDim;
  /* open the file  and check handel */
  FILE * fileHandle_p = fopen(filename_p, "w");
  if (fileHandle_p==NULL) {
         Finley_ErrorCode=IO_ERROR;
         sprintf(Finley_ErrorMsg,"File %s could not be opened for writing.",filename_p);
         return;
  }
  /* positions */
  fprintf(fileHandle_p, "object 1 class array type float rank 1 shape %d items %d data follows\n",
                                                                   nDim, mesh_p->Nodes->reducedNumNodes);
  for (i = 0; i < mesh_p->Nodes->numNodes; i++) {
    if (mesh_p->Nodes->toReduced[i]>=0) {
       fprintf(fileHandle_p, "%g", mesh_p->Nodes->Coordinates[INDEX2(0, i, nDim)]);
       for (j = 1; j < nDim; j++) fprintf(fileHandle_p, " %g",mesh_p->Nodes->Coordinates[INDEX2(j, i, nDim)]);
       fprintf(fileHandle_p, "\n");
    }
  } 
  /* connections */
  /* get a pointer to the relevant mesh component */
  if (isEmpty(data_p)) {
      elements=mesh_p->Elements;
  } else {
      switch(getFunctionSpaceType(data_p)) {
       case(FINLEY_DEGREES_OF_FREEDOM):
          nodetype=FINLEY_DEGREES_OF_FREEDOM;
          isCellCentered=FALSE;
          elements=mesh_p->Elements;
          break;
       case(FINLEY_REDUCED_DEGREES_OF_FREEDOM):
          nodetype=FINLEY_REDUCED_DEGREES_OF_FREEDOM;
          isCellCentered=FALSE;
          break;
          elements=mesh_p->Elements;
       case(FINLEY_NODES):
          nodetype=FINLEY_NODES;
          isCellCentered=FALSE;
          elements=mesh_p->Elements;
          break;
       case(FINLEY_ELEMENTS):
          isCellCentered=TRUE;
          elements=mesh_p->Elements;
          break;
       case(FINLEY_FACE_ELEMENTS):
          isCellCentered=TRUE;
          elements=mesh_p->FaceElements;
          break;
       case(FINLEY_POINTS):
          isCellCentered=TRUE;
          elements=mesh_p->Points;
          break;
       case(FINLEY_CONTACT_ELEMENTS_1):
       case(FINLEY_CONTACT_ELEMENTS_2):
          isCellCentered=TRUE;
          elements=mesh_p->ContactElements;
          break;
       default:
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"Finley does not know anything about function space type %d",getFunctionSpaceType(data_p));
          return;
     }
  }
  /* if no element table is present jump over the connection table */
  if (elements!=NULL) {
       ElementTypeId TypeId = elements->ReferenceElement->Type->TypeId;
       if (TypeId==Line2 || TypeId==Line3 || TypeId==Line4 ) {
          numDXNodesPerElement=2;
          resortIndex=resort[0];
          strcpy(elemTypeStr, "lines");
       } else if (TypeId==Tri3 || TypeId==Tri6 || TypeId==Tri9 || TypeId==Tri10 ) {
          numDXNodesPerElement = 3;
          resortIndex=resort[1];
          strcpy(elemTypeStr, "triangles");
       } else if (TypeId==Rec4 || TypeId==Rec8 || TypeId==Rec9 || TypeId==Rec12 || TypeId==Rec16 ) {
          numDXNodesPerElement = 4;
          resortIndex=resort[3];
          strcpy(elemTypeStr, "quads");
        } else if (TypeId==Tet4 || TypeId==Tet10 || TypeId==Tet16 ) {
          numDXNodesPerElement = 4;
          resortIndex=resort[2];
          strcpy(elemTypeStr, "tetrahedra");
        } else if (TypeId==Hex8 || TypeId==Hex20 || TypeId==Hex32 ) {
          numDXNodesPerElement = 8;
          resortIndex=resort[4];
          strcpy(elemTypeStr, "cubes");
        } else {
          Finley_ErrorCode=VALUE_ERROR;
          sprintf(Finley_ErrorMsg, "Element type %s is not supported by DX",elements->ReferenceElement->Type->Name);
          return;
        } 
        int NN=elements->ReferenceElement->Type->numNodes;
        fprintf(fileHandle_p, "object 2 class array type int rank 1 shape %d items %d data follows\n",numDXNodesPerElement, elements->numElements);
        for (i = 0; i < elements->numElements; i++) {
          fprintf(fileHandle_p,"%d",mesh_p->Nodes->toReduced[elements->Nodes[INDEX2(resortIndex[0], i, NN)]]);
          for (j = 1; j < numDXNodesPerElement; j++) {
             fprintf(fileHandle_p," %d",mesh_p->Nodes->toReduced[elements->Nodes[INDEX2(resortIndex[j], i, NN)]]);
          }
          fprintf(fileHandle_p, "\n");
        } 
        fprintf(fileHandle_p, "attribute \"element type\" string \"%s\"\n",elemTypeStr);
        fprintf(fileHandle_p, "attribute \"ref\" string \"positions\"\n");

  }
  /* data */
  if (!isEmpty(data_p)) {
      int rank=getDataPointRank(data_p);
      int nComp=getDataPointSize(data_p);
      fprintf(fileHandle_p, "object 3 class array type float rank %d ", rank);
      if (0 < rank) {
         fprintf(fileHandle_p, "shape ");
         for (i = 0; i < rank; i++) fprintf(fileHandle_p, "%d ", getDataPointShape(data_p,i));
      }
      if (isCellCentered) {
          int numPointsPerSample=elements->ReferenceElement->numQuadNodes;
          if (numPointsPerSample>0) {
             fprintf(fileHandle_p, "items %d data follows\n", elements->numElements);
             for (i=0;i<elements->numElements;i++) {
                 values=getSampleData(data_p,i);
                 for (k=0;k<nComp;k++) {
                     rtmp=0.;
                     for (j=0;j<numPointsPerSample;j++) rtmp+=values[INDEX2(k,j,nComp)];
                     fprintf(fileHandle_p, " %g", rtmp/numPointsPerSample);
                 }
	         fprintf(fileHandle_p, "\n");
             }
             fprintf(fileHandle_p, "attribute \"dep\" string \"connections\"\n");
         }
      } else {
          fprintf(fileHandle_p, "items %d data follows\n", mesh_p->Nodes->reducedNumNodes);
          for (i=0;i<mesh_p->Nodes->numNodes;i++) {
              if (mesh_p->Nodes->toReduced[i]>=0) {
                 switch (nodetype) {
                    case(FINLEY_DEGREES_OF_FREEDOM):
                       values=getSampleData(data_p,mesh_p->Nodes->degreeOfFreedom[i]);
                       break;
                    case(FINLEY_REDUCED_DEGREES_OF_FREEDOM):
                       values=getSampleData(data_p,mesh_p->Nodes->reducedDegreeOfFreedom[i]);
                       break;
                    case(FINLEY_NODES):
                       values=getSampleData(data_p,i);
                       break;
                 }
                 for (k=0;k<nComp;k++) fprintf(fileHandle_p, " %g", values[k]);
	         fprintf(fileHandle_p, "\n");
              }
          }
          fprintf(fileHandle_p, "attribute \"dep\" string \"positions\"\n");
      }
  }

  /* and finish it up */
  fprintf(fileHandle_p, "object 4 class field\n");
  fprintf(fileHandle_p, "component \"positions\" value 1\n");
  if (elements!=NULL) fprintf(fileHandle_p, "component \"connections\" value 2\n");
  if (!isEmpty(data_p)) fprintf(fileHandle_p, "component \"data\" value 3\n");
  fprintf(fileHandle_p, "end\n");
  /* close the file */
  fclose(fileHandle_p);
  return;
}

/*
 * $Log$
 * Revision 1.4  2005/07/08 04:07:55  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.5  2005/06/29 02:34:53  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1.2.4  2005/03/03 05:01:12  gross
 * bug in saveDX fixed
 *
 * Revision 1.1.1.1.2.3  2005/02/17 23:43:06  cochrane
 * Fixed error throwing bug.  Default case of switch statement should have ended
 * with return instead of break, hence errors weren't being thrown (but they now
 * should be).
 *
 * Revision 1.1.1.1.2.2  2005/02/17 05:53:26  gross
 * some bug in saveDX fixed: in fact the bug was in
 * DataC/getDataPointShape
 *
 * Revision 1.1.1.1.2.1  2005/02/17 03:23:01  gross
 * some performance improvements in MVM
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.1  2004/07/27 08:27:11  gross
 * Finley: saveDX added: now it is possible to write data on boundary and contact elements
 *
 */
