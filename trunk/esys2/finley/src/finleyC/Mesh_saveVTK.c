/* $Id$ */

/**************************************************************/

/*   writes data and mesh in a vtk file */

/**************************************************************/

/*   Copyrights by ACcESS, Australia 2004 */
/*   Author: Paul Cochrane, cochrane@esscc.uq.edu.au */

/**************************************************************/

#include "Finley.h"
#include "Common.h"
#include "Mesh.h"
#include "escript/Data/DataC.h"

void Finley_Mesh_saveVTK(const char * filename_p, Finley_Mesh *mesh_p, escriptDataC* data_p) {
  /* if there is no mesh we just return */
  if (mesh_p==NULL) return;
  /* some tables needed for reordering */
  int resort[6][9]={
                    {0,1},             /* line */
                    {0,1,2},           /* triangle */
                    {0,1,2,3},         /* tetrahedron */
                    {0,3,1,2},         /* quadrilateral */
                    {3,0,7,4,2,1,6,5}, /* hexahedron */
                   };
  Finley_ElementFile* elements=NULL;
  char elemTypeStr[32];
  int i,j,k,numVTKNodesPerElement,*resortIndex,isCellCentered,nodetype;
  double* values,rtmp;
  int nDim = mesh_p->Nodes->numDim;
  /* open the file  and check handle */
  FILE * fileHandle_p = fopen(filename_p, "w");
  if (fileHandle_p==NULL) {
         Finley_ErrorCode=IO_ERROR;
         sprintf(Finley_ErrorMsg,"File %s could not be opened for writing.",filename_p);
         return;
  }
  /* xml header */
  fprintf(fileHandle_p, "<?xml version=\"1.0\"?>\n");
  fprintf(fileHandle_p, 
	  "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n");

  /* finley uses an unstructured mesh, so UnstructuredGrid *should* work */
  fprintf(fileHandle_p, "<UnstructuredGrid>\n");

  /* is there only one "piece" to the data?? */
  /* fprintf(fileHandle_p, 
	  "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",);
	  */

  /* now for the point data */
  /* fprintf(fileHandle_p, 
	  "<PointData Scalars=\"%s\" Vectors=\"%s\">\n",);
	  */

  fprintf(fileHandle_p, "</PointData>\n");

  /* now for the cell data */
  /* fprintf(fileHandle_p,
	  "<CellData Scalars=\"%s\" Vectors=\"%s\">\n",);
	  */

  fprintf(fileHandle_p, "</CellData>\n");

  /* now for the points */
  fprintf(fileHandle_p, "<Points>\n");
  fprintf(fileHandle_p, "<DataArray NumberOfComponents=\"%d\">\n");
  fprintf(fileHandle_p, "</DataArray>\n");
  fprintf(fileHandle_p, "</Points>\n");

  /* now for the cells */
  fprintf(fileHandle_p, "<Cells>\n");
  /* fprintf(fileHandle_p, 
	  "<DataArray type=\"%s\" Name=\"%s\" format=\"ascii\">\n",);
	  */
  fprintf(fileHandle_p,
	  "</DataArray>\n");
  fprintf(fileHandle_p, "</Cells>\n");

  /* finish off the piece */
  fprintf(fileHandle_p, "</Piece>\n");

  /* positions */
  fprintf(fileHandle_p, "object 1 class array type float rank 1 shape %d items %d data follows\n", nDim, mesh_p->Nodes->reducedNumNodes);
  for (i = 0; i < mesh_p->Nodes->numNodes; i++) {
    if (mesh_p->Nodes->toReduced[i]>=0) {
       fprintf(fileHandle_p, "%e", mesh_p->Nodes->Coordinates[INDEX2(0, i, nDim)]);
       for (j = 1; j < nDim; j++) fprintf(fileHandle_p, " %f",mesh_p->Nodes->Coordinates[INDEX2(j, i, nDim)]);
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
       case(FINLEY_REDUCED_DEGREES_OF_FREEDOM):
          nodetype=FINLEY_REDUCED_DEGREES_OF_FREEDOM;
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
          break;
     }
  }
  /* if no element table is present jump over the connection table */
  if (elements!=NULL) {
       ElementTypeId TypeId = elements->ReferenceElement->Type->TypeId;
       if (TypeId==Line2 || TypeId==Line3 || TypeId==Line4 ) {
          numVTKNodesPerElement=2;
          resortIndex=resort[0];
          strcpy(elemTypeStr, "lines");
       } else if (TypeId==Tri3 || TypeId==Tri6 || TypeId==Tri9 || TypeId==Tri10 ) {
          numVTKNodesPerElement = 3;
          resortIndex=resort[1];
          strcpy(elemTypeStr, "triangles");
       } else if (TypeId==Rec4 || TypeId==Rec8 || TypeId==Rec9 || TypeId==Rec12 || TypeId==Rec16 ) {
          numVTKNodesPerElement = 4;
          resortIndex=resort[3];
          strcpy(elemTypeStr, "quads");
        } else if (TypeId==Tet4 || TypeId==Tet10 || TypeId==Tet16 ) {
          numVTKNodesPerElement = 4;
          resortIndex=resort[2];
          strcpy(elemTypeStr, "tetrahedra");
        } else if (TypeId==Hex8 || TypeId==Hex20 || TypeId==Hex32 ) {
          numVTKNodesPerElement = 8;
          resortIndex=resort[4];
          strcpy(elemTypeStr, "cubes");
        } else {
          Finley_ErrorCode=VALUE_ERROR;
          sprintf(Finley_ErrorMsg, "Element type %s is not supported by VTK",elements->ReferenceElement->Type->Name);
          return;
        } 
        int NN=elements->ReferenceElement->Type->numNodes;
        fprintf(fileHandle_p, "object 2 class array type int rank 1 shape %d items %d data follows\n",numVTKNodesPerElement, elements->numElements);
        for (i = 0; i < elements->numElements; i++) {
          fprintf(fileHandle_p,"%d",mesh_p->Nodes->toReduced[mesh_p->Elements->Nodes[INDEX2(resortIndex[0], i, NN)]]);
          for (j = 1; j < numVTKNodesPerElement; j++) {
             fprintf(fileHandle_p," %d",mesh_p->Nodes->toReduced[mesh_p->Elements->Nodes[INDEX2(resortIndex[j], i, NN)]]);
          }
          fprintf(fileHandle_p, "\n");
        } 
        fprintf(fileHandle_p, "attribute \"element type\" string \"%s\"\n",elemTypeStr);
        fprintf(fileHandle_p, "attribute \"ref\" string \"positions\"\n");

  }
  /* data */
  if (!isEmpty(data_p)) {
      int rank=getDataPointRank(data_p);
      int* shape=getDataPointShape(data_p);
      int nComp=getDataPointSize(data_p);
      fprintf(fileHandle_p, "object 3 class array type float rank %d ", rank);
      if (0 < rank) {
         fprintf(fileHandle_p, "shape ");
         for (i = 0; i < rank; i++) fprintf(fileHandle_p, "%d ", shape[i]);
      }
      if (isCellCentered) {
          int numPointsPerSample=elements->ReferenceElement->numQuadNodes;
          if (numPointsPerSample) {
             fprintf(fileHandle_p, "items %d data follows\n", elements->numElements);
             for (i=0;i<elements->numElements;i++) {
                 values=getSampleData(data_p,i);
                 for (k=0;k<nComp;k++) {
                     rtmp=0.;
                     for (j=0;j<numPointsPerSample;j++) rtmp+=values[INDEX2(k,j,nComp)];
                     fprintf(fileHandle_p, " %f", rtmp/numPointsPerSample);
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
              }
              for (k=0;k<nComp;k++) fprintf(fileHandle_p, " %f", values[k]);
	      fprintf(fileHandle_p, "\n");
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

  fprintf(fileHandle_p, "</UnstructuredGrid>\n");
  /* write the xml footer */
  fprintf(fileHandle_p, "</VTKFile>\n");
  /* close the file */
  fclose(fileHandle_p);
  return;
}

/*
 * $Log$
 * Revision 1.2  2005/02/14 04:14:42  jgs
 * *** empty log message ***
 *
 * Revision 1.1.2.2  2005/02/10 01:34:22  cochrane
 * Quick fix to make sure that saveVTK compiles so that finley is still buildable.  Apologies to those this has affected.
 *
 * Revision 1.1.2.1  2005/02/09 06:53:15  cochrane
 * Initial import to repository.  This is the file to implement saving finley/escript meshes out to vtk formatted files.  It is basically just a hack of the opendx equivalent, with a lot of the opendx stuff still in the file, so it doesn't actually work just yet, but it probably needs to be added to the cvs.
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.1  2004/07/27 08:27:11  gross
 * Finley: saveDX added: now it is possible to write data on boundary and contact elements
 *
 */
