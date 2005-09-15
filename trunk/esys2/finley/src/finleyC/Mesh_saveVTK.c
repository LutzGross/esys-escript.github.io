/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2003,2004,2005 -  All Rights Reserved              *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/


/**************************************************************/

/*   writes data and mesh in a vtk file */

/**************************************************************/

/*   Author: Paul Cochrane, cochrane@esscc.uq.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"
#include "vtkCellType.h"  /* copied from vtk source directory !!! */

/**************************************************************/

void Finley_Mesh_saveVTK(const char * filename_p, Finley_Mesh *mesh_p, escriptDataC* data_p) {
  char error_msg[LenErrorMsg_MAX];
  /* if there is no mesh we just return */
  if (mesh_p==NULL) return;
  Finley_ElementFile* elements=NULL;
  char elemTypeStr[32];
  int i, j, k, numVTKNodesPerElement, isCellCentered=FALSE, nodetype=FINLEY_DEGREES_OF_FREEDOM;
  index_t j2;
  double* values, rtmp;
  int nDim = mesh_p->Nodes->numDim;
  int numPoints=0;
  
  if (nDim==1) {
      Finley_setError(VALUE_ERROR,"saveVTK: 1-dimensional domains are not supported.");
      return;
  }
  /* get a pointer to the relevant mesh component */
  if (isEmpty(data_p)) {
    numPoints = mesh_p->Nodes->numNodes;
    elements=mesh_p->Elements;
  } else {
    switch(getFunctionSpaceType(data_p)) {
    case(FINLEY_DEGREES_OF_FREEDOM):
      numPoints = mesh_p->Nodes->numNodes;
      nodetype = FINLEY_DEGREES_OF_FREEDOM;
      isCellCentered = FALSE;
      elements = mesh_p->Elements;
      break;
    case(FINLEY_REDUCED_DEGREES_OF_FREEDOM):
      numPoints = mesh_p->Nodes->reducedNumNodes;
      nodetype =FINLEY_REDUCED_DEGREES_OF_FREEDOM;
      isCellCentered = FALSE;
      elements = mesh_p->Elements;
      break;
    case(FINLEY_NODES):
      numPoints = mesh_p->Nodes->numNodes;
      nodetype=FINLEY_NODES;
      isCellCentered=FALSE;
      elements=mesh_p->Elements;
      break;
    case(FINLEY_ELEMENTS):
      numPoints = mesh_p->Nodes->numNodes;
      nodetype=FINLEY_NODES;
      isCellCentered=TRUE;
      elements=mesh_p->Elements;
      break;
    case(FINLEY_FACE_ELEMENTS):
      numPoints = mesh_p->Nodes->numNodes;
      nodetype=FINLEY_NODES;
      isCellCentered=TRUE;
      elements=mesh_p->FaceElements;
      break;
    case(FINLEY_POINTS):
      numPoints = mesh_p->Nodes->numNodes;
      nodetype=FINLEY_NODES;
      isCellCentered=TRUE;
      elements=mesh_p->Points;
      break;
    case(FINLEY_CONTACT_ELEMENTS_1):
    case(FINLEY_CONTACT_ELEMENTS_2):
      numPoints = mesh_p->Nodes->numNodes;
      nodetype=FINLEY_NODES;
      isCellCentered=TRUE;
      elements=mesh_p->ContactElements;
      break;
    default:
      sprintf(error_msg,"saveVTK: Finley does not know anything about function space type %d",getFunctionSpaceType(data_p));
      Finley_setError(TYPE_ERROR,error_msg);
      return;
    }
  }

  /* the number of points */

  /* the number of cells */
  if (elements == NULL) {
    Finley_setError(VALUE_ERROR,"saveVTK: elements object is NULL; cannot proceed");
    return;
  }
  int numCells = elements->numElements;   
  /* open the file and check handle */
  FILE * fileHandle_p = fopen(filename_p, "w");
  if (fileHandle_p==NULL) {
    sprintf(error_msg, "saveVTK: File %s could not be opened for writing.", filename_p);
    Finley_setError(IO_ERROR,error_msg);
    return;
  }
  /* xml header */
  fprintf(fileHandle_p, "<?xml version=\"1.0\"?>\n");
  fprintf(fileHandle_p, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n");

  /* finley uses an unstructured mesh, so UnstructuredGrid *should* work */
  fprintf(fileHandle_p, "<UnstructuredGrid>\n");

  /* is there only one "piece" to the data?? */
  fprintf(fileHandle_p, "<Piece "
	      "NumberOfPoints=\"%d\" "
	      "NumberOfCells=\"%d\">\n",
	      numPoints, numCells);
  /* now for the points; equivalent to positions section in saveDX() */
  /* "The points element explicitly defines coordinates for each point
   * individually.  It contains one DataArray element describing an array
   * with three components per value, each specifying the coordinates of one
   * point" - from Vtk User's Guide
   */
  fprintf(fileHandle_p, "<Points>\n");
  /* 
   * the reason for this if statement is explained in the long comment below
   */
  if (nDim < 3) {
    fprintf(fileHandle_p, "<DataArray "
	    "NumberOfComponents=\"3\" "
	    "type=\"Float32\" "
	    "format=\"ascii\">\n");
  } else {
    fprintf(fileHandle_p, "<DataArray "
	    "NumberOfComponents=\"%d\" "
	    "type=\"Float32\" "
	    "format=\"ascii\">\n", 
	    nDim);
  }
  /* vtk/mayavi doesn't like 2D data, it likes 3D data with a degenerate
     * third dimension to handle 2D data (like a sheet of paper).  So, if
     * nDim is 2, we have to append zeros to the array to get this third
     * dimension, and keep the visualisers happy.
     * Indeed, if nDim is less than 3, must pad all empty dimensions, so
     * that the total number of dims is 3.
  */
  if (nodetype==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
     for (i = 0; i < mesh_p->Nodes->numNodes; i++) {
       if (mesh_p->Nodes->toReduced[i]>=0) {
          for (j = 0; j < nDim; j++) 
            fprintf(fileHandle_p, " %e",mesh_p->Nodes->Coordinates[INDEX2(j, i, nDim)]);
          for (k=0; k<3-nDim; k++) fprintf(fileHandle_p, " %e",0.);
          fprintf(fileHandle_p, "\n");
       }
     } 
  } else {
     for (i = 0; i < mesh_p->Nodes->numNodes; i++) {
       for (j = 0; j < nDim; j++) 
         fprintf(fileHandle_p, " %e",mesh_p->Nodes->Coordinates[INDEX2(j, i, nDim)]);
       for (k=0; k<3-nDim; k++) fprintf(fileHandle_p, " %e",0.);
       fprintf(fileHandle_p, "\n");
     }
  } 
  fprintf(fileHandle_p, "</DataArray>\n");
  fprintf(fileHandle_p, "</Points>\n");

  /* connections */
  /* now for the cells */
  /* "The Cells element defines cells explicitly by specifying point
   * connectivity and cell types.  It contains three DataArray elements.  The
   * first array specifies the point connectivity.  All cells' point lists
   * are concatenated together.  The second array specifies th eoffset into
   * the connectivity array for the end of each cell.  The third array
   * specifies the type of each cell.
   */
  /* if no element table is present jump over the connection table */
  if (elements!=NULL) {
    int cellType;
    ElementTypeId TypeId;

    if (nodetype==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
      TypeId = elements->LinearReferenceElement->Type->TypeId;
    } else {
      TypeId = elements->ReferenceElement->Type->TypeId;
    }

    switch(TypeId) {
    case Point1:
      cellType = VTK_VERTEX;
      break;
    case Line2:
      cellType = VTK_LINE;
      break;
    case Line3:
      cellType = VTK_QUADRATIC_EDGE;
      break;
    case Tri3:
      cellType = VTK_TRIANGLE;
      break;
    case Tri6:
      cellType = VTK_QUADRATIC_TRIANGLE;
      break;
    case Rec4:
      cellType = VTK_QUAD;
      break;
    case Rec8:
      cellType = VTK_QUADRATIC_QUAD;
      break;
    case Tet4:
      cellType = VTK_TETRA;
      break;
    case Tet10:
      cellType = VTK_QUADRATIC_TETRA;
      break;
    case Hex8:
      cellType = VTK_HEXAHEDRON;
      break;
    case Hex20:
      cellType = VTK_QUADRATIC_HEXAHEDRON;
      break;
    case Line2Face:
      cellType = VTK_VERTEX;
      break;
    case Line3Face:
      cellType = VTK_VERTEX;
      break;
    case Tri3Face:
      cellType = VTK_LINE;
      break;
    case Tri6Face:
      cellType = VTK_QUADRATIC_EDGE;
      break;
    case Rec4Face:
      cellType = VTK_LINE;
      break;
    case Rec8Face:
      cellType = VTK_QUADRATIC_EDGE;
      break;
    case Tet4Face:
      cellType = VTK_TRIANGLE;
      break;
    case Tet10Face:
      cellType = VTK_QUADRATIC_TRIANGLE;
      break;
    case Hex8Face:
      cellType = VTK_QUAD;
      break;
    case Hex20Face:
      cellType = VTK_QUADRATIC_QUAD;
      break;
    case Point1_Contact:
      cellType = VTK_VERTEX;
      break;
    case Line2_Contact:
      cellType = VTK_LINE;
      break;
    case Line3_Contact:
      cellType = VTK_QUADRATIC_EDGE;
      break;
    case Tri3_Contact:
      cellType = VTK_LINE;
      break;
    case Tri6_Contact:
      cellType = VTK_QUADRATIC_TRIANGLE;
      break;
    case Rec4_Contact:
      cellType = VTK_QUAD;
      break;
    case Rec8_Contact:
      cellType = VTK_QUADRATIC_QUAD;
      break;
    case Line2Face_Contact:
      cellType = VTK_VERTEX;
      break;
    case Line3Face_Contact:
      cellType = VTK_VERTEX;
      break;
    case Tri3Face_Contact:
      cellType = VTK_LINE;
      break;
    case Tri6Face_Contact:
      cellType = VTK_QUADRATIC_EDGE;
      break;
    case Rec4Face_Contact:
      cellType = VTK_LINE;
      break;
    case Rec8Face_Contact:
      cellType = VTK_QUADRATIC_EDGE;
      break;
    case Tet4Face_Contact:
      cellType = VTK_TRIANGLE;
      break;
    case Tet10Face_Contact:
      cellType = VTK_QUADRATIC_TRIANGLE;
      break;
    case Hex8Face_Contact:
      cellType = VTK_QUAD;
      break;
    case Hex20Face_Contact:
      cellType = VTK_QUADRATIC_QUAD;
      break;
    default: 
      sprintf(error_msg, "saveVTK: Element type %s is not supported by VTK",elements->ReferenceElement->Type->Name);
      Finley_setError(VALUE_ERROR,error_msg);
      return;
    } 

    switch(cellType) {
    case VTK_VERTEX:
      numVTKNodesPerElement = 1;
      strcpy(elemTypeStr, "VTK_VERTEX");
      break;
    case VTK_LINE:
      numVTKNodesPerElement = 2;
      strcpy(elemTypeStr, "VTK_LINE");
      break;
    case VTK_TRIANGLE:
      numVTKNodesPerElement = 3;
      strcpy(elemTypeStr, "VTK_TRIANGLE");
      break;
    case VTK_QUAD:
      numVTKNodesPerElement = 4;
      strcpy(elemTypeStr, "VTK_QUAD");
      break;
    case VTK_TETRA:
      numVTKNodesPerElement = 4;
      strcpy(elemTypeStr, "VTK_TETRA");
      break;
    case VTK_HEXAHEDRON:
      numVTKNodesPerElement = 8;
      strcpy(elemTypeStr, "VTK_HEXAHEDRON");
      break;
    case VTK_QUADRATIC_EDGE:
      numVTKNodesPerElement = 3;
      strcpy(elemTypeStr, "VTK_QUADRATIC_EDGE");
      break;
    case VTK_QUADRATIC_TRIANGLE:
      numVTKNodesPerElement = 6;
      strcpy(elemTypeStr, "VTK_QUADRATIC_TRIANGLE");
      break;
    case VTK_QUADRATIC_QUAD:
      numVTKNodesPerElement = 8;
      strcpy(elemTypeStr, "VTK_QUADRATIC_QUAD");
      break;
    case VTK_QUADRATIC_TETRA:
      numVTKNodesPerElement = 10;
      strcpy(elemTypeStr, "VTK_QUADRATIC_TETRA");
      break;
    case VTK_QUADRATIC_HEXAHEDRON:
      numVTKNodesPerElement = 20;
      strcpy(elemTypeStr, "VTK_QUADRATIC_HEXAHEDRON");
      break;
    default:
      sprintf(error_msg,"saveVTK: Cell type %d is not supported by VTK", cellType);
      Finley_setError(VALUE_ERROR,error_msg);
      return;
    }

    /* write out the DataArray element for the connectivity */
    int NN = elements->ReferenceElement->Type->numNodes;
    fprintf(fileHandle_p, "<Cells>\n");
    fprintf(fileHandle_p, "<DataArray "
	    "Name=\"connectivity\" "
	    "type=\"Int32\" "
	    "format=\"ascii\">\n");
    if (nodetype==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
       for (i = 0; i < numCells; i++) {
          for (j = 0; j < numVTKNodesPerElement; j++) {
               j2=elements->ReferenceElement->Type->linearNodes[j];
               fprintf(fileHandle_p,"%d ",mesh_p->Nodes->toReduced[elements->Nodes[INDEX2(j2, i, NN)]]); 
          }
          fprintf(fileHandle_p, "\n");
       }
    } else if (VTK_QUADRATIC_HEXAHEDRON==cellType) {
          for (i = 0; i < numCells; i++) {
            fprintf(fileHandle_p,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", 
                               elements->Nodes[INDEX2(0, i, NN)],
                               elements->Nodes[INDEX2(1, i, NN)],
                               elements->Nodes[INDEX2(2, i, NN)],
                               elements->Nodes[INDEX2(3, i, NN)],
                               elements->Nodes[INDEX2(4, i, NN)],
                               elements->Nodes[INDEX2(5, i, NN)],
                               elements->Nodes[INDEX2(6, i, NN)],
                               elements->Nodes[INDEX2(7, i, NN)],
                               elements->Nodes[INDEX2(8, i, NN)],
                               elements->Nodes[INDEX2(9, i, NN)],
                               elements->Nodes[INDEX2(10, i, NN)],
                               elements->Nodes[INDEX2(11, i, NN)],
                               elements->Nodes[INDEX2(16, i, NN)],
                               elements->Nodes[INDEX2(17, i, NN)],
                               elements->Nodes[INDEX2(18, i, NN)],
                               elements->Nodes[INDEX2(19, i, NN)],
                               elements->Nodes[INDEX2(12, i, NN)],
                               elements->Nodes[INDEX2(13, i, NN)],
                               elements->Nodes[INDEX2(14, i, NN)],
                               elements->Nodes[INDEX2(15, i, NN)]);
          }
     } else if (numVTKNodesPerElement!=NN) {
          for (i = 0; i < numCells; i++) {
            for (j = 0; j < numVTKNodesPerElement; j++) {
                 j2=elements->ReferenceElement->Type->geoNodes[j];
                 fprintf(fileHandle_p,"%d ", elements->Nodes[INDEX2(j2, i, NN)]);
            }
            fprintf(fileHandle_p, "\n");
          }
     } else {
          for (i = 0; i < numCells; i++) {
            for (j = 0; j < numVTKNodesPerElement; j++) fprintf(fileHandle_p,"%d ", elements->Nodes[INDEX2(j, i, NN)]);
            fprintf(fileHandle_p, "\n");
          }
     } 
    fprintf(fileHandle_p, "</DataArray>\n");

    /* write out the DataArray element for the offsets */
    fprintf(fileHandle_p, "<DataArray "
	    "Name=\"offsets\" "
	    "type=\"Int32\" "
	    "format=\"ascii\">\n");
    for (i=numVTKNodesPerElement; i<=numCells*numVTKNodesPerElement; i+=numVTKNodesPerElement) 
        fprintf(fileHandle_p, "%d\n", i);
    fprintf(fileHandle_p, "</DataArray>\n");

    /* write out the DataArray element for the types */
    fprintf(fileHandle_p, "<DataArray "
	    "Name=\"types\" "
	    "type=\"UInt8\" "
	    "format=\"ascii\">\n");
    for (i=0; i<numCells; i++) 
        fprintf(fileHandle_p, "%d\n", cellType);
    fprintf(fileHandle_p, "</DataArray>\n");

    /* finish off the <Cells> element */
    fprintf(fileHandle_p, "</Cells>\n");
  }

  /* data */
  if (!isEmpty(data_p)) {
    int rank = getDataPointRank(data_p);
    int nComp = getDataPointSize(data_p);
    /* barf if rank is greater than two */
    if (rank > 2) {
      sprintf(error_msg, "saveVTK: Vtk can't handle objects with rank greater than 2. object rank = %d\n", rank);
      Finley_setError(VALUE_ERROR,error_msg);
      return;
    }
    /* if the rank == 0:   --> scalar data
     * if the rank == 1:   --> vector data
     * if the rank == 2:   --> tensor data
     */
    char dataNameStr[31], dataTypeStr[63];
    int nCompReqd=1;   /* the number of components required by vtk */
    if (rank == 0) {
      strcpy(dataNameStr, "scalar");
      sprintf(dataTypeStr, "Scalars=\"%s\"", dataNameStr);
      nCompReqd = 1;
    }
    else if (rank == 1) {
      strcpy(dataNameStr, "vector");
      sprintf(dataTypeStr, "Vectors=\"%s\"", dataNameStr);
      nCompReqd = 3;
    }
    else if (rank == 2) {
      strcpy(dataNameStr, "tensor");
      sprintf(dataTypeStr, "Tensors=\"%s\"", dataNameStr);
      nCompReqd = 9;
    }
    /* if have cell centred data then use CellData element, 
     * if have node centred data, then use PointData element
     */
    if (isCellCentered) {
      /* now for the cell data */
      fprintf(fileHandle_p, "<CellData %s>\n", dataTypeStr);
      fprintf(fileHandle_p, 
	      "<DataArray "
	      "Name=\"%s\" "
	      "type=\"Float32\" "
	      "NumberOfComponents=\"%d\" "
	      "format=\"ascii\">\n",
	      dataNameStr, nCompReqd);
      int numPointsPerSample = elements->ReferenceElement->numQuadNodes;
      if (numPointsPerSample) {
	int shape = getDataPointShape(data_p, 0);
	if (shape > 3) {
	    sprintf(error_msg,"saveVTK: shape should be 1, 2, or 3; I got %d\n", shape);
            Finley_setError(VALUE_ERROR,error_msg);
	    return;
	}
	for (i=0; i<numCells; i++) {
	  values = getSampleData(data_p, i);
	  double sampleAvg[nComp];
	  for (k=0; k<nComp; k++) {
	    /* averaging over the number of points in the sample */
	    rtmp = 0.;
	    for (j=0; j<numPointsPerSample; j++) 
	      rtmp += values[INDEX2(k,j,nComp)];
	    sampleAvg[k] = rtmp/numPointsPerSample;
	  }
	  /* if the number of required components is more than the number
	   * of actual components, pad with zeros
	   */
	  /* probably only need to get shape of first element */
	  /* write the data different ways for scalar, vector and tensor */
	  if (nCompReqd == 1) {
	    fprintf(fileHandle_p, " %e", sampleAvg[0]);
	  } else if (nCompReqd == 3) {
	    /* write out the data */
	    for (int m=0; m<shape; m++) {
	      fprintf(fileHandle_p, " %e", sampleAvg[m]);
	    }
	    /* pad with zeros */
	    for (int m=0; m<nCompReqd-shape; m++) 
	      fprintf(fileHandle_p, " %e", 0.);
	  } else if (nCompReqd == 9) { 
	    /* tensor data, so have a 3x3 matrix to output as a row 
	     * of 9 data points */
	    int count = 0;
	    for (int m=0; m<shape; m++) {
	      for (int n=0; n<shape; n++) {
		fprintf(fileHandle_p, " %e", sampleAvg[count]);
		count++;
	      }
	      for (int n=0; n<3-shape; n++) 
	        fprintf(fileHandle_p, " %e", 0.);
	    }
	    for (int m=0; m<3-shape; m++) {
	        for (int n=0; n<3; n++) 
	           fprintf(fileHandle_p, " %e", 0.);
            }
	  }
	  fprintf(fileHandle_p, "\n");
	}
      }
      fprintf(fileHandle_p, "</DataArray>\n");
      fprintf(fileHandle_p, "</CellData>\n");
    } else {
      /* now for the point data */
      fprintf(fileHandle_p, "<PointData %s>\n", dataTypeStr);
      fprintf(fileHandle_p, "<DataArray "
	      "Name=\"%s\" "
	      "type=\"Float32\" "
	      "NumberOfComponents=\"%d\" "
	      "format=\"ascii\">\n",
	      dataNameStr, nCompReqd);
	/* write out the data */
	/* if the number of required components is more than the number
	 * of actual components, pad with zeros
       */
      bool_t do_write=TRUE;
      int shape = getDataPointShape(data_p, 0);
      if (shape > 3) {
	  sprintf(error_msg,"shape should be 1, 2, or 3; I got %d\n", shape);
          Finley_setError(VALUE_ERROR,error_msg);
	  return;
      }
      for (i=0; i<mesh_p->Nodes->numNodes; i++) {
	switch (nodetype) {
	case(FINLEY_DEGREES_OF_FREEDOM): 
	  values = getSampleData(data_p,mesh_p->Nodes->degreeOfFreedom[i]);
	  break;
	case(FINLEY_REDUCED_DEGREES_OF_FREEDOM):
	  if (mesh_p->Nodes->toReduced[i]>=0) {
            do_write=TRUE;
	    values = getSampleData(data_p,mesh_p->Nodes->reducedDegreeOfFreedom[i]);
	  } else {
            do_write=FALSE;
          }
	  break;
	case(FINLEY_NODES):
	  values = getSampleData(data_p,i);
	  break;
	}
	/* write the data different ways for scalar, vector and tensor */
        if (do_write) {
	    if (nCompReqd == 1) {
	      fprintf(fileHandle_p, " %e", values[0]);
	    }
	    else if (nCompReqd == 3) {
	      /* write out the data */
	      for (int m=0; m<shape; m++) 
	        fprintf(fileHandle_p, " %e", values[m]);
	      /* pad with zeros */
	      for (int m=0; m<nCompReqd-shape; m++) 
	          fprintf(fileHandle_p, " %e", 0.);
	    }
	    else if (nCompReqd == 9) { 
	      /* tensor data, so have a 3x3 matrix to output as a row 
	       * of 9 data points */
	      int count = 0;
	      for (int m=0; m<shape; m++) {
	        for (int n=0; n<shape; n++) {
	          fprintf(fileHandle_p, " %e", values[count]);
	          count++;
	        }
	        for (int n=0; n<3-shape; n++) 
	          fprintf(fileHandle_p, " %e", 0.);
	      }
	      for (int m=0; m<3-shape; m++)  {
	          for (int n=0; n<3; n++) 
	              fprintf(fileHandle_p, " %e", 0.);
              }
	    }
	    fprintf(fileHandle_p, "\n");
         }
      }
      fprintf(fileHandle_p, "</DataArray>\n");
      fprintf(fileHandle_p, "</PointData>\n");
    }
  }

  /* finish off the piece */
  fprintf(fileHandle_p, "</Piece>\n");

  fprintf(fileHandle_p, "</UnstructuredGrid>\n");
  /* write the xml footer */
  fprintf(fileHandle_p, "</VTKFile>\n");
  /* close the file */
  fclose(fileHandle_p);
  return;
}

/*
 * Revision 1.6  2005/08/12 01:45:43  jgs
 * erge of development branch dev-02 back to main trunk on 2005-08-12
 *
 * Revision 1.5.2.4  2005/09/09 08:15:17  gross
 * bugs in vtk and dx writer fixed
 *
 * Revision 1.5.2.3  2005/09/08 08:28:39  gross
 * some cleanup in savevtk
 *
 * Revision 1.5.2.2  2005/09/07 06:26:20  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.5.2.1  2005/08/10 06:14:37  gross
 * QUADRATIC HEXAHEDRON elements fixed
 *
 * Revision 1.5  2005/07/08 04:07:55  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.4  2005/05/06 04:26:15  jgs
 * Merge of development branch back to main trunk on 2005-05-06
 *
 * Revision 1.1.2.7  2005/06/29 02:34:54  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.2.6  2005/05/06 01:17:19  cochrane
 * Fixed incorrect reporting of number of components in PointData arrays for
 * vector data.
 *
 * Revision 1.1.2.5  2005/05/05 05:38:44  cochrane
 * Improved formatting of VTK file output.
 *
 * Revision 1.1.2.4  2005/02/22 10:03:54  cochrane
 * Implementation of writing of vtk xml files from finley.  This function will
 * require more testing, but on the cases that I have tried (and with the help
 * of Lutz and mayavi), it looks like it's producing the correct output.  Testing
 * with more realistic data would be good.  I'm at least confident enough to
 * commit my changes.
 *
 * Revision 1.1.2.3  2005/02/17 05:53:26  gross
 * some bug in saveDX fixed: in fact the bug was in
 * DataC/getDataPointShape
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
