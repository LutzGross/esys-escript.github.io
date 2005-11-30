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

void Finley_Mesh_saveVTK(const char * filename_p, Finley_Mesh *mesh_p, const dim_t num_data,char* *names_p, escriptDataC* *data_pp) {
  char error_msg[LenErrorMsg_MAX];
  /* if there is no mesh we just return */
  if (mesh_p==NULL) return;

  int i, j, k, numVTKNodesPerElement,i_data;
  index_t j2;
  double* values, rtmp;
  
  /* open the file and check handle */

  FILE * fileHandle_p = fopen(filename_p, "w");
  if (fileHandle_p==NULL) {
    sprintf(error_msg, "saveVTK: File %s could not be opened for writing.", filename_p);
    Finley_setError(IO_ERROR,error_msg);
    return;
  }
  /* find the mesh type to be written */
  int nodetype=FINLEY_DEGREES_OF_FREEDOM;
  int elementtype=FINLEY_UNKNOWN;
  bool_t isCellCentered[num_data],write_celldata=FALSE,write_pointdata=FALSE;
  for (i_data=0;i_data<num_data;++i_data) {
     if (! isEmpty(data_pp[i_data])) {
        switch(getFunctionSpaceType(data_pp[i_data])) {
           case FINLEY_DEGREES_OF_FREEDOM:
             nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS) {
                 elementtype=FINLEY_ELEMENTS;
             } else {
                 Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
                 fclose(fileHandle_p);
                 return;
             }
             isCellCentered[i_data]=FALSE;
             break;
           case FINLEY_REDUCED_DEGREES_OF_FREEDOM:
             nodetype = FINLEY_REDUCED_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS) {
                 elementtype=FINLEY_ELEMENTS;
             } else {
                 Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
                 fclose(fileHandle_p);
                 return;
             }
             isCellCentered[i_data]=FALSE;
             break;
           case FINLEY_NODES:
             nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS) {
                 elementtype=FINLEY_ELEMENTS;
             } else {
                 Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
                 fclose(fileHandle_p);
                 return;
             }
             isCellCentered[i_data]=FALSE;
             break;
           case FINLEY_ELEMENTS:
             nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS) {
                 elementtype=FINLEY_ELEMENTS;
             } else {
                 Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
                 fclose(fileHandle_p);
                 return;
             }
             isCellCentered[i_data]=TRUE;
             break;
           case FINLEY_FACE_ELEMENTS:
             nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_FACE_ELEMENTS) {
                 elementtype=FINLEY_FACE_ELEMENTS;
             } else {
                 Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
                 fclose(fileHandle_p);
                 return;
             }
             isCellCentered[i_data]=TRUE;
             break;
           case FINLEY_POINTS:
             nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_POINTS) {
                 elementtype=FINLEY_POINTS;
             } else {
                 Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
                 fclose(fileHandle_p);
                 return;
             }
             isCellCentered[i_data]=TRUE;
             break;
           case FINLEY_CONTACT_ELEMENTS_1:
             nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_CONTACT_ELEMENTS_1) {
                 elementtype=FINLEY_CONTACT_ELEMENTS_1;
             } else {
                 Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
                 fclose(fileHandle_p);
                 return;
             }
             isCellCentered[i_data]=TRUE;
             break;
           case FINLEY_CONTACT_ELEMENTS_2:
             nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_CONTACT_ELEMENTS_1) {
                 elementtype=FINLEY_CONTACT_ELEMENTS_1;
             } else {
                 Finley_setError(TYPE_ERROR,"saveVTK: cannot write given data in single file.");
                 fclose(fileHandle_p);
                 return;
             } 
             isCellCentered[i_data]=TRUE;
             break;
           default:
             sprintf(error_msg,"saveVTK: Finley does not know anything about function space type %d",getFunctionSpaceType(data_pp[i_data]));
             Finley_setError(TYPE_ERROR,error_msg);
             fclose(fileHandle_p);
             return;
        }
        if (isCellCentered[i_data]) {
           write_celldata=TRUE;
        } else {
           write_pointdata=TRUE;
        }
     }
  }
  /* select nomber of points and the mesh component */
  int numPoints = mesh_p->Nodes->numNodes;
  if (nodetype==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
      numPoints = mesh_p->Nodes->reducedNumNodes;
  } else {
      numPoints = mesh_p->Nodes->numNodes;
  }
  if (elementtype==FINLEY_UNKNOWN) elementtype=FINLEY_ELEMENTS;
  Finley_ElementFile* elements=NULL;
  switch(elementtype) {
    case FINLEY_ELEMENTS:
      elements=mesh_p->Elements;
      break;
    case FINLEY_FACE_ELEMENTS:
      elements=mesh_p->FaceElements;
      break;
    case FINLEY_POINTS:
      elements=mesh_p->Points;
      break;
    case FINLEY_CONTACT_ELEMENTS_1:
      elements=mesh_p->ContactElements;
      break;
  }
  if (elements==NULL) {
     Finley_setError(SYSTEM_ERROR,"saveVTK: undefined element file");
     fclose(fileHandle_p);
     return;
  }
  /* map finley element type to VTK element type */
  int numCells = elements->numElements;   
  int cellType;
  ElementTypeId TypeId;
  char elemTypeStr[32];
  if (nodetype==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
     TypeId = elements->LinearReferenceElement->Type->TypeId;
  } else {
     TypeId = elements->ReferenceElement->Type->TypeId;
  }

  switch(TypeId) {
    case Point1:
    case Line2Face:
    case Line3Face:
    case Point1_Contact:
    case Line2Face_Contact:
    case Line3Face_Contact:
      cellType = VTK_VERTEX;
      numVTKNodesPerElement = 1;
      strcpy(elemTypeStr, "VTK_VERTEX");
      break;

    case Line2:
    case Tri3Face:
    case Rec4Face:
    case Line2_Contact:
    case Tri3_Contact:
    case Tri3Face_Contact:
    case Rec4Face_Contact:
      cellType = VTK_LINE;
      numVTKNodesPerElement = 2;
      strcpy(elemTypeStr, "VTK_LINE");
      break;

    case Tri3:
    case Tet4Face:
    case Tet4Face_Contact:
      cellType = VTK_TRIANGLE;
      numVTKNodesPerElement = 3;
      strcpy(elemTypeStr, "VTK_TRIANGLE");
      break;

    case Rec4:
    case Hex8Face:
    case Rec4_Contact:
    case Hex8Face_Contact:
      cellType = VTK_QUAD;
      numVTKNodesPerElement = 4;
      strcpy(elemTypeStr, "VTK_QUAD");
      break;

    case Tet4:
      cellType = VTK_TETRA;
      numVTKNodesPerElement = 4;
      strcpy(elemTypeStr, "VTK_TETRA");
      break;

    case Hex8:
      cellType = VTK_HEXAHEDRON;
      numVTKNodesPerElement = 8;
      strcpy(elemTypeStr, "VTK_HEXAHEDRON");
      break;

    case Line3:
    case Tri6Face:
    case Rec8Face:
    case Line3_Contact:
    case Tri6Face_Contact:
    case Rec8Face_Contact:
      cellType = VTK_QUADRATIC_EDGE;
      numVTKNodesPerElement = 3;
      strcpy(elemTypeStr, "VTK_QUADRATIC_EDGE");
      break;

    case Tri6:
    case Tet10Face:
    case Tri6_Contact:
    case Tet10Face_Contact:
      cellType = VTK_QUADRATIC_TRIANGLE;
      numVTKNodesPerElement = 6;
      strcpy(elemTypeStr, "VTK_QUADRATIC_TRIANGLE");
      break;

    case Rec8:
    case Hex20Face:
    case Rec8_Contact:
    case Hex20Face_Contact:
      cellType = VTK_QUADRATIC_QUAD;
      numVTKNodesPerElement = 8;
      strcpy(elemTypeStr, "VTK_QUADRATIC_QUAD");
      break;

    case Tet10:
      cellType = VTK_QUADRATIC_TETRA;
      numVTKNodesPerElement = 10;
      strcpy(elemTypeStr, "VTK_QUADRATIC_TETRA");
      break;

    case Hex20:
      cellType = VTK_QUADRATIC_HEXAHEDRON;
      numVTKNodesPerElement = 20;
      strcpy(elemTypeStr, "VTK_QUADRATIC_HEXAHEDRON");
      break;

    default: 
      sprintf(error_msg, "saveVTK: Element type %s is not supported by VTK",elements->ReferenceElement->Type->Name);
      Finley_setError(VALUE_ERROR,error_msg);
      fclose(fileHandle_p);
      return;
  } 
  /* xml header */
  fprintf(fileHandle_p, "<?xml version=\"1.0\"?>\n");
  fprintf(fileHandle_p, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n");

  /* finley uses an unstructured mesh, so UnstructuredGrid *should* work */
  fprintf(fileHandle_p, "<UnstructuredGrid>\n");

  /* is there only one "piece" to the data?? */
  fprintf(fileHandle_p, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",numPoints, numCells);
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
  int nDim = mesh_p->Nodes->numDim;
  fprintf(fileHandle_p, "<DataArray NumberOfComponents=\"%d\" type=\"Float32\" format=\"ascii\">\n",MAX(3,nDim));
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
          for (j = 0; j < nDim; j++) fprintf(fileHandle_p, " %e",(float) (mesh_p->Nodes->Coordinates[INDEX2(j, i, nDim)]));
          for (k=0; k<3-nDim; k++) fprintf(fileHandle_p, " %e",(float) 0.);
          fprintf(fileHandle_p, "\n");
       }
     } 
  } else {
     for (i = 0; i < mesh_p->Nodes->numNodes; i++) {
       for (j = 0; j < nDim; j++) fprintf(fileHandle_p, " %e",(float) (mesh_p->Nodes->Coordinates[INDEX2(j, i, nDim)]));
       for (k=0; k<3-nDim; k++) fprintf(fileHandle_p, " %e",(float) 0.);
       fprintf(fileHandle_p, "\n");
     }
  } 
  fprintf(fileHandle_p, "</DataArray>\n");
  fprintf(fileHandle_p, "</Points>\n");

  /* write out the DataArray element for the connectivity */

  int NN = elements->ReferenceElement->Type->numNodes;
  fprintf(fileHandle_p, "<Cells>\n");
  fprintf(fileHandle_p, "<DataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\">\n");

  if (nodetype==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
     for (i = 0; i < numCells; i++) {
        for (j = 0; j < numVTKNodesPerElement; j++) 
             fprintf(fileHandle_p,"%d ",mesh_p->Nodes->toReduced[elements->Nodes[INDEX2(elements->ReferenceElement->Type->linearNodes[j], i, NN)]]); 
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
        for (j = 0; j < numVTKNodesPerElement; j++) fprintf(fileHandle_p,"%d ", elements->Nodes[INDEX2(elements->ReferenceElement->Type->geoNodes[j], i, NN)]);
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
  fprintf(fileHandle_p, "<DataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\">\n");
  for (i=numVTKNodesPerElement; i<=numCells*numVTKNodesPerElement; i+=numVTKNodesPerElement) fprintf(fileHandle_p, "%d\n", i);
  fprintf(fileHandle_p, "</DataArray>\n");

  /* write out the DataArray element for the types */
  fprintf(fileHandle_p, "<DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\">\n");
  for (i=0; i<numCells; i++) fprintf(fileHandle_p, "%d\n", cellType);
  fprintf(fileHandle_p, "</DataArray>\n");

  /* finish off the <Cells> element */
  fprintf(fileHandle_p, "</Cells>\n");

  /* cell data */
  if (write_celldata) {
       /* mark the active data arrays */
       bool_t set_scalar=FALSE,set_vector=FALSE, set_tensor=FALSE;
       fprintf(fileHandle_p, "<CellData");
       for (i_data =0 ;i_data<num_data;++i_data) {
            if (! isEmpty(data_pp[i_data]) && isCellCentered[i_data]) {
                /* if the rank == 0:   --> scalar data
                 * if the rank == 1:   --> vector data
                 * if the rank == 2:   --> tensor data
                 */
                switch(getDataPointRank(data_pp[i_data])) {
                   case 0:
                       if (! set_scalar) {
                             fprintf(fileHandle_p," Scalars=\"%s\"",names_p[i_data]);
                             set_scalar=TRUE;
                       }
                       break;
                   case 1:
                       if (! set_vector) {
                             fprintf(fileHandle_p," Vectors=\"%s\"",names_p[i_data]);
                             set_vector=TRUE;
                       }
                       break;
                   case 2:
                       if (! set_tensor) {
                             fprintf(fileHandle_p," Tensors=\"%s\"",names_p[i_data]);
                             set_tensor=TRUE;
                       }
                       break;
                   default:
                       sprintf(error_msg, "saveVTK: data %s: Vtk can't handle objects with rank greater than 2.",names_p[i_data]);
                       Finley_setError(VALUE_ERROR,error_msg);
                       fclose(fileHandle_p);
                       return;
                }
            }
       }
       fprintf(fileHandle_p, ">\n");
       /* write the arrays */
       for (i_data =0 ;i_data<num_data;++i_data) {
          if (! isEmpty(data_pp[i_data]) && isCellCentered[i_data]) {
             int numPointsPerSample = elements->ReferenceElement->numQuadNodes;
             int rank = getDataPointRank(data_pp[i_data]);
             int nComp = getDataPointSize(data_pp[i_data]);
             int nCompReqd=1;   /* the number of components required by vtk */
             int shape=0;
             if (rank == 0) {
                nCompReqd = 1;
             } else if (rank == 1) {
                 shape=getDataPointShape(data_pp[i_data], 0);
                 if  (shape>3) {
                     Finley_setError(VALUE_ERROR, "saveVTK: rank 1 object must have less then 4 components");
                     fclose(fileHandle_p);
                     return;
                 }
                 nCompReqd = 3;
             } else {
                 shape=getDataPointShape(data_pp[i_data], 0);
                 if  (shape>3 || shape != getDataPointShape(data_pp[i_data], 1)) {
                     Finley_setError(VALUE_ERROR, "saveVTK: rank 2 object must have less then 4x4 components and must have a square shape");
                     fclose(fileHandle_p);
                     return;
                 }
                 nCompReqd = 9;
             }
             fprintf(fileHandle_p, "<DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"%d\" format=\"ascii\">\n",names_p[i_data], nCompReqd);
   
   	     double sampleAvg[nComp];
   	     for (i=0; i<numCells; i++) {
   	        values = getSampleData(data_pp[i_data], i);
   	        /* averaging over the number of points in the sample */
   	        for (k=0; k<nComp; k++) {
   	           rtmp = 0.;
   	           for (j=0; j<numPointsPerSample; j++) rtmp += values[INDEX2(k,j,nComp)];
   	           sampleAvg[k] = rtmp/numPointsPerSample;
   	         }
   	         /* if the number of required components is more than the number
   	          * of actual components, pad with zeros
   	          */
   	         /* probably only need to get shape of first element */
   	         /* write the data different ways for scalar, vector and tensor */
   	         if (nCompReqd == 1) {
   	           fprintf(fileHandle_p, " %e", (float) sampleAvg[0]);
   	         } else if (nCompReqd == 3) {
   	           /* write out the data */
   	           for (int m=0; m<shape; m++) fprintf(fileHandle_p, " %e", (float) sampleAvg[m]);
   	           for (int m=0; m<nCompReqd-shape; m++) fprintf(fileHandle_p, " %e", (float) 0.);
   	         } else if (nCompReqd == 9) { 
   	           /* tensor data, so have a 3x3 matrix to output as a row 
   	            * of 9 data points */
   	            int count = 0;
   	            for (int m=0; m<shape; m++) {
   	              for (int n=0; n<shape; n++) {
   	                 fprintf(fileHandle_p, " %e", (float) sampleAvg[count]);
   	                 count++;
   	              }
   	              for (int n=0; n<3-shape; n++) fprintf(fileHandle_p, " %e", (float) 0.);
   	            }
   	            for (int m=0; m<3-shape; m++) 
   	               for (int n=0; n<3; n++) fprintf(fileHandle_p, " %e", (float) 0.);
   	            }
   	          fprintf(fileHandle_p, "\n");
             }
             fprintf(fileHandle_p, "</DataArray>\n");
         }
       }
       fprintf(fileHandle_p, "</CellData>\n");
  }
  /* point data */
  if (write_pointdata) {
       /* mark the active data arrays */
       bool_t set_scalar=FALSE,set_vector=FALSE, set_tensor=FALSE;
       fprintf(fileHandle_p, "<PointData");
       for (i_data =0 ;i_data<num_data;++i_data) {
            if (! isEmpty(data_pp[i_data]) && !isCellCentered[i_data]) {
                /* if the rank == 0:   --> scalar data
                 * if the rank == 1:   --> vector data
                 * if the rank == 2:   --> tensor data
                 */
                switch(getDataPointRank(data_pp[i_data])) {
                   case 0:
                       if (! set_scalar) {
                             fprintf(fileHandle_p," Scalars=\"%s\"",names_p[i_data]);
                             set_scalar=TRUE;
                       }
                       break;
                   case 1:
                       if (! set_vector) {
                             fprintf(fileHandle_p," Vectors=\"%s\"",names_p[i_data]);
                             set_vector=TRUE;
                       }
                       break;
                   case 2:
                       if (! set_tensor) {
                             fprintf(fileHandle_p," Tensors=\"%s\"",names_p[i_data]);
                             set_tensor=TRUE;
                       }
                       break;
                   default:
                       sprintf(error_msg, "saveVTK: data %s: Vtk can't handle objects with rank greater than 2.",names_p[i_data]);
                       Finley_setError(VALUE_ERROR,error_msg);
                       fclose(fileHandle_p);
                       return;
                }
            }
       }
       fprintf(fileHandle_p, ">\n");
       /* write the arrays */
       for (i_data =0 ;i_data<num_data;++i_data) {
          if (! isEmpty(data_pp[i_data]) && !isCellCentered[i_data]) {
             int numPointsPerSample = elements->ReferenceElement->numQuadNodes;
             int rank = getDataPointRank(data_pp[i_data]);
             int nComp = getDataPointSize(data_pp[i_data]);
             int nCompReqd=1;   /* the number of components required by vtk */
             int shape=0;
             if (rank == 0) {
                nCompReqd = 1;
             } else if (rank == 1) {
                 shape=getDataPointShape(data_pp[i_data], 0);
                 if  (shape>3) {
                     Finley_setError(VALUE_ERROR, "saveVTK: rank 1 object must have less then 4 components");
                     fclose(fileHandle_p);
                     return;
                 }
                 nCompReqd = 3;
             } else {
                 shape=getDataPointShape(data_pp[i_data], 0);
                 if  (shape>3 || shape != getDataPointShape(data_pp[i_data], 1)) {
                     Finley_setError(VALUE_ERROR, "saveVTK: rank 2 object must have less then 4x4 components and must have a square shape");
                     fclose(fileHandle_p);
                     return;
                 }
                 nCompReqd = 9;
             }
             fprintf(fileHandle_p, "<DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"%d\" format=\"ascii\">\n",names_p[i_data], nCompReqd);
   	     /* write out the data */
   	     /* if the number of required components is more than the number
   	      * of actual components, pad with zeros
               */
             bool_t do_write=TRUE;
             for (i=0; i<mesh_p->Nodes->numNodes; i++) {
               if (nodetype==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
   	            if (mesh_p->Nodes->toReduced[i]>=0) {
                       switch(getFunctionSpaceType(data_pp[i_data])) {
                          case FINLEY_DEGREES_OF_FREEDOM:
                             values = getSampleData(data_pp[i_data],mesh_p->Nodes->degreeOfFreedom[i]);
                             break;
                          case FINLEY_REDUCED_DEGREES_OF_FREEDOM:
                             values = getSampleData(data_pp[i_data],mesh_p->Nodes->reducedDegreeOfFreedom[i]);
                             break;
                          case FINLEY_NODES:
   	                     values = getSampleData(data_pp[i_data],i);
                             break;
                       }
                       do_write=TRUE;
                    } else {
                       do_write=FALSE;
                    }
               } else {
                    do_write=TRUE;
                    switch(getFunctionSpaceType(data_pp[i_data])) {
                       case FINLEY_DEGREES_OF_FREEDOM:
   	                  values = getSampleData(data_pp[i_data],mesh_p->Nodes->degreeOfFreedom[i]);
                          break;
                       case FINLEY_NODES:
   	                  values = getSampleData(data_pp[i_data],i);
                          break;
                    }
               }
   	       /* write the data different ways for scalar, vector and tensor */
               if (do_write) {
   	          if (nCompReqd == 1) {
   	            fprintf(fileHandle_p, " %e", (float) values[0]);
   	          } else if (nCompReqd == 3) {
   	            for (int m=0; m<shape; m++) fprintf(fileHandle_p, " %e", (float) values[m]);
   	            for (int m=0; m<nCompReqd-shape; m++) fprintf(fileHandle_p, " %e", (float) 0.);
   	          } else if (nCompReqd == 9) { 
   	            /* tensor data, so have a 3x3 matrix to output as a row 
   	             * of 9 data points */
   	            int count = 0;
   	            for (int m=0; m<shape; m++) {
   	              for (int n=0; n<shape; n++) {
   	                fprintf(fileHandle_p, " %e", (float) values[count]);
   	                count++;
   	              }
   	              for (int n=0; n<3-shape; n++) fprintf(fileHandle_p, " %e", (float) 0.);
   	            }
   	            for (int m=0; m<3-shape; m++)  
   	                for (int n=0; n<3; n++) fprintf(fileHandle_p, " %e", (float) 0.);
   	          }
      	          fprintf(fileHandle_p, "\n");
               }
             }
             fprintf(fileHandle_p, "</DataArray>\n");
           }
       }
       fprintf(fileHandle_p, "</PointData>\n");
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
