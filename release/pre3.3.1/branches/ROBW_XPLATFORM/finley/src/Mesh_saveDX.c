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

/*   writes data and mesh in an opendx file */

/**************************************************************/

/*   Author: Lutz Gross, gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Finley_Mesh_saveDX(const char * filename_p, Finley_Mesh *mesh_p, const dim_t num_data,char* *names_p,escriptDataC* *data_pp) {
  char error_msg[LenErrorMsg_MAX];
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
  int i,j,k,i_data;
  /* open the file  and check handel */

  FILE * fileHandle_p = fopen(filename_p, "w");
  if (fileHandle_p==NULL) {
    sprintf(error_msg,"File %s could not be opened for writing.",filename_p);
    Finley_setError(IO_ERROR,error_msg);
    return;
  }
  /* find the mesh type to be written */
  int nodetype=FINLEY_DEGREES_OF_FREEDOM;
  int elementtype=FINLEY_UNKNOWN;
  bool_t isCellCentered[num_data];
  for (i_data=0;i_data<num_data;++i_data) {
     if (! isEmpty(data_pp[i_data])) {
        switch(getFunctionSpaceType(data_pp[i_data])) {
           case FINLEY_DEGREES_OF_FREEDOM:
             nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS) {
                 elementtype=FINLEY_ELEMENTS;
             } else {
                 Finley_setError(TYPE_ERROR,"saveDX: cannot write given data in single file.");
                 return;
             }
             isCellCentered[i_data]=FALSE;
             break;
           case FINLEY_REDUCED_DEGREES_OF_FREEDOM:
             nodetype = FINLEY_REDUCED_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS) {
                 elementtype=FINLEY_ELEMENTS;
             } else {
                 Finley_setError(TYPE_ERROR,"saveDX: cannot write given data in single file.");
                 return;
             }
             isCellCentered[i_data]=FALSE;
             break;
           case FINLEY_NODES:
             nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS) {
                 elementtype=FINLEY_ELEMENTS;
             } else {
                 Finley_setError(TYPE_ERROR,"saveDX: cannot write given data in single file.");
                 return;
             }
             isCellCentered[i_data]=FALSE;
             break;
           case FINLEY_ELEMENTS:
             nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_ELEMENTS) {
                 elementtype=FINLEY_ELEMENTS;
             } else {
                 Finley_setError(TYPE_ERROR,"saveDX: cannot write given data in single file.");
                 return;
             }
             isCellCentered[i_data]=TRUE;
             break;
           case FINLEY_FACE_ELEMENTS:
             nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_FACE_ELEMENTS) {
                 elementtype=FINLEY_FACE_ELEMENTS;
             } else {
                 Finley_setError(TYPE_ERROR,"saveDX: cannot write given data in single file.");
                 return;
             }
             isCellCentered[i_data]=TRUE;
             break;
           case FINLEY_POINTS:
             nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_POINTS) {
                 elementtype=FINLEY_POINTS;
             } else {
                 Finley_setError(TYPE_ERROR,"saveDX: cannot write given data in single file.");
                 return;
             }
             isCellCentered[i_data]=TRUE;
             break;
           case FINLEY_CONTACT_ELEMENTS_1:
             nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_CONTACT_ELEMENTS_1) {
                 elementtype=FINLEY_CONTACT_ELEMENTS_1;
             } else {
                 Finley_setError(TYPE_ERROR,"saveDX: cannot write given data in single file.");
                 return;
             }
             isCellCentered[i_data]=TRUE;
             break;
           case FINLEY_CONTACT_ELEMENTS_2:
             nodetype = (nodetype == FINLEY_REDUCED_DEGREES_OF_FREEDOM) ? FINLEY_REDUCED_DEGREES_OF_FREEDOM : FINLEY_DEGREES_OF_FREEDOM;
             if (elementtype==FINLEY_UNKNOWN || elementtype==FINLEY_CONTACT_ELEMENTS_1) {
                 elementtype=FINLEY_CONTACT_ELEMENTS_1;
             } else {
                 Finley_setError(TYPE_ERROR,"saveDX: cannot write given data in single file.");
                 return;
             }
             isCellCentered[i_data]=TRUE;
             break;
           default:
             sprintf(error_msg,"saveDX: Finley does not know anything about function space type %d",getFunctionSpaceType(data_pp[i_data]));
             Finley_setError(TYPE_ERROR,error_msg);
             return;
        }
     }
  }
  /* select number of points and the mesh component */
  int numPoints = mesh_p->Nodes->numNodes;
  int nDim = mesh_p->Nodes->numDim; 
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
     Finley_setError(SYSTEM_ERROR,"saveDX: undefined element file");
     return;
  }

  /* map finley element type to DX element type */
  ElementTypeId TypeId = elements->ReferenceElement->Type->TypeId;
  int *resortIndex=NULL;
  int numDXNodesPerElement=0;
  int numCells = elements->numElements;
  char elemTypeStr[32];
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
     sprintf(error_msg,"saveDX: Element type %s is not supported by DX",elements->ReferenceElement->Type->Name);
     Finley_setError(VALUE_ERROR,error_msg);
     return;
   } 

  /* positions */
  fprintf(fileHandle_p, "object 1 class array type float rank 1 shape %d items %d data follows\n",nDim, mesh_p->Nodes->reducedNumNodes);
  for (i = 0; i < mesh_p->Nodes->numNodes; i++) {
    if (mesh_p->Nodes->toReduced[i]>=0) {
       for (j = 0; j < nDim; j++) fprintf(fileHandle_p, " %g",mesh_p->Nodes->Coordinates[INDEX2(j, i, nDim)]);
       fprintf(fileHandle_p, "\n");
    }
  } 
  /* connection table */
  int NN=elements->ReferenceElement->Type->numNodes;
  fprintf(fileHandle_p, "object 2 class array type int rank 1 shape %d items %d data follows\n",numDXNodesPerElement, numCells);
  for (i = 0; i < numCells; i++) {
      for (j = 0; j < numDXNodesPerElement; j++) fprintf(fileHandle_p," %d",mesh_p->Nodes->toReduced[elements->Nodes[INDEX2(resortIndex[j], i, NN)]]);
      fprintf(fileHandle_p, "\n");
  } 
  fprintf(fileHandle_p, "attribute \"element type\" string \"%s\"\n",elemTypeStr);
  fprintf(fileHandle_p, "attribute \"ref\" string \"positions\"\n");

  /* data */
  int object_count=2;
  for (i_data =0 ;i_data<num_data;++i_data) {
      if (! isEmpty(data_pp[i_data])) {
         object_count++;
         int rank=getDataPointRank(data_pp[i_data]);
         int nComp=getDataPointSize(data_pp[i_data]);
         double* values,rtmp;
         fprintf(fileHandle_p, "object %d class array type float rank %d ",object_count,rank);
         if (0 < rank) {
            fprintf(fileHandle_p, "shape ");
            for (i = 0; i < rank; i++) fprintf(fileHandle_p, "%d ", getDataPointShape(data_pp[i_data],i));
         }
         if (isCellCentered[i_data]) {
             int numPointsPerSample=elements->ReferenceElement->numQuadNodes;
             if (numPointsPerSample>0) {
                fprintf(fileHandle_p, "items %d data follows\n", numCells);
                for (i=0;i<elements->numElements;i++) {
                    values=getSampleData(data_pp[i_data],i);
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
                    switch (getFunctionSpaceType(data_pp[i_data])) {
                       case FINLEY_DEGREES_OF_FREEDOM:
                          values=getSampleData(data_pp[i_data],mesh_p->Nodes->degreeOfFreedom[i]);
                          break;
                       case FINLEY_REDUCED_DEGREES_OF_FREEDOM:
                          values=getSampleData(data_pp[i_data],mesh_p->Nodes->reducedDegreeOfFreedom[i]);
                          break;
                       case FINLEY_NODES:
                          values=getSampleData(data_pp[i_data],i);
                          break;
                    }
                    for (k=0;k<nComp;k++) fprintf(fileHandle_p, " %g", values[k]);
	            fprintf(fileHandle_p, "\n");
                 }
             }
             fprintf(fileHandle_p, "attribute \"dep\" string \"positions\"\n");
         }
     }
  }

  /* and finish it up */
  if (num_data==0) {
     fprintf(fileHandle_p, "object %d class field\n",object_count+1);
     fprintf(fileHandle_p, "component \"positions\" value 1\n");
     fprintf(fileHandle_p, "component \"connections\" value 2\n");
  } else {
     object_count=2;
     for (i_data=0; i_data<num_data;++i_data) {
         if (! isEmpty(data_pp[i_data])) {
            object_count++;
            fprintf(fileHandle_p, "object \"%s\" class field\n",names_p[i_data]);
            fprintf(fileHandle_p, "component \"positions\" value 1\n");
            fprintf(fileHandle_p, "component \"connections\" value 2\n");
            fprintf(fileHandle_p, "component \"data\" value %d\n",object_count);
         }
     }
  }
  fprintf(fileHandle_p, "end\n");
  /* close the file */
  fclose(fileHandle_p);
  return;
}
