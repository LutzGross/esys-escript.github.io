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

/*   Finley: write Mesh */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

/*  writes the mesh to the external file fname unsing the Finley file format: */

void Finley_Mesh_write(Finley_Mesh *in,char* fname) {
  char error_msg[LenErrorMsg_MAX];
  FILE *f;
  int NN,i,j,numDim;
  Finley_TagMap* tag_map=in->TagMap;

  if (in->MPIInfo->size >1 ) {
    Finley_setError(IO_ERROR,"Mesh_write: only single processor runs are supported.");
    return;

  }
  /* open file */
  f=fopen(fname,"w");
  if (f==NULL) {
    sprintf(error_msg,"Mesh_write: Opening file %s for writing failed.",fname);
    Finley_setError(IO_ERROR,error_msg);
    return;
  }

  /* write header */

  fprintf(f,"%s\n",in->Name);
  
  /*  write nodes: */
  
  if (in->Nodes!=NULL) {
    numDim=Finley_Mesh_getDim(in);
    fprintf(f,"%1dD-Nodes %d\n", numDim, in->Nodes->numNodes);
    for (i=0;i<in->Nodes->numNodes;i++) {
      fprintf(f,"%d %d %d",in->Nodes->Id[i],in->Nodes->globalDegreesOfFreedom[i],in->Nodes->Tag[i]);
      for (j=0;j<numDim;j++) fprintf(f," %20.15e",in->Nodes->Coordinates[INDEX2(j,i,numDim)]);
      fprintf(f,"\n");
    }
  } else {
    fprintf(f,"0D-Nodes 0\n");
  }
  
  /*  write elements: */

  if (in->Elements!=NULL) {
    fprintf(f, "%s %d\n",in->Elements->ReferenceElement->Type->Name,in->Elements->numElements);
    NN=in->Elements->numNodes;
    for (i=0;i<in->Elements->numElements;i++) {
      fprintf(f,"%d %d",in->Elements->Id[i],in->Elements->Tag[i]);
      for (j=0;j<NN;j++) fprintf(f," %d",in->Nodes->Id[in->Elements->Nodes[INDEX2(j,i,NN)]]);
      fprintf(f,"\n");
    }
  } else {
    fprintf(f,"Tet4 0\n");
  }

  /*  write face elements: */
  if (in->FaceElements!=NULL) {
    fprintf(f, "%s %d\n", in->FaceElements->ReferenceElement->Type->Name,in->FaceElements->numElements);
    NN=in->FaceElements->numNodes;
    for (i=0;i<in->FaceElements->numElements;i++) {
      fprintf(f,"%d %d",in->FaceElements->Id[i],in->FaceElements->Tag[i]);
      for (j=0;j<NN;j++) fprintf(f," %d",in->Nodes->Id[in->FaceElements->Nodes[INDEX2(j,i,NN)]]);
      fprintf(f,"\n");
    }
  } else {
    fprintf(f,"Tri3 0\n");
  }

  /*  write Contact elements : */
  if (in->ContactElements!=NULL) {
    fprintf(f, "%s %d\n",in->ContactElements->ReferenceElement->Type->Name,in->ContactElements->numElements);
    NN=in->ContactElements->numNodes;
    for (i=0;i<in->ContactElements->numElements;i++) {
      fprintf(f,"%d %d",in->ContactElements->Id[i],in->ContactElements->Tag[i]);
      for (j=0;j<NN;j++) fprintf(f," %d",in->Nodes->Id[in->ContactElements->Nodes[INDEX2(j,i,NN)]]);
      fprintf(f,"\n");
    }
  } else {
    fprintf(f,"Tri3_Contact 0\n");
  }
  
  /*  write points: */
  if (in->Points!=NULL) {
    fprintf(f, "%s %d\n",in->Points->ReferenceElement->Type->Name,in->Points->numElements);
    for (i=0;i<in->Points->numElements;i++) {
      fprintf(f,"%d %d %d\n",in->Points->Id[i],in->Points->Tag[i],in->Nodes->Id[in->Points->Nodes[INDEX2(0,i,1)]]);
    }
  } else {
    fprintf(f,"Point1 0\n");
  }

  /*  write tags:*/
  if (tag_map) {
     fprintf(f,"Tags\n");
     while (tag_map) {
        fprintf(f,"%s %d\n",tag_map->name,tag_map->tag_key);
        tag_map=tag_map->next;
     }
  }
  fclose(f);
  #ifdef Finley_TRACE
  printf("mesh %s has been written to file %s\n",in->Name,fname);
  #endif
}

/*
* $Log$
* Revision 1.2  2005/09/15 03:44:23  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.1.1.1.6.1  2005/09/07 06:26:20  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/

