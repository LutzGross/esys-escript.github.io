/**************************************************************/

/*   Finley: Mesh */

/*   takes nodes elements, etc of in2 and copies them into in1 */
/*   Ids of in2 are shifted by the maximum Id of in1 */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003/04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Finley.h"
#include "Mesh.h"
#include "Util.h"

/**************************************************************/

static double  Finley_Mesh_lockingGridSize=0;

Finley_Mesh* Finley_Mesh_merge(int numMsh, Finley_Mesh** msh) {
  Finley_Mesh* out=NULL;
  int i;
  char newName[LenString_MAX];
  if (numMsh==0) {
     Finley_ErrorCode=VALUE_ERROR;
     sprintf(Finley_ErrorMsg,"Empty mesh list");
  } else {
    int order=msh[0]->order;
    int numDim=msh[0]->Nodes->numDim;
    ElementTypeId elementTypeId=NoType;
    ElementTypeId faceElementTypeId=NoType;
    ElementTypeId pointTypeId=NoType;
    ElementTypeId contactTypeId=NoType;
    strcpy(newName,"");
    int numNodes=0;
    int numElements=0;
    int numFaceElements=0;
    int numContactElements=0;
    int numPoints=0;
    for (i=0;i<numMsh;i++) {
       /* check if all mesh have the same type and dimensions */
       order=MAX(order,msh[i]->order);
       numNodes+=msh[i]->Nodes->numNodes;
       if (numDim!=msh[i]->Nodes->numDim) {
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"Spatial dimensions of meshes don't match.");
       }

       if (msh[i]->Elements!=NULL) {
          numElements+=msh[i]->Elements->numElements;
	  if (elementTypeId==NoType) {
             elementTypeId=msh[i]->Elements->ReferenceElement->Type->TypeId;
	  } else {
             if (elementTypeId!=msh[i]->Elements->ReferenceElement->Type->TypeId ) {
               Finley_ErrorCode=TYPE_ERROR;
               sprintf(Finley_ErrorMsg,"element types of meshes don't match.");
             }
          }
       }

       if (msh[i]->FaceElements!=NULL) {
          numFaceElements+=msh[i]->FaceElements->numElements;
	  if (faceElementTypeId==NoType) {
             faceElementTypeId=msh[i]->FaceElements->ReferenceElement->Type->TypeId;
	  } else {
             if (faceElementTypeId!=msh[i]->FaceElements->ReferenceElement->Type->TypeId ) {
               Finley_ErrorCode=TYPE_ERROR;
               sprintf(Finley_ErrorMsg,"face element types of meshes don't match.");
             }
          }
       }

       if (msh[i]->ContactElements!=NULL) {
          numContactElements+=msh[i]->ContactElements->numElements;
	  if (contactTypeId==NoType) {
             contactTypeId=msh[i]->ContactElements->ReferenceElement->Type->TypeId;
	  } else {
             if (contactTypeId!=msh[i]->ContactElements->ReferenceElement->Type->TypeId ) {
               Finley_ErrorCode=TYPE_ERROR;
               sprintf(Finley_ErrorMsg,"contact element types of meshes don't match.");
             }
          }
       }

       if (msh[i]->Points!=NULL) {
          numPoints+=msh[i]->Points->numElements;
	  if (pointTypeId==NoType) {
             pointTypeId=msh[i]->Points->ReferenceElement->Type->TypeId;
	  } else {
             if (pointTypeId!=msh[i]->Points->ReferenceElement->Type->TypeId ) {
               Finley_ErrorCode=TYPE_ERROR;
               sprintf(Finley_ErrorMsg,"point element types of meshes don't match.");
             }
          }
       }

       strncat(newName,"+",LenString_MAX-strlen(newName));
       strncat(newName,msh[i]->Name,LenString_MAX-strlen(newName)-1);
    }

    /* allocate */

    if (Finley_ErrorCode == NO_ERROR) out=Finley_Mesh_alloc(newName,numDim,order);

    out->Elements=Finley_ElementFile_alloc(elementTypeId,out->order);
    out->FaceElements=Finley_ElementFile_alloc(faceElementTypeId,out->order);
    out->Points=Finley_ElementFile_alloc(pointTypeId,out->order);
    out->ContactElements=Finley_ElementFile_alloc(contactTypeId,out->order);

    /* allocate new tables */

    if (Finley_ErrorCode == NO_ERROR) {
        Finley_NodeFile_allocTable(out->Nodes,numNodes);
        Finley_ElementFile_allocTable(out->Elements,numElements);
        Finley_ElementFile_allocTable(out->FaceElements,numFaceElements);
        Finley_ElementFile_allocTable(out->ContactElements,numContactElements);
        Finley_ElementFile_allocTable(out->Points,numPoints);
    }

    /* copy tables :*/

    if (Finley_ErrorCode == NO_ERROR) {

        int numNodes=0;
        int numElements=0;
        int numFaceElements=0;
        int numContactElements=0;
        int numPoints=0;
        int maxNodeID=0;
        int maxDOF=0;
        int maxElementID=0;
        int maxElementID2=0;

        for (i=0;i<numMsh;i++) {

           Finley_NodeFile_copyTable(numNodes,out->Nodes,maxNodeID,maxDOF,msh[i]->Nodes); 
           Finley_ElementFile_copyTable(numElements,out->Elements,numNodes,maxElementID,msh[i]->Elements); 
           Finley_ElementFile_copyTable(numFaceElements,out->FaceElements,numNodes,maxElementID,msh[i]->FaceElements); 
           Finley_ElementFile_copyTable(numContactElements,out->ContactElements,numNodes,maxElementID,msh[i]->ContactElements); 
           Finley_ElementFile_copyTable(numPoints,out->Points,numNodes,maxElementID,msh[i]->Points); 

           numNodes=+msh[i]->Nodes->numNodes;
           numElements=+msh[i]->Elements->numElements;
           numFaceElements=+msh[i]->FaceElements->numElements;
           numContactElements=+msh[i]->ContactElements->numElements;
           numPoints=+msh[i]->Points->numElements;

           if (msh[i]->Nodes->numNodes>0) 
              maxNodeID+=Finley_Util_getMaxInt(1,msh[i]->Nodes->numNodes,msh[i]->Nodes->Id)+1;
              maxDOF+=Finley_Util_getMaxInt(1,msh[i]->Nodes->numNodes,msh[i]->Nodes->degreeOfFreedom)+1;
           maxElementID2=0;
           if (msh[i]->Elements->numElements>0) 
              maxElementID2=MAX(maxElementID2,Finley_Util_getMaxInt(1,msh[i]->Elements->numElements,msh[i]->Elements->Id));
           if (msh[i]->FaceElements->numElements>0) 
              maxElementID2=MAX(maxElementID2,Finley_Util_getMaxInt(1,msh[i]->FaceElements->numElements,msh[i]->FaceElements->Id));
           if (msh[i]->ContactElements->numElements>0) 
              maxElementID2=MAX(maxElementID2,Finley_Util_getMaxInt(1,msh[i]->ContactElements->numElements,msh[i]->ContactElements->Id));
           if (msh[i]->Points->numElements) 
              maxElementID2=MAX(maxElementID2,Finley_Util_getMaxInt(1,msh[i]->Points->numElements,msh[i]->Points->Id));
           maxElementID+=maxElementID2+1;
        }
    }
    /* all done  */

    if (Finley_ErrorCode != NO_ERROR) {
       Finley_Mesh_dealloc(out);
    } else {
       Finley_Mesh_prepare(out);
       #ifdef Finley_TRACE
       printf("%d meshes merged.\n",numMsh);
       #endif
    }
  }
  return out;
}

/*
* $Log$
* Revision 1.1  2004/10/26 06:53:57  jgs
* Initial revision
*
* Revision 1.2  2004/07/30 04:37:06  gross
* escript and finley are linking now and RecMeshTest.py has been passed
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/

