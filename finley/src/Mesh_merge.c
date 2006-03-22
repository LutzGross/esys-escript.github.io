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

/*   Finley: Mesh */

/*   takes nodes elements, etc of in2 and copies them into in1 */
/*   Ids of in2 are shifted by the maximum Id of in1 */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "Mesh.h"
#include "Util.h"

/**************************************************************/

static double  Finley_Mesh_lockingGridSize=0;

Finley_Mesh* Finley_Mesh_merge(dim_t numMsh, Finley_Mesh** msh) {
  Finley_Mesh* out=NULL;
  dim_t i;
  char newName[LenString_MAX];
  if (numMsh==0) {
     Finley_setError(VALUE_ERROR,"__FILE__: Empty mesh list");
  } else {
    index_t order=msh[0]->order;
    dim_t numDim=msh[0]->Nodes->numDim;
    ElementTypeId elementTypeId=NoType;
    ElementTypeId faceElementTypeId=NoType;
    ElementTypeId pointTypeId=NoType;
    ElementTypeId contactTypeId=NoType;
    strcpy(newName,"");
    dim_t numNodes=0;
    dim_t numElements=0;
    dim_t numFaceElements=0;
    dim_t numContactElements=0;
    dim_t numPoints=0;
    for (i=0;i<numMsh;i++) {
       /* check if all mesh have the same type and dimensions */
       order=MAX(order,msh[i]->order);
       numNodes+=msh[i]->Nodes->numNodes;
       if (numDim!=msh[i]->Nodes->numDim) {
          Finley_setError(TYPE_ERROR,"__FILE__: Spatial dimensions of meshes don't match.");
       }

       if (msh[i]->Elements!=NULL) {
          numElements+=msh[i]->Elements->numElements;
	  if (elementTypeId==NoType) {
             elementTypeId=msh[i]->Elements->ReferenceElement->Type->TypeId;
	  } else {
             if (elementTypeId!=msh[i]->Elements->ReferenceElement->Type->TypeId ) {
               Finley_setError(TYPE_ERROR,"__FILE__: element types of meshes don't match.");
             }
          }
       }

       if (msh[i]->FaceElements!=NULL) {
          numFaceElements+=msh[i]->FaceElements->numElements;
	  if (faceElementTypeId==NoType) {
             faceElementTypeId=msh[i]->FaceElements->ReferenceElement->Type->TypeId;
	  } else {
             if (faceElementTypeId!=msh[i]->FaceElements->ReferenceElement->Type->TypeId ) {
               Finley_setError(TYPE_ERROR,"__FILE__: face element types of meshes don't match.");
             }
          }
       }

       if (msh[i]->ContactElements!=NULL) {
          numContactElements+=msh[i]->ContactElements->numElements;
	  if (contactTypeId==NoType) {
             contactTypeId=msh[i]->ContactElements->ReferenceElement->Type->TypeId;
	  } else {
             if (contactTypeId!=msh[i]->ContactElements->ReferenceElement->Type->TypeId ) {
               Finley_setError(TYPE_ERROR,"__FILE__: contact element types of meshes don't match.");
             }
          }
       }

       if (msh[i]->Points!=NULL) {
          numPoints+=msh[i]->Points->numElements;
	  if (pointTypeId==NoType) {
             pointTypeId=msh[i]->Points->ReferenceElement->Type->TypeId;
	  } else {
             if (pointTypeId!=msh[i]->Points->ReferenceElement->Type->TypeId ) {
               Finley_setError(TYPE_ERROR,"__FILE__: point element types of meshes don't match.");
             }
          }
       }

       strncat(newName,"+",LenString_MAX-strlen(newName));
       strncat(newName,msh[i]->Name,LenString_MAX-strlen(newName)-1);
    }

    /* allocate */

    if (Finley_noError()) out=Finley_Mesh_alloc(newName,numDim,order);

    out->Elements=Finley_ElementFile_alloc(elementTypeId,out->order);
    out->FaceElements=Finley_ElementFile_alloc(faceElementTypeId,out->order);
    out->Points=Finley_ElementFile_alloc(pointTypeId,out->order);
    out->ContactElements=Finley_ElementFile_alloc(contactTypeId,out->order);

    /* allocate new tables */

    if (Finley_noError()) {
        Finley_NodeFile_allocTable(out->Nodes,numNodes);
        Finley_ElementFile_allocTable(out->Elements,numElements);
        Finley_ElementFile_allocTable(out->FaceElements,numFaceElements);
        Finley_ElementFile_allocTable(out->ContactElements,numContactElements);
        Finley_ElementFile_allocTable(out->Points,numPoints);
    }

    /* copy tables :*/

    if (Finley_noError()) {

        dim_t numNodes=0;
        dim_t numElements=0;
        dim_t numFaceElements=0;
        dim_t numContactElements=0;
        dim_t numPoints=0;
        index_t maxNodeID=0;
        index_t maxDOF=0;
        index_t maxElementID=0;
        index_t maxElementID2=0;

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

    if (! Finley_noError()) {
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
* Revision 1.3  2005/09/15 03:44:22  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.2.2.1  2005/09/07 06:26:19  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.2  2005/07/08 04:07:53  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.1.1.1.2.1  2005/06/29 02:34:52  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.2  2004/07/30 04:37:06  gross
* escript and finley are linking now and RecMeshTest.py has been passed
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/

