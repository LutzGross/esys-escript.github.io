
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*	 Dudley: Mesh */

/*	 takes nodes elements, etc of in2 and copies them into in1 */
/*	 Ids of in2 are shifted by the maximum Id of in1 */

/**************************************************************/

#include "Mesh.h"
#include "Util.h"

/**************************************************************/

/*static double	 Dudley_Mesh_lockingGridSize=0;*/

Dudley_Mesh* Dudley_Mesh_merge(dim_t numMsh, Dudley_Mesh** msh) {
  Paso_MPIInfo *mpi_info=NULL;
  Dudley_ReferenceElementSet *refPoints=NULL, *refContactElements=NULL, *refFaceElements=NULL, *refElements=NULL;
  Dudley_Mesh* out=NULL;
  dim_t numNodes=0;
  dim_t numElements=0;
  dim_t numFaceElements=0;
  dim_t numContactElements=0;
  dim_t numPoints=0;
  dim_t i;
  index_t order, reduced_order;
  dim_t numDim;
  ElementTypeId elementTypeId=NoRef ;
  ElementTypeId faceElementTypeId=NoRef ;
  ElementTypeId pointTypeId=NoRef ;
  ElementTypeId contactTypeId=NoRef ;
  index_t maxNodeID=0;
  index_t maxDOF=0;
  index_t maxElementID=0;
  index_t maxElementID2=0;
  char newName[LenString_MAX];
  if (numMsh==0) {
	 Dudley_setError(VALUE_ERROR,"Dudley_Mesh_merge: Empty mesh list");
  } else {
	for	 (i=0;i<numMsh;i++) {
		 if (msh[i]->MPIInfo->size > 1) {
			  Dudley_setError(TYPE_ERROR,"Dudley_Mesh_merge: more than processor is not supported yet.");
			  return NULL;
		 }
	}
	order=msh[0]->integrationOrder;
	reduced_order=msh[0]->reducedIntegrationOrder;
	numDim=msh[0]->Nodes->numDim;
	mpi_info=msh[0]->MPIInfo;
	strcpy(newName,"");
	for (i=0;i<numMsh;i++) {
	   /* check if all mesh have the same type and dimensions */
	   order=MAX(order,msh[i]->integrationOrder);
	   reduced_order=MIN(reduced_order,msh[i]->reducedIntegrationOrder);
	   numNodes+=msh[i]->Nodes->numNodes;
	   if (mpi_info->comm!=msh[i]->MPIInfo->comm) {
		  Dudley_setError(TYPE_ERROR,"Dudley_Mesh_merge: MPI communicators of meshes don't match.");
	   }
	   if (numDim!=msh[i]->Nodes->numDim) {
		  Dudley_setError(TYPE_ERROR,"Dudley_Mesh_merge: Spatial dimensions of meshes don't match.");
	   }

	   if (msh[i]->Elements!=NULL) {
		  numElements+=msh[i]->Elements->numElements;
		  if (elementTypeId==NoRef ) {
			 elementTypeId=msh[i]->Elements->referenceElementSet->referenceElement->Type->TypeId;
		  } else {
			 if (elementTypeId!=msh[i]->Elements->referenceElementSet->referenceElement->Type->TypeId ) {
			   Dudley_setError(TYPE_ERROR,"Dudley_Mesh_merge: element types of meshes don't match.");
			 }
		  }
	   }

	   if (msh[i]->FaceElements!=NULL) {
		  numFaceElements+=msh[i]->FaceElements->numElements;
		  if (faceElementTypeId==NoRef ) {
			 faceElementTypeId=msh[i]->FaceElements->referenceElementSet->referenceElement->Type->TypeId;
		  } else {
			 if (faceElementTypeId!=msh[i]->FaceElements->referenceElementSet->referenceElement->Type->TypeId ) {
			   Dudley_setError(TYPE_ERROR,"Dudley_Mesh_merge: face element types of meshes don't match.");
			 }
		  }
	   }

	   if (msh[i]->ContactElements!=NULL) {
		  numContactElements+=msh[i]->ContactElements->numElements;
		  if (contactTypeId==NoRef ) {
			 contactTypeId=msh[i]->ContactElements->referenceElementSet->referenceElement->Type->TypeId;
		  } else {
			 if (contactTypeId!=msh[i]->ContactElements->referenceElementSet->referenceElement->Type->TypeId ) {
			   Dudley_setError(TYPE_ERROR,"Dudley_Mesh_merge: contact element types of meshes don't match.");
			 }
		  }
	   }

	   if (msh[i]->Points!=NULL) {
		  numPoints+=msh[i]->Points->numElements;
		  if (pointTypeId==NoRef ) {
			 pointTypeId=msh[i]->Points->referenceElementSet->referenceElement->Type->TypeId;
		  } else {
			 if (pointTypeId!=msh[i]->Points->referenceElementSet->referenceElement->Type->TypeId ) {
			   Dudley_setError(TYPE_ERROR,"Dudley_Mesh_merge: point element types of meshes don't match.");
			 }
		  }
	   }

	   strncat(newName,"+",LenString_MAX-strlen(newName));
	   strncat(newName,msh[i]->Name,LenString_MAX-strlen(newName)-1);
	}

	if (mpi_info->size >1 ) {
		Dudley_setError(TYPE_ERROR,"Dudley_Mesh_merge: only single processor runs are supported.");
	}
	/* allocate */

	if (Dudley_noError()) {
	  out=Dudley_Mesh_alloc(newName,numDim,mpi_info);
    }
	if (Dudley_noError()) {
		refElements= Dudley_ReferenceElementSet_alloc(elementTypeId,order,reduced_order);
	  	refFaceElements=Dudley_ReferenceElementSet_alloc(faceElementTypeId, order,reduced_order);
		refContactElements=Dudley_ReferenceElementSet_alloc(contactTypeId, order,reduced_order);
		refPoints=Dudley_ReferenceElementSet_alloc(pointTypeId, order,reduced_order);
	}
	if (Dudley_noError()) {
	  out->Elements=Dudley_ElementFile_alloc(refElements,mpi_info);
	  out->FaceElements=Dudley_ElementFile_alloc(refFaceElements,mpi_info);
	  out->Points=Dudley_ElementFile_alloc(refPoints,mpi_info);
	  out->ContactElements=Dudley_ElementFile_alloc(refContactElements,mpi_info);

	}
	/* allocate new tables */
	if (Dudley_noError()) {
		Dudley_NodeFile_allocTable(out->Nodes,numNodes);
		Dudley_ElementFile_allocTable(out->Elements,numElements);
		Dudley_ElementFile_allocTable(out->FaceElements,numFaceElements);
		Dudley_ElementFile_allocTable(out->ContactElements,numContactElements);
		Dudley_ElementFile_allocTable(out->Points,numPoints);		
	}

	/* copy tables :*/

	if (Dudley_noError()) {
	   numNodes=0;
	   numElements=0;
	   numFaceElements=0;
	   numContactElements=0;
	   numPoints=0;

		for (i=0;i<numMsh;i++) {

		   Dudley_NodeFile_copyTable(numNodes,out->Nodes,maxNodeID,maxDOF,msh[i]->Nodes); 
		   Dudley_ElementFile_copyTable(numElements,out->Elements,numNodes,maxElementID,msh[i]->Elements); 
		   Dudley_ElementFile_copyTable(numFaceElements,out->FaceElements,numNodes,maxElementID,msh[i]->FaceElements); 
		   Dudley_ElementFile_copyTable(numContactElements,out->ContactElements,numNodes,maxElementID,msh[i]->ContactElements); 
		   Dudley_ElementFile_copyTable(numPoints,out->Points,numNodes,maxElementID,msh[i]->Points); 

		   numNodes=+msh[i]->Nodes->numNodes;
		   numElements=+msh[i]->Elements->numElements;
		   numFaceElements=+msh[i]->FaceElements->numElements;
		   numContactElements=+msh[i]->ContactElements->numElements;
		   numPoints=+msh[i]->Points->numElements;

		   if (msh[i]->Nodes->numNodes>0) 
			  maxNodeID+=Dudley_Util_getMaxInt(1,msh[i]->Nodes->numNodes,msh[i]->Nodes->Id)+1;
			  maxDOF+=Dudley_Util_getMaxInt(1,msh[i]->Nodes->numNodes,msh[i]->Nodes->globalDegreesOfFreedom)+1;
		   maxElementID2=0;
		   if (msh[i]->Elements->numElements>0) 
			  maxElementID2=MAX(maxElementID2,Dudley_Util_getMaxInt(1,msh[i]->Elements->numElements,msh[i]->Elements->Id));
		   if (msh[i]->FaceElements->numElements>0) 
			  maxElementID2=MAX(maxElementID2,Dudley_Util_getMaxInt(1,msh[i]->FaceElements->numElements,msh[i]->FaceElements->Id));
		   if (msh[i]->ContactElements->numElements>0) 
			  maxElementID2=MAX(maxElementID2,Dudley_Util_getMaxInt(1,msh[i]->ContactElements->numElements,msh[i]->ContactElements->Id));
		   if (msh[i]->Points->numElements) 
			  maxElementID2=MAX(maxElementID2,Dudley_Util_getMaxInt(1,msh[i]->Points->numElements,msh[i]->Points->Id));
		   maxElementID+=maxElementID2+1;
		}
	}
	/* all done	 */
	Dudley_ReferenceElementSet_dealloc(refPoints);
	Dudley_ReferenceElementSet_dealloc(refContactElements);
	Dudley_ReferenceElementSet_dealloc(refFaceElements);
	Dudley_ReferenceElementSet_dealloc(refElements);
	if (! Dudley_noError()) {
	   Dudley_Mesh_free(out);
	} else {
	   Dudley_Mesh_prepare(out, FALSE);
	}
  }
  return out;
}
