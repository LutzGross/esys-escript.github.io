
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

/*   Finley: ElementFile */

/*   mark the used nodes with offeset: */

/**************************************************************/

#include "ElementFile.h"

/**************************************************************/

void Finley_ElementFile_markNodes(index_t* mask,index_t offset,dim_t numNodes,Finley_ElementFile* in,bool_t useLinear) {
   dim_t i,NN,NN2,e;
   index_t *lin_nodes;
   Finley_ReferenceElement* refElement=NULL;
   
   if (in!=NULL) {
	     refElement=Finley_ReferenceElementSet_borrowReferenceElement(in->referenceElementSet, FALSE);	   
		 NN2=in->numNodes;
		 if (useLinear) {
			 NN=refElement->numLinearNodes;
			 lin_nodes=refElement->Type->linearNodes;
			 #pragma omp parallel for private(e,i) schedule(static)
			 for (e=0;e<in->numElements;e++) {
				 for (i=0;i<NN;i++) {
					   mask[in->Nodes[INDEX2(lin_nodes[i],e,NN2)]-offset]=1;
				 }
			 }
		 } else {
			 NN=refElement->Type->numNodes;
			 #pragma omp parallel for private(e,i) schedule(static)
			 for (e=0;e<in->numElements;e++) {
				 for (i=0;i<NN;i++) {
					 mask[in->Nodes[INDEX2(i,e,NN2)]-offset]=1;
				 }
			 }
		 }
   	}
}

void Finley_ElementFile_markDOFsConnectedToRange(index_t* mask,index_t offset,index_t marker,index_t firstDOF,index_t lastDOF,index_t *dofIndex,Finley_ElementFile*in ,bool_t useLinear) 
{
   dim_t i,NN,NN2,e,j;
   index_t color,*lin_nodes;
   Finley_ReferenceElement* refElement=NULL;
   register index_t k;
   
   if (in!=NULL) {
	     refElement=Finley_ReferenceElementSet_borrowReferenceElement(in->referenceElementSet, FALSE);	   
		 NN2=in->numNodes;
		 if (useLinear) {
			 NN=refElement->numLinearNodes;
			 lin_nodes=refElement->Type->linearNodes;
			 for (color=in->minColor;color<=in->maxColor;color++) {
				 #pragma omp parallel for private(e,i,j,k) schedule(static)
				 for (e=0;e<in->numElements;e++) {
					 if (in->Color[e]==color) {
						 for (i=0;i<NN;i++) {
							 k=dofIndex[in->Nodes[INDEX2(lin_nodes[i],e,NN2)]];
							 if ( (firstDOF<=k) && (k<lastDOF) ) {
								 for (j=0;j<NN;j++) mask[dofIndex[in->Nodes[INDEX2(lin_nodes[j],e,NN2)]]-offset]=marker;
								 break;
							 }
						 }
					 }
				 }
			 }
		 } else {
			 NN=refElement->Type->numNodes;
			 for (color=in->minColor;color<=in->maxColor;color++) {
				 #pragma omp parallel for private(e,i,j,k) schedule(static)
				 for (e=0;e<in->numElements;e++) {
					 if (in->Color[e]==color) {
						 for (i=0;i<NN;i++) {
							 k=dofIndex[in->Nodes[INDEX2(i,e,NN2)]];
							 if ( (firstDOF<=k) && (k<lastDOF) ) {
								 for (j=0;j<NN;j++) mask[dofIndex[in->Nodes[INDEX2(j,e,NN2)]]-offset]=marker;
								 break;
							 }
						 }
					 }
				 }
			 }
		 }
   	}	   
}

