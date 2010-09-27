
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

/*   Dudley: ElementFile */

/*   mark the used nodes with offeset: */

/**************************************************************/

#include "ElementFile.h"

/**************************************************************/

void Dudley_ElementFile_markNodes(index_t* mask,index_t offset,dim_t numNodes,Dudley_ElementFile* in,bool_t useLinear) {
    dim_t i,NN,e;
//    Dudley_ReferenceElement* refElement=NULL;
   
    if (in!=NULL)
    {
//	refElement=Dudley_ReferenceElementSet_borrowReferenceElement(in->referenceElementSet, FALSE);	   
	NN=in->numNodes;
	#pragma omp parallel for private(e,i) schedule(static)
	for (e=0;e<in->numElements;e++)
	{
	    for (i=0;i<NN;i++)
	    {
		mask[in->Nodes[INDEX2(i,e,NN)]-offset]=1;
	    }
	}
    }
}

void Dudley_ElementFile_markDOFsConnectedToRange(index_t* mask,index_t offset,index_t marker,index_t firstDOF,index_t lastDOF,index_t *dofIndex,Dudley_ElementFile*in ,bool_t useLinear) 
{
    dim_t i,NN,e,j;
    index_t color;
//    Dudley_ReferenceElement* refElement=NULL;
    register index_t k;
   
    if (in!=NULL)
    {
//	refElement=Dudley_ReferenceElementSet_borrowReferenceElement(in->referenceElementSet, FALSE);	   
	NN=in->numNodes;
	for (color=in->minColor;color<=in->maxColor;color++)
	{
		#pragma omp parallel for private(e,i,j,k) schedule(static)
		for (e=0;e<in->numElements;e++)
		{
			if (in->Color[e]==color)
			{
				for (i=0;i<NN;i++)
				{
					k=dofIndex[in->Nodes[INDEX2(i,e,NN)]];
					if ( (firstDOF<=k) && (k<lastDOF) )
					{
						 for (j=0;j<NN;j++) mask[dofIndex[in->Nodes[INDEX2(j,e,NN)]]-offset]=marker;
						 break;
					}
				}
			}
		}
	}
    }	   
}

