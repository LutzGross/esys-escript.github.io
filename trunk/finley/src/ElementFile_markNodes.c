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

/*   Finley: ElementFile */

/*   mark the used nodes with offeset: */

/**************************************************************/

/*  Copyrights by ACcESS Australia 2003,2004,2005 */
/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "ElementFile.h"

/**************************************************************/

void Finley_ElementFile_markNodes(index_t* mask,index_t offset,Finley_ElementFile* in,bool_t useLinear) {
   dim_t i,NN,NN2,e;
   index_t color,*lin_node;
   if (in!=NULL) {
     index_t id[in->ReferenceElement->Type->numNodes];
     for (i=0;i<in->ReferenceElement->Type->numNodes;i++) id[i]=i;
     if (useLinear) {
        NN=in->LinearReferenceElement->Type->numNodes;
        lin_node=in->ReferenceElement->Type->linearNodes;
     } else {
        NN=in->ReferenceElement->Type->numNodes;
        lin_node=id;
     }
     NN2=in->ReferenceElement->Type->numNodes;
/* knock out the OpenMP directives for the MPI version, get the OpenMP running
   after the MPI is working satisfactorily */
#ifndef PASO_MPI
     if ((in->maxColor-in->minColor+1)*NN<in->numElements) {
        #pragma omp parallel private(color)
        {
           for (color=in->minColor;color<=in->maxColor;color++) {
             #pragma omp for private(e,i) schedule(static)
             for (e=0;e<in->numElements;e++) {
               if (in->Color[e]==color) {
                  for (i=0;i<NN;i++) 
                    mask[in->Nodes[INDEX2(lin_node[i],e,NN2)]-offset]=1;
               }
             }
           }
           #pragma omp barrier
        }
      } else {
        #pragma omp parallel for private(e,i) schedule(static)
#endif
        for (e=0;e<in->numElements;e++) {
           for (i=0;i<NN;i++) 
             mask[in->Nodes[INDEX2(lin_node[i],e,NN2)]-offset]=1;
        }
#ifndef PASO_MPI
      }
#endif
   }
}

#ifdef PASO_MPI

/* these have been trimmed down from the colors-based scheme for simplicity of testing
   in the MPI scheme. */
void Finley_ElementFile_markBoundaryElementNodes(index_t* mask,index_t offset,Finley_ElementFile* in,bool_t useLinear) {
   dim_t i,NN,NN2,e;
   index_t color,*lin_node;
   if (in!=NULL&& in->numElements>0) {
     index_t id[in->ReferenceElement->Type->numNodes];
     for (i=0;i<in->ReferenceElement->Type->numNodes;i++) id[i]=i;
     if (useLinear) {
        NN=in->LinearReferenceElement->Type->numNodes;
        lin_node=in->ReferenceElement->Type->linearNodes;
     } else {
        NN=in->ReferenceElement->Type->numNodes;
        lin_node=id;
     }
     NN2=in->ReferenceElement->Type->numNodes;
     for (e=0;e<in->elementDistribution->numLocal;e++) 
			 if( in->Id[e]>=in->elementDistribution->numInternal )
				 for (i=0;i<NN;i++) 
					 mask[in->Nodes[INDEX2(lin_node[i],e,NN2)]-offset]=1;
   }
}

void Finley_ElementFile_markInternalElementNodes(index_t* mask,index_t offset,Finley_ElementFile* in,bool_t useLinear) {
   dim_t i,NN,NN2,e;
   index_t color,*lin_node;
   if (in!=NULL && in->numElements>0) {
     index_t id[in->ReferenceElement->Type->numNodes];
     for (i=0;i<in->ReferenceElement->Type->numNodes;i++) id[i]=i;
     if (useLinear) {
        NN=in->LinearReferenceElement->Type->numNodes;
        lin_node=in->ReferenceElement->Type->linearNodes;
     } else {
        NN=in->ReferenceElement->Type->numNodes;
        lin_node=id;
     }
     NN2=in->ReferenceElement->Type->numNodes;

     for (e=0;e<in->elementDistribution->numLocal;e++)
			 if( in->Id[e]<in->elementDistribution->numInternal )
				 for (i=0;i<NN;i++) 
					 mask[in->Nodes[INDEX2(lin_node[i],e,NN2)]-offset]=1;
   }
}

void Finley_ElementFile_markBoundaryElementDOF(index_t* mask,index_t offset,index_t *degreeOfFreedom, Finley_ElementFile* in,bool_t useLinear, index_t startElement) {
	dim_t i,NN,NN2,e;
	index_t color,*lin_node;
	if (in!=NULL&& in->numElements>0) {

		index_t id[in->ReferenceElement->Type->numNodes];
		for (i=0;i<in->ReferenceElement->Type->numNodes;i++) id[i]=i;
			if (useLinear) {
				NN=in->LinearReferenceElement->Type->numNodes;
				lin_node=in->ReferenceElement->Type->linearNodes;
			} else {
				NN=in->ReferenceElement->Type->numNodes;
				lin_node=id;
			}
			NN2=in->ReferenceElement->Type->numNodes;
			for (e=0;e<in->elementDistribution->numLocal;e++) 
	 			if( in->Id[e]-startElement>=in->elementDistribution->numInternal )
					for (i=0;i<NN;i++) 
						mask[degreeOfFreedom[in->Nodes[INDEX2(lin_node[i],e,NN2)]]-offset]=1;
	}

}

void Finley_ElementFile_markInternalElementDOF(index_t* mask,index_t offset,index_t *degreeOfFreedom, Finley_ElementFile* in,bool_t useLinear, index_t startElement) {
   dim_t i,NN,NN2,e;
   index_t color,*lin_node;
   if (in!=NULL && in->numElements>0) {
     index_t id[in->ReferenceElement->Type->numNodes];
     for (i=0;i<in->ReferenceElement->Type->numNodes;i++) id[i]=i;
     if (useLinear) {
        NN=in->LinearReferenceElement->Type->numNodes;
        lin_node=in->ReferenceElement->Type->linearNodes;
     } else {
        NN=in->ReferenceElement->Type->numNodes;
        lin_node=id;
     }
     NN2=in->ReferenceElement->Type->numNodes;

     for (e=0;e<in->elementDistribution->numLocal;e++) 
			 if( in->Id[e]-startElement<in->elementDistribution->numInternal )
				 for (i=0;i<NN;i++) 
					 mask[degreeOfFreedom[in->Nodes[INDEX2(lin_node[i],e,NN2)]]-offset]=1;
   }
}
#endif

/* 
* $Log$
* Revision 1.4  2005/09/15 03:44:22  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.3.2.1  2005/09/07 06:26:18  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.3  2005/07/22 03:53:08  jgs
* Merge of development branch back to main trunk on 2005-07-22
*
* Revision 1.2  2005/07/08 04:07:50  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.1.1.1.2.2  2005/07/18 10:34:54  gross
* some informance improvements when reading meshes
*
* Revision 1.1.1.1.2.1  2005/06/29 02:34:49  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/
