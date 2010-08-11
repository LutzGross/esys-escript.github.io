
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

/*	  assemblage routines: calculates the normal vector at quadrature points on face elements */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/


void Dudley_Assemble_setNormal(Dudley_NodeFile* nodes, Dudley_ElementFile* elements, escriptDataC* normal) {
  double *local_X=NULL, *dVdv=NULL,*normal_array;
  index_t sign,node_offset;
  Dudley_ReferenceElement* reference_element=NULL;
  dim_t e,q, NN, NS, numDim, numQuad, numDim_local;
  if (nodes==NULL || elements==NULL) return;
  Dudley_resetError();
  reference_element=Dudley_ReferenceElementSet_borrowReferenceElement(elements->referenceElementSet, Dudley_Assemble_reducedIntegrationOrder(normal));
  NN=elements->numNodes;
  numDim=nodes->numDim;

  numQuad=reference_element->Parametrization->numQuadNodes;
  numDim_local=reference_element->Parametrization->Type->numDim;
  NS=reference_element->Parametrization->Type->numShapes;


  
  /* set some parameter */

  node_offset=reference_element->Type->offsets[0];
  sign=1;
  /* check the dimensions of normal */
  if (! (numDim==numDim_local || numDim-1==numDim_local)) {
	   Dudley_setError(TYPE_ERROR,"Dudley_Assemble_setNormal: Cannot calculate normal vector");
  } else if (! isDataPointShapeEqual(normal,1,&(numDim))) {
	   Dudley_setError(TYPE_ERROR,"Dudley_Assemble_setNormal: illegal number of samples of normal Data object");
  } else if (! numSamplesEqual(normal,numQuad,elements->numElements)) {
	   Dudley_setError(TYPE_ERROR,"Dudley_Assemble_setNormal: illegal number of samples of normal Data object");
  } else if (! isDataPointShapeEqual(normal,1,&(numDim))) {
	   Dudley_setError(TYPE_ERROR,"Dudley_Assemble_setNormal: illegal point data shape of normal Data object");
  }	 else if (!isExpanded(normal)) {
	   Dudley_setError(TYPE_ERROR,"Dudley_Assemble_setNormal: expanded Data object is expected for normal.");
  }
   
  /* now we can start */
	if (Dudley_noError()) {
	  requireWrite(normal);
		  #pragma omp parallel private(local_X,dVdv)
		  {
			 local_X=dVdv=NULL;
			 /* allocation of work arrays */
			 local_X=THREAD_MEMALLOC(NS*numDim,double); 
			 dVdv=THREAD_MEMALLOC(numQuad*numDim*numDim_local,double); 
			 if (!(Dudley_checkPtr(local_X) || Dudley_checkPtr(dVdv) ) ) {
					   /* open the element loop */
					   #pragma omp for private(e,q,normal_array) schedule(static)
					   for(e=0;e<elements->numElements;e++) {
						  /* gather local coordinates of nodes into local_X: */
						  Dudley_Util_Gather_double(NS,&(elements->Nodes[INDEX2(node_offset,e,NN)]),numDim,nodes->Coordinates,local_X);
						  /*  calculate dVdv(i,j,q)=local_X(i,n)*DSDv(n,j,q) */
						  Dudley_Util_SmallMatMult(numDim,numDim_local*numQuad,dVdv,NS,local_X,reference_element->Parametrization->dSdv);
						  /* get normalized vector:	 */
						  normal_array=getSampleDataRW(normal,e);
						  Dudley_NormalVector(numQuad,numDim,numDim_local,dVdv,normal_array);
						  for (q=0;q<numQuad*numDim;q++) normal_array[q]*=sign;
					   }
				 }
				 THREAD_MEMFREE(local_X);
				 THREAD_MEMFREE(dVdv);
			 }
	}
}

