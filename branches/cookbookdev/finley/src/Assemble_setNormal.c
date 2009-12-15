
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
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


void Finley_Assemble_setNormal(Finley_NodeFile* nodes, Finley_ElementFile* elements, escriptDataC* normal) {
  double *local_X=NULL, *dVdv=NULL,*normal_array;
  index_t sign,node_offset;
  Finley_ReferenceElement* reference_element=NULL;
  dim_t e,q, NN, NS, numDim, numQuad, numDim_local;
  if (nodes==NULL || elements==NULL) return;
  Finley_resetError();
  reference_element=Finley_ReferenceElementSet_borrowReferenceElement(elements->referenceElementSet, Finley_Assemble_reducedIntegrationOrder(normal));
  NN=elements->numNodes;
  numDim=nodes->numDim;

  numQuad=reference_element->Parametrization->numQuadNodes;
  numDim_local=reference_element->Parametrization->Type->numDim;
  NS=reference_element->Parametrization->Type->numShapes;


  
  /* set some parameter */

  if (getFunctionSpaceType(normal)==FINLEY_CONTACT_ELEMENTS_2) {
	  node_offset=reference_element->Type->offsets[1];
	  sign=-1;
  } else {
	  node_offset=reference_element->Type->offsets[0];
	  sign=1;
  }
  /* check the dimensions of normal */
  if (! (numDim==numDim_local || numDim-1==numDim_local)) {
	   Finley_setError(TYPE_ERROR,"Finley_Assemble_setNormal: Cannot calculate normal vector");
  } else if (! isDataPointShapeEqual(normal,1,&(numDim))) {
	   Finley_setError(TYPE_ERROR,"Finley_Assemble_setNormal: illegal number of samples of normal Data object");
  } else if (! numSamplesEqual(normal,numQuad,elements->numElements)) {
	   Finley_setError(TYPE_ERROR,"Finley_Assemble_setNormal: illegal number of samples of normal Data object");
  } else if (! isDataPointShapeEqual(normal,1,&(numDim))) {
	   Finley_setError(TYPE_ERROR,"Finley_Assemble_setNormal: illegal point data shape of normal Data object");
  }	 else if (!isExpanded(normal)) {
	   Finley_setError(TYPE_ERROR,"Finley_Assemble_setNormal: expanded Data object is expected for normal.");
  }
   
  /* now we can start */
	if (Finley_noError()) {
	  requireWrite(normal);
		  #pragma omp parallel private(local_X,dVdv)
		  {
			 local_X=dVdv=NULL;
			 /* allocation of work arrays */
			 local_X=THREAD_MEMALLOC(NS*numDim,double); 
			 dVdv=THREAD_MEMALLOC(numQuad*numDim*numDim_local,double); 
			 if (!(Finley_checkPtr(local_X) || Finley_checkPtr(dVdv) ) ) {
					   /* open the element loop */
					   #pragma omp for private(e,q,normal_array) schedule(static)
					   for(e=0;e<elements->numElements;e++) {
						  /* gather local coordinates of nodes into local_X: */
						  Finley_Util_Gather_double(NS,&(elements->Nodes[INDEX2(node_offset,e,NN)]),numDim,nodes->Coordinates,local_X);
						  /*  calculate dVdv(i,j,q)=local_X(i,n)*DSDv(n,j,q) */
						  Finley_Util_SmallMatMult(numDim,numDim_local*numQuad,dVdv,NS,local_X,reference_element->Parametrization->dSdv);
						  /* get normalized vector:	 */
						  normal_array=getSampleDataRW(normal,e);
						  Finley_NormalVector(numQuad,numDim,numDim_local,dVdv,normal_array);
						  for (q=0;q<numQuad*numDim;q++) normal_array[q]*=sign;
					   }
				 }
				 THREAD_MEMFREE(local_X);
				 THREAD_MEMFREE(dVdv);
			 }
	}
}

