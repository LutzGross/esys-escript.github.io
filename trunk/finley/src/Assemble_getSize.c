
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*    assemblage routines: */

/*    calculates the minimum distance between two vertices of elements and assigns the value to each  */
/*    quadrature point in element_size */


/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/
void Finley_Assemble_getSize(Finley_NodeFile* nodes, Finley_ElementFile* elements, escriptDataC* element_size) {

  Finley_ShapeFunction *shape=NULL;
  Finley_ReferenceElement*  refElement=NULL;
  double *local_X=NULL,*element_size_array;
  dim_t e,n0,n1,q,i, NVertices, NN, NS, numQuad, numDim;
  index_t node_offset;
  double d,diff,max_diff, f;
  Finley_resetError();

  if (nodes==NULL || elements==NULL) return;
  refElement=Finley_ReferenceElementSet_borrowReferenceElement(elements->referenceElementSet,Finley_Assemble_reducedIntegrationOrder(element_size));
  shape= refElement->Parametrization;
  numDim=nodes->numDim;
  numQuad=shape->numQuadNodes;
  NN=elements->numNodes;
  NS=refElement->Parametrization->Type->numShapes;
  NVertices=refElement->Parametrization->Type->numVertices;
  
  
  /* set a few more parameters */

  if (getFunctionSpaceType(element_size)==FINLEY_CONTACT_ELEMENTS_2) {
      node_offset=refElement->Type->offsets[1];
  } else {
      node_offset=refElement->Type->offsets[0];
  }
  f=pow(0.5, pow((double)(refElement->Type->numSubElements), 1./(double)(numDim))-1);

  /* check the dimensions of element_size */

  if (! numSamplesEqual(element_size,numQuad,elements->numElements)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_getSize: illegal number of samples of element size Data object");
  } else if (! isDataPointShapeEqual(element_size,0,&(numDim))) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_getSize: illegal data point shape of element size Data object");
  }  else if (!isExpanded(element_size)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_getSize: expanded Data object is expected for element size.");
  }
  /* now we can start: */

  if (Finley_noError()) {
    requireWrite(element_size);
        #pragma omp parallel private(local_X)
        {
			/* allocation of work arrays */
			local_X=THREAD_MEMALLOC(NN*numDim,double);
			if (! Finley_checkPtr(local_X) ) {
				/* open the element loop */
				#pragma omp for private(e,max_diff,diff,n0,n1,d,q,i,element_size_array) schedule(static)
				for(e=0;e<elements->numElements;e++) {
					/* gather local coordinates of nodes into local_X(numDim,NN): */
					Finley_Util_Gather_double(NS,&(elements->Nodes[INDEX2(node_offset,e,NN)]),numDim,nodes->Coordinates,local_X);
					/* calculate minimal differences */
					max_diff=0.;
					for (n0=0;n0<NVertices;n0++) {
						for (n1=n0+1;n1<NVertices;n1++) {
							diff=0;
							for (i=0;i<numDim;i++) {
								d=local_X[INDEX2(i,n0,numDim)]-local_X[INDEX2(i,n1,numDim)];
								diff+=d*d;
							}
							max_diff=MAX(max_diff,diff);
						}
					}
					max_diff=sqrt(max_diff)*f;
					/* set all values to max_diff */
					element_size_array=getSampleDataRW(element_size,e);
					for (q=0;q<numQuad;q++) element_size_array[q]=max_diff;
				}
			}
			THREAD_MEMFREE(local_X);
       } /* end of parallel region */
  }
  return;
}
