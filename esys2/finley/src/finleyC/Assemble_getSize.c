/* $Id$ */

/**************************************************************/

/*    assemblage routines: */

/*    calculates the minimum distance between two vertices of elements and assigns the value to each  */
/*    quadrature point in element_size                                                                         */


/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004 */
/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "escript/Data/DataC.h"
#include "Finley.h"
#include "Assemble.h"
#include "NodeFile.h"
#include "ElementFile.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/
void Finley_Assemble_getSize(Finley_NodeFile* nodes, Finley_ElementFile* elements, escriptDataC* element_size) {

  double *local_X=NULL,*element_size_array;
  int e,n0,n1,q,i,node_offset;
  double d,diff,min_diff;

  if (nodes==NULL || elements==NULL) return;
  int NVertices=elements->ReferenceElement->Type->numVertices;
  int NN=elements->ReferenceElement->Type->numNodes;
  int NS=elements->ReferenceElement->Type->numShapes;
  int numQuad=elements->ReferenceElement->numQuadNodes;
  int numDim=nodes->numDim;

  /* set a few more parameters */

  if (getFunctionSpaceType(element_size)==FINLEY_CONTACT_ELEMENTS_2) {
      node_offset=NN-NS;
  } else {
      node_offset=0;
  }

  /* check the dimensions of element_size */

  if (numDim!=elements->ReferenceElement->Type->numDim) {
     Finley_ErrorCode=TYPE_ERROR;
     sprintf(Finley_ErrorMsg,"Gradient: Spatial and element dimension must match.");
  } else if (! numSamplesEqual(element_size,numQuad,elements->numElements)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"illegal number of samples of element size Data object");
  } else if (! isDataPointShapeEqual(element_size,0,&(numDim))) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"illegal data point shape of element size Data object");
  }  else if (!isExpanded(element_size)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"expanded Data object is expected for element size.");
  }
  /* now we can start: */

  if (Finley_ErrorCode==NO_ERROR) {
        #pragma omp parallel private(local_X)
        {
	   /* allocation of work arrays */
	   local_X=(double*) THREAD_MEMALLOC(NN*numDim*sizeof(double));
	   if (! Finley_checkPtr(local_X) ) {
	     /* open the element loop */
             #pragma omp for private(e,min_diff,diff,n0,n1,d,q,i,element_size_array)
	     for(e=0;e<elements->numElements;e++) {
	       /* gather local coordinates of nodes into local_X(numDim,NN): */
	       Finley_Util_Gather_double(NS,&(elements->Nodes[INDEX2(node_offset,e,NN)]),numDim,nodes->Coordinates,local_X);
	       /* calculate minimal differences */
	       min_diff=-1;
	       for (n0=0;n0<NVertices;n0++) {
	         for (n1=n0+1;n1<NVertices;n1++) {
		   diff=0;
		   for (i=0;i<numDim;i++) {
		     d=local_X[INDEX2(i,n0,numDim)]-local_X[INDEX2(i,n1,numDim)];
		     diff+=d*d;
		   }
		   if (min_diff<0) {
		     min_diff=diff;
		   } else {
		     min_diff=MIN(min_diff,diff);
		   }
	         }
	       }
	       min_diff=sqrt(MAX(min_diff,0));
	       /* set all values to min_diff */
               element_size_array=getSampleData(element_size,e);
	       for (q=0;q<numQuad;q++) element_size_array[q]=min_diff;
	     }
	   }
	   THREAD_MEMFREE(local_X);
        }
  }
  return;
}
/*
 * $Log$
 * Revision 1.3  2004/12/15 03:48:44  jgs
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.3  2004/07/30 04:37:06  gross
 * escript and finley are linking now and RecMeshTest.py has been passed
 *
 * Revision 1.2  2004/07/21 05:00:54  gross
 * name changes in DataC
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
