/* $Id$ */

/**************************************************************/

/*    assemblage routines: calculates the normal vector at quadrature points on face elements */

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


void Finley_Assemble_setNormal(Finley_NodeFile* nodes, Finley_ElementFile* elements, escriptDataC* normal) {
  double *local_X=NULL, *dVdv=NULL,*normal_array;
  int e,sign,q,node_offset;
  if (nodes==NULL || elements==NULL) return;
  int NN=elements->ReferenceElement->Type->numNodes;
  int NS=elements->ReferenceElement->Type->numShapes;
  int numDim=nodes->numDim;
  int numQuad=elements->ReferenceElement->numQuadNodes;
  int numDim_local=elements->ReferenceElement->Type->numDim;

  /* set some parameter */

  if (getFunctionSpaceType(normal)==FINLEY_CONTACT_ELEMENTS_2) {
      node_offset=NN-NS;
      sign=-1;
  } else {
      node_offset=0;
      sign=1;
  }
  /* check the dimensions of normal */
  if (numDim==numDim_local || numDim==numDim_local-1) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"Cannot calculate normal vector");
  } else if (! isDataPointShapeEqual(normal,1,&(numDim))) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"illegal number of samples of normal Data object");
  } else if (! numSamplesEqual(normal,numQuad,elements->numElements)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"illegal number of samples of normal Data object");
  } else if (! isDataPointShapeEqual(normal,1,&(numDim))) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"illegal point data shape of normal Data object");
  }  else if (!isExpanded(normal)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"expanded Data object is expected for normal.");
  }
   
  /* now we can start */
    if (Finley_ErrorCode==NO_ERROR) {
          #pragma omp parallel private(local_X,dVdv)
          {
             local_X=dVdv=NULL;
             /* allocation of work arrays */
             local_X=(double*) THREAD_MEMALLOC(NS*numDim*sizeof(double)); 
             dVdv=(double*) THREAD_MEMALLOC(numQuad*numDim*numDim_local*sizeof(double)); 
             if (!(Finley_checkPtr(local_X) || Finley_checkPtr(dVdv) ) ) {
                       /* open the element loop */
                       #pragma omp for private(e,q,normal_array) schedule(static)
                       for(e=0;e<elements->numElements;e++) {
                          /* gather local coordinates of nodes into local_X: */
                          Finley_Util_Gather_double(NS,&(elements->Nodes[INDEX2(node_offset,e,NN)]),numDim,nodes->Coordinates,local_X);
                          /*  calculate dVdv(i,j,q)=local_X(i,n)*DSDv(n,j,q) */
                          Finley_Util_SmallMatMult(numDim,numDim_local*numQuad,dVdv,NS,local_X,elements->ReferenceElement->dSdv);
                          /* get normalized vector:  */
                          normal_array=getSampleData(normal,e);
                          Finley_NormalVector(numQuad,numDim,numDim_local,dVdv,normal_array);
			  for (q=0;q<numQuad*numDim;q++) normal_array[q]*=sign;
                       }
                 }
                 THREAD_MEMFREE(local_X);
                 THREAD_MEMFREE(dVdv);
             }
    }
}
/*
 * $Log$
 * Revision 1.3  2004/12/15 03:48:45  jgs
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
