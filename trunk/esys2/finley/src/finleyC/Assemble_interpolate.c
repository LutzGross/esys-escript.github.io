/* $Id$ */

/**************************************************************/

/*    assemblage routines: interpolates nodal data in a data array onto elements (=integration points) */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004 */
/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"
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


void Finley_Assemble_interpolate(Finley_NodeFile *nodes, Finley_ElementFile* elements,escriptDataC* data,escriptDataC* interpolated_data) {
  double* local_data=NULL,*S=NULL,*data_array; 
  int dof_offset,q,*resort_nodes,type,i,NS_DOF,NN_DOF,numNodes,e;
  #define NODES 0
  #define DOF 1
  #define REDUCED_DOF 2
  if (nodes==NULL || elements==NULL) return;
  int NN=elements->ReferenceElement->Type->numNodes;
  int NS=elements->ReferenceElement->Type->numShapes;
  int numComps=getDataPointSize(data);
  int data_type=getFunctionSpaceType(data);
  int numQuad=elements->ReferenceElement->numQuadNodes;
  int id[NN];
  for (i=0;i<NN;i++) id[i]=i;

  /* set some parameter */

  if (data_type==FINLEY_NODES) {
       type=NODES;
       resort_nodes=id;
       NN_DOF=elements->ReferenceElement->Type->numNodes;
       NS_DOF=elements->ReferenceElement->Type->numShapes;
       S=elements->ReferenceElement->S;
       numNodes=nodes->numNodes;
  } else if (data_type==FINLEY_DEGREES_OF_FREEDOM) {
       type=DOF;
       resort_nodes=id;
       NN_DOF=elements->ReferenceElement->Type->numNodes;
       NS_DOF=elements->ReferenceElement->Type->numShapes;
       S=elements->ReferenceElement->S;
       numNodes=nodes->numDegreesOfFreedom;
  } else if (data_type==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
       type=REDUCED_DOF;
       resort_nodes=elements->ReferenceElement->Type->linearNodes;
       NN_DOF=elements->LinearReferenceElement->Type->numNodes;
       NS_DOF=elements->LinearReferenceElement->Type->numShapes;
       S=elements->LinearReferenceElement->S;
       numNodes=nodes->reducedNumDegreesOfFreedom;
   } else {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"Cannot interpolate data");
  }

  if (getFunctionSpaceType(interpolated_data)==FINLEY_CONTACT_ELEMENTS_2) {
       dof_offset=NN_DOF-NS_DOF;
  } else {
       dof_offset=0;
  }

  /* check the dimensions of interpolated_data and data */

  if (! numSamplesEqual(interpolated_data,numQuad,elements->numElements)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"illegal number of samples of output Data object");
  } else if (! numSamplesEqual(data,1,numNodes)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"illegal number of samples of input Data object");
  } else if (numComps!=getDataPointSize(interpolated_data)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"number of components of input and interpolated Data do not match.");
  }  else if (!isExpanded(interpolated_data)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"expanded Data object is expected for output data.");
  }

  /* now we can start */

  if (Finley_ErrorCode==NO_ERROR) {
       #pragma omp parallel private(local_data)
       {
          local_data=NULL; 
          /* allocation of work arrays */
          local_data=(double*) THREAD_MEMALLOC(NS*numComps*sizeof(double)); 
          if (! Finley_checkPtr(local_data)) {

	    /* open the element loop */

            #pragma omp for private(e,q,i,data_array) schedule(static)
	    for(e=0;e<elements->numElements;e++) {
   	      /* gather local data into local_data(numComps,NS_DOF): */
              switch (type) {
                 case NODES:
                        for (q=0;q<NS_DOF;q++) {
                           i=elements->Nodes[INDEX2(resort_nodes[dof_offset+q],e,NN)];
                           data_array=getSampleData(data,i);
                           Finley_copyDouble(numComps,data_array,local_data+q*numComps);
                        break;
                 case DOF:
                        for (q=0;q<NS_DOF;q++) {
                           i=elements->Nodes[INDEX2(resort_nodes[dof_offset+q],e,NN)];
                           data_array=getSampleData(data,nodes->degreeOfFreedom[i]);
                           Finley_copyDouble(numComps,data_array,local_data+q*numComps);
                        }
                        break;
                 case REDUCED_DOF:
                        for (q=0;q<NS_DOF;q++) {
                           i=elements->Nodes[INDEX2(resort_nodes[dof_offset+q],e,NN)];
                           data_array=getSampleData(data,nodes->reducedDegreeOfFreedom[i]);
                           Finley_copyDouble(numComps,data_array,local_data+q*numComps);
                        }
                        break;
              }

	      /*  calculate interpolated_data=local_data*S */

	      Finley_Util_SmallMatMult(numComps,numQuad,getSampleData(interpolated_data,e),NS_DOF,local_data,S);

	    } /* end of element loop */

          }
	  THREAD_MEMFREE(local_data);
        }
      }
  }
  #undef NODES 
  #undef DOF 
  #undef REDUCED_DOF 
}
/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.2  2004/07/21 05:00:54  gross
 * name changes in DataC
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
