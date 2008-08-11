
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/*    assemblage routines: interpolates nodal data in a data array onto elements (=integration points) */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/


void Finley_Assemble_interpolate(Finley_NodeFile *nodes, Finley_ElementFile* elements,escriptDataC* data,escriptDataC* interpolated_data) {
  double* local_data=NULL,*S=NULL,*data_array; 
  index_t dof_offset, NN, NS;
  bool_t reduced_integration=FALSE;
  dim_t q,i,NS_DOF,NN_DOF,numNodes,e, numQuad;
  Finley_RefElement* reference_element=NULL;
  dim_t numComps=getDataPointSize(data);
  index_t id[MAX_numNodes], *resort_nodes, *map;
  type_t data_type=getFunctionSpaceType(data);
  type_t type;
  size_t numComps_size;
  Finley_resetError();
  #define NODES 0
  #define REDUCED_NODES 3
  #define DOF 1
  #define REDUCED_DOF 2
  if (nodes==NULL || elements==NULL) return;
  
  NN=elements->numNodes;
  NS=elements->ReferenceElement->Type->numShapes;
  reduced_integration = Finley_Assemble_reducedIntegrationOrder(interpolated_data);
  for (i=0;i<NN;i++) id[i]=i;

  /* set some parameter */

  if (data_type==FINLEY_NODES) {
       type=NODES;
       resort_nodes=id;
       if (reduced_integration) {
          reference_element=elements->ReferenceElementReducedOrder;
       } else {
          reference_element=elements->ReferenceElement;
       } 
       numNodes=Finley_NodeFile_getNumNodes(nodes);
       map=Finley_NodeFile_borrowTargetNodes(nodes);
  } else if (data_type==FINLEY_REDUCED_NODES) {
       type=REDUCED_NODES;
       resort_nodes=elements->ReferenceElement->Type->linearNodes;
       if (reduced_integration) {
           reference_element=elements->LinearReferenceElementReducedOrder;
       } else {
           reference_element=elements->LinearReferenceElement;
       } 
       numNodes=Finley_NodeFile_getNumReducedNodes(nodes);
       map=Finley_NodeFile_borrowTargetReducedNodes(nodes);
  } else if (data_type==FINLEY_DEGREES_OF_FREEDOM) {
       if (elements->MPIInfo->size>1) {
          Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: for more than one processor DEGREES_OF_FREEDOM data are not accepted as input.");
          return;
       }
       type=DOF;
       resort_nodes=id;
       if (reduced_integration) {
           reference_element=elements->ReferenceElementReducedOrder;
       } else {
           reference_element=elements->ReferenceElement;
       }
       numNodes=Finley_NodeFile_getNumDegreesOfFreedom(nodes);
       map=Finley_NodeFile_borrowTargetDegreesOfFreedom(nodes);
  } else if (data_type==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
       if (elements->MPIInfo->size>1) {
          Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: for more than one processor REDUCED_DEGREES_OF_FREEDOM data are not accepted as input.");
          return;
       }
       type=REDUCED_DOF;
       resort_nodes=elements->ReferenceElement->Type->linearNodes;
       if (reduced_integration) {
           reference_element=elements->LinearReferenceElementReducedOrder;
       } else {
           reference_element=elements->LinearReferenceElement;
       }
       numNodes=Finley_NodeFile_getNumReducedDegreesOfFreedom(nodes);
       map=Finley_NodeFile_borrowTargetReducedDegreesOfFreedom(nodes);
   } else {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: Cannot interpolate data");
  }
  NN_DOF=reference_element->Type->numNodes;
  NS_DOF=reference_element->Type->numShapes;
  S=reference_element->S;
  numQuad=reference_element->numQuadNodes;
  if (getFunctionSpaceType(interpolated_data)==FINLEY_CONTACT_ELEMENTS_2) {
       dof_offset=NN_DOF-NS_DOF;
  } else {
       dof_offset=0;
  }

  /* check the dimensions of interpolated_data and data */

  if (! numSamplesEqual(interpolated_data,numQuad,elements->numElements)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: illegal number of samples of output Data object");
  } else if (! numSamplesEqual(data,1,numNodes)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: illegal number of samples of input Data object");
  } else if (numComps!=getDataPointSize(interpolated_data)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: number of components of input and interpolated Data do not match.");
  }  else if (!isExpanded(interpolated_data)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: expanded Data object is expected for output data.");
  }

  /* now we can start */

  if (Finley_noError()) {
       #pragma omp parallel private(local_data, numComps_size)
       {
          local_data=NULL; 
          /* allocation of work arrays */
          local_data=THREAD_MEMALLOC(NS*numComps,double); 
          if (! Finley_checkPtr(local_data)) {

            numComps_size=(size_t)numComps*sizeof(double);
	    /* open the element loop */

            #pragma omp for private(e,q,i,data_array) schedule(static)
	    for(e=0;e<elements->numElements;e++) {
              for (q=0;q<NS_DOF;q++) {
                      i=elements->Nodes[INDEX2(resort_nodes[dof_offset+q],e,NN)];
                      data_array=getSampleData(data,map[i]);
                      memcpy(local_data+q*numComps, data_array, numComps_size);
              }
	      /*  calculate interpolated_data=local_data*S */

	      Finley_Util_SmallMatMult(numComps,numQuad,getSampleData(interpolated_data,e),NS_DOF,local_data,S);

	    } /* end of element loop */

          }
	  THREAD_MEMFREE(local_data);
     } /* end of parallel region */
  }
  #undef NODES 
  #undef REDUCED_NODES 
  #undef DOF 
  #undef REDUCED_DOF 
}
