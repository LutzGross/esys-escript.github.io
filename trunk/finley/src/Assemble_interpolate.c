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

/*    assemblage routines: interpolates nodal data in a data array onto elements (=integration points) */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/


void Finley_Assemble_interpolate(Finley_NodeFile *nodes, Finley_ElementFile* elements,escriptDataC* data,escriptDataC* interpolated_data) {
  double* local_data=NULL,*S=NULL,*data_array; 
  index_t dof_offset,*resort_nodesi, NN, NS;
  bool_t reduced_integration=FALSE;
  dim_t q,i,NS_DOF,NN_DOF,numNodes,e, numQuad;
  Finley_RefElement* reference_element=NULL;
  dim_t numComps=getDataPointSize(data);
  index_t id[MAX_numNodes], *resort_nodes;
  type_t data_type=getFunctionSpaceType(data);
  type_t type;
  Finley_resetError();
  #define NODES 0
  #define REDUCED_NODES 3
  #define DOF 1
  #define REDUCED_DOF 2
  if (nodes==NULL || elements==NULL) return;
  NN=elements->ReferenceElement->Type->numNodes;
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
       numNodes=nodes->numNodes;
  } else if (data_type==FINLEY_REDUCED_NODES) {
       type=REDUCED_NODES;
       if (reduced_integration) {
           reference_element=elements->LinearReferenceElementReducedOrder;
       } else {
           reference_element=elements->LinearReferenceElement;
       } 
       /* TODO */
       Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: input from reduced nodes is not supported yet.");
  } else if (data_type==FINLEY_DEGREES_OF_FREEDOM) {
       type=DOF;
       resort_nodes=id;
       if (reduced_integration) {
           reference_element=elements->ReferenceElementReducedOrder;
       } else {
           reference_element=elements->ReferenceElement;
       }
       numNodes=nodes->numDegreesOfFreedom;
  } else if (data_type==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
       type=REDUCED_DOF;
       resort_nodes=elements->ReferenceElement->Type->linearNodes;
       if (reduced_integration) {
           reference_element=elements->LinearReferenceElementReducedOrder;
       } else {
           reference_element=elements->LinearReferenceElement;
       }
       numNodes=nodes->reducedNumDegreesOfFreedom;
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
       #pragma omp parallel private(local_data)
       {
          local_data=NULL; 
          /* allocation of work arrays */
          local_data=THREAD_MEMALLOC(NS*numComps,double); 
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
                        }
                        break;
                 case REDUCED_NODES:
                        for (q=0;q<NS_DOF;q++) {
                           i=elements->Nodes[INDEX2(resort_nodes[dof_offset+q],e,NN)];
                           data_array=getSampleData(data,i); /* TODO */
                           Finley_copyDouble(numComps,data_array,local_data+q*numComps);
                        }
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
     } /* end of parallel region */
  }
  #undef NODES 
  #undef REDUCED_NODES 
  #undef DOF 
  #undef REDUCED_DOF 
}
/*
 * $Log$
 * Revision 1.6  2005/09/15 03:44:21  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.5.2.1  2005/09/07 06:26:18  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.5  2005/07/08 04:07:48  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.4  2004/12/15 07:08:32  jgs
 * *** empty log message ***
 * Revision 1.1.1.1.2.3  2005/06/29 02:34:48  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1.2.2  2004/11/24 01:37:12  gross
 * some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
 *
 *
 *
 */
