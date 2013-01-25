
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

/*	  assemblage routines: interpolates nodal data in a data array onto elements (=integration points) */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/


void Finley_Assemble_interpolate(Finley_NodeFile *nodes, Finley_ElementFile* elements,escriptDataC* data,escriptDataC* interpolated_data) {
  __const double *data_array;
  double* local_data=NULL; 
  index_t dof_offset;
  bool_t reduced_integration=FALSE;
  dim_t q, i, NS_DOF, NN, numNodes=0, e, numQuad=0, numSub=0, isub, numShapesTotal2;
  Finley_ReferenceElement* reference_element=NULL;
  Finley_ShapeFunction *basis=NULL;
  dim_t numComps=getDataPointSize(data);
  index_t *resort_nodes=NULL, *map=NULL;
  type_t data_type=getFunctionSpaceType(data);
  size_t numComps_size;
  Finley_resetError();
  if (nodes==NULL || elements==NULL) return;
  reduced_integration = Finley_Assemble_reducedIntegrationOrder(interpolated_data);
  reference_element= Finley_ReferenceElementSet_borrowReferenceElement(elements->referenceElementSet, reduced_integration);
  NN=elements->numNodes;
  
  /* set some parameter */

  if (data_type==FINLEY_NODES) {
	   resort_nodes=reference_element->Type->subElementNodes;
	   numSub=reference_element->Type->numSubElements;
	   basis=reference_element->BasisFunctions;
	   numNodes=Finley_NodeFile_getNumNodes(nodes);
	   map=Finley_NodeFile_borrowTargetNodes(nodes);
  	   if (getFunctionSpaceType(interpolated_data)==FINLEY_CONTACT_ELEMENTS_2) {
	            dof_offset=reference_element->Type->offsets[1];
           } else {
	            dof_offset=reference_element->Type->offsets[0];
           }
  } else if (data_type==FINLEY_REDUCED_NODES) {
	   numSub=1;
	   resort_nodes=reference_element->Type->linearNodes;
	   basis=reference_element->LinearBasisFunctions;
	   numNodes=Finley_NodeFile_getNumReducedNodes(nodes);
	   map=Finley_NodeFile_borrowTargetReducedNodes(nodes);
  	   if (getFunctionSpaceType(interpolated_data)==FINLEY_CONTACT_ELEMENTS_2) {
	            dof_offset=reference_element->LinearType->offsets[1];
           } else {
	            dof_offset=reference_element->LinearType->offsets[0];
           }
  } else if (data_type==FINLEY_DEGREES_OF_FREEDOM) {
	   if (elements->MPIInfo->size>1) {
		  Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: for more than one processor DEGREES_OF_FREEDOM data are not accepted as input.");
		  return;
	   }
	   numSub=reference_element->Type->numSubElements;
	   resort_nodes=reference_element->Type->subElementNodes;
	   basis=reference_element->BasisFunctions;	
	   numNodes=Finley_NodeFile_getNumDegreesOfFreedom(nodes);
	   map=Finley_NodeFile_borrowTargetDegreesOfFreedom(nodes);
  	   if (getFunctionSpaceType(interpolated_data)==FINLEY_CONTACT_ELEMENTS_2) {
	            dof_offset=reference_element->Type->offsets[1];
           } else {
	            dof_offset=reference_element->Type->offsets[0];
           }
  } else if (data_type==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
	   if (elements->MPIInfo->size>1) {
		  Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: for more than one processor REDUCED_DEGREES_OF_FREEDOM data are not accepted as input.");
		  return;
	   }
	   numSub=1;
	   resort_nodes=reference_element->Type->linearNodes;
	   basis=reference_element->LinearBasisFunctions;
	   numNodes=Finley_NodeFile_getNumReducedDegreesOfFreedom(nodes);
	   map=Finley_NodeFile_borrowTargetReducedDegreesOfFreedom(nodes);
  	   if (getFunctionSpaceType(interpolated_data)==FINLEY_CONTACT_ELEMENTS_2) {
	            dof_offset=reference_element->LinearType->offsets[1];
           } else {
	            dof_offset=reference_element->LinearType->offsets[0];
           }
   } else {
	   Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: Cannot interpolate data");
	   return;
  }

  numQuad=basis->numQuadNodes;
  numShapesTotal2=basis->Type->numShapes * reference_element->Type->numSides;
  NS_DOF=basis->Type->numShapes;
  
  /* check the dimensions of interpolated_data and data */

  if (! numSamplesEqual(interpolated_data,numQuad*numSub,elements->numElements)) {
	   Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: illegal number of samples of output Data object");
  } else if (! numSamplesEqual(data,1,numNodes)) {
	   Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: illegal number of samples of input Data object");
  } else if (numComps!=getDataPointSize(interpolated_data)) {
	   Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: number of components of input and interpolated Data do not match.");
  }	 else if (!isExpanded(interpolated_data)) {
	   Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: expanded Data object is expected for output data.");
  }

  /* now we can start */

  if (Finley_noError()) {
	   requireWrite(interpolated_data);
	   #pragma omp parallel private(local_data, numComps_size, isub)
	   {
		  local_data=NULL; 
		  /* allocation of work arrays */
		  local_data=THREAD_MEMALLOC(NS_DOF*numComps*numSub,double); 
		  if (! Finley_checkPtr(local_data)) {
			  numComps_size=(size_t)numComps*sizeof(double);
			  /* open the element loop */
			  #pragma omp for private(e,q,i,data_array,isub) schedule(static)
			  for(e=0;e<elements->numElements;e++) {
				  for (isub=0; isub<numSub; isub++) {
					  for (q=0;q<NS_DOF;q++) {											  
						  i=elements->Nodes[INDEX2(resort_nodes[INDEX2(dof_offset+q,isub,numShapesTotal2)],e,NN)];
						  data_array=getSampleDataRO(data,map[i]);
						  memcpy(&(local_data[INDEX3(0,q,isub, numComps,NS_DOF)]), data_array, numComps_size);
					  }
				  }
				  /*  calculate interpolated_data=local_data*S */
				  Finley_Util_SmallMatSetMult1(numSub,numComps,numQuad,getSampleDataRW(interpolated_data,e),NS_DOF,local_data,basis->S);
			  } /* end of element loop */
		  }
		   THREAD_MEMFREE(local_data);
	 	} /* end of parallel region */
  }
}
