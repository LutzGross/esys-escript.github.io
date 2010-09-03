
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


void Dudley_Assemble_interpolate(Dudley_NodeFile *nodes, Dudley_ElementFile* elements,escriptDataC* data,escriptDataC* interpolated_data) {
  __const double *data_array;
  double* local_data=NULL; 
  bool_t reduced_integration=FALSE;
  dim_t q, i, NS_DOF, NN, numNodes=0, e, numQuad=0, numShapesTotal2;
  Dudley_ReferenceElement* reference_element=NULL;
  Dudley_ShapeFunction *basis=NULL;
  dim_t numComps=getDataPointSize(data);
  index_t *map=NULL;
  type_t data_type=getFunctionSpaceType(data);
  type_t type;
  size_t numComps_size;
  Dudley_resetError();
  #define NODES 0
  #define REDUCED_NODES 3
  #define DOF 1
  #define REDUCED_DOF 2
  if (nodes==NULL || elements==NULL) return;
  reduced_integration = Dudley_Assemble_reducedIntegrationOrder(interpolated_data);
  reference_element= Dudley_ReferenceElementSet_borrowReferenceElement(elements->referenceElementSet, reduced_integration);
  NN=elements->numNodes;
  
  /* set some parameter */

  if (data_type==DUDLEY_NODES) {
	   type=NODES;
	   basis=reference_element->BasisFunctions;
	   numNodes=Dudley_NodeFile_getNumNodes(nodes);
	   map=Dudley_NodeFile_borrowTargetNodes(nodes);
  } else if (data_type==DUDLEY_REDUCED_NODES) {
	   type=REDUCED_NODES;
	   basis=reference_element->BasisFunctions;
	   numNodes=Dudley_NodeFile_getNumReducedNodes(nodes);
	   map=Dudley_NodeFile_borrowTargetReducedNodes(nodes);
  } else if (data_type==DUDLEY_DEGREES_OF_FREEDOM) {
	   if (elements->MPIInfo->size>1) {
		  Dudley_setError(TYPE_ERROR,"Dudley_Assemble_interpolate: for more than one processor DEGREES_OF_FREEDOM data are not accepted as input.");
		  return;
	   }
	   type=DOF;
	   basis=reference_element->BasisFunctions;	
	   numNodes=Dudley_NodeFile_getNumDegreesOfFreedom(nodes);
	   map=Dudley_NodeFile_borrowTargetDegreesOfFreedom(nodes);
  } else if (data_type==DUDLEY_REDUCED_DEGREES_OF_FREEDOM) {
	   if (elements->MPIInfo->size>1) {
		  Dudley_setError(TYPE_ERROR,"Dudley_Assemble_interpolate: for more than one processor REDUCED_DEGREES_OF_FREEDOM data are not accepted as input.");
		  return;
	   }
	   type=REDUCED_DOF;
	   basis=reference_element->BasisFunctions;
	   numNodes=Dudley_NodeFile_getNumReducedDegreesOfFreedom(nodes);
	   map=Dudley_NodeFile_borrowTargetReducedDegreesOfFreedom(nodes);
   } else {
	   Dudley_setError(TYPE_ERROR,"Dudley_Assemble_interpolate: Cannot interpolate data");
	   return;
  }

  numQuad=basis->numQuadNodes;
  numShapesTotal2=basis->Type->numShapes;
  NS_DOF=basis->Type->numShapes;
  
  /* check the dimensions of interpolated_data and data */

  if (! numSamplesEqual(interpolated_data,numQuad,elements->numElements)) {
	   Dudley_setError(TYPE_ERROR,"Dudley_Assemble_interpolate: illegal number of samples of output Data object");
  } else if (! numSamplesEqual(data,1,numNodes)) {
	   Dudley_setError(TYPE_ERROR,"Dudley_Assemble_interpolate: illegal number of samples of input Data object");
  } else if (numComps!=getDataPointSize(interpolated_data)) {
	   Dudley_setError(TYPE_ERROR,"Dudley_Assemble_interpolate: number of components of input and interpolated Data do not match.");
  }	 else if (!isExpanded(interpolated_data)) {
	   Dudley_setError(TYPE_ERROR,"Dudley_Assemble_interpolate: expanded Data object is expected for output data.");
  }

  /* now we can start */

  if (Dudley_noError())
  {
	requireWrite(interpolated_data);
	#pragma omp parallel private(local_data, numComps_size, isub)
	{
	    local_data=NULL; 
	    /* allocation of work arrays */
	    local_data=THREAD_MEMALLOC(NS_DOF*numComps,double); 
	    if (! Dudley_checkPtr(local_data))
	    {
		numComps_size=(size_t)numComps*sizeof(double);
		/* open the element loop */
		#pragma omp for private(e,q,i,data_array,isub) schedule(static)
		for(e=0;e<elements->numElements;e++)
		{
			for (q=0;q<NS_DOF;q++)
			{
				    i=elements->Nodes[INDEX2(q,e,NN)];
				    data_array=getSampleDataRO(data,map[i]);
				    memcpy(&(local_data[INDEX3(0,q,0, numComps,NS_DOF)]), data_array, numComps_size);
			}
			/*  calculate interpolated_data=local_data*S */
			Dudley_Util_SmallMatSetMult1(1, numComps, numQuad, getSampleDataRW(interpolated_data,e),
				 NS_DOF,local_data,basis->S);
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
