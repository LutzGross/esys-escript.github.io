
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

/*    assemblage jacobians: calculate the gradient of nodal data at quadrature points */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif
/*****************************************************************/


void Finley_Assemble_gradient(Finley_NodeFile* nodes, Finley_ElementFile* elements,
                              escriptDataC* grad_data,escriptDataC* data) {

  Finley_ReferenceElement*  refElement=NULL;
  size_t localGradSize=0;
  register dim_t e,q,l,s,n;
  register __const double *data_array;
  register double *grad_data_e;
  dim_t numNodes=0, numShapes=0, numShapesTotal=0, numComps, NN=0, numDim=0, numShapesTotal2=0, numQuad=0, numSub=0, isub=0;
  type_t data_type=getFunctionSpaceType(data);
  bool_t reducedShapefunction=FALSE, reducedIntegrationOrder=FALSE;
  index_t s_offset=0,  *nodes_selector=NULL;
  Finley_ElementFile_Jacobeans* jac=NULL;
  type_t grad_data_type=getFunctionSpaceType(grad_data);
  
  Finley_resetError();
  if (nodes==NULL || elements==NULL) return;
  numComps=getDataPointSize(data);
  NN=elements->numNodes;
  reducedIntegrationOrder=Finley_Assemble_reducedIntegrationOrder(grad_data);

  if (data_type==FINLEY_NODES) {
       reducedShapefunction=FALSE;
       numNodes=nodes->nodesMapping->numTargets;
  } else if (data_type==FINLEY_REDUCED_NODES) { 
       reducedShapefunction=TRUE;
       numNodes=nodes->reducedNodesMapping->numTargets;
  } else if (data_type==FINLEY_DEGREES_OF_FREEDOM) {
       if (elements->MPIInfo->size>1) {
          Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: for more than one processor DEGREES_OF_FREEDOM data are not accepted as input.");
          return;
       }
       reducedShapefunction=FALSE;
       numNodes=nodes->degreesOfFreedomMapping->numTargets;
  } else if (data_type==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
       if (elements->MPIInfo->size>1) {
          Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: for more than one processor REDUCED_DEGREES_OF_FREEDOM data are not accepted as input.");
          return;
       }
       reducedShapefunction=TRUE;
       numNodes=nodes->reducedDegreesOfFreedomMapping->numTargets;
  } else {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: Cannot calculate gradient of data because of unsuitable input data representation.");
  }

  jac=Finley_ElementFile_borrowJacobeans(elements,nodes,reducedShapefunction,reducedIntegrationOrder);
  refElement=Finley_ReferenceElementSet_borrowReferenceElement(elements->referenceElementSet, reducedIntegrationOrder);
  
  if (Finley_noError()) {

	  numDim=jac->numDim;
          numShapes=jac->BasisFunctions->Type->numShapes;
	  numShapesTotal=jac->numShapesTotal;
	  numSub=jac->numSub;
	  numQuad=jac->numQuadTotal/numSub;
      	  if (grad_data_type==FINLEY_CONTACT_ELEMENTS_2 || grad_data_type== FINLEY_REDUCED_CONTACT_ELEMENTS_2)  {
       	 	   s_offset=jac->offsets[1];
       	 	   s_offset=jac->offsets[1];
      	  } else {
       	  	  s_offset=jac->offsets[0];
       	          s_offset=jac->offsets[0];
      	  }
	  localGradSize=sizeof(double)*numDim*numQuad*numSub*numComps;
	  if ( (data_type==FINLEY_REDUCED_NODES) || (FINLEY_REDUCED_DEGREES_OF_FREEDOM==data_type) )  {
		  nodes_selector=refElement->Type->linearNodes;
		  numShapesTotal2=refElement->LinearBasisFunctions->Type->numShapes * refElement->Type->numSides;
	  } else { 
		  nodes_selector=refElement->Type->subElementNodes;
		  numShapesTotal2=refElement->BasisFunctions->Type->numShapes * refElement->Type->numSides;
	  }
      /* check the dimensions of data */

      if (! numSamplesEqual(grad_data,numQuad*numSub,elements->numElements)) {
           Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: illegal number of samples in gradient Data object");
      } else if (! numSamplesEqual(data,1,numNodes)) {
           Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: illegal number of samples of input Data object");
      } else if (numDim*numComps!=getDataPointSize(grad_data)) {
           Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: illegal number of components in gradient data object.");
      }  else if (!isExpanded(grad_data)) {
           Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: expanded Data object is expected for output data.");
      } else if (! (s_offset+numShapes <= numShapesTotal)) {
           Finley_setError(SYSTEM_ERROR,"Finley_Assemble_gradient: nodes per element is inconsistent with number of jacobians.");
      } else if (! (s_offset+numShapes <= numShapesTotal)) {
           Finley_setError(SYSTEM_ERROR,"Finley_Assemble_gradient: offset test failed.");
      }
  }
  /* now we can start */

  if (Finley_noError()) {
      requireWrite(grad_data);
      #pragma omp parallel private(e,q,l,s,n,data_array,grad_data_e, isub)
      {

         if (data_type==FINLEY_NODES) {
            if (numDim==1) {
                #define DIM 1
                #pragma omp for schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
					for (isub=0; isub<numSub; isub++) {
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
							data_array=getSampleDataRO(data,n);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
                                      grad_data_e[INDEX4(l,0,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
								}
							}
						}
					}
                }
                #undef DIM
            } else if (numDim==2) {
                #define DIM 2
                #pragma omp for schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
					for (isub=0; isub<numSub; isub++) {
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
							data_array=getSampleDataRO(data,n);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
									grad_data_e[INDEX4(l,1,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
/*printf("data_array of l=%d = %e\n",l,data_array[l]); */
								}
							}
                       }
                    }
                }
                #undef DIM
            } else if (numDim==3) {
                #define DIM 3
                #pragma omp for private(e,grad_data_e,s,n,data_array,q,l) schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e); 
                    memset(grad_data_e,0, localGradSize);
					for (isub=0; isub<numSub; isub++) {
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
							data_array=getSampleDataRO(data,n);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
									grad_data_e[INDEX4(l,1,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
									grad_data_e[INDEX4(l,2,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,2,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
								}
							}
						}
                    }
                }
                #undef DIM
            }
         } else if (data_type==FINLEY_REDUCED_NODES) {
            if (numDim==1) {
                #define DIM 1
                #pragma omp for schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
					for (isub=0; isub<numSub; isub++) {
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
							data_array=getSampleDataRO(data,nodes->reducedNodesMapping->target[n]);            
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {								
									grad_data_e[INDEX4(l,0,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
								}
							}
						}
                    }
                }
                #undef DIM
            } else if (numDim==2) {
                #define DIM 2
                #pragma omp for schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
					for (isub=0; isub<numSub; isub++) {
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
							data_array=getSampleDataRO(data,nodes->reducedNodesMapping->target[n]);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
									grad_data_e[INDEX4(l,1,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
								}
							}
						}
                    }
                }
                #undef DIM
            } else if (numDim==3) {
                #define DIM 3
                #pragma omp for schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
					for (isub=0; isub<numSub; isub++) {
                    	for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
							data_array=getSampleDataRO(data,nodes->reducedNodesMapping->target[n]);
							for (q=0;q<numQuad;q++) {	
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
									grad_data_e[INDEX4(l,1,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
									grad_data_e[INDEX4(l,2,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,2,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
								}
							}
						}
					}
                }
                #undef DIM
            }
         } else if (data_type==FINLEY_DEGREES_OF_FREEDOM) {

            if (numDim==1) {
                #define DIM 1
                #pragma omp for schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
					for (isub=0; isub<numSub; isub++) {
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
							data_array=getSampleDataRO(data,nodes->degreesOfFreedomMapping->target[n]);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
								}
							}
						}
                    }
                }
                #undef DIM
            } else if (numDim==2) {
                #define DIM 2
                #pragma omp for schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
					for (isub=0; isub<numSub; isub++) {
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
							data_array=getSampleDataRO(data,nodes->degreesOfFreedomMapping->target[n]);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
                           	        grad_data_e[INDEX4(l,1,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
								}
							}
						}
					}
				}
                #undef DIM
            } else if (numDim==3) {
                #define DIM 3
                #pragma omp for schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
					for (isub=0; isub<numSub; isub++) {
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
							data_array=getSampleDataRO(data,nodes->degreesOfFreedomMapping->target[n]);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
									grad_data_e[INDEX4(l,1,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
									grad_data_e[INDEX4(l,2,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,2,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
								}
							}
						}
					}
				}
                #undef DIM
            }
         } else if (data_type==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            if (numDim==1) {
                #define DIM 1
                #pragma omp for schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
					for (isub=0; isub<numSub; isub++) {
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
							data_array=getSampleDataRO(data,nodes->reducedDegreesOfFreedomMapping->target[n]);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
								}
							}
						}
					}
				}
                #undef DIM
            } else if (numDim==2) {
                #define DIM 2
                #pragma omp for schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
					for (isub=0; isub<numSub; isub++) {
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
							data_array=getSampleDataRO(data,nodes->reducedDegreesOfFreedomMapping->target[n]);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
									grad_data_e[INDEX4(l,1,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
								}
							}
						}
                    }
                }
                #undef DIM

            } else if (numDim==3) {
                #define DIM 3
                #pragma omp for schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
					for (isub=0; isub<numSub; isub++) {
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(nodes_selector[INDEX2(s_offset+s,isub,numShapesTotal2)],e, NN)];
							data_array=getSampleDataRO(data,nodes->reducedDegreesOfFreedomMapping->target[n]);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
									grad_data_e[INDEX4(l,1,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
									grad_data_e[INDEX4(l,2,q,isub, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,2,q,isub,e, numShapesTotal,DIM,numQuad,numSub)];
								}
							}
						}
                    }
                }
                #undef DIM
            }
         }
      } /* end parallel region */
  }
}
