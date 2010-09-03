
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

/*    assemblage jacobiens: calculate the gradient of nodal data at quadrature points */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif
/*****************************************************************/


void Dudley_Assemble_gradient(Dudley_NodeFile* nodes, Dudley_ElementFile* elements,
                              escriptDataC* grad_data,escriptDataC* data) {

  Dudley_ReferenceElement*  refElement=NULL;
  size_t localGradSize=0;
  register dim_t e,q,l,s,n;
  register __const double *data_array;
  register double *grad_data_e;
  dim_t numNodes=0, numShapes=0, numShapesTotal=0, numComps, NN=0, numDim=0, numShapesTotal2=0, numQuad=0;
  type_t data_type=getFunctionSpaceType(data);
  bool_t reducedShapefunction=FALSE, reducedIntegrationOrder=FALSE;
  index_t s_offset=0;
  Dudley_ElementFile_Jacobeans* jac=NULL;
  
  Dudley_resetError();
  if (nodes==NULL || elements==NULL) return;
  numComps=getDataPointSize(data);
  NN=elements->numNodes;
  reducedIntegrationOrder=Dudley_Assemble_reducedIntegrationOrder(grad_data);

  if (data_type==DUDLEY_NODES) {
       reducedShapefunction=FALSE;
       numNodes=nodes->nodesMapping->numTargets;
  } else if (data_type==DUDLEY_REDUCED_NODES) { 
       reducedShapefunction=TRUE;
       numNodes=nodes->reducedNodesMapping->numTargets;
  } else if (data_type==DUDLEY_DEGREES_OF_FREEDOM) {
       if (elements->MPIInfo->size>1) {
          Dudley_setError(TYPE_ERROR,"Dudley_Assemble_gradient: for more than one processor DEGREES_OF_FREEDOM data are not accepted as input.");
          return;
       }
       reducedShapefunction=FALSE;
       numNodes=nodes->degreesOfFreedomMapping->numTargets;
  } else if (data_type==DUDLEY_REDUCED_DEGREES_OF_FREEDOM) {
       if (elements->MPIInfo->size>1) {
          Dudley_setError(TYPE_ERROR,"Dudley_Assemble_gradient: for more than one processor REDUCED_DEGREES_OF_FREEDOM data are not accepted as input.");
          return;
       }
       reducedShapefunction=TRUE;
       numNodes=nodes->reducedDegreesOfFreedomMapping->numTargets;
  } else {
       Dudley_setError(TYPE_ERROR,"Dudley_Assemble_gradient: Cannot calculate gradient of data because of unsuitable input data representation.");
  }

  jac=Dudley_ElementFile_borrowJacobeans(elements,nodes,reducedShapefunction,reducedIntegrationOrder);
  refElement=Dudley_ReferenceElementSet_borrowReferenceElement(elements->referenceElementSet, reducedIntegrationOrder);
  
  if (Dudley_noError()) {

	  numDim=jac->numDim;
          numShapes=jac->BasisFunctions->Type->numShapes;
	  numShapesTotal=jac->numShapesTotal;
	  numQuad=jac->numQuadTotal;
       	  s_offset=jac->offsets[0];
       	  s_offset=jac->offsets[0];
	  localGradSize=sizeof(double)*numDim*numQuad*numComps;
	  if ( (data_type==DUDLEY_REDUCED_NODES) || (DUDLEY_REDUCED_DEGREES_OF_FREEDOM==data_type) )  {
		  numShapesTotal2=refElement->LinearBasisFunctions->Type->numShapes;
	  } else { 
		  numShapesTotal2=refElement->BasisFunctions->Type->numShapes;
	  }
      /* check the dimensions of data */

      if (! numSamplesEqual(grad_data,numQuad,elements->numElements)) {
           Dudley_setError(TYPE_ERROR,"Dudley_Assemble_gradient: illegal number of samples in gradient Data object");
      } else if (! numSamplesEqual(data,1,numNodes)) {
           Dudley_setError(TYPE_ERROR,"Dudley_Assemble_gradient: illegal number of samples of input Data object");
      } else if (numDim*numComps!=getDataPointSize(grad_data)) {
           Dudley_setError(TYPE_ERROR,"Dudley_Assemble_gradient: illegal number of components in gradient data object.");
      }  else if (!isExpanded(grad_data)) {
           Dudley_setError(TYPE_ERROR,"Dudley_Assemble_gradient: expanded Data object is expected for output data.");
      } else if (! (s_offset+numShapes <= numShapesTotal)) {
           Dudley_setError(SYSTEM_ERROR,"Dudley_Assemble_gradient: nodes per element is inconsistent with number of jacobeans.");
      } else if (! (s_offset+numShapes <= numShapesTotal)) {
           Dudley_setError(SYSTEM_ERROR,"Dudley_Assemble_gradient: offset test failed.");
      }
  }
  /* now we can start */

  if (Dudley_noError()) {
      requireWrite(grad_data);
      #pragma omp parallel private(e,q,l,s,n,data_array,grad_data_e)
      {

         if (data_type==DUDLEY_NODES) {
            if (numDim==1) {
                #define DIM 1
                #pragma omp for schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(s,e, NN)];
							data_array=getSampleDataRO(data,n);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
                                      grad_data_e[INDEX4(l,0,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,0,e, numShapesTotal,DIM,numQuad,1)];
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
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(s,e, NN)];
							data_array=getSampleDataRO(data,n);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,0,e, numShapesTotal,DIM,numQuad,1)];
									grad_data_e[INDEX4(l,1,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,0,e, numShapesTotal,DIM,numQuad,1)];
/*printf("data_array of l=%d = %e\n",l,data_array[l]); */
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
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(s,e, NN)];
							data_array=getSampleDataRO(data,n);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,0,e, numShapesTotal,DIM,numQuad,1)];
									grad_data_e[INDEX4(l,1,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,0,e, numShapesTotal,DIM,numQuad,1)];
									grad_data_e[INDEX4(l,2,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,2,q,0,e, numShapesTotal,DIM,numQuad,1)];
								}
							}
						}
                }
                #undef DIM
            }
         } else if (data_type==DUDLEY_REDUCED_NODES) {
            if (numDim==1) {
                #define DIM 1
                #pragma omp for schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(s,e, NN)];
							data_array=getSampleDataRO(data,nodes->reducedNodesMapping->target[n]);            
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {								
									grad_data_e[INDEX4(l,0,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,0,e, numShapesTotal,DIM,numQuad,1)];
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
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(s,e, NN)];
							data_array=getSampleDataRO(data,nodes->reducedNodesMapping->target[n]);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,0,e, numShapesTotal,DIM,numQuad,1)];
									grad_data_e[INDEX4(l,1,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,0,e, numShapesTotal,DIM,numQuad,1)];
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
                    	for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(s,e, NN)];
							data_array=getSampleDataRO(data,nodes->reducedNodesMapping->target[n]);
							for (q=0;q<numQuad;q++) {	
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,0,e, numShapesTotal,DIM,numQuad,1)];
									grad_data_e[INDEX4(l,1,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,0,e, numShapesTotal,DIM,numQuad,1)];
									grad_data_e[INDEX4(l,2,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,2,q,0,e, numShapesTotal,DIM,numQuad,1)];
								}
							}
						}
                }
                #undef DIM
            }
         } else if (data_type==DUDLEY_DEGREES_OF_FREEDOM) {

            if (numDim==1) {
                #define DIM 1
                #pragma omp for schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(s,e, NN)];
							data_array=getSampleDataRO(data,nodes->degreesOfFreedomMapping->target[n]);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,0,e, numShapesTotal,DIM,numQuad,1)];
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
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(s,e, NN)];
							data_array=getSampleDataRO(data,nodes->degreesOfFreedomMapping->target[n]);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,0,e, numShapesTotal,DIM,numQuad,1)];
                           	        grad_data_e[INDEX4(l,1,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,0,e, numShapesTotal,DIM,numQuad,1)];
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
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(s,e, NN)];
							data_array=getSampleDataRO(data,nodes->degreesOfFreedomMapping->target[n]);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,0,e, numShapesTotal,DIM,numQuad,1)];
									grad_data_e[INDEX4(l,1,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,0,e, numShapesTotal,DIM,numQuad,1)];
									grad_data_e[INDEX4(l,2,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,2,q,0,e, numShapesTotal,DIM,numQuad,1)];
								}
							}
					}
				}
                #undef DIM
            }
         } else if (data_type==DUDLEY_REDUCED_DEGREES_OF_FREEDOM) {
            if (numDim==1) {
                #define DIM 1
                #pragma omp for schedule(static)
				for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleDataRW(grad_data,e);
                    memset(grad_data_e,0, localGradSize);
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(s,e, NN)];
							data_array=getSampleDataRO(data,nodes->reducedDegreesOfFreedomMapping->target[n]);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,0,e, numShapesTotal,DIM,numQuad,1)];
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
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(s,e, NN)];
							data_array=getSampleDataRO(data,nodes->reducedDegreesOfFreedomMapping->target[n]);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,0,e, numShapesTotal,DIM,numQuad,1)];
									grad_data_e[INDEX4(l,1,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,0,e, numShapesTotal,DIM,numQuad,1)];
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
						for (s=0;s<numShapes;s++) {
							n=elements->Nodes[INDEX2(s,e, NN)];
							data_array=getSampleDataRO(data,nodes->reducedDegreesOfFreedomMapping->target[n]);
							for (q=0;q<numQuad;q++) {
								#pragma ivdep
								for (l=0;l<numComps;l++) {
									grad_data_e[INDEX4(l,0,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,0,q,0,e, numShapesTotal,DIM,numQuad,1)];
									grad_data_e[INDEX4(l,1,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,1,q,0,e, numShapesTotal,DIM,numQuad,1)];
									grad_data_e[INDEX4(l,2,q,0, numComps,DIM,numQuad)]+=data_array[l]*jac->DSDX[INDEX5(s_offset+s,2,q,0,e, numShapesTotal,DIM,numQuad,1)];
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
