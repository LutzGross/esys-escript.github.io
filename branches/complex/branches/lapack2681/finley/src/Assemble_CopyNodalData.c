
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*    assemblage routines: copies data between different types nodal representation   */

/**************************************************************/

#include "Util.h"
#include "Assemble.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void Finley_Assemble_CopyNodalData(Finley_NodeFile* nodes,escriptDataC* out,escriptDataC* in) {
    dim_t n,k,l, mpiSize;
    dim_t numComps=getDataPointSize(out);
    Paso_Coupler *coupler=NULL;
    type_t in_data_type=getFunctionSpaceType(in);
    type_t out_data_type=getFunctionSpaceType(out);
    index_t upperBound;
    double* recv_buffer;
    size_t numComps_size=0;
    Finley_resetError();
    if (nodes==NULL) return;
    mpiSize=nodes->MPIInfo->size;

    /* check out and in */
    if (numComps!=getDataPointSize(in)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: number of components of input and output Data do not match.");
    } else if (!isExpanded(out)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: expanded Data object is expected for output data.");
    }

    /* more sophisticated test needed for overlapping node/DOF counts */
    if (in_data_type == FINLEY_NODES) {
        if (! numSamplesEqual(in,1,Finley_NodeFile_getNumNodes(nodes))) {
               Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of input Data object");
       }
    } else if (in_data_type == FINLEY_REDUCED_NODES) {
       if (! numSamplesEqual(in,1,Finley_NodeFile_getNumReducedNodes(nodes))) {
               Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of input Data object");
       }
    } else if (in_data_type == FINLEY_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(in,1,Finley_NodeFile_getNumDegreesOfFreedom(nodes))) {
               Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of input Data object");
       }
       if ( (((out_data_type == FINLEY_NODES) || (out_data_type == FINLEY_DEGREES_OF_FREEDOM)) && !isExpanded(in) && (mpiSize>1))) {

               Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: FINLEY_DEGREES_OF_FREEDOM to FINLEY_NODES or FINLEY_DEGREES_OF_FREEDOM requires expanded input data on more than one processor.");
       }
    } else if (in_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(in,1,Finley_NodeFile_getNumReducedDegreesOfFreedom(nodes))) {
            Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of input Data object");
       }
       if ( (out_data_type == FINLEY_DEGREES_OF_FREEDOM) && !isExpanded(in) && (mpiSize>1)) {

               Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: FINLEY_REDUCED_DEGREES_OF_FREEDOM to FINLEY_DEGREES_OF_FREEDOM requires expanded input data on more than one processor.");
       }
    } else {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal function space type for target object");
    }
    
    if (out_data_type == FINLEY_NODES) {
        if (! numSamplesEqual(out,1,Finley_NodeFile_getNumNodes(nodes))) {
               Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of output Data object");
       }
    } else if (out_data_type == FINLEY_REDUCED_NODES) {
       if (! numSamplesEqual(out,1,Finley_NodeFile_getNumReducedNodes(nodes))) {
               Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of output Data object");
       }
    } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(out,1,Finley_NodeFile_getNumDegreesOfFreedom(nodes))) {
               Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of output Data object");
       }
    } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(out,1,Finley_NodeFile_getNumReducedDegreesOfFreedom(nodes))) {
               Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal number of samples of output Data object");
       }
    } else {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: illegal function space type for source object");
    }

    /* now we can start */
 
    if (Finley_noError()) {
        /*********************** FINLEY_NODES **************************************************/
        numComps_size=(size_t)numComps*sizeof(double);
        if (in_data_type == FINLEY_NODES) {
	   requireWrite(out);
           if  (out_data_type == FINLEY_NODES) {
	      void* buff=allocSampleBuffer(in);
	      #pragma omp parallel private(n)
	      {

        	#pragma omp parallel for private(n) schedule(static)
              	for (n=0;n<nodes->nodesMapping->numNodes;n++) {
                   memcpy(getSampleDataRWFast(out,n), getSampleDataROFast(in,n,buff), numComps_size);
              	}
	      }
	      freeSampleBuffer(buff);
           } else if (out_data_type == FINLEY_REDUCED_NODES) {
	      void* buff=allocSampleBuffer(in);
	      #pragma omp parallel private(n)
	      {
        	#pragma omp for schedule(static)
                for (n=0;n<nodes->reducedNodesMapping->numTargets;n++) {
                   memcpy(getSampleDataRWFast(out,n),
                          getSampleDataROFast(in,nodes->reducedNodesMapping->map[n],buff),
                          numComps_size);
		}
              }
	      freeSampleBuffer(buff);
           } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
	      int nComps = Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution);
	      void* buff=allocSampleBuffer(in);
	      #pragma omp parallel private(n)
	      {
        	#pragma omp for schedule(static)
              	for (n=0;n<nComps;n++) {
                   memcpy(getSampleDataRWFast(out,n),
                          getSampleDataROFast(in,nodes->degreesOfFreedomMapping->map[n],buff),
                          numComps_size);
		}
              }
	      freeSampleBuffer(buff);
           } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
	      int nComps = Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
	      void* buff=allocSampleBuffer(in);
	      #pragma omp parallel private(n)
	      {
        	#pragma omp for schedule(static)
              	for (n=0;n<nComps;n++) {
                   memcpy(getSampleDataRWFast(out,n),
                          getSampleDataROFast(in,nodes->reducedDegreesOfFreedomMapping->map[n],buff),
                          numComps_size);
		}
              }
	      freeSampleBuffer(buff);
           }
        /*********************** FINLEY_REDUCED_NODES **************************************************/
        } else if (in_data_type == FINLEY_REDUCED_NODES) {
	   requireWrite(out);
           if  (out_data_type == FINLEY_NODES) {
             Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: cannot copy from reduced nodes to nodes.");

           } else if (out_data_type == FINLEY_REDUCED_NODES) {
	      void* buff=allocSampleBuffer(in);
	      #pragma omp parallel private(n)
	      {
        	#pragma omp for schedule(static)
              	for (n=0;n<nodes->reducedNodesMapping->numNodes;n++) {
                   memcpy(getSampleDataRWFast(out,n),getSampleDataROFast(in,n,buff),numComps_size);
		}
              }
	      freeSampleBuffer(buff);
           } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
             Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: cannot copy from reduced nodes to degrees of freedom.");
           } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
	      int nComps = Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
	      void* buff=allocSampleBuffer(in);
	      #pragma omp parallel private(n,k)
	      {
        	#pragma omp for schedule(static)
              	for (n=0;n<nComps;n++) {
                   k=nodes->reducedDegreesOfFreedomMapping->map[n];
                   memcpy(getSampleDataRWFast(out,n),
                          getSampleDataROFast(in,nodes->reducedNodesMapping->target[k],buff),
                          numComps_size);
		}
              }
	      freeSampleBuffer(buff);
           }

        /*********************** FINLEY_DEGREES_OF_FREEDOM **************************************************/
        } else if (in_data_type == FINLEY_DEGREES_OF_FREEDOM) {
	    requireWrite(out);
            if  (out_data_type == FINLEY_NODES) {
               coupler=Paso_Coupler_alloc(nodes->degreesOfFreedomConnector,numComps);
               if (Paso_noError()) {
		    void* buff=allocSampleBuffer(in);
		    /* It is not immediately clear whether coupler can be trusted with constant data so I'll assume RW */
		    /* Also, it holds pointers so it might not be safe to use on lazy data anyway?*/
		    requireWrite(in);
                    Paso_Coupler_startCollect(coupler,getDataRW(in));
                    recv_buffer=Paso_Coupler_finishCollect(coupler);
                    upperBound=Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution);
		    #pragma omp parallel private(n,k)
		    {
			#pragma omp for schedule(static)
                    	for (n=0;n<nodes->numNodes;n++) {
                          k=nodes->degreesOfFreedomMapping->target[n];
                          if (k < upperBound) {
                                memcpy(getSampleDataRWFast(out,n),
                                       getSampleDataROFast(in,k,buff),
                                       numComps_size);
                           } else {
                                memcpy(getSampleDataRWFast(out,n),
                                       &recv_buffer[(k-upperBound)*numComps],
                                       numComps_size);
                           }
                    	}
		    }
		    freeSampleBuffer(buff);
               }
               Paso_Coupler_free(coupler);
            } else if  (out_data_type == FINLEY_REDUCED_NODES) {
               coupler=Paso_Coupler_alloc(nodes->degreesOfFreedomConnector,numComps);
               if (Paso_noError()) {
		    void* buff=allocSampleBuffer(in);
		    requireWrite(in);	/* See comment above about coupler and const */
                    Paso_Coupler_startCollect(coupler,getDataRW(in));
                    recv_buffer=Paso_Coupler_finishCollect(coupler);
                    upperBound=Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution);
		    requireWrite(out);

		    #pragma omp parallel private(n,k,l)
		    {
			#pragma omp for schedule(static)
                    	for (n=0;n<nodes->reducedNodesMapping->numTargets;n++) {
                          l=nodes->reducedNodesMapping->map[n];
                          k=nodes->degreesOfFreedomMapping->target[l];
                          if (k < upperBound) {
                                memcpy(getSampleDataRWFast(out,n),
                                       getSampleDataROFast(in,k,buff),
                                       numComps_size);
                           } else {
                                memcpy(getSampleDataRWFast(out,n),
                                       &recv_buffer[(k-upperBound)*numComps],
                                       numComps_size);
                           }
                    	}
		    }
		    freeSampleBuffer(buff);
               }
               Paso_Coupler_free(coupler);
           } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
	      void* buff=allocSampleBuffer(in);
	      int nComps = Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution);
	      requireWrite(out);
	      #pragma omp parallel private(n)
	      {
		#pragma omp for schedule(static)
              	for (n=0;n<nComps;n++) {
                    memcpy(getSampleDataRWFast(out,n),getSampleDataROFast(in,n,buff),numComps_size);
              	}
	      }
	      freeSampleBuffer(buff);
           } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
	      void* buff=allocSampleBuffer(in);
	      int nComps = Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
	      requireWrite(out);
	      #pragma omp parallel private(n,k)
	      {
		#pragma omp for schedule(static)
              	for (n=0;n<nComps;n++) {
                   k=nodes->reducedDegreesOfFreedomMapping->map[n];
                   memcpy(getSampleDataRWFast(out,n),
                          getSampleDataROFast(in,nodes->degreesOfFreedomMapping->target[k],buff),
                          numComps_size);
              	}
	      }
	      freeSampleBuffer(buff);
           }

        /*********************** FINLEY_REDUCED_DEGREES_OF_FREEDOM **************************************************/
        } else if (in_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {

           if  (out_data_type == FINLEY_NODES) {
             Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: cannot copy from reduced degrees of freedom to nodes.");
           } else if (out_data_type == FINLEY_REDUCED_NODES) {
               coupler=Paso_Coupler_alloc(nodes->reducedDegreesOfFreedomConnector,numComps);
               if (Paso_noError()) {
		    void* buff=allocSampleBuffer(in);
                    upperBound=Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
		    requireWrite(in);			/* See comment about coupler and const */
                    Paso_Coupler_startCollect(coupler,getDataRW(in));
                    recv_buffer=Paso_Coupler_finishCollect(coupler);
		    requireWrite(out);
		    #pragma omp parallel private(n,k,l)
		    {
			#pragma omp for schedule(static)
                    	for (n=0;n<nodes->reducedNodesMapping->numTargets;n++) {
                          l=nodes->reducedNodesMapping->map[n];
                          k=nodes->reducedDegreesOfFreedomMapping->target[l];
                          if (k < upperBound) {
                                memcpy(getSampleDataRWFast(out,n),
                                       getSampleDataROFast(in,k,buff),
                                       numComps_size);
                           } else {
                                memcpy(getSampleDataRWFast(out,n),
                                       &recv_buffer[(k-upperBound)*numComps],
                                       numComps_size);
                           }
                    	}
		    }
		    freeSampleBuffer(buff);
               }
               Paso_Coupler_free(coupler);
           } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
	      int nComps = Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
	      void* buff=allocSampleBuffer(in);
	      requireWrite(out);
	      #pragma omp parallel private(n)
	      {
		#pragma omp for schedule(static)
              	for (n=0;n<nComps;n++) {
                    memcpy(getSampleDataRWFast(out,n),getSampleDataROFast(in,n,buff),numComps_size);
              	}
	      }
           } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM ) {
             Finley_setError(TYPE_ERROR,"Finley_Assemble_CopyNodalData: cannot copy from reduced degrees of freedom to degrees of freedom.");
           }
        }
   }
   return;
}

