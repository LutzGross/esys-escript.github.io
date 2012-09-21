
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/************************************************************************************/

/*    assemblage routines: integrates data on quadrature points   */

/************************************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/************************************************************************************/

void Finley_Assemble_integrate(Finley_NodeFile* nodes, Finley_ElementFile* elements,escriptDataC* data,double* out) {
/*    type_t data_type=getFunctionSpaceType(data);*/
    dim_t numQuadTotal;
    dim_t numComps=getDataPointSize(data);
    Finley_ElementFile_Jacobeans* jac=NULL;
    Esys_MPI_rank my_mpi_rank;
    
    Finley_resetError();
    if (nodes==NULL || elements==NULL) return;
    my_mpi_rank = nodes->MPIInfo->rank;
    jac=Finley_ElementFile_borrowJacobeans(elements,nodes,FALSE,Finley_Assemble_reducedIntegrationOrder(data));
    if (Finley_noError()) {
		numQuadTotal=jac->numQuadTotal;
        /* check the shape of the data  */
        if (! numSamplesEqual(data,numQuadTotal,elements->numElements)) {
           Finley_setError(TYPE_ERROR,"Finley_Assemble_integrate: illegal number of samples of integrant kernel Data object");
        }
        /* now we can start */

        if (Finley_noError()) {
            dim_t q,e,i;
        	__const double *data_array=NULL;
            double *out_local=NULL, rtmp;
            for (q=0;q<numComps;q++) out[q]=0;
            #pragma omp parallel private(q,i,rtmp,data_array,out_local)
            {
                out_local=THREAD_MEMALLOC(numComps,double);
                if (! Finley_checkPtr(out_local) ) {
                   /* initialize local result */

				   for (i=0;i<numComps;i++) out_local[i]=0;

                   /* open the element loop */
   
                   if (isExpanded(data)) {
                       #pragma omp for private(e) schedule(static)
                       for(e=0;e<elements->numElements;e++) {
                          if (elements->Owner[e] == my_mpi_rank) {
                            data_array=getSampleDataRO(data,e);
                            for (q=0;q<numQuadTotal;q++) {
                                  for (i=0;i<numComps;i++) out_local[i]+=data_array[INDEX2(i,q,numComps)]*jac->volume[INDEX2(q,e,numQuadTotal)];
                            }
                         }
                       }
                   } else {
                      #pragma omp for private(e) schedule(static)
                      for(e=0;e<elements->numElements;e++) {
                          if (elements->Owner[e] == my_mpi_rank) {
                           data_array=getSampleDataRO(data,e);
                           rtmp=0.;
                           for (q=0;q<numQuadTotal;q++) rtmp+=jac->volume[INDEX2(q,e,numQuadTotal)];
                           for (i=0;i<numComps;i++) out_local[i]+=data_array[i]*rtmp;
                          }
                      }
                   }
                   /* add local results to global result */
                   #pragma omp critical
                   for (i=0;i<numComps;i++) out[i]+=out_local[i];
                }
                THREAD_MEMFREE(out_local);
            }
		}
	}
}

