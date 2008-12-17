
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*    assemblage routines: integrates data on quadrature points   */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/

void Finley_Assemble_integrate(Finley_NodeFile* nodes, Finley_ElementFile* elements,escriptDataC* data,double* out) {
/*    type_t data_type=getFunctionSpaceType(data);*/
    dim_t numComps=getDataPointSize(data);
    Finley_ElementFile_Jacobeans* jac=NULL;
    Paso_MPI_rank my_mpi_rank;
    
    Finley_resetError();
    if (nodes==NULL || elements==NULL) return;
    my_mpi_rank = nodes->MPIInfo->rank;
    /* set some parameter */
    jac=Finley_ElementFile_borrowJacobeans(elements,nodes,FALSE,Finley_Assemble_reducedIntegrationOrder(data));
    if (Finley_noError()) {

        /* check the shape of the data  */
        if (! numSamplesEqual(data,jac->ReferenceElement->numQuadNodes,elements->numElements)) {
           Finley_setError(TYPE_ERROR,"Finley_Assemble_integrate: illegal number of samples of integrant kernel Data object");
        }
        /* now we can start */

        if (Finley_noError()) {
            dim_t q,e,i;
            double *out_local=NULL, rtmp,*data_array=NULL;
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
                            data_array=getSampleData(data,e);
                            for (q=0;q<jac->ReferenceElement->numQuadNodes;q++) {
                                  for (i=0;i<numComps;i++) out_local[i]+=data_array[INDEX2(i,q,numComps)]*jac->volume[INDEX2(q,e,jac->ReferenceElement->numQuadNodes)];
                            }
                         }
                       }
                   } else {
                      #pragma omp for private(e) schedule(static)
                      for(e=0;e<elements->numElements;e++) {
                          if (elements->Owner[e] == my_mpi_rank) {
                           data_array=getSampleData(data,e);
                           rtmp=0.;
                           for (q=0;q<jac->ReferenceElement->numQuadNodes;q++) rtmp+=jac->volume[INDEX2(q,e,jac->ReferenceElement->numQuadNodes)];
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
 * Revision 1.1.1.1.2.2  2005/06/29 02:34:48  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1.2.1  2004/11/24 01:37:12  gross
 * some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
 *
 *
 *
 */
