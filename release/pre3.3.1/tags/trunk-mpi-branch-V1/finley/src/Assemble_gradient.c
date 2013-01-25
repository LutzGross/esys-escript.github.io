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

/*    assemblage rjacines: calculate the gradient of nodal data at quadrature points */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004,2005 */
/*   author: gross@access.edu.au */
/*   version: $Id$ */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif
/*****************************************************************/


void Finley_Assemble_gradient(Finley_NodeFile* nodes, Finley_ElementFile* elements,
                              escriptDataC* grad_data,escriptDataC* data) {

  dim_t numNodes, numShapes, numLocalNodes, numComps, NN;
  type_t data_type=getFunctionSpaceType(data);
  type_t grad_data_type=getFunctionSpaceType(grad_data);
  bool_t reducedShapefunction=FALSE, reducedIntegrationOrder=FALSE;
  index_t dof_offset, s_offset;
  Finley_ElementFile_Jacobeans* jac=NULL;
  Finley_resetError();
  if (nodes==NULL || elements==NULL) return;
  numComps=getDataPointSize(data);
  NN=elements->ReferenceElement->Type->numNodes;
  reducedIntegrationOrder=Finley_Assemble_reducedIntegrationOrder(grad_data);

  if (data_type==FINLEY_NODES) {
       reducedShapefunction=FALSE;
       numNodes=nodes->numNodes;
  } else if (data_type==FINLEY_REDUCED_NODES) { /* TODO */
       reducedShapefunction=FALSE;
       Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: FINLEY_REDUCED_NODES is not supported yet.");
  }
  /* lock these two options jac for the MPI version */
  #ifndef PASO_MPI
  else if (data_type==FINLEY_DEGREES_OF_FREEDOM) {
       reducedShapefunction=FALSE;
       numNodes=nodes->numDegreesOfFreedom;
  } else if (data_type==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
       reducedShapefunction=TRUE;
       numNodes=nodes->reducedNumDegreesOfFreedom;
  } 
  #endif
  else {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: Cannot calculate gradient of data because of unsuitable input data representation.");
  }

  jac=Finley_ElementFile_borrowJacobeans(elements,nodes,reducedShapefunction,reducedIntegrationOrder);
  if (Finley_noError()) {

      numShapes=jac->ReferenceElement->Type->numShapes;
      numLocalNodes=jac->ReferenceElement->Type->numNodes;
      if (grad_data_type==FINLEY_CONTACT_ELEMENTS_2 || grad_data_type== FINLEY_REDUCED_CONTACT_ELEMENTS_2)  {
       dof_offset=numShapes;
       s_offset=jac->ReferenceElement->Type->numShapes;
      } else {
       dof_offset=0;
       s_offset=0;
      }

      /* check the dimensions of data */

      if (! numSamplesEqual(grad_data,jac->ReferenceElement->numQuadNodes,elements->numElements)) {
           Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: illegal number of samples in gradient Data object");
      } else if (! numSamplesEqual(data,1,numNodes)) {
           Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: illegal number of samples of input Data object");
      } else if (jac->numDim*numComps!=getDataPointSize(grad_data)) {
           Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: illegal number of components in gradient data object.");
      }  else if (!isExpanded(grad_data)) {
           Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: expanded Data object is expected for output data.");
      } else if (! (dof_offset+jac->ReferenceElement->Type->numShapes <= numLocalNodes)) {
           Finley_setError(SYSTEM_ERROR,"Finley_Assemble_gradient: nodes per element is inconsistent with number of jacobeans.");
      } else if (! (s_offset+jac->ReferenceElement->Type->numShapes <= jac->ReferenceElement->Type->numNodes)) {
           Finley_setError(SYSTEM_ERROR,"Finley_Assemble_gradient: offset test failed.");
      }
  }
  /* now we can start */

  if (Finley_noError()) {
      register dim_t e,q,l,s,n;
      register double* data_array,  *grad_data_e;
      #pragma omp parallel private(e,q,l,s,n,data_array,grad_data_e)
      {
         if (data_type==FINLEY_NODES) {
            if (jac->numDim==1) {
                #define DIM 1
                #pragma omp for schedule(static)
   	        for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleData(grad_data,e);
                    for (q=0;q<DIM*(jac->ReferenceElement->numQuadNodes)*numComps; q++) grad_data_e[q]=0;
                    for (s=0;s<jac->ReferenceElement->Type->numShapes;s++) {
                       n=elements->Nodes[INDEX2(dof_offset+s,e,NN)];
                       data_array=getSampleData(data,n);
                       for (q=0;q<jac->ReferenceElement->numQuadNodes;q++) {
                           for (l=0;l<numComps;l++) {
                                grad_data_e[INDEX3(l,0,q,numComps,DIM)]+=data_array[l]* 
                                      jac->DSDX[INDEX4(s_offset+s,0,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                           }
                       }
                    }
                }

                #undef DIM
            } else if (jac->numDim==2) {
                #define DIM 2
                #pragma omp for schedule(static)
   	        for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleData(grad_data,e);
                    for (q=0;q<DIM*(jac->ReferenceElement->numQuadNodes)*numComps; q++) grad_data_e[q]=0;
                    for (s=0;s<jac->ReferenceElement->Type->numShapes;s++) {
                       n=elements->Nodes[INDEX2(dof_offset+s,e,NN)];
                       data_array=getSampleData(data,n);
                       for (q=0;q<jac->ReferenceElement->numQuadNodes;q++) {
                           for (l=0;l<numComps;l++) {
                               grad_data_e[INDEX3(l,0,q,numComps,DIM)]+=data_array[l]* 
                                        jac->DSDX[INDEX4(s_offset+s,0,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                               grad_data_e[INDEX3(l,1,q,numComps,DIM)]+=data_array[l]* 
                                        jac->DSDX[INDEX4(s_offset+s,1,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                           }
                       }
                    }
                }
                #undef DIM
            } else if (jac->numDim==3) {
                #define DIM 3
                #pragma omp for private(e,grad_data_e,s,n,data_array,q,l) schedule(static)
   	        for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleData(grad_data,e); 
                    for (q=0;q<DIM*(jac->ReferenceElement->numQuadNodes)*numComps; q++) grad_data_e[q]=0;
                    for (s=0;s<jac->ReferenceElement->Type->numShapes;s++) {
                       n=elements->Nodes[INDEX2(dof_offset+s,e,NN)];
                       data_array=getSampleData(data,n);
                       for (q=0;q<jac->ReferenceElement->numQuadNodes;q++) {
                           for (l=0;l<numComps;l++) {
                               grad_data_e[INDEX3(l,0,q,numComps,DIM)]+=data_array[l]*
                                    jac->DSDX[INDEX4(s_offset+s,0,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                               grad_data_e[INDEX3(l,1,q,numComps,DIM)]+=data_array[l]*
                                    jac->DSDX[INDEX4(s_offset+s,1,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                               grad_data_e[INDEX3(l,2,q,numComps,DIM)]+=data_array[l]*
                                    jac->DSDX[INDEX4(s_offset+s,2,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                           }
                       }
                    }
                }
                #undef DIM
            }

         } else if (data_type==FINLEY_DEGREES_OF_FREEDOM) {
            if (jac->numDim==1) {
                #define DIM 1
                #pragma omp for schedule(static)
   	        for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleData(grad_data,e);
                    for (q=0;q<DIM*(jac->ReferenceElement->numQuadNodes)*numComps; q++) grad_data_e[q]=0;
                    for (s=0;s<jac->ReferenceElement->Type->numShapes;s++) {
                       n=elements->Nodes[INDEX2(dof_offset+s,e,NN)];
                       data_array=getSampleData(data,nodes->degreeOfFreedom[n]);
                       for (q=0;q<jac->ReferenceElement->numQuadNodes;q++) {
                           for (l=0;l<numComps;l++) {
                               grad_data_e[INDEX3(l,0,q,numComps,DIM)]+=data_array[l]* 
                                        jac->DSDX[INDEX4(s_offset+s,0,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                           }
                       }
                    }
                }

                #undef DIM
            } else if (jac->numDim==2) {
                #define DIM 2
                #pragma omp for schedule(static)
   	        for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleData(grad_data,e);
                    for (q=0;q<DIM*(jac->ReferenceElement->numQuadNodes)*numComps; q++) grad_data_e[q]=0;
                    for (s=0;s<jac->ReferenceElement->Type->numShapes;s++) {
                       n=elements->Nodes[INDEX2(dof_offset+s,e,NN)];
                       data_array=getSampleData(data,nodes->degreeOfFreedom[n]);
                       for (q=0;q<jac->ReferenceElement->numQuadNodes;q++) {
                           for (l=0;l<numComps;l++) {
                                   grad_data_e[INDEX3(l,0,q,numComps,DIM)]+=data_array[l]* 
                                        jac->DSDX[INDEX4(s_offset+s,0,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                                   grad_data_e[INDEX3(l,1,q,numComps,DIM)]+=data_array[l]* 
                                        jac->DSDX[INDEX4(s_offset+s,1,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                           }
                       }
                    }
                }
                #undef DIM
            } else if (jac->numDim==3) {
                #define DIM 3
                #pragma omp for schedule(static)
   	        for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleData(grad_data,e);
                    for (q=0;q<DIM*(jac->ReferenceElement->numQuadNodes)*numComps; q++) grad_data_e[q]=0;
                    for (s=0;s<jac->ReferenceElement->Type->numShapes;s++) {
                       n=elements->Nodes[INDEX2(dof_offset+s,e,NN)];
                       data_array=getSampleData(data,nodes->degreeOfFreedom[n]);
                       for (q=0;q<jac->ReferenceElement->numQuadNodes;q++) {
                           for (l=0;l<numComps;l++) {
                               grad_data_e[INDEX3(l,0,q,numComps,DIM)]+=data_array[l]* 
                                    jac->DSDX[INDEX4(s_offset+s,0,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                               grad_data_e[INDEX3(l,1,q,numComps,DIM)]+=data_array[l]* 
                                    jac->DSDX[INDEX4(s_offset+s,1,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                               grad_data_e[INDEX3(l,2,q,numComps,DIM)]+=data_array[l]* 
                                    jac->DSDX[INDEX4(s_offset+s,2,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                           }
                       }
                    }
                }
                #undef DIM
            }
         } else if (data_type==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            if (jac->numDim==1) {
                #define DIM 1
                #pragma omp for schedule(static)
   	        for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleData(grad_data,e);
                    for (q=0;q<DIM*(jac->ReferenceElement->numQuadNodes)*numComps; q++) grad_data_e[q]=0;
                    for (s=0;s<jac->ReferenceElement->Type->numShapes;s++) {
                       n=elements->Nodes[INDEX2(elements->ReferenceElement->Type->linearNodes[dof_offset+s],e,NN)];
                       data_array=getSampleData(data,nodes->reducedDegreeOfFreedom[n]);
                       for (q=0;q<jac->ReferenceElement->numQuadNodes;q++) {
                           for (l=0;l<numComps;l++) {
                               grad_data_e[INDEX3(l,0,q,numComps,DIM)]+=data_array[l]* 
                                        jac->DSDX[INDEX4(s_offset+s,0,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                           }
                       }
                    }
                }

                #undef DIM
            } else if (jac->numDim==2) {
                #define DIM 2
                #pragma omp for schedule(static)
   	        for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleData(grad_data,e);
                    for (q=0;q<DIM*(jac->ReferenceElement->numQuadNodes)*numComps; q++) grad_data_e[q]=0;
                    for (s=0;s<jac->ReferenceElement->Type->numShapes;s++) {
                       n=elements->Nodes[INDEX2(elements->ReferenceElement->Type->linearNodes[dof_offset+s],e,NN)];
                       data_array=getSampleData(data,nodes->reducedDegreeOfFreedom[n]);
                       for (q=0;q<jac->ReferenceElement->numQuadNodes;q++) {
                           for (l=0;l<numComps;l++) {
                               grad_data_e[INDEX3(l,0,q,numComps,DIM)]+=data_array[l]* 
                                     jac->DSDX[INDEX4(s_offset+s,0,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                               grad_data_e[INDEX3(l,1,q,numComps,DIM)]+=data_array[l]* 
                                     jac->DSDX[INDEX4(s_offset+s,1,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                           }
                       }
                    }
                }

                #undef DIM

            } else if (jac->numDim==3) {
                #define DIM 3
                #pragma omp for schedule(static)
   	        for (e=0;e<elements->numElements;e++) {
                    grad_data_e=getSampleData(grad_data,e);
                    for (q=0;q<DIM*(jac->ReferenceElement->numQuadNodes)*numComps; q++) grad_data_e[q]=0;
                    for (s=0;s<jac->ReferenceElement->Type->numShapes;s++) {
                       n=elements->Nodes[INDEX2(elements->ReferenceElement->Type->linearNodes[dof_offset+s],e,NN)];
                       data_array=getSampleData(data,nodes->reducedDegreeOfFreedom[n]);
                       for (q=0;q<jac->ReferenceElement->numQuadNodes;q++) {
                           for (l=0;l<numComps;l++) {
                               grad_data_e[INDEX3(l,0,q,numComps,DIM)]+=data_array[l]* 
                                        jac->DSDX[INDEX4(s_offset+s,0,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                               grad_data_e[INDEX3(l,1,q,numComps,DIM)]+=data_array[l]* 
                                        jac->DSDX[INDEX4(s_offset+s,1,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
                               grad_data_e[INDEX3(l,2,q,numComps,DIM)]+=data_array[l]* 
                                        jac->DSDX[INDEX4(s_offset+s,2,q,e,jac->ReferenceElement->Type->numNodes,DIM,jac->ReferenceElement->numQuadNodes)];
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
/*
 * $Log$
 * Revision 1.6  2005/09/15 03:44:21  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.5.2.1  2005/09/07 06:26:17  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.5  2005/07/08 04:07:47  jgs
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
