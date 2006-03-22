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

/*    assemblage routines: calculate the gradient of nodal data at quadrature points */

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


#define NODES 0
#define DOF 1
#define REDUCED_DOF 2

void Finley_Assemble_gradient(Finley_NodeFile* nodes, Finley_ElementFile* elements,
             escriptDataC* grad_data,escriptDataC* data) {

  double *local_X=NULL, *local_data=NULL, *dVdv=NULL, *dvdV=NULL, *Vol=NULL, *d_datadv=NULL, *gradS=NULL,*data_array;
  index_t node_offset,*resort_nodes=FALSE,dof_offset;
  dim_t numNodes=0,e,i,q,NS_DOF=0,NN_DOF=0;
  type_t type=DOF;
  if (nodes==NULL || elements==NULL) return;
  dim_t NN=elements->ReferenceElement->Type->numNodes;
  dim_t NS=elements->ReferenceElement->Type->numShapes;
  index_t id[NN];
  dim_t numDim=nodes->numDim;
  type_t data_type=getFunctionSpaceType(data);
  dim_t numComps=getDataPointSize(data);
  dim_t numQuad=elements->ReferenceElement->numQuadNodes;
  for (i=0;i<NN;i++) id[i]=i;
  Finley_resetError();

    /* set some parameter */

  if (data_type==FINLEY_NODES) {
       type=NODES;
       resort_nodes=id;
       NN_DOF=elements->ReferenceElement->Type->numNodes;
       NS_DOF=elements->ReferenceElement->Type->numShapes;
       gradS=elements->ReferenceElement->dSdv;
       numNodes=nodes->numNodes;
  } else if (data_type==FINLEY_DEGREES_OF_FREEDOM) {
       type=DOF;
       resort_nodes=id;
       NN_DOF=elements->ReferenceElement->Type->numNodes;
       NS_DOF=elements->ReferenceElement->Type->numShapes;
       gradS=elements->ReferenceElement->dSdv;
       numNodes=nodes->numDegreesOfFreedom;
  } else if (data_type==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
       type=REDUCED_DOF;
       resort_nodes=elements->ReferenceElement->Type->linearNodes;
       NN_DOF=elements->LinearReferenceElement->Type->numNodes;
       NS_DOF=elements->LinearReferenceElement->Type->numShapes;
       gradS=elements->LinearReferenceElement->dSdv;
       numNodes=nodes->reducedNumDegreesOfFreedom;
   } else {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: Cannot calculate gradient of data");
  }
  if (getFunctionSpaceType(grad_data)==FINLEY_CONTACT_ELEMENTS_2) {
       node_offset=NN-NS;
       dof_offset=NN_DOF-NS_DOF;
  } else {
       node_offset=0;
       dof_offset=0;
  }
                                                                                                                                                         
  /* check the dimensions of interpolated_data and data */

  if (numDim!=elements->ReferenceElement->Type->numDim) {
     Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: Spatial and element dimension must match.");
  } else if (! numSamplesEqual(grad_data,numQuad,elements->numElements)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: illegal number of samples in gradient Data object");
  } else if (! numSamplesEqual(data,1,numNodes)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: illegal number of samples of input Data object");
  } else if (numDim*numComps!=getDataPointSize(grad_data)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: illegal number of components in gradient data object.");
  }  else if (!isExpanded(grad_data)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_gradient: expanded Data object is expected for output data.");
  }
  
  /* now we can start */

  if (Finley_noError()) {
          #pragma omp parallel private(local_X,local_data,dvdV,dVdv,Vol,d_datadv)
          {
               local_X=local_data=dVdv=dvdV=Vol=d_datadv=NULL;
   	      /* allocation of work arrays */
   	      local_X=THREAD_MEMALLOC(NS*numDim,double);
   	      local_data=THREAD_MEMALLOC(NS*numComps,double);
   	      dVdv=THREAD_MEMALLOC(numQuad*numDim*numDim,double);
   	      dvdV=THREAD_MEMALLOC(numQuad*numDim*numDim,double);
   	      Vol=THREAD_MEMALLOC(numQuad,double);
   	      d_datadv=THREAD_MEMALLOC(numQuad*numComps*numDim,double);
   	      if (!(Finley_checkPtr(local_X) || Finley_checkPtr(dVdv) || Finley_checkPtr(dvdV) || Finley_checkPtr(Vol) || Finley_checkPtr(d_datadv) || Finley_checkPtr(local_data) ))  {
   	        /* open the element loop */
                #pragma omp for private(e,i,q,data_array) schedule(static)
   	        for(e=0;e<elements->numElements;e++) {
   	          /* gather local coordinates of nodes into local_X: */
   	          Finley_Util_Gather_double(NS,&(elements->Nodes[INDEX2(node_offset,e,NN)]),numDim,nodes->Coordinates,local_X);
   	          /*  calculate dVdv(i,j,q)=local_X(i,n)*DSDv(n,j,q) */
   	          Finley_Util_SmallMatMult(numDim,numDim*numQuad,dVdv,NS,local_X,elements->ReferenceElement->dSdv);
   	          /*  dvdV=invert(dVdv) */
   	          Finley_Util_InvertSmallMat(numQuad,numDim,dVdv,dvdV,Vol);
   	          /* gather local data into local_data(numComps,NS_DOF): */
                  switch (type) {
                     case NODES:
                        for (q=0;q<NS_DOF;q++) {
                           i=elements->Nodes[INDEX2(resort_nodes[dof_offset+q],e,NN)];
                           data_array=getSampleData(data,i);
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
   	          /*  calculate grad_data(l,i,q)=local_data(l,n)* DSDV(n,i,q) */
   	          // Finley_Util_SmallMatMult(numQuad,numComps,numQuad*numDim,getSampleData(grad_data,e),NS_DOF,local_data,dSdV);

   	          /*  calculate d_datadv(l,i,q)=local_data(l,n)*DSDv(n,i,q) */
   	          Finley_Util_SmallMatMult(numComps,numDim*numQuad,d_datadv,NS_DOF,local_data,gradS);
   	          /*  calculate grad_data(l,i,q)=d_datadv(l,k,q)*dvdV(k,i,q) */
   	          Finley_Util_SmallMatSetMult(numQuad,numComps,numDim,getSampleData(grad_data,e),numDim,d_datadv,dvdV);
   	        } /* for */
   	      }
   	      THREAD_MEMFREE(local_X);
   	      THREAD_MEMFREE(dVdv);
   	      THREAD_MEMFREE(dvdV);
   	      THREAD_MEMFREE(Vol);
   	      THREAD_MEMFREE(local_data);
   	      THREAD_MEMFREE(d_datadv);
      }
  }
}
#undef NODES
#undef DOF
#undef REDUCED_DOF
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
