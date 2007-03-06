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

/*    assemblage routines: prepares the assemble parameter set */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004,2005 */
/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Assemble.h"

/**************************************************************/

void Assemble_getAssembleParameters(Finley_NodeFile* nodes,Finley_ElementFile* elements,Paso_SystemMatrix* S, 
                                        escriptDataC* F, bool_t reducedIntegrationOrder, Assemble_Parameters *parm) {
  Finley_resetError();

  dim_t i;
  for (i=0;i<MAX_numNodes;i++) parm->id[i]=i;

  if (!isEmpty(F) && !isExpanded(F) ) {
      Finley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: Right hand side is not expanded.");
      return;
  }
  /*  check the dimensions of S and F */
  if (S!=NULL && !isEmpty(F)) {
    printf("ksteube S->myNumRows=%d S->row_block_size=%d S->logical_row_block_size=%d\n", S->myNumRows, S->row_block_size, S->logical_row_block_size);
    if (! numSamplesEqual(F,1,(S->myNumRows*S->row_block_size)/S->logical_row_block_size)) {
      Finley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: number of rows of matrix and length of right hand side don't match.");
      return;
    }
  }
  /* get the number of equations and components */
  if (S==NULL) {
     if (isEmpty(F)) {
        parm->numEqu=1;
        parm->numComp=1;
     } else {
        parm->numEqu=getDataPointSize(F);
        parm->numComp=parm->numEqu;
     }
  } else {
     if (isEmpty(F)) {
        parm->numEqu=S->logical_row_block_size;
        parm->numComp=S->logical_col_block_size;
     } else {
        if ( getDataPointSize(F)!=S->logical_row_block_size) {
          Finley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: matrix row block size and number of components of right hand side don't match.");
          return;
        }
        parm->numEqu=S->logical_row_block_size;
        parm->numComp=S->logical_col_block_size;
     }
  }

  parm->col_DOF=nodes->degreeOfFreedom;
  parm->row_DOF=nodes->degreeOfFreedom;
  /* get the information for the labeling of the degrees of freedom from matrix */
  /* printf("ksteube in Assemble_getAssembleParameters cpu=%d NR=%d BS=%d DOF=%d\n", elements->elementDistribution->MPIInfo->rank, S->myNumRows, S->row_block_size, parm->numEqu*nodes->degreeOfFreedomDistribution->numLocal); */
  if (S!=NULL) {
      /* Make sure # rows in matrix == num DOF for one of: full or reduced (use numLocalDOF for MPI) */
#ifndef PASO_MPI
      if (S->myNumRows*S->row_block_size==parm->numEqu*nodes->numDegreesOfFreedom) {
           parm->row_DOF_UpperBound = nodes->numDegreesOfFreedom;
#else
      if (S->myNumRows*S->row_block_size==parm->numEqu*nodes->degreeOfFreedomDistribution->numLocal) {
           parm->row_DOF_UpperBound = nodes->degreeOfFreedomDistribution->numLocal;
#endif
           parm->row_DOF=nodes->degreeOfFreedom;
           parm->row_node=&(parm->id[0]);
           parm->row_jac=Finley_ElementFile_borrowJacobeans(elements,nodes,FALSE,reducedIntegrationOrder);
      } 
#ifndef PASO_MPI
      else if (S->myNumRows*S->row_block_size==parm->numEqu*nodes->reducedNumDegreesOfFreedom) {
           parm->row_DOF_UpperBound = nodes->reducedNumDegreesOfFreedom;
#else
      else if (S->myNumRows*S->row_block_size==parm->numEqu*nodes->reducedDegreeOfFreedomDistribution->numLocal) {
           parm->row_DOF_UpperBound = nodes->reducedDegreeOfFreedomDistribution->numLocal;
#endif
           parm->row_DOF=nodes->reducedDegreeOfFreedom;
           parm->row_jac=Finley_ElementFile_borrowJacobeans(elements,nodes,TRUE,reducedIntegrationOrder);
           parm->row_node=parm->row_jac->ReferenceElement->Type->linearNodes;
      } else {
           Finley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: number of rows in matrix does not match the number of degrees of freedom in mesh");
      }
      /* Make sure # cols in matrix == num DOF for one of: full or reduced (use numLocalDOF for MPI) */
#ifndef PASO_MPI      
      if (S->myNumCols*S->col_block_size==parm->numComp*nodes->numDegreesOfFreedom) {
           parm->col_DOF_UpperBound = nodes->numDegreesOfFreedom;
#else
      if (S->myNumCols*S->col_block_size==parm->numComp*nodes->degreeOfFreedomDistribution->numLocal) {
           parm->col_DOF_UpperBound = nodes->degreeOfFreedomDistribution->numLocal;
#endif
           parm->col_jac=Finley_ElementFile_borrowJacobeans(elements,nodes,FALSE,reducedIntegrationOrder);
           parm->col_DOF=nodes->degreeOfFreedom;
           parm->col_node=&(parm->id[0]);
      } 
#ifndef PASO_MPI
      else if (S->myNumCols*S->col_block_size==parm->numComp*nodes->reducedNumDegreesOfFreedom) {
           parm->col_DOF_UpperBound = nodes->reducedNumDegreesOfFreedom;
#else
      else if (S->myNumCols*S->col_block_size==parm->numComp*nodes->reducedDegreeOfFreedomDistribution->numLocal) {
           parm->col_DOF_UpperBound = nodes->reducedDegreeOfFreedomDistribution->numLocal;
#endif
           parm->col_DOF=nodes->reducedDegreeOfFreedom;
           parm->col_jac=Finley_ElementFile_borrowJacobeans(elements,nodes,TRUE,reducedIntegrationOrder);
           parm->col_node=parm->row_jac->ReferenceElement->Type->linearNodes;
      } else {
           Finley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: number of columns in matrix does not match the number of degrees of freedom in mesh");
      }
  }
  if (! Finley_noError()) return;
  /* get the information from right hand side */
  if (!isEmpty(F)) {
#ifndef PASO_MPI
      if (numSamplesEqual(F,1,nodes->numDegreesOfFreedom)) {
           parm->row_DOF_UpperBound=nodes->numDegreesOfFreedom;
#else
      if (numSamplesEqual(F,1,nodes->degreeOfFreedomDistribution->numLocal)) {
           parm->row_DOF_UpperBound = nodes->degreeOfFreedomDistribution->numLocal;
#endif
           parm->row_DOF=nodes->degreeOfFreedom;
           parm->row_jac=Finley_ElementFile_borrowJacobeans(elements,nodes,FALSE,reducedIntegrationOrder);
           parm->row_node=&(parm->id[0]);
      } 
#ifndef PASO_MPI
      else if (numSamplesEqual(F,1,nodes->reducedNumDegreesOfFreedom)) {
           parm->row_DOF_UpperBound=nodes->reducedNumDegreesOfFreedom;
#else
      else if (numSamplesEqual(F,1,nodes->reducedDegreeOfFreedomDistribution->numLocal)) {
           parm->row_DOF_UpperBound = nodes->reducedDegreeOfFreedomDistribution->numLocal;
#endif
           parm->row_DOF=nodes->reducedDegreeOfFreedom;
           parm->row_jac=Finley_ElementFile_borrowJacobeans(elements,nodes,TRUE,reducedIntegrationOrder);
           parm->row_node=parm->row_jac->ReferenceElement->Type->linearNodes;
      } else {
           Finley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: length of RHS vector does not match the number of degrees of freedom in mesh");
      }
      if (S==NULL) {
           parm->col_DOF_UpperBound=parm->row_DOF_UpperBound;
           parm->col_DOF=parm->row_DOF;
           parm->col_node=parm->row_node;
           parm->col_jac=parm->row_jac;
      } 
  }
  if (! Finley_noError()) return;
  if (parm->row_jac!=parm->col_jac) {
      Finley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: assemblage cannot handle different shape functions for rows and columns (yet).");
  }
  if (! Finley_noError()) return;
  parm->NN=elements->ReferenceElement->Type->numNodes;
  parm->numDim=parm->row_jac->numDim;
  parm->numQuad=parm->row_jac->ReferenceElement->numQuadNodes;
  parm->row_NN=parm->row_jac->ReferenceElement->Type->numNodes;
  parm->row_NS=parm->row_jac->ReferenceElement->Type->numShapes;
  parm->col_NN=parm->col_jac->ReferenceElement->Type->numNodes;
  parm->col_NS=parm->col_jac->ReferenceElement->Type->numShapes;
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
 * Revision 1.1.1.1.2.4  2005/07/01 07:02:13  gross
 * some bug with OPENMP fixed
 *
 * Revision 1.1.1.1.2.3  2005/06/29 02:34:48  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1.2.2  2004/11/12 06:58:18  gross
 * a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
 *
 * Revision 1.1.1.1.2.1  2004/10/28 22:59:24  gross
 * finley's RecTest.py is running now: problem in SystemMatrixAdapater fixed
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.2  2004/07/21 05:00:54  gross
 * name changes in DataC
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */

