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
  dim_t i;
  Finley_resetError();

  for (i=0;i<MAX_numNodes;i++) parm->id[i]=i;

  if (!isEmpty(F) && !isExpanded(F) ) {
      Finley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: Right hand side is not expanded.");
      return;
  }
  /*  check the dimensions of S and F */
  if (S!=NULL && !isEmpty(F)) {
    if (! numSamplesEqual(F,1,(Paso_Distribution_getMyNumComponents(S->row_distribution)*S->row_block_size)/S->logical_row_block_size)) {
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

  parm->col_DOF=nodes->degreesOfFreedomMapping->target;
  parm->row_DOF=nodes->degreesOfFreedomMapping->target;
  /* get the information for the labeling of the degrees of freedom from matrix */
  if (S!=NULL) {
      /* Make sure # rows in matrix == num DOF for one of: full or reduced (use numLocalDOF for MPI) */
      if ( Paso_Distribution_getMyNumComponents(S->row_distribution)*S->row_block_size==parm->numEqu* Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution)) {
           parm->row_DOF_UpperBound =  Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution);
           parm->row_DOF=nodes->degreesOfFreedomMapping->target;

           parm->row_node=&(parm->id[0]);
           parm->row_jac=Finley_ElementFile_borrowJacobeans(elements,nodes,FALSE,reducedIntegrationOrder);
      } 
      else if ( Paso_Distribution_getMyNumComponents(S->row_distribution)*S->row_block_size==parm->numEqu* Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution)) {
           parm->row_DOF_UpperBound =  Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
           parm->row_DOF=nodes->reducedDegreesOfFreedomMapping->target;
           parm->row_jac=Finley_ElementFile_borrowJacobeans(elements,nodes,TRUE,reducedIntegrationOrder);
           parm->row_node=parm->row_jac->ReferenceElement->Type->linearNodes;
      } else {
           Finley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: number of rows in matrix does not match the number of degrees of freedom in mesh");
      }
      /* Make sure # cols in matrix == num DOF for one of: full or reduced (use numLocalDOF for MPI) */
      if ( Paso_Distribution_getMyNumComponents(S->col_distribution)*S->col_block_size==parm->numComp* Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution)) {
           parm->col_DOF_UpperBound =  Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution);
           parm->col_jac=Finley_ElementFile_borrowJacobeans(elements,nodes,FALSE,reducedIntegrationOrder);
           parm->col_DOF=nodes->degreesOfFreedomMapping->target;
           parm->col_node=&(parm->id[0]);
      } 
      else if ( Paso_Distribution_getMyNumComponents(S->col_distribution)*S->col_block_size==parm->numComp* Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution)) {
           parm->col_DOF_UpperBound =  Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
           parm->col_DOF=nodes->reducedDegreesOfFreedomMapping->target;
           parm->col_jac=Finley_ElementFile_borrowJacobeans(elements,nodes,TRUE,reducedIntegrationOrder);
           parm->col_node=parm->row_jac->ReferenceElement->Type->linearNodes;
      } else {
           Finley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: number of columns in matrix does not match the number of degrees of freedom in mesh");
      }
  }
  if (! Finley_noError()) return;
  /* get the information from right hand side */
  if (!isEmpty(F)) {
      if (numSamplesEqual(F,1, Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution))) {
           parm->row_DOF_UpperBound =  Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution);
           parm->row_DOF=nodes->degreesOfFreedomMapping->target;
           parm->row_jac=Finley_ElementFile_borrowJacobeans(elements,nodes,FALSE,reducedIntegrationOrder);
           parm->row_node=&(parm->id[0]);
      } 
      else if (numSamplesEqual(F,1, Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution))) {
           parm->row_DOF_UpperBound =  Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
           parm->row_DOF=nodes->reducedDegreesOfFreedomMapping->target;
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
