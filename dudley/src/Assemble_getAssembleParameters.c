
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

/*    assemblage routines: prepares the assemble parameter set */

/**************************************************************/

#include "Assemble.h"

/**************************************************************/

void Assemble_getAssembleParameters(Dudley_NodeFile* nodes,Dudley_ElementFile* elements,Paso_SystemMatrix* S, 
                                        escriptDataC* F, bool_t reducedIntegrationOrder, Assemble_Parameters *parm) {
  Dudley_resetError();
  
  if (!isEmpty(F) && !isExpanded(F) ) {
      Dudley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: Right hand side is not expanded.");
      return;
  }
  /*  check the dimensions of S and F */
  if (S!=NULL && !isEmpty(F)) {
    if (! numSamplesEqual(F,1,(Paso_Distribution_getMyNumComponents(S->row_distribution)*S->row_block_size)/S->logical_row_block_size)) {
      Dudley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: number of rows of matrix and length of right hand side don't match.");
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
          Dudley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: matrix row block size and number of components of right hand side don't match.");
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
           parm->row_jac=Dudley_ElementFile_borrowJacobeans(elements,nodes,FALSE,reducedIntegrationOrder);
      } 
      else if ( Paso_Distribution_getMyNumComponents(S->row_distribution)*S->row_block_size==parm->numEqu* Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution)) {
           parm->row_DOF_UpperBound =  Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
           parm->row_DOF=nodes->reducedDegreesOfFreedomMapping->target;
	   parm->row_jac=Dudley_ElementFile_borrowJacobeans(elements,nodes,TRUE,reducedIntegrationOrder);
      } else {
           Dudley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: number of rows in matrix does not match the number of degrees of freedom in mesh");
      }
      /* Make sure # cols in matrix == num DOF for one of: full or reduced (use numLocalDOF for MPI) */
      if ( Paso_Distribution_getMyNumComponents(S->col_distribution)*S->col_block_size==parm->numComp* Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution)) {
           parm->col_DOF_UpperBound =  Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution);
		   parm->col_DOF=nodes->degreesOfFreedomMapping->target;
		   parm->col_jac=Dudley_ElementFile_borrowJacobeans(elements,nodes,FALSE,reducedIntegrationOrder);
      } else if ( Paso_Distribution_getMyNumComponents(S->col_distribution)*S->col_block_size==parm->numComp* Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution)) {
           parm->col_DOF_UpperBound =  Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
           parm->col_DOF=nodes->reducedDegreesOfFreedomMapping->target;
           parm->col_jac=Dudley_ElementFile_borrowJacobeans(elements,nodes,TRUE,reducedIntegrationOrder);
      } else {
           Dudley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: number of columns in matrix does not match the number of degrees of freedom in mesh");
      }
  }
  if (! Dudley_noError()) return;
  /* get the information from right hand side */
  if (!isEmpty(F)) {
      if (numSamplesEqual(F,1, Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution))) {
           parm->row_DOF_UpperBound =  Paso_Distribution_getMyNumComponents(nodes->degreesOfFreedomDistribution);
           parm->row_DOF=nodes->degreesOfFreedomMapping->target;
           parm->row_jac=Dudley_ElementFile_borrowJacobeans(elements,nodes,FALSE,reducedIntegrationOrder);
      } 
      else if (numSamplesEqual(F,1, Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution))) {
           parm->row_DOF_UpperBound =  Paso_Distribution_getMyNumComponents(nodes->reducedDegreesOfFreedomDistribution);
           parm->row_DOF=nodes->reducedDegreesOfFreedomMapping->target;
	   parm->row_jac=Dudley_ElementFile_borrowJacobeans(elements,nodes,TRUE,reducedIntegrationOrder);
      } else {
           Dudley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: length of RHS vector does not match the number of degrees of freedom in mesh");
      }
      if (S==NULL) {
           parm->col_DOF_UpperBound=parm->row_DOF_UpperBound;
           parm->col_DOF=parm->row_DOF;
           parm->col_jac=parm->row_jac;
      } 
  }

  if ( parm->row_jac->numDim !=parm->col_jac->numDim) {
      Dudley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: spacial dimension for row and column shape function must match.");
  }
  
  if (elements->numNodes < parm->row_jac->numShapesTotal) {
      Dudley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: too many nodes are expected by row.");  
  }
  if (elements->numNodes < parm->row_jac->numShapesTotal) {
      Dudley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: too many nodes are expected by col.");  
  }
  if ( parm->row_jac->numElements !=elements->numElements) {
      Dudley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: number of elements for row is wrong.");
  }
  if ( parm->col_jac->numElements !=elements->numElements) {
      Dudley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: number of elements for column is wrong.");
  }
  if ( parm->row_jac->numQuadTotal !=parm->col_jac->numQuadTotal) {
      Dudley_setError(TYPE_ERROR,"Assemble_getAssembleParameters: number of quadrature points for row and column shape functions must match.");
  }
  
  parm->numQuadTotal=parm->row_jac->numQuadTotal; 
  parm->NN=elements->numNodes;
  parm->numElements=elements->numElements;
  parm->numDim=parm->row_jac->numDim;
  parm->col_node=parm->col_jac->node_selection;
  parm->row_node=parm->row_jac->node_selection;
  parm->row_numShapesTotal=parm->row_jac->numShapesTotal;
  parm->row_numShapes=parm->row_jac->BasisFunctions->Type->numShapes;
  parm->col_numShapesTotal=parm->col_jac->numShapesTotal;
  parm->col_numShapes=parm->col_jac->BasisFunctions->Type->numShapes;
	

}
