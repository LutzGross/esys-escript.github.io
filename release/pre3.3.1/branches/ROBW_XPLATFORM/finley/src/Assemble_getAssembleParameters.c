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
                                        escriptDataC* F,Assemble_Parameters *parm) {
  Finley_resetError();
  dim_t i;
  parm->NN=elements->ReferenceElement->Type->numNodes;
  for (i=0;i<parm->NN;i++) parm->id[i]=i;
  parm->NS=elements->ReferenceElement->Type->numShapes;
  parm->referenceElement=elements->ReferenceElement;
  parm->numQuad=parm->referenceElement->numQuadNodes;
  parm->numDim=nodes->numDim;
  parm->numElementDim=parm->referenceElement->Type->numDim;

  if (!isEmpty(F) && !isExpanded(F) ) {
      Finley_setError(TYPE_ERROR,"__FILE__: Right hand side is not expanded.");
      return;
  }
  /*  check the dimensions of S and F */
  if (S!=NULL && !isEmpty(F)) {
    if ( getDataPointSize(F)!=S->logical_row_block_size) {
      Finley_setError(TYPE_ERROR,"__FILE__: matrix row block size and number of components of right hand side don't match.");
      return;
    }

    if (! numSamplesEqual(F,1,(S->num_rows*S->row_block_size)/S->logical_row_block_size)) {
      Finley_setError(TYPE_ERROR,"__FILE__: number of rows of matrix and length of right hand side don't match.");
      return;
    }
  }
  /* get the number of equations and components */
  if (S!=NULL) {
    parm->numEqu=S->logical_row_block_size;
    parm->numComp=S->logical_col_block_size;
  } else {
    parm->numEqu=1;
    parm->numComp=1;
  }
  if (!isEmpty(F)) parm->numEqu=getDataPointSize(F);
  
  parm->label_col=nodes->degreeOfFreedom;
  parm->label_row=nodes->degreeOfFreedom;
  parm->referenceElement_row=elements->ReferenceElement;
  parm->referenceElement_col=elements->ReferenceElement;
  /* get the information for the labeling of the degrees of freedom */
  if (S!=NULL) {
      if (S->num_rows*S->row_block_size==parm->numEqu*nodes->numDegreesOfFreedom) {
           parm->label_row=nodes->degreeOfFreedom;
           parm->row_node=&(parm->id[0]);
           parm->referenceElement_row=elements->ReferenceElement;
      } else if (S->num_rows*S->row_block_size==parm->numEqu*nodes->reducedNumDegreesOfFreedom) {
           parm->label_row=nodes->reducedDegreeOfFreedom;
           parm->row_node=parm->referenceElement->Type->linearNodes;
           parm->referenceElement_row=elements->LinearReferenceElement;
      } else {
           Finley_setError(TYPE_ERROR,"__FILE__: number of rows in matrix does not match the number of degrees of freedom in mesh");
           return;
      }
      if (S->num_cols*S->col_block_size==parm->numComp*nodes->numDegreesOfFreedom) {
           parm->label_col=nodes->degreeOfFreedom;
           parm->col_node=&(parm->id[0]);
           parm->referenceElement_col=elements->ReferenceElement;
      } else if (S->num_cols*S->col_block_size==parm->numComp*nodes->reducedNumDegreesOfFreedom) {
           parm->label_col=nodes->reducedDegreeOfFreedom;
           parm->col_node=parm->referenceElement->Type->linearNodes;
           parm->referenceElement_col=elements->LinearReferenceElement;
      } else {
           Finley_setError(TYPE_ERROR,"__FILE__: number of columns in matrix does not match the number of degrees of freedom in mesh");
           return;
      }
  }
  if (!isEmpty(F)) {
      if (numSamplesEqual(F,1,nodes->numDegreesOfFreedom)) {
           parm->row_node=&(parm->id[0]);
           parm->label_row=nodes->degreeOfFreedom;
           parm->referenceElement_row=elements->ReferenceElement;
      } else if (numSamplesEqual(F,1,nodes->reducedNumDegreesOfFreedom)) {
           parm->label_row=nodes->reducedDegreeOfFreedom;
           parm->row_node=parm->referenceElement->Type->linearNodes;
           parm->referenceElement_row=elements->LinearReferenceElement;
      } else {
           Finley_setError(TYPE_ERROR,"__FILE__: length of RHS vector does not match the number of degrees of freedom in mesh");
           return;
      }
      if (S==NULL) {
           parm->label_col=parm->label_row;
           parm->col_node=parm->row_node;
           parm->referenceElement_col=parm->referenceElement_row;
      } 
  }
  if (parm->referenceElement_row!=parm->referenceElement_col) {
      Finley_setError(TYPE_ERROR,"__FILE__: assemblage cannot handel different shape functions for rows and columns (yet).");
      return;
  }
  parm->NN_row=parm->referenceElement_row->Type->numNodes;
  parm->NS_row=parm->referenceElement_row->Type->numShapes;
  parm->NN_col=parm->referenceElement_col->Type->numNodes;
  parm->NS_col=parm->referenceElement_col->Type->numShapes;
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

