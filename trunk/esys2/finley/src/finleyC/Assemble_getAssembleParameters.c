/* $Id$ */

/**************************************************************/

/*    assemblage routines: prepares the assemble parameter set */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004 */
/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "escript/Data/DataC.h"
#include "Finley.h"
#include "Assemble.h"
#include "NodeFile.h"
#include "ElementFile.h"

/**************************************************************/

void Assemble_getAssembleParameters(Finley_NodeFile* nodes,Finley_ElementFile* elements,Finley_SystemMatrix* S, 
                                        escriptDataC* F,Assemble_Parameters *parm) {
  int i;
  parm->NN=elements->ReferenceElement->Type->numNodes;
  for (i=0;i<parm->NN;i++) parm->id[i]=i;
  parm->NS=elements->ReferenceElement->Type->numShapes;
  parm->referenceElement=elements->ReferenceElement;
  parm->numQuad=parm->referenceElement->numQuadNodes;
  parm->numDim=nodes->numDim;
  parm->numElementDim=parm->referenceElement->Type->numDim;

  if (!isExpanded(F) ) {
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"Right hand side is not expanded.");
      return;
  }
  /*  check the dimensions of S and F */
  if (S!=NULL && !isEmpty(F)) {
    if ( getDataPointSize(F)!=S->total_row_block_size) {
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"matrix rows and number of components of right hand side don't match.");
      return;
    }
    if (! numSamplesEqual(F,1,S->num_rows*S->row_block_size)) {
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"number of rows of matrix and length of right hand side don't match.");
      return;
    }
  }
  /* get the number of equations and components */
  if (S!=NULL) {
    parm->numEqu=S->total_row_block_size;
    parm->numComp=S->total_col_block_size;
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
           Finley_ErrorCode=TYPE_ERROR;
           sprintf(Finley_ErrorMsg,"number of rows in matrix does not match the number of degrees of freedom in mesh");
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
           Finley_ErrorCode=TYPE_ERROR;
           sprintf(Finley_ErrorMsg,"number of columns in matrix does not match the number of degrees of freedom in mesh");
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
           Finley_ErrorCode=TYPE_ERROR;
           sprintf(Finley_ErrorMsg,"length of RHS vector does not match the number of degrees of freedom in mesh");
           return;
      }
      if (S==NULL) {
           parm->label_col=parm->label_row;
           parm->col_node=parm->row_node;
           parm->referenceElement_col=parm->referenceElement_row;
      } 
  }
  if (parm->referenceElement_row!=parm->referenceElement_col) {
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"assemblage cannot handel different shape functions for rows and columns (yet).");
      return;
  }
  parm->NN_row=parm->referenceElement_row->Type->numNodes;
  parm->NS_row=parm->referenceElement_row->Type->numShapes;
  parm->NN_col=parm->referenceElement_col->Type->numNodes;
  parm->NS_col=parm->referenceElement_col->Type->numShapes;
}

/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.2  2004/07/21 05:00:54  gross
 * name changes in DataC
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */

