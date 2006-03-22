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

/* Finley: Converting an element list into a matrix shape     */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "IndexList.h"

/**************************************************************/
/* inserts the contributions from the element matrices of elements
   into the row index col. If symmetric is set, only the upper
   triangle of the matrix is stored. */

void Finley_IndexList_insertElements(Finley_IndexList* index_list, Finley_ElementFile* elements,
                                       bool_t reduce_row_order, index_t* row_Label,
                                       bool_t reduce_col_order, index_t* col_Label) {
  index_t color;
  dim_t e,kr,kc,NN_row,NN_col,i,icol,irow;

  if (elements!=NULL) {
    dim_t NN=elements->ReferenceElement->Type->numNodes;
    index_t id[NN],*row_node,*col_node;
    for (i=0;i<NN;i++) id[i]=i;
    if (reduce_col_order) {
       col_node=elements->ReferenceElement->Type->linearNodes;
       NN_col=elements->LinearReferenceElement->Type->numNodes;
    } else {
       col_node=id;
       NN_col=elements->ReferenceElement->Type->numNodes;
    }
    if (reduce_row_order) {
       row_node=elements->ReferenceElement->Type->linearNodes;
       NN_row=elements->LinearReferenceElement->Type->numNodes;
    } else {
       row_node=id;
       NN_row=elements->ReferenceElement->Type->numNodes;
    }
    for (color=elements->minColor;color<=elements->maxColor;color++) {
        #pragma omp for private(e,irow,kr,kc,icol) schedule(static)
        for (e=0;e<elements->numElements;e++) {
            if (elements->Color[e]==color) {
                for (kr=0;kr<NN_row;kr++) {
                  irow=row_Label[elements->Nodes[INDEX2(row_node[kr],e,NN)]];
                  for (kc=0;kc<NN_col;kc++) {
                       icol=col_Label[elements->Nodes[INDEX2(col_node[kc],e,NN)]];
                       Finley_IndexList_insertIndex(&(index_list[irow]),icol);
                  }
                }
            }
        }
      }
  }
  return;
}

/* inserts row index row into the Finley_IndexList in if it does not exist */

void Finley_IndexList_insertIndex(Finley_IndexList* in, index_t index) {
  dim_t i;
  /* is index in in? */
  for (i=0;i<in->n;i++) {
    if (in->index[i]==index)  return;
  }
  /* index could not be found */
  if (in->n==INDEXLIST_LENGTH) {
     /* if in->index is full check the extension */
     if (in->extension==NULL) {
        in->extension=TMPMEMALLOC(1,Finley_IndexList);
        if (Finley_checkPtr(in->extension)) return;
        in->extension->n=0;
        in->extension->extension=NULL;
     }
     Finley_IndexList_insertIndex(in->extension,index);
  } else {
     /* insert index into in->index*/
     in->index[in->n]=index;
     in->n++;
  }
}

/* counts the number of row indices in the Finley_IndexList in */

dim_t Finley_IndexList_count(Finley_IndexList* in) {
  if (in==NULL) {
     return 0;
  } else {
     return (in->n)+Finley_IndexList_count(in->extension);
  }
}

/* count the number of row indices in the Finley_IndexList in */

void Finley_IndexList_toArray(Finley_IndexList* in, index_t* array) {
  dim_t i;
  if (in!=NULL) {
    for (i=0;i<in->n;i++) array[i]=in->index[i];
    Finley_IndexList_toArray(in->extension,&(array[in->n]));
  }
}

/* deallocates the Finley_IndexList in by recursive calls */

void Finley_IndexList_free(Finley_IndexList* in) {
  if (in!=NULL) {
    Finley_IndexList_free(in->extension);
    TMPMEMFREE(in);
  }
}

/*
 * $Log$
 * Revision 1.6  2005/09/15 03:44:22  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.5.2.1  2005/09/07 06:26:18  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.5  2005/07/08 04:07:51  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.4  2004/12/15 07:08:32  jgs
 * *** empty log message ***
 * Revision 1.1.1.1.2.3  2005/06/29 02:34:50  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1.2.2  2004/11/24 01:37:13  gross
 * some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
 *
 *
 *
 */
