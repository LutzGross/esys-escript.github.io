
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

/* Dudley: Converting an element list into a matrix shape     */

/**************************************************************/

#include "IndexList.h"

/* Translate from distributed/local array indices to global indices */

/**************************************************************/
/* inserts the contributions from the element matrices of elements
   into the row index col. If symmetric is set, only the upper
   triangle of the matrix is stored. */

void Dudley_IndexList_insertElements(Dudley_IndexList* index_list, Dudley_ElementFile* elements,
                                       bool_t reduce_row_order, index_t* row_map,
                                       bool_t reduce_col_order, index_t* col_map)
{
  /* index_list is an array of linked lists. Each entry is a row (DOF) and contains the indices to the non-zero columns */
  index_t color;
  Dudley_ReferenceElement*refElement; 
  dim_t e, kr, kc, NN_row, NN_col, icol, irow, NN, *row_node=NULL, *col_node=NULL;
  if (elements!=NULL)
  {
    NN=elements->numNodes;
    refElement= Dudley_ReferenceElementSet_borrowReferenceElement(elements->referenceElementSet, FALSE);
    if (reduce_col_order)
    {
	col_node=refElement->Type->linearNodes;
	NN_col=(refElement->LinearBasisFunctions->Type->numShapes);
    }
    else
    {
	col_node=refElement->Type->subElementNodes;
	NN_col=(refElement->BasisFunctions->Type->numShapes);
    }
    if (reduce_row_order)
    {
	row_node=refElement->Type->linearNodes;
	NN_row=(refElement->LinearBasisFunctions->Type->numShapes);
    } else {
	row_node=refElement->Type->subElementNodes;
	NN_row=(refElement->BasisFunctions->Type->numShapes) ;
    }

    for (color=elements->minColor;color<=elements->maxColor;color++)
    {
	#pragma omp for private(e,irow,kr,kc,icol,isub) schedule(static)
	for (e=0;e<elements->numElements;e++) 
	{
		if (elements->Color[e]==color)
		{
			for (kr=0;kr<NN_row;kr++)
			{
				irow=row_map[elements->Nodes[INDEX2(row_node[INDEX2(kr,0,NN_row)],e,NN)]];
				for (kc=0;kc<NN_col;kc++)
				{
				icol=col_map[elements->Nodes[INDEX2(col_node[INDEX2(kc,0,NN_col)],e,NN)]];
				Dudley_IndexList_insertIndex(&(index_list[irow]),icol);
				}
			}
		}
		}
	}
  }
  return;
}


void Dudley_IndexList_insertElementsWithRowRange(Dudley_IndexList* index_list, index_t firstRow, index_t lastRow,
                                                 Dudley_ElementFile* elements, index_t* row_map, index_t* col_map)
{
/* this does not resolve macro elements */
	index_t color;
  dim_t e,kr,kc,icol,irow, NN;
  if (elements!=NULL) {
    NN=elements->numNodes;
    for (color=elements->minColor;color<=elements->maxColor;color++) {
           #pragma omp for private(e,irow,kr,kc,icol) schedule(static)
           for (e=0;e<elements->numElements;e++) {
               if (elements->Color[e]==color) {
                   for (kr=0;kr<NN;kr++) {
                     irow=row_map[elements->Nodes[INDEX2(kr,e,NN)]];
                     if ((firstRow<=irow) && (irow < lastRow)) {
                          irow-=firstRow;
                          for (kc=0;kc<NN;kc++) {
                              icol=col_map[elements->Nodes[INDEX2(kc,e,NN)]];
                              Dudley_IndexList_insertIndex(&(index_list[irow]),icol);
                          }
                      }
                  }
               }
           }
    }
  }
}
void Dudley_IndexList_insertElementsWithRowRangeNoMainDiagonal(Dudley_IndexList* index_list, index_t firstRow, index_t lastRow,
                                                              Dudley_ElementFile* elements, index_t* row_map, index_t* col_map)
{
  /* this does not resolve macro elements */
  index_t color;
  dim_t e,kr,kc,icol,irow, NN,irow_loc;
  if (elements!=NULL) {
    NN=elements->numNodes;
    for (color=elements->minColor;color<=elements->maxColor;color++) {
           #pragma omp for private(e,irow,kr,kc,icol,irow_loc) schedule(static)
           for (e=0;e<elements->numElements;e++) {
               if (elements->Color[e]==color) {
                   for (kr=0;kr<NN;kr++) {
                     irow=row_map[elements->Nodes[INDEX2(kr,e,NN)]];
                     if ((firstRow<=irow) && (irow < lastRow)) {
                          irow_loc=irow-firstRow;
                          for (kc=0;kc<NN;kc++) {
                              icol=col_map[elements->Nodes[INDEX2(kc,e,NN)]];
                              if (icol != irow) Dudley_IndexList_insertIndex(&(index_list[irow_loc]),icol);
                          }
                      }
                  }
               }
           }
    }
  }
}

/* inserts row index row into the Dudley_IndexList in if it does not exist */

void Dudley_IndexList_insertIndex(Dudley_IndexList* in, index_t index) {
  dim_t i;
  /* is index in in? */
  for (i=0;i<in->n;i++) {
    if (in->index[i]==index)  return;
  }
  /* index could not be found */
  if (in->n==INDEXLIST_LENGTH) {
     /* if in->index is full check the extension */
     if (in->extension==NULL) {
        in->extension=TMPMEMALLOC(1,Dudley_IndexList);
        if (Dudley_checkPtr(in->extension)) return;
        in->extension->n=0;
        in->extension->extension=NULL;
     }
     Dudley_IndexList_insertIndex(in->extension,index);
  } else {
     /* insert index into in->index*/
     in->index[in->n]=index;
     in->n++;
  }
}

/* counts the number of row indices in the Dudley_IndexList in */

dim_t Dudley_IndexList_count(Dudley_IndexList* in, index_t range_min,index_t range_max) {
  dim_t i;
  dim_t out=0;
  register index_t itmp;
  if (in==NULL) {
     return 0;
  } else {
    for (i=0;i<in->n;i++) {
          itmp=in->index[i];
          if ((itmp>=range_min) && (range_max>itmp)) ++out;
    }
     return out+Dudley_IndexList_count(in->extension, range_min,range_max);
  }
}

/* count the number of row indices in the Dudley_IndexList in */

void Dudley_IndexList_toArray(Dudley_IndexList* in, index_t* array, index_t range_min,index_t range_max, index_t index_offset) {
  dim_t i, ptr;
  register index_t itmp;
  if (in!=NULL) {
    ptr=0;
    for (i=0;i<in->n;i++) {
          itmp=in->index[i];
          if ((itmp>=range_min) && (range_max>itmp)) {
             array[ptr]=itmp+index_offset;
             ptr++;
          }

    }
    Dudley_IndexList_toArray(in->extension,&(array[ptr]), range_min, range_max, index_offset);
  }
}

/* deallocates the Dudley_IndexList in by recursive calls */

void Dudley_IndexList_free(Dudley_IndexList* in) {
  if (in!=NULL) {
    Dudley_IndexList_free(in->extension);
    TMPMEMFREE(in);
  }
}

/* creates a Paso_pattern from a range of indices */
Paso_Pattern* Dudley_IndexList_createPattern(dim_t n0, dim_t n,Dudley_IndexList* index_list,index_t range_min,index_t range_max,index_t index_offset)
{
   dim_t *ptr=NULL;
   register dim_t s,i,itmp;
   index_t *index=NULL;
   Paso_Pattern* out=NULL;

   ptr=MEMALLOC(n+1-n0,index_t);
   if (! Dudley_checkPtr(ptr) ) {
       /* get the number of connections per row */
       #pragma omp parallel for schedule(static) private(i)
       for(i=n0;i<n;++i) {
              ptr[i-n0]=Dudley_IndexList_count(&index_list[i],range_min,range_max);
       }
       /* accumulate ptr */
       s=0;
       for(i=n0;i<n;++i) {
               itmp=ptr[i-n0];
               ptr[i-n0]=s;
               s+=itmp;
       }
       ptr[n-n0]=s;
       /* fill index */
       index=MEMALLOC(ptr[n-n0],index_t);
       if (! Dudley_checkPtr(index)) {
              #pragma omp parallel for schedule(static)
              for(i=n0;i<n;++i) {
                  Dudley_IndexList_toArray(&index_list[i],&index[ptr[i-n0]],range_min,range_max,index_offset);
              }
              out=Paso_Pattern_alloc(PATTERN_FORMAT_DEFAULT,n-n0,range_max+index_offset,ptr,index);
       }
  }
  if (! Dudley_noError()) {
        MEMFREE(ptr);
        MEMFREE(index);
        Paso_Pattern_free(out);
  }
  return out;
}
