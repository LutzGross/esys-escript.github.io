/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrix and SystemVector */

/*  adds the matrix array[row_block_size,col_block_size,NN,NN] onto the matrix in. */
/* the rows/columns are given by */
/*  ir+row_block_size*Nodes_row[Nodes[jr]] (ir=0:row_block_size; jr=0:NN_row). */
/*  the routine has to be called from a parallel region                        */

/*  This routine assumes that in->row_block_size=in->col_block_size=1, i.e. */
/*  array is fully packed. */
/* TODO: the case in->row_block_size!=1  */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/*  Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "System.h"

/**************************************************************/

void  Finley_SystemMatrix_add(Finley_SystemMatrix* in,int NN_row,maybelong* Nodes_row, int row_block_size, 
                                                      int NN_col,maybelong* Nodes_col,int col_block_size, double* array) {
  int kr,jr,jc,kc,ir,ic;
  maybelong irow,icol,k;
  switch(in->type) {
     case CSR:
          for (kr=0;kr<NN_row;kr++) {
            jr=Nodes_row[kr];
            for (ir=0;ir<row_block_size;ir++) {
               irow=ir+row_block_size*jr;
               for (kc=0;kc<NN_col;kc++) {
	         jc=Nodes_col[kc];
	         for (ic=0;ic<col_block_size;ic++) {
	            icol=ic+col_block_size*jc+INDEX_OFFSET;
	            for(k=in->ptr[irow]-PTR_OFFSET;k<in->ptr[irow+1]-PTR_OFFSET;k++) {
	               if (in->index[k]==icol) {
                          #pragma omp atomic
		          in->val[k]+=array[INDEX4(ir,ic,kr,kc,row_block_size,col_block_size,NN_row)];
                          break;
                       }
	            }
	         }
              } 
          }
      }
      break;
   case CSC:
        for (kr=0;kr<NN_row;kr++) {
          jr=Nodes_row[kr];
          for (ir=0;ir<row_block_size;ir++) {
             irow=ir+row_block_size*jr;
             for (kc=0;kc<NN_col;kc++) {
	        jc=Nodes_col[kc];
	        for (ic=0;ic<col_block_size;ic++) {
	           icol=ic+col_block_size*jc+INDEX_OFFSET;
	           for(k=in->ptr[icol]-PTR_OFFSET;k<in->ptr[icol+1]-PTR_OFFSET;k++) {
	              if (in->index[k]==irow) {
                         #pragma omp atomic
		         in->val[k]+=array[INDEX4(ir,ic,kr,kc,row_block_size,col_block_size,NN_row)];
                         break;
                      }
                   }
                }
	     }
	  }
        }
	break;
   default:
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"Unknown matrix type.");
  } /* switch in->type */
}
/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
