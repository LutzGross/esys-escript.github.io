/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrix and SystemVector */

/*  adds the matrix array[Equa,Sol,NN,NN] onto the matrix in. */
/* the rows/columns are given by */
/*  i_Equa+Equa*Nodes_Equa[Nodes[j_Equa]] (i_Equa=0:Equa; j_Equa=0:NN_Equa). */
/*  the routine has to be called from a parallel region                        */

/*  This routine assumes that in->Equa=in->Sol=1, i.e. */
/*  array is fully packed. */
/* TODO: the case in->Equa!=1  */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/*  Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "System.h"

/**************************************************************/

void  Finley_SystemMatrix_add(Finley_SystemMatrix* in,int NN_Equa,maybelong* Nodes_Equa, int num_Equa, 
                                                      int NN_Sol,maybelong* Nodes_Sol, int num_Sol, double* array) {
  int k_Equa,j_Equa,j_Sol,k_Sol,i_Equa,i_Sol,l_col,l_row,ic,ir,index,k,iptr;
  int row_block_size=in->row_block_size;
  int col_block_size=in->col_block_size;
  int block_size=in->block_size;
  int num_subblocks_Equa=num_Equa/row_block_size;
  int num_subblocks_Sol=num_Sol/col_block_size;

  if (in->type==CSR) {
          for (k_Equa=0;k_Equa<NN_Equa;++k_Equa) {
            j_Equa=Nodes_Equa[k_Equa];
            for (l_row=0;l_row<num_subblocks_Equa;++l_row) {
               iptr=j_Equa*num_subblocks_Equa+l_row;
               for (k_Sol=0;k_Sol<NN_Sol;++k_Sol) {
                 j_Sol=Nodes_Sol[k_Sol];
                 for (l_col=0;l_col<num_subblocks_Sol;++l_col) {
                    index=j_Sol*num_subblocks_Sol+INDEX_OFFSET+l_col;
	            for (k=in->pattern->ptr[iptr]-PTR_OFFSET;k<in->pattern->ptr[iptr+1]-PTR_OFFSET;++k) {
	                if (in->pattern->index[k]==index) {
                          for (ic=0;ic<col_block_size;++ic) {
                                i_Sol=ic+col_block_size*l_col;
                                for (ir=0;ir<row_block_size;++ir) {
                                   i_Equa=ir+row_block_size*l_row;
		                   in->val[k*block_size+ir+row_block_size*ic]+=
                                           array[INDEX4(i_Equa,i_Sol,k_Equa,k_Sol,num_Equa,num_Sol,NN_Equa)];
                                }
                          }
                          break;
                        }
                    }
                 }
               }
             }
          }
     } else {
         for (k_Sol=0;k_Sol<NN_Sol;k_Sol++) {
            j_Sol=Nodes_Sol[k_Sol];
            for (l_col=0;l_col<num_subblocks_Sol;++l_col) {
               iptr=j_Sol*num_subblocks_Sol+l_col;
               for (k_Equa=0;k_Equa<NN_Equa;k_Equa++) {
                 j_Equa=Nodes_Equa[k_Equa];
                 for (l_row=0;l_row<num_subblocks_Equa;++l_row) {
                    index=j_Equa*num_subblocks_Equa+INDEX_OFFSET+l_row;
	            for (k=in->pattern->ptr[iptr]-PTR_OFFSET;k<in->pattern->ptr[iptr+1]-PTR_OFFSET;++k) {
	                if (in->pattern->index[k]==index) {
                          for (ic=0;ic<col_block_size;++ic) {
                                i_Sol=ic+num_subblocks_Sol*l_col;
                                for (ir=0;ir<row_block_size;++ir) {
                                   i_Equa=ir+num_subblocks_Equa*l_row;
		                   in->val[k*block_size+ir+row_block_size*ic]+=
                                           array[INDEX4(i_Equa,i_Sol,k_Equa,k_Sol,num_Equa,num_Sol,NN_Equa)];
                                }
                          }
                          break;
                        }
                    }
                 }
               }
            }
         }
   }
}
/*
 * $Log$
 * Revision 1.4  2004/12/15 07:08:33  jgs
 * *** empty log message ***
 *
 *
 *
 */
