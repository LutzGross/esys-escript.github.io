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

/* Finley: SystemMatrix and SystemVector */

/*  adds the matrix array[Equa,Sol,NN,NN] onto the matrix in. */
/* the rows/columns are given by */
/*  i_Equa+Equa*Nodes_Equa[Nodes[j_Equa]] (i_Equa=0:Equa; j_Equa=0:NN_Equa). */
/*  the routine has to be called from a parallel region                        */

/*  This routine assumes that in->Equa=in->Sol=1, i.e. */
/*  array is fully packed. */
/* TODO: the case in->Equa!=1  */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "Assemble.h"
#include "IndexList.h"

/**************************************************************/

void  Finley_Assemble_addToSystemMatrix(Paso_SystemMatrix* in,dim_t NN_Equa,index_t* Nodes_Equa, dim_t num_Equa, 
                                        dim_t NN_Sol,index_t* Nodes_Sol, dim_t num_Sol, double* array) {
  index_t index_offset=(in->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  dim_t k_Equa,j_Equa,j_Sol,k_Sol,i_Equa,i_Sol,l_col,l_row,ic,ir,index,k,i_row, i_col;
  index_t *mainBlock_ptr, *mainBlock_index, *coupleBlock_ptr, *coupleBlock_index;
  double *mainBlock_val, *coupleBlock_val;
  dim_t row_block_size=in->row_block_size;
  dim_t col_block_size=in->col_block_size;
  dim_t block_size=in->block_size;
  dim_t num_subblocks_Equa=num_Equa/row_block_size;
  dim_t num_subblocks_Sol=num_Sol/col_block_size;
  dim_t numMyCols=in->pattern->mainPattern->numInput;
  dim_t numMyRows=in->pattern->mainPattern->numOutput;


  if (in->type & MATRIX_FORMAT_CSC) {
         /* MATRIX_FORMAT_CSC does not support MPI !!!!! */
         mainBlock_ptr=in->mainBlock->pattern->ptr;
         mainBlock_index=in->mainBlock->pattern->index;
         mainBlock_val=in->mainBlock->val;
         for (k_Sol=0;k_Sol<NN_Sol;k_Sol++) {
            j_Sol=Nodes_Sol[k_Sol];
            for (l_col=0;l_col<num_subblocks_Sol;++l_col) {
               i_col=j_Sol*num_subblocks_Sol+l_col;
               for (k_Equa=0;k_Equa<NN_Equa;k_Equa++) {
                 j_Equa=Nodes_Equa[k_Equa];
                 for (l_row=0;l_row<num_subblocks_Equa;++l_row) {
                    i_row=j_Equa*num_subblocks_Equa+index_offset+l_row;
	            for (k=mainBlock_ptr[i_col]-index_offset; k<mainBlock_ptr[i_col+1]-index_offset;++k) {
	                if (mainBlock_index[k] == i_row) {
                          for (ic=0;ic<col_block_size;++ic) {
                                i_Sol=ic+col_block_size*l_col;
                                for (ir=0;ir<row_block_size;++ir) {
                                   i_Equa=ir+row_block_size*l_row;
		                   mainBlock_val[k*block_size+ir+row_block_size*ic]+=
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
   } else if (in->type & MATRIX_FORMAT_TRILINOS_CRS) {
       /* this needs to be modified */
       #ifdef TRILINOS
          for (k_Equa=0;k_Equa<NN_Equa;++k_Equa) { /* Down columns of array */
            j_Equa=Nodes_Equa[k_Equa];
	    if (j_Equa < in->mainBlock->pattern->output_node_distribution->numLocal) {
              for (k_Sol=0;k_Sol<NN_Sol;++k_Sol) { /* Across rows of array */
                j_Sol=Nodes_Sol[k_Sol];
                for (l_row=0;l_row<num_subblocks_Equa;++l_row) {
                  irow=j_Equa*row_block_size+l_row;
                  for (l_col=0;l_col<col_block_size;++l_col) {
                     icol=j_Sol*col_block_size+index_offset+l_col;
		     // irow is local and icol is global
		     Trilinos_SumIntoMyValues(in->trilinos_data, irow, icol, array[INDEX4(l_row,l_col,k_Equa,k_Sol,num_Equa,num_Sol,NN_Equa)]);
	          }
                }
              }
            }
          }
       #endif
   } else {

         mainBlock_ptr=in->mainBlock->pattern->ptr;
         mainBlock_index=in->mainBlock->pattern->index;
         mainBlock_val=in->mainBlock->val;
         coupleBlock_ptr=in->coupleBlock->pattern->ptr;
         coupleBlock_index=in->coupleBlock->pattern->index;
         coupleBlock_val=in->coupleBlock->val;

         for (k_Equa=0;k_Equa<NN_Equa;++k_Equa) { /* Down columns of array */
                j_Equa=Nodes_Equa[k_Equa];
                for (l_row=0;l_row<num_subblocks_Equa;++l_row) {
                   i_row=j_Equa*num_subblocks_Equa+l_row;
                   /* only look at the matrix rows stored on this processor */
                   if (i_row < numMyRows) {
                      for (k_Sol=0;k_Sol<NN_Sol;++k_Sol) { /* Across rows of array */
                        j_Sol=Nodes_Sol[k_Sol];
                        for (l_col=0;l_col<num_subblocks_Sol;++l_col) {
                           /* only look at the matrix rows stored on this processor */
                           i_col=j_Sol*num_subblocks_Sol+index_offset+l_col;
                           if (i_col < numMyCols + index_offset ) {
	                       for (k=mainBlock_ptr[i_row]-index_offset;k<mainBlock_ptr[i_row+1]-index_offset;++k) {
	                           if (mainBlock_index[k]==i_col) {
                                     /* Entry array(k_Sol, j_Equa) is a block (row_block_size x col_block_size) */
                                     for (ic=0;ic<col_block_size;++ic) { 
                                           i_Sol=ic+col_block_size*l_col;
                                           for (ir=0;ir<row_block_size;++ir) {
                                              i_Equa=ir+row_block_size*l_row;
		                              mainBlock_val[k*block_size+ir+row_block_size*ic]+=
                                                      array[INDEX4(i_Equa,i_Sol,k_Equa,k_Sol,num_Equa,num_Sol,NN_Equa)];
                                           }
                                     }
                                     break;
                                   }
                               }
                           } else {
	                       for (k=coupleBlock_ptr[i_row]-index_offset;k<coupleBlock_ptr[i_row+1]-index_offset;++k) {
	                           if (coupleBlock_index[k] == i_col-numMyCols) {
                                     /* Entry array(k_Sol, j_Equa) is a block (row_block_size x col_block_size) */
                                     for (ic=0;ic<col_block_size;++ic) { 
                                           i_Sol=ic+col_block_size*l_col;
                                           for (ir=0;ir<row_block_size;++ir) {
                                              i_Equa=ir+row_block_size*l_row;
		                              coupleBlock_val[k*block_size+ir+row_block_size*ic]+=
                                                      array[INDEX4(i_Equa,i_Sol,k_Equa,k_Sol,num_Equa,num_Sol,NN_Equa)];
                                           }
                                     }
                                     break;
                                   }
                               }
                           }
                        }
                      }
                    } /* end i_row check */
              }
        }
   }
}
