
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/****************************************************************************

  Finley: SystemMatrix

  Adds the matrix array[Equa,Sol,NN,NN] to the matrix in.
  The rows/columns are given by
    i_Equa+Equa*Nodes_Equa[Nodes[j_Equa]] (i_Equa=0:Equa; j_Equa=0:NN_Equa).

  The routine has to be called from a parallel region.
  This routine assumes that in->Equa=in->Sol=1, i.e. array is fully packed.
  TODO: the case in->Equa!=1
  WARNING: MATRIX_FORMAT_CSC does not support MPI!!

*****************************************************************************/

#include "Assemble.h"

namespace finley {

void Assemble_addToSystemMatrix_CSC(Paso_SystemMatrix* in,
        const int NN_Equa, const int* Nodes_Equa, const int num_Equa,
        const int NN_Sol, const int* Nodes_Sol, const int num_Sol,
        const double* array);

void Assemble_addToSystemMatrix_Trilinos(Paso_SystemMatrix* in,
        const int NN_Equa, const int* Nodes_Equa, const int num_Equa,
        const int NN_Sol, const int* Nodes_Sol, const int num_Sol,
        const double* array);

void Assemble_addToSystemMatrix_CSR(Paso_SystemMatrix* in,
        const int NN_Equa, const int* Nodes_Equa, const int num_Equa,
        const int NN_Sol, const int* Nodes_Sol, const int num_Sol,
        const double* array);

void Assemble_addToSystemMatrix(Paso_SystemMatrix* in,
        const int NN_Equa, const int* Nodes_Equa, const int num_Equa,
        const int NN_Sol, const int* Nodes_Sol, const int num_Sol,
        const double* array)
{
    // call the right function depending on storage type
    if (in->type & MATRIX_FORMAT_CSC) {
        Assemble_addToSystemMatrix_CSC(in, NN_Equa, Nodes_Equa,
                                  num_Equa, NN_Sol, Nodes_Sol, num_Sol, array);
    } else if (in->type & MATRIX_FORMAT_TRILINOS_CRS) {
        Assemble_addToSystemMatrix_Trilinos(in, NN_Equa, Nodes_Equa,
                                  num_Equa, NN_Sol, Nodes_Sol, num_Sol, array);
    } else { // type == CSR
        Assemble_addToSystemMatrix_CSR(in, NN_Equa, Nodes_Equa,
                                  num_Equa, NN_Sol, Nodes_Sol, num_Sol, array);
    }
}

void Assemble_addToSystemMatrix_CSC(Paso_SystemMatrix* in, const int NN_Equa,
                                    const int* Nodes_Equa, const int num_Equa,
                                    const int NN_Sol, const int* Nodes_Sol,
                                    const int num_Sol, const double* array)
{
    const int index_offset=(in->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    const int row_block_size=in->row_block_size;
    const int col_block_size=in->col_block_size;
    const int block_size=in->block_size;
    const int num_subblocks_Equa=num_Equa/row_block_size;
    const int num_subblocks_Sol=num_Sol/col_block_size;
    const int numMyCols=in->pattern->mainPattern->numInput;
    const int numMyRows=in->pattern->mainPattern->numOutput;
    const int *mainBlock_ptr=in->mainBlock->pattern->ptr;
    const int *mainBlock_index=in->mainBlock->pattern->index;
    double *mainBlock_val=in->mainBlock->val;
    const int *col_coupleBlock_ptr=in->col_coupleBlock->pattern->ptr;
    const int *col_coupleBlock_index=in->col_coupleBlock->pattern->index;
    double *col_coupleBlock_val=in->col_coupleBlock->val;
    //const int *row_coupleBlock_ptr=in->row_coupleBlock->pattern->ptr;
    const int *row_coupleBlock_index=in->row_coupleBlock->pattern->index;
    double *row_coupleBlock_val=in->row_coupleBlock->val;

    for (int k_Sol=0; k_Sol<NN_Sol; ++k_Sol) {
        // Down columns of array
        const int j_Sol=Nodes_Sol[k_Sol];
        for (int l_col=0; l_col<num_subblocks_Sol; ++l_col) {
            const int i_col=j_Sol*num_subblocks_Sol+l_col;
            if (i_col < numMyCols) {
                for (int k_Equa=0;k_Equa<NN_Equa;++k_Equa) {
                    // Across cols of array
                    const int j_Equa=Nodes_Equa[k_Equa];
                    for (int l_row=0; l_row<num_subblocks_Equa; ++l_row) {
                        const int i_row=j_Equa*num_subblocks_Equa+index_offset+l_row;
                        if (i_row < numMyRows + index_offset ) {
                            for (int k=mainBlock_ptr[i_col]-index_offset; k<mainBlock_ptr[i_col+1]-index_offset; ++k) {
                                if (mainBlock_index[k]==i_row) {
                                    // Entry array(k_Equa, j_Sol) is a block (col_block_size x col_block_size)
                                    for (int ic=0; ic<col_block_size; ++ic) {
                                        const int i_Sol=ic+col_block_size*l_col;
                                        for (int ir=0; ir<row_block_size; ++ir) {
                                            const int i_Equa=ir+row_block_size*l_row;
                                            mainBlock_val[k*block_size+ir+row_block_size*ic]+=
                                                    array[INDEX4(i_Equa,i_Sol,k_Equa,k_Sol,num_Equa,num_Sol,NN_Equa)];
                                        }
                                    }
                                    break;
                                }
                            }
                        } else {
                            for (int k=col_coupleBlock_ptr[i_col]-index_offset; k<col_coupleBlock_ptr[i_col+1]-index_offset; ++k) {
                                if (row_coupleBlock_index[k] == i_row-numMyRows) {
                                    for (int ic=0; ic<col_block_size; ++ic) {
                                        const int i_Sol=ic+col_block_size*l_col;
                                        for (int ir=0; ir<row_block_size; ++ir) {
                                            const int i_Equa=ir+row_block_size*l_row;
                                            row_coupleBlock_val[k*block_size+ir+row_block_size*ic]+=
                                                array[INDEX4(i_Equa,i_Sol,k_Equa,k_Sol,num_Equa,num_Sol,NN_Equa)];
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            } else { // i_col >= numMyCols
                for (int k_Equa=0;k_Equa<NN_Equa;++k_Equa) {
                    // Across rows of array
                    const int j_Equa=Nodes_Equa[k_Equa];
                    for (int l_row=0; l_row<num_subblocks_Equa; ++l_row) {
                        const int i_row=j_Equa*num_subblocks_Equa+index_offset+l_row;
                        if (i_row < numMyRows + index_offset ) {
                            for (int k=col_coupleBlock_ptr[i_col-numMyCols]-index_offset; k<col_coupleBlock_ptr[i_col-numMyCols+1]-index_offset; ++k) {
                                if (col_coupleBlock_index[k] == i_row) {
                                    for (int ic=0; ic<col_block_size; ++ic) {
                                        const int i_Sol=ic+col_block_size*l_col;
                                        for (int ir=0; ir<row_block_size; ++ir) {
                                            const int i_Equa=ir+row_block_size*l_row;
                                            col_coupleBlock_val[k*block_size+ir+row_block_size*ic]+=
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
}

void Assemble_addToSystemMatrix_Trilinos(Paso_SystemMatrix* in,
                                         const int NN_Equa,
                                         const int* Nodes_Equa,
                                         const int num_Equa,
                                         const int NN_Sol,
                                         const int* Nodes_Sol,
                                         const int num_Sol,
                                         const double* array)
{
    // FIXME: this needs to be modified
#ifdef TRILINOS
    const int index_offset=(in->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    const int row_block_size=in->row_block_size;
    const int col_block_size=in->col_block_size;
    const int num_subblocks_Equa=num_Equa/row_block_size;
    for (int k_Equa=0;k_Equa<NN_Equa;++k_Equa) { // Down columns of array
        const int j_Equa=Nodes_Equa[k_Equa];
        if (j_Equa < in->mainBlock->pattern->output_node_distribution->numLocal) {
            for (int k_Sol=0; k_Sol<NN_Sol; ++k_Sol) {
                // Across rows of array
                const int j_Sol=Nodes_Sol[k_Sol];
                for (int l_row=0; l_row<num_subblocks_Equa; ++l_row) {
                    const int irow=j_Equa*row_block_size+l_row;
                    for (int l_col=0; l_col<col_block_size; ++l_col) {
                        const int icol=j_Sol*col_block_size+index_offset+l_col;
                        // irow is local and icol is global
                        Trilinos_SumIntoMyValues(in->trilinos_data, irow, icol, array[INDEX4(l_row,l_col,k_Equa,k_Sol,num_Equa,num_Sol,NN_Equa)]);
                    }
                }
            }
        }
    }
#endif
}

void Assemble_addToSystemMatrix_CSR(Paso_SystemMatrix* in, const int NN_Equa,
                                    const int* Nodes_Equa, const int num_Equa,
                                    const int NN_Sol, const int* Nodes_Sol,
                                    const int num_Sol, const double* array)
{
    const int index_offset=(in->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    const int row_block_size=in->row_block_size;
    const int col_block_size=in->col_block_size;
    const int block_size=in->block_size;
    const int num_subblocks_Equa=num_Equa/row_block_size;
    const int num_subblocks_Sol=num_Sol/col_block_size;
    const int numMyCols=in->pattern->mainPattern->numInput;
    const int numMyRows=in->pattern->mainPattern->numOutput;
    const int *mainBlock_ptr=in->mainBlock->pattern->ptr;
    const int *mainBlock_index=in->mainBlock->pattern->index;
    double *mainBlock_val=in->mainBlock->val;
    const int *col_coupleBlock_ptr=in->col_coupleBlock->pattern->ptr;
    const int *col_coupleBlock_index=in->col_coupleBlock->pattern->index;
    double *col_coupleBlock_val=in->col_coupleBlock->val;
    const int *row_coupleBlock_ptr=in->row_coupleBlock->pattern->ptr;
    const int *row_coupleBlock_index=in->row_coupleBlock->pattern->index;
    double *row_coupleBlock_val=in->row_coupleBlock->val;

    for (int k_Equa=0; k_Equa<NN_Equa; ++k_Equa) {
        // Down columns of array
        const int j_Equa=Nodes_Equa[k_Equa];
        for (int l_row=0; l_row<num_subblocks_Equa; ++l_row) {
            const int i_row=j_Equa*num_subblocks_Equa+l_row;
            // only look at the matrix rows stored on this processor
            if (i_row < numMyRows) {
                for (int k_Sol=0; k_Sol<NN_Sol; ++k_Sol) {
                    // Across rows of array
                    const int j_Sol=Nodes_Sol[k_Sol];
                    for (int l_col=0; l_col<num_subblocks_Sol; ++l_col) {
                        // only look at the matrix rows stored on this processor
                        const int i_col=j_Sol*num_subblocks_Sol+index_offset+l_col;
                        if (i_col < numMyCols + index_offset ) {
                            for (int k=mainBlock_ptr[i_row]-index_offset; k<mainBlock_ptr[i_row+1]-index_offset; ++k) {
                                if (mainBlock_index[k]==i_col) {
                                    // Entry array(k_Sol, j_Equa) is a block
                                    // (row_block_size x col_block_size)
                                    for (int ic=0; ic<col_block_size; ++ic) {
                                        const int i_Sol=ic+col_block_size*l_col;
                                        for (int ir=0; ir<row_block_size; ++ir) {
                                            const int i_Equa=ir+row_block_size*l_row;
                                            mainBlock_val[k*block_size+ir+row_block_size*ic]+=
                                                  array[INDEX4(i_Equa,i_Sol,k_Equa,k_Sol,num_Equa,num_Sol,NN_Equa)];
                                        }
                                    }
                                    break;
                                }
                            }
                        } else {
                            for (int k=col_coupleBlock_ptr[i_row]-index_offset; k<col_coupleBlock_ptr[i_row+1]-index_offset; ++k) {
                                if (col_coupleBlock_index[k] == i_col-numMyCols) {
                                    // Entry array(k_Sol, j_Equa) is a block
                                    // (row_block_size x col_block_size)
                                    for (int ic=0; ic<col_block_size; ++ic) {
                                        const int i_Sol=ic+col_block_size*l_col;
                                        for (int ir=0; ir<row_block_size; ++ir) {
                                            const int i_Equa=ir+row_block_size*l_row;
                                            col_coupleBlock_val[k*block_size+ir+row_block_size*ic]+=
                                                  array[INDEX4(i_Equa,i_Sol,k_Equa,k_Sol,num_Equa,num_Sol,NN_Equa)];
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            } else { // i_row >= numMyRows
                for (int k_Sol=0; k_Sol<NN_Sol; ++k_Sol) {
                    // Across rows of array
                    const int j_Sol=Nodes_Sol[k_Sol];
                    for (int l_col=0; l_col<num_subblocks_Sol; ++l_col) {
                        const int i_col=j_Sol*num_subblocks_Sol+index_offset+l_col;
                        if (i_col < numMyCols + index_offset ) {
                            for (int k=row_coupleBlock_ptr[i_row-numMyRows]-index_offset; k<row_coupleBlock_ptr[i_row-numMyRows+1]-index_offset; ++k) {
                                if (row_coupleBlock_index[k] == i_col) {
                                    // Entry array(k_Sol, j_Equa) is a block
                                    // (row_block_size x col_block_size)
                                    for (int ic=0; ic<col_block_size; ++ic) {
                                        const int i_Sol=ic+col_block_size*l_col;
                                        for (int ir=0; ir<row_block_size; ++ir) {
                                            const int i_Equa=ir+row_block_size*l_row;
                                            row_coupleBlock_val[k*block_size+ir+row_block_size*ic]+=
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
}

} // namespace finley

