
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
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
  WARNING: MATRIX_FORMAT_CSC is not supported under MPI!

*****************************************************************************/

#include "Assemble.h"
#include <paso/SystemMatrix.h>

#ifdef ESYS_HAVE_TRILINOS
#include <trilinoswrap/TrilinosMatrixAdapter.h>

using esys_trilinos::TrilinosMatrixAdapter;
#endif

namespace finley {

void Assemble_addToSystemMatrix_CSC(paso::SystemMatrix* in, int NN_Equa,
                                    const index_t* Nodes_Equa, int num_Equa,
                                    int NN_Sol, const index_t* Nodes_Sol,
                                    int num_Sol, const double* array);

#ifdef ESYS_HAVE_TRILINOS
void Assemble_addToSystemMatrix_Trilinos(TrilinosMatrixAdapter* in,
                        int NN_Equa, const index_t* Nodes_Equa, int num_Equa,
                        int NN_Sol, const index_t* Nodes_Sol, int num_Sol,
                        const double* array);
#endif

void Assemble_addToSystemMatrix_CSR(paso::SystemMatrix* in, int NN_Equa,
                                    const index_t* Nodes_Equa, int num_Equa,
                                    int NN_Sol, const index_t* Nodes_Sol,
                                    int num_Sol, const double* array);

void Assemble_addToSystemMatrix(escript::ASM_ptr in, int NN_Equa,
                                const index_t* Nodes_Equa, int num_Equa,
                                int NN_Sol, const index_t* Nodes_Sol,
                                int num_Sol, const double* array)
{
    using paso::SystemMatrix;
    SystemMatrix* pmat(dynamic_cast<SystemMatrix*>(in.get()));

    if (pmat) {
        // call the right function depending on storage type
        if (pmat->type & MATRIX_FORMAT_CSC) {
            Assemble_addToSystemMatrix_CSC(pmat, NN_Equa, Nodes_Equa,
                                           num_Equa, NN_Sol, Nodes_Sol,
                                           num_Sol, array);
        } else { // type == CSR
            Assemble_addToSystemMatrix_CSR(pmat, NN_Equa, Nodes_Equa,
                                           num_Equa, NN_Sol, Nodes_Sol,
                                           num_Sol, array);
        }
        return;
    }
#ifdef ESYS_HAVE_TRILINOS
    TrilinosMatrixAdapter* tmat(dynamic_cast<TrilinosMatrixAdapter*>(in.get()));
    if (tmat) {
        Assemble_addToSystemMatrix_Trilinos(tmat, NN_Equa, Nodes_Equa,
                                            num_Equa, NN_Sol, Nodes_Sol,
                                            num_Sol, array);
        return;
    }
#endif
    throw FinleyException("Assemble_addToSystemMatrix: unknown system matrix type.");
}

void Assemble_addToSystemMatrix_CSC(paso::SystemMatrix* in, int NN_Equa,
                                    const index_t* Nodes_Equa, int num_Equa,
                                    int NN_Sol, const index_t* Nodes_Sol,
                                    int num_Sol, const double* array)
{
    const int index_offset=(in->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    const int row_block_size=in->row_block_size;
    const int col_block_size=in->col_block_size;
    const int block_size=in->block_size;
    const int num_subblocks_Equa=num_Equa/row_block_size;
    const int num_subblocks_Sol=num_Sol/col_block_size;
    const dim_t numMyCols=in->pattern->mainPattern->numInput;
    const dim_t numMyRows=in->pattern->mainPattern->numOutput;
    const index_t *mainBlock_ptr=in->mainBlock->pattern->ptr;
    const index_t *mainBlock_index=in->mainBlock->pattern->index;
    double *mainBlock_val=in->mainBlock->val;
    const index_t *col_coupleBlock_ptr=in->col_coupleBlock->pattern->ptr;
    const index_t *col_coupleBlock_index=in->col_coupleBlock->pattern->index;
    double *col_coupleBlock_val=in->col_coupleBlock->val;
    //const index_t *row_coupleBlock_ptr=in->row_coupleBlock->pattern->ptr;
    const index_t *row_coupleBlock_index=in->row_coupleBlock->pattern->index;
    double *row_coupleBlock_val=in->row_coupleBlock->val;

    for (int k_Sol=0; k_Sol<NN_Sol; ++k_Sol) {
        // Down columns of array
        const index_t j_Sol=Nodes_Sol[k_Sol];
        for (int l_col=0; l_col<num_subblocks_Sol; ++l_col) {
            const index_t i_col=j_Sol*num_subblocks_Sol+l_col;
            if (i_col < numMyCols) {
                for (int k_Equa=0;k_Equa<NN_Equa;++k_Equa) {
                    // Across cols of array
                    const index_t j_Equa=Nodes_Equa[k_Equa];
                    for (int l_row=0; l_row<num_subblocks_Equa; ++l_row) {
                        const index_t i_row=j_Equa*num_subblocks_Equa+index_offset+l_row;
                        if (i_row < numMyRows + index_offset ) {
                            for (index_t k=mainBlock_ptr[i_col]-index_offset; k<mainBlock_ptr[i_col+1]-index_offset; ++k) {
                                if (mainBlock_index[k]==i_row) {
                                    // Entry array(k_Equa, j_Sol) is a block (col_block_size x col_block_size)
                                    for (int ic=0; ic<col_block_size; ++ic) {
                                        const int i_Sol=ic+col_block_size*l_col;
                                        for (int ir=0; ir<row_block_size; ++ir) {
                                            const index_t i_Equa=ir+row_block_size*l_row;
                                            mainBlock_val[k*block_size+ir+row_block_size*ic]+=
                                                    array[INDEX4(i_Equa,i_Sol,k_Equa,k_Sol,num_Equa,num_Sol,NN_Equa)];
                                        }
                                    }
                                    break;
                                }
                            }
                        } else {
                            for (index_t k=col_coupleBlock_ptr[i_col]-index_offset; k<col_coupleBlock_ptr[i_col+1]-index_offset; ++k) {
                                if (row_coupleBlock_index[k] == i_row-numMyRows) {
                                    for (int ic=0; ic<col_block_size; ++ic) {
                                        const int i_Sol=ic+col_block_size*l_col;
                                        for (int ir=0; ir<row_block_size; ++ir) {
                                            const index_t i_Equa=ir+row_block_size*l_row;
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
                    const index_t j_Equa=Nodes_Equa[k_Equa];
                    for (int l_row=0; l_row<num_subblocks_Equa; ++l_row) {
                        const index_t i_row=j_Equa*num_subblocks_Equa+index_offset+l_row;
                        if (i_row < numMyRows + index_offset ) {
                            for (index_t k=col_coupleBlock_ptr[i_col-numMyCols]-index_offset; k<col_coupleBlock_ptr[i_col-numMyCols+1]-index_offset; ++k) {
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

void Assemble_addToSystemMatrix_CSR(paso::SystemMatrix* in, int NN_Equa,
                                    const index_t* Nodes_Equa, int num_Equa,
                                    int NN_Sol, const index_t* Nodes_Sol,
                                    int num_Sol, const double* array)
{
    const int index_offset=(in->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    const int row_block_size=in->row_block_size;
    const int col_block_size=in->col_block_size;
    const int block_size=in->block_size;
    const int num_subblocks_Equa=num_Equa/row_block_size;
    const int num_subblocks_Sol=num_Sol/col_block_size;
    const dim_t numMyCols=in->pattern->mainPattern->numInput;
    const dim_t numMyRows=in->pattern->mainPattern->numOutput;
    const index_t *mainBlock_ptr=in->mainBlock->pattern->ptr;
    const index_t *mainBlock_index=in->mainBlock->pattern->index;
    double *mainBlock_val=in->mainBlock->val;
    const index_t *col_coupleBlock_ptr=in->col_coupleBlock->pattern->ptr;
    const index_t *col_coupleBlock_index=in->col_coupleBlock->pattern->index;
    double *col_coupleBlock_val=in->col_coupleBlock->val;
    const index_t *row_coupleBlock_ptr=in->row_coupleBlock->pattern->ptr;
    const index_t *row_coupleBlock_index=in->row_coupleBlock->pattern->index;
    double *row_coupleBlock_val=in->row_coupleBlock->val;

    for (int k_Equa=0; k_Equa<NN_Equa; ++k_Equa) {
        // Down columns of array
        const index_t j_Equa=Nodes_Equa[k_Equa];
        for (int l_row=0; l_row<num_subblocks_Equa; ++l_row) {
            const index_t i_row=j_Equa*num_subblocks_Equa+l_row;
            // only look at the matrix rows stored on this processor
            if (i_row < numMyRows) {
                for (int k_Sol=0; k_Sol<NN_Sol; ++k_Sol) {
                    // Across rows of array
                    const index_t j_Sol=Nodes_Sol[k_Sol];
                    for (int l_col=0; l_col<num_subblocks_Sol; ++l_col) {
                        // only look at the matrix rows stored on this processor
                        const index_t i_col=j_Sol*num_subblocks_Sol+index_offset+l_col;
                        if (i_col < numMyCols + index_offset ) {
                            for (index_t k=mainBlock_ptr[i_row]-index_offset; k<mainBlock_ptr[i_row+1]-index_offset; ++k) {
                                if (mainBlock_index[k]==i_col) {
                                    // Entry array(k_Sol, j_Equa) is a block
                                    // (row_block_size x col_block_size)
                                    for (int ic=0; ic<col_block_size; ++ic) {
                                        const int i_Sol=ic+col_block_size*l_col;
                                        for (int ir=0; ir<row_block_size; ++ir) {
                                            const index_t i_Equa=ir+row_block_size*l_row;
                                            mainBlock_val[k*block_size+ir+row_block_size*ic]+=
                                                  array[INDEX4(i_Equa,i_Sol,k_Equa,k_Sol,num_Equa,num_Sol,NN_Equa)];
                                        }
                                    }
                                    break;
                                }
                            }
                        } else {
                            for (index_t k=col_coupleBlock_ptr[i_row]-index_offset; k<col_coupleBlock_ptr[i_row+1]-index_offset; ++k) {
                                if (col_coupleBlock_index[k] == i_col-numMyCols) {
                                    // Entry array(k_Sol, j_Equa) is a block
                                    // (row_block_size x col_block_size)
                                    for (int ic=0; ic<col_block_size; ++ic) {
                                        const index_t i_Sol=ic+col_block_size*l_col;
                                        for (int ir=0; ir<row_block_size; ++ir) {
                                            const index_t i_Equa=ir+row_block_size*l_row;
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
                    const index_t j_Sol=Nodes_Sol[k_Sol];
                    for (int l_col=0; l_col<num_subblocks_Sol; ++l_col) {
                        const index_t i_col=j_Sol*num_subblocks_Sol+index_offset+l_col;
                        if (i_col < numMyCols + index_offset ) {
                            for (index_t k=row_coupleBlock_ptr[i_row-numMyRows]-index_offset; k<row_coupleBlock_ptr[i_row-numMyRows+1]-index_offset; ++k) {
                                if (row_coupleBlock_index[k] == i_col) {
                                    // Entry array(k_Sol, j_Equa) is a block
                                    // (row_block_size x col_block_size)
                                    for (int ic=0; ic<col_block_size; ++ic) {
                                        const int i_Sol=ic+col_block_size*l_col;
                                        for (int ir=0; ir<row_block_size; ++ir) {
                                            const index_t i_Equa=ir+row_block_size*l_row;
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

#ifdef ESYS_HAVE_TRILINOS
void Assemble_addToSystemMatrix_Trilinos(TrilinosMatrixAdapter* in,
                                         int NN_Equa, const index_t* Nodes_Equa,
                                         int num_Equa, int NN_Sol,
                                         const index_t* Nodes_Sol, int num_Sol,
                                         const double* array)
{
    std::vector<index_t> rowIdx(Nodes_Equa, Nodes_Equa+NN_Equa);
    //std::vector<index_t> colIdx(Nodes_Sol, Nodes_Sol+NN_Sol);
    std::vector<double> arr(array, array+(NN_Equa*NN_Sol*num_Sol*num_Equa));
    in->add(rowIdx, arr);
}
#endif

} // namespace finley

