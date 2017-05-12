
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
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

#include "Assemble.h"

#ifdef ESYS_HAVE_PASO
#include <paso/SystemMatrix.h>
#endif

#ifdef ESYS_HAVE_TRILINOS
#include <trilinoswrap/TrilinosMatrixAdapter.h>

using esys_trilinos::TrilinosMatrixAdapter;
#endif

namespace dudley {

using escript::DataTypes::real_t;
using escript::DataTypes::cplx_t;

#ifdef ESYS_HAVE_PASO
static void addToSystemMatrixPasoCSC(paso::SystemMatrix* S,
                                     const std::vector<index_t>& Nodes,
                                     int numEq,
                                     const std::vector<double>& array);

static void addToSystemMatrixPasoCSR(paso::SystemMatrix* S,
                                     const std::vector<index_t>& Nodes,
                                     int numEq,
                                     const std::vector<double>& array);
#endif

template<>
void Assemble_addToSystemMatrix<real_t>(escript::AbstractSystemMatrix* S,
                                        const std::vector<index_t>& Nodes,
                                        int numEq,
                                        const std::vector<real_t>& array)
{
#ifdef ESYS_HAVE_PASO
    paso::SystemMatrix* pmat = dynamic_cast<paso::SystemMatrix*>(S);
    if (pmat) {
        // call the right function depending on storage type
        if (pmat->type & MATRIX_FORMAT_CSC) {
            addToSystemMatrixPasoCSC(pmat, Nodes, numEq, array);
        } else { // type == CSR
            addToSystemMatrixPasoCSR(pmat, Nodes, numEq, array);
        }
        return;
    }
#endif
#ifdef ESYS_HAVE_TRILINOS
    TrilinosMatrixAdapter* tmat = dynamic_cast<TrilinosMatrixAdapter*>(S);
    if (tmat) {
        tmat->add(Nodes, array);
        return;
    }
#endif
    throw DudleyException("Assemble_addToSystemMatrix: unsupported system "
                          "matrix type.");
}

template<>
void Assemble_addToSystemMatrix<cplx_t>(escript::AbstractSystemMatrix* S,
                                        const std::vector<index_t>& Nodes,
                                        int numEq,
                                        const std::vector<cplx_t>& array)
{
#ifdef ESYS_HAVE_TRILINOS
    TrilinosMatrixAdapter* tmat = dynamic_cast<TrilinosMatrixAdapter*>(S);
    if (tmat) {
        tmat->add(Nodes, array);
        return;
    }
#endif
    throw DudleyException("addToSystemMatrix: only Trilinos matrices support "
                          "complex-valued assembly!");
}

#ifdef ESYS_HAVE_PASO
void addToSystemMatrixPasoCSC(paso::SystemMatrix* in,
                              const std::vector<index_t>& Nodes,
                              int numEq, const std::vector<double>& array)
{
    const int index_offset = (in->type & MATRIX_FORMAT_OFFSET1 ? 1 : 0);
    const int row_block_size = in->row_block_size;
    const int col_block_size = in->col_block_size;
    const int block_size = in->block_size;
    const int num_subblocks_Eq = numEq / row_block_size;
    const int num_subblocks_Sol = numEq / col_block_size;
    const dim_t numMyCols = in->pattern->mainPattern->numInput;
    const dim_t numMyRows = in->pattern->mainPattern->numOutput;
    const int NN = Nodes.size();

    const index_t* mainBlock_ptr = in->mainBlock->pattern->ptr;
    const index_t* mainBlock_index = in->mainBlock->pattern->index;
    double* mainBlock_val = in->mainBlock->val;
    const index_t* col_coupleBlock_ptr = in->col_coupleBlock->pattern->ptr;
    const index_t* col_coupleBlock_index = in->col_coupleBlock->pattern->index;
    double* col_coupleBlock_val = in->col_coupleBlock->val;
    const index_t* row_coupleBlock_index = in->row_coupleBlock->pattern->index;
    double* row_coupleBlock_val = in->row_coupleBlock->val;

    for (int k_Sol = 0; k_Sol < NN; ++k_Sol) { // Down columns of array
        const index_t j_Sol = Nodes[k_Sol];
        for (int l_col = 0; l_col < num_subblocks_Sol; ++l_col) {
            const index_t i_col = j_Sol * num_subblocks_Sol + l_col;
            if (i_col < numMyCols) {
                for (int k_Eq = 0; k_Eq < NN; ++k_Eq) {
                    // Across cols of array
                    const index_t j_Eq = Nodes[k_Eq];
                    for (int l_row = 0; l_row < num_subblocks_Eq; ++l_row) {
                        const index_t i_row = j_Eq * num_subblocks_Eq + index_offset + l_row;
                        if (i_row < numMyRows + index_offset) {
                            for (index_t k = mainBlock_ptr[i_col]-index_offset;
                                 k < mainBlock_ptr[i_col + 1]-index_offset; ++k) {
                                if (mainBlock_index[k] == i_row) {
                                    // Entry array(k_Eq, j_Sol) is a block
                                    // (col_block_size x col_block_size)
                                    for (int ic = 0; ic < col_block_size; ++ic) {
                                        const int i_Sol = ic + col_block_size * l_col;
                                        for (int ir = 0; ir < row_block_size; ++ir) {
                                            const int i_Eq = ir + row_block_size * l_row;
                                            mainBlock_val[k*block_size + ir + row_block_size*ic] +=
                                                array[INDEX4
                                                  (i_Eq, i_Sol, k_Eq, k_Sol, numEq, numEq, NN)];
                                        }
                                    }
                                    break;
                                }
                            }
                        } else {
                            for (index_t k = col_coupleBlock_ptr[i_col]-index_offset;
                                 k < col_coupleBlock_ptr[i_col + 1]-index_offset; ++k) {
                                if (row_coupleBlock_index[k] == i_row - numMyRows) {
                                    for (int ic = 0; ic < col_block_size; ++ic) {
                                        const int i_Sol = ic + col_block_size * l_col;
                                        for (int ir = 0; ir < row_block_size; ++ir) {
                                            const int i_Eq = ir + row_block_size * l_row;
                                            row_coupleBlock_val[k*block_size + ir + row_block_size*ic] +=
                                                array[INDEX4
                                                  (i_Eq, i_Sol, k_Eq, k_Sol, numEq, numEq, NN)];
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            } else { // i_col >= numMyCols
                for (int k_Eq = 0; k_Eq < NN; ++k_Eq) {
                    // Across rows of array
                    const index_t j_Eq = Nodes[k_Eq];
                    for (int l_row = 0; l_row < num_subblocks_Eq; ++l_row) {
                        const index_t i_row = j_Eq * num_subblocks_Eq + index_offset + l_row;
                        if (i_row < numMyRows + index_offset) {
                            for (index_t k = col_coupleBlock_ptr[i_col-numMyCols]-index_offset;
                                 k < col_coupleBlock_ptr[i_col - numMyCols + 1] - index_offset; ++k) {
                                if (col_coupleBlock_index[k] == i_row) {
                                    for (int ic = 0; ic < col_block_size; ++ic) {
                                        const int i_Sol = ic + col_block_size * l_col;
                                        for (int ir = 0; ir < row_block_size; ++ir) {
                                            const int i_Eq = ir + row_block_size * l_row;
                                            col_coupleBlock_val[k*block_size + ir + row_block_size*ic] +=
                                                array[INDEX4
                                                  (i_Eq, i_Sol, k_Eq, k_Sol, numEq, numEq, NN)];
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

void addToSystemMatrixPasoCSR(paso::SystemMatrix* in,
                              const std::vector<index_t>& Nodes,
                              int numEq, const std::vector<double>& array)
{
    const int index_offset = (in->type & MATRIX_FORMAT_OFFSET1 ? 1 : 0);
    const int row_block_size = in->row_block_size;
    const int col_block_size = in->col_block_size;
    const int block_size = in->block_size;
    const int num_subblocks_Eq = numEq / row_block_size;
    const int num_subblocks_Sol = numEq / col_block_size;
    const dim_t numMyCols = in->pattern->mainPattern->numInput;
    const dim_t numMyRows = in->pattern->mainPattern->numOutput;
    const int NN = Nodes.size();

    const index_t* mainBlock_ptr = in->mainBlock->pattern->ptr;
    const index_t* mainBlock_index = in->mainBlock->pattern->index;
    double* mainBlock_val = in->mainBlock->val;
    const index_t* col_coupleBlock_ptr = in->col_coupleBlock->pattern->ptr;
    const index_t* col_coupleBlock_index = in->col_coupleBlock->pattern->index;
    double* col_coupleBlock_val = in->col_coupleBlock->val;
    const index_t* row_coupleBlock_ptr = in->row_coupleBlock->pattern->ptr;
    const index_t* row_coupleBlock_index = in->row_coupleBlock->pattern->index;
    double* row_coupleBlock_val = in->row_coupleBlock->val;

    for (int k_Eq = 0; k_Eq < NN; ++k_Eq) { // Down columns of array
        const index_t j_Eq = Nodes[k_Eq];
        for (int l_row = 0; l_row < num_subblocks_Eq; ++l_row) {
        const index_t i_row = j_Eq * num_subblocks_Eq + l_row;
        // only look at the matrix rows stored on this processor
        if (i_row < numMyRows) {
            for (int k_Sol = 0; k_Sol < NN; ++k_Sol) { // Across rows of array
            const index_t j_Sol = Nodes[k_Sol];
            for (int l_col = 0; l_col < num_subblocks_Sol; ++l_col) {
                // only look at the matrix rows stored on this processor
                const index_t i_col = j_Sol * num_subblocks_Sol + index_offset + l_col;
                if (i_col < numMyCols + index_offset) {
                    for (index_t k = mainBlock_ptr[i_row] - index_offset;
                         k < mainBlock_ptr[i_row + 1] - index_offset; ++k) {
                        if (mainBlock_index[k] == i_col) {
                            // Entry array(k_Sol, j_Eq) is a block
                            // (row_block_size x col_block_size)
                            for (int ic = 0; ic < col_block_size; ++ic) {
                                const int i_Sol = ic + col_block_size * l_col;
                                for (int ir = 0; ir < row_block_size; ++ir) {
                                    const int i_Eq = ir + row_block_size * l_row;
                                    mainBlock_val[k*block_size + ir + row_block_size*ic] +=
                                        array[INDEX4
                                          (i_Eq, i_Sol, k_Eq, k_Sol, numEq, numEq, NN)];
                                }
                            }
                            break;
                        }
                    }
                } else {
                    for (index_t k = col_coupleBlock_ptr[i_row] - index_offset;
                         k < col_coupleBlock_ptr[i_row + 1] - index_offset; ++k) {
                        if (col_coupleBlock_index[k] == i_col - numMyCols) {
                            // Entry array(k_Sol, j_Eq) is a block
                            // (row_block_size x col_block_size)
                            for (int ic = 0; ic < col_block_size; ++ic) {
                                const int i_Sol = ic + col_block_size * l_col;
                                for (int ir = 0; ir < row_block_size; ++ir) {
                                    const int i_Eq = ir+row_block_size*l_row;
                                    col_coupleBlock_val[k*block_size + ir + row_block_size*ic] +=
                                        array[INDEX4
                                          (i_Eq, i_Sol, k_Eq, k_Sol, numEq, numEq, NN)];
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
            }
        } else {
            for (int k_Sol = 0; k_Sol < NN; ++k_Sol) { // Across rows of array
                const index_t j_Sol = Nodes[k_Sol];
                for (int l_col = 0; l_col < num_subblocks_Sol; ++l_col) {
                    const index_t i_col = j_Sol * num_subblocks_Sol + index_offset + l_col;
                    if (i_col < numMyCols + index_offset) {
                        for (index_t k = row_coupleBlock_ptr[i_row - numMyRows] - index_offset;
                             k < row_coupleBlock_ptr[i_row - numMyRows + 1] - index_offset; ++k) {
                            if (row_coupleBlock_index[k] == i_col) {
                                // Entry array(k_Sol, j_Eq) is a block
                                // (row_block_size x col_block_size)
                                for (int ic = 0; ic < col_block_size; ++ic) {
                                    const int i_Sol = ic + col_block_size * l_col;
                                    for (int ir = 0; ir < row_block_size; ++ir) {
                                        const int i_Eq = ir + row_block_size * l_row;
                                        row_coupleBlock_val[k*block_size + ir + row_block_size*ic] +=
                                            array[INDEX4
                                              (i_Eq, i_Sol, k_Eq, k_Sol, numEq, numEq, NN)];
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
#endif // ESYS_HAVE_PASO

} // namespace dudley

