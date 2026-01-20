
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "Assemble.h"

#ifdef ESYS_HAVE_PASO
#include <paso/SystemMatrix.h>
#endif

#ifdef ESYS_HAVE_TRILINOS
#include <trilinoswrap/TrilinosMatrixAdapter.h>

using esys_trilinos::TrilinosMatrixAdapter;
#endif

namespace finley {

using escript::DataTypes::real_t;
using escript::DataTypes::cplx_t;

#ifdef ESYS_HAVE_PASO
static void addToSystemMatrixPasoCSC(paso::SystemMatrix<double>* S, int NN_Equa,
                                     const index_t* Nodes_Equa, int num_Equa,
                                     int NN_Sol, const index_t* Nodes_Sol,
                                     int num_Sol, const real_t* array);

template <typename T>
static void addToSystemMatrixPasoCSR(paso::SystemMatrix<T>* S, int NN_Equa,
                                    const index_t* Nodes_Equa, int num_Equa,
                                    int NN_Sol, const index_t* Nodes_Sol,
                                    int num_Sol, const T* array);
#endif

template<>
void Assemble_addToSystemMatrix<real_t>(escript::ASM_ptr S, int NN_Equa,
                                    const index_t* Nodes_Equa, int num_Equa,
                                    int NN_Sol, const index_t* Nodes_Sol,
                                    int num_Sol, const real_t* array)
{
#ifdef ESYS_HAVE_PASO
    paso::SystemMatrix<real_t>* pmat = dynamic_cast<paso::SystemMatrix<real_t>*>(S.get());
    if (pmat) {
        // call the right function depending on storage type
        if (pmat->type & MATRIX_FORMAT_CSC) {
            addToSystemMatrixPasoCSC(pmat, NN_Equa, Nodes_Equa,
                                     num_Equa, NN_Sol, Nodes_Sol,
                                     num_Sol, array);
        } else { // type == CSR
            addToSystemMatrixPasoCSR(pmat, NN_Equa, Nodes_Equa,
                                     num_Equa, NN_Sol, Nodes_Sol,
                                     num_Sol, array);
        }
        return;
    }
#endif
#ifdef ESYS_HAVE_TRILINOS
    TrilinosMatrixAdapter* tmat(dynamic_cast<TrilinosMatrixAdapter*>(S.get()));
    if (tmat) {
        IndexVector rowIdx(Nodes_Equa, Nodes_Equa+NN_Equa);
        //IndexVector colIdx(Nodes_Sol, Nodes_Sol+NN_Sol);
        std::vector<real_t> arr(array, array+(NN_Equa*NN_Sol*num_Sol*num_Equa));
        tmat->add(rowIdx, arr);
        return;
    }
#endif
    throw FinleyException("Assemble_addToSystemMatrix: unknown system "
                          "matrix type.");
}

template<>
void Assemble_addToSystemMatrix<cplx_t>(escript::ASM_ptr S, int NN_Equa,
                                    const index_t* Nodes_Equa, int num_Equa,
                                    int NN_Sol, const index_t* Nodes_Sol,
                                    int num_Sol, const cplx_t* array)
{
#ifdef ESYS_HAVE_MUMPS
    paso::SystemMatrix<cplx_t>* pmat = dynamic_cast<paso::SystemMatrix<cplx_t>*>(S.get());
    if (pmat) {
        if (pmat->type & (MATRIX_FORMAT_OFFSET1 + MATRIX_FORMAT_BLK1)) {
            addToSystemMatrixPasoCSR(pmat, NN_Equa, Nodes_Equa,
                                     num_Equa, NN_Sol, Nodes_Sol,
                                     num_Sol, array);
        } else {
            throw FinleyException("addToSystemMatrix: MUMPS requires CSR format with "
                                  "index offset 1 and block size 1.");
        }
        return;
    }
#endif
#ifdef ESYS_HAVE_TRILINOS
    TrilinosMatrixAdapter* tmat = dynamic_cast<TrilinosMatrixAdapter*>(S.get());
    if (tmat) {
        IndexVector rowIdx(Nodes_Equa, Nodes_Equa+NN_Equa);
        //IndexVector colIdx(Nodes_Sol, Nodes_Sol+NN_Sol);
        std::vector<cplx_t> arr(array, array+(NN_Equa*NN_Sol*num_Sol*num_Equa));
        tmat->add(rowIdx, arr);
        return;
    }
#endif
    throw FinleyException("addToSystemMatrix: only Trilinos matrices support "
                          "complex-valued assembly!");
}

#ifdef ESYS_HAVE_PASO
void addToSystemMatrixPasoCSC(paso::SystemMatrix<double>* in, int NN_Equa,
                              const index_t* Nodes_Equa, int num_Equa,
                              int NN_Sol, const index_t* Nodes_Sol,
                              int num_Sol, const real_t* array)
{
    const int index_offset = (in->type & MATRIX_FORMAT_OFFSET1 ? 1 : 0);
    const int row_block_size = in->row_block_size;
    const int col_block_size = in->col_block_size;
    const int block_size = in->block_size;
    const int num_subblocks_Equa = num_Equa/row_block_size;
    const int num_subblocks_Sol = num_Sol/col_block_size;
    const dim_t numMyCols = in->pattern->mainPattern->numInput;
    const dim_t numMyRows = in->pattern->mainPattern->numOutput;

    const index_t* mainBlock_ptr = in->mainBlock->pattern->ptr;
    const index_t* mainBlock_index = in->mainBlock->pattern->index;
    real_t* mainBlock_val = in->mainBlock->val;
    const index_t* col_coupleBlock_ptr = in->col_coupleBlock->pattern->ptr;
    const index_t* col_coupleBlock_index = in->col_coupleBlock->pattern->index;
    real_t* col_coupleBlock_val = in->col_coupleBlock->val;
    //const index_t* row_coupleBlock_ptr = in->row_coupleBlock->pattern->ptr;
    const index_t* row_coupleBlock_index = in->row_coupleBlock->pattern->index;
    real_t* row_coupleBlock_val = in->row_coupleBlock->val;

    for (int k_Sol = 0; k_Sol < NN_Sol; ++k_Sol) {
        // Down columns of array
        const index_t j_Sol = Nodes_Sol[k_Sol];
        for (int l_col = 0; l_col < num_subblocks_Sol; ++l_col) {
            const index_t i_col = j_Sol * num_subblocks_Sol + l_col;
            if (i_col < numMyCols) {
                for (int k_Equa = 0; k_Equa < NN_Equa; ++k_Equa) {
                    // Across cols of array
                    const index_t j_Equa = Nodes_Equa[k_Equa];
                    for (int l_row = 0; l_row < num_subblocks_Equa; ++l_row) {
                        const index_t i_row = j_Equa*num_subblocks_Equa+index_offset+l_row;
                        if (i_row < numMyRows + index_offset ) {
                            for (index_t k = mainBlock_ptr[i_col]-index_offset;
                                 k < mainBlock_ptr[i_col + 1]-index_offset; ++k) {
                                if (mainBlock_index[k] == i_row) {
                                    // Entry array(k_Equa, j_Sol) is a block
                                    // (col_block_size x col_block_size)
                                    for (int ic = 0; ic < col_block_size; ++ic) {
                                        const int i_Sol = ic + col_block_size * l_col;
                                        for (int ir = 0; ir < row_block_size; ++ir) {
                                            const int i_Eq = ir + row_block_size * l_row;
                                            mainBlock_val[k*block_size + ir + row_block_size*ic] +=
                                                array[INDEX4
                                                  (i_Eq, i_Sol, k_Equa, k_Sol, num_Equa, num_Sol, NN_Equa)];
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
                                                  (i_Eq, i_Sol, k_Equa, k_Sol, num_Equa, num_Sol, NN_Equa)];
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            } else { // i_col >= numMyCols
                for (int k_Equa = 0; k_Equa < NN_Equa; ++k_Equa) {
                    // Across rows of array
                    const index_t j_Equa = Nodes_Equa[k_Equa];
                    for (int l_row = 0; l_row < num_subblocks_Equa; ++l_row) {
                        const index_t i_row = j_Equa * num_subblocks_Equa + index_offset + l_row;
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
                                                  (i_Eq, i_Sol, k_Equa, k_Sol, num_Equa, num_Sol, NN_Equa)];
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

template <typename T>
void addToSystemMatrixPasoCSR(paso::SystemMatrix<T>* in, int NN_Equa,
                              const index_t* Nodes_Equa, int num_Equa,
                              int NN_Sol, const index_t* Nodes_Sol,
                              int num_Sol, const T* array)
{
    const int index_offset = (in->type & MATRIX_FORMAT_OFFSET1 ? 1 : 0);
    const int row_block_size = in->row_block_size;
    const int col_block_size = in->col_block_size;
    const int block_size = in->block_size;
    const int num_subblocks_Equa = num_Equa / row_block_size;
    const int num_subblocks_Sol = num_Sol / col_block_size;
    const dim_t numMyCols = in->pattern->mainPattern->numInput;
    const dim_t numMyRows = in->pattern->mainPattern->numOutput;

    const index_t* mainBlock_ptr = in->mainBlock->pattern->ptr;
    const index_t* mainBlock_index = in->mainBlock->pattern->index;
    T* mainBlock_val = in->mainBlock->val;
    const index_t* col_coupleBlock_ptr = in->col_coupleBlock->pattern->ptr;
    const index_t* col_coupleBlock_index = in->col_coupleBlock->pattern->index;
    T* col_coupleBlock_val = in->col_coupleBlock->val;
    const index_t* row_coupleBlock_ptr = in->row_coupleBlock->pattern->ptr;
    const index_t* row_coupleBlock_index = in->row_coupleBlock->pattern->index;
    T* row_coupleBlock_val = in->row_coupleBlock->val;

    for (int k_Equa = 0; k_Equa < NN_Equa; ++k_Equa) {
        // Down columns of array
        const index_t j_Equa = Nodes_Equa[k_Equa];
        for (int l_row = 0; l_row<num_subblocks_Equa; ++l_row) {
            const index_t i_row = j_Equa*num_subblocks_Equa+l_row;
            // only look at the matrix rows stored on this processor
            if (i_row < numMyRows) {
                for (int k_Sol=0; k_Sol<NN_Sol; ++k_Sol) {
                    // Across rows of array
                    const index_t j_Sol=Nodes_Sol[k_Sol];
                    for (int l_col=0; l_col<num_subblocks_Sol; ++l_col) {
                        // only look at the matrix rows stored on this processor
                        const index_t i_col = j_Sol * num_subblocks_Sol + index_offset + l_col;
                        if (i_col < numMyCols + index_offset) {
                            for (index_t k = mainBlock_ptr[i_row] - index_offset;
                                 k < mainBlock_ptr[i_row + 1] - index_offset; ++k) {
                                if (mainBlock_index[k] == i_col) {
                                    // Entry array(k_Sol, j_Equa) is a block
                                    // (row_block_size x col_block_size)
                                    for (int ic = 0; ic < col_block_size; ++ic) {
                                        const int i_Sol = ic + col_block_size * l_col;
                                        for (int ir = 0; ir < row_block_size; ++ir) {
                                            const int i_Eq = ir + row_block_size * l_row;
                                            mainBlock_val[k*block_size + ir + row_block_size*ic]+=
                                                  array[INDEX4
                                                    (i_Eq, i_Sol, k_Equa, k_Sol, num_Equa, num_Sol, NN_Equa)];
                                        }
                                    }
                                    break;
                                }
                            }
                        } else {
                            for (index_t k = col_coupleBlock_ptr[i_row] - index_offset;
                                 k < col_coupleBlock_ptr[i_row + 1] - index_offset; ++k) {
                                if (col_coupleBlock_index[k] == i_col - numMyCols) {
                                    // Entry array(k_Sol, j_Equa) is a block
                                    // (row_block_size x col_block_size)
                                    for (int ic = 0; ic < col_block_size; ++ic) {
                                        const int i_Sol = ic + col_block_size * l_col;
                                        for (int ir = 0; ir < row_block_size; ++ir) {
                                            const int i_Eq = ir+row_block_size*l_row;
                                            col_coupleBlock_val[k*block_size + ir + row_block_size*ic]+=
                                                  array[INDEX4
                                                    (i_Eq, i_Sol, k_Equa, k_Sol, num_Equa, num_Sol, NN_Equa)];
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            } else { // i_row >= numMyRows
                for (int k_Sol = 0; k_Sol < NN_Sol; ++k_Sol) {
                    // Across rows of array
                    const index_t j_Sol = Nodes_Sol[k_Sol];
                    for (int l_col = 0; l_col < num_subblocks_Sol; ++l_col) {
                        const index_t i_col = j_Sol * num_subblocks_Sol + index_offset + l_col;
                        if (i_col < numMyCols + index_offset) {
                            for (index_t k = row_coupleBlock_ptr[i_row - numMyRows] - index_offset;
                                 k < row_coupleBlock_ptr[i_row - numMyRows + 1] - index_offset; ++k) {
                                if (row_coupleBlock_index[k] == i_col) {
                                    // Entry array(k_Sol, j_Equa) is a block
                                    // (row_block_size x col_block_size)
                                    for (int ic = 0; ic < col_block_size; ++ic) {
                                        const int i_Sol = ic + col_block_size * l_col;
                                        for (int ir = 0; ir < row_block_size; ++ir) {
                                            const int i_Eq = ir + row_block_size * l_row;
                                            row_coupleBlock_val[k*block_size + ir + row_block_size*ic]+=
                                                array[INDEX4
                                                  (i_Eq, i_Sol, k_Equa, k_Sol, num_Equa, num_Sol, NN_Equa)];
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

} // namespace finley

