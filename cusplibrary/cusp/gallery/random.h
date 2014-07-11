/*
 *  Copyright 2008-2009 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

/*! \file random.h
 *  \brief Random matrix generators
 */

#pragma once

#include <cusp/detail/config.h>

#include <cusp/cds_matrix.h>
#include <cusp/coo_matrix.h>
#include <thrust/unique.h>
#include <thrust/sort.h>

#include <stdlib.h> // XXX remove when we switch RNGs

namespace cusp
{
namespace gallery
{
/*! \addtogroup gallery Matrix Gallery
 *  \ingroup gallery
 *  \{
 */

// TODO use thrust RNGs, add seed parameter defaulting to num_rows ^ num_cols ^ num_samples
// TODO document
template <class MatrixType>
void random(size_t num_rows, size_t num_cols, size_t num_samples, MatrixType& output)
{
    typedef typename MatrixType::index_type IndexType;
    typedef typename MatrixType::value_type ValueType;

    cusp::coo_matrix<IndexType,ValueType,cusp::host_memory> coo(num_rows, num_cols, num_samples);

    srand(num_rows ^ num_cols ^ num_samples);

    for(size_t n = 0; n < num_samples; n++)
    {
        coo.row_indices[n]    = rand() % num_rows;
        coo.column_indices[n] = rand() % num_cols;
        coo.values[n]         = ValueType(1);
    }

    // sort indices by (row,column)
    coo.sort_by_row_and_column();

    size_t num_entries = thrust::unique(thrust::make_zip_iterator(thrust::make_tuple(coo.row_indices.begin(), coo.column_indices.begin())),
                                        thrust::make_zip_iterator(thrust::make_tuple(coo.row_indices.end(),   coo.column_indices.end())))
                         - thrust::make_zip_iterator(thrust::make_tuple(coo.row_indices.begin(), coo.column_indices.begin()));

    coo.resize(num_rows, num_cols, num_entries);
    
    output = coo;
}

template <class MatrixType>
void randomblock(size_t num_rows, size_t num_diagonals, size_t block_size, MatrixType& output)
{
    typedef typename MatrixType::index_type IndexType;
    typedef typename MatrixType::value_type ValueType;

    if (num_rows % block_size != 0)
        throw cusp::invalid_input_exception("number of rows must be a multiple of block size!");

    cusp::cds_matrix<IndexType,ValueType,cusp::host_memory> cds(num_rows, 0, num_diagonals, block_size);

    srand(num_rows ^ num_diagonals);

    // instead of entirely random, let's try and even out the number of
    // subdiagonals and superdiagonals
    const size_t max_offset = num_rows/32 - 1;
    for(size_t n = 0; n < num_diagonals/2; n++)
    {
        const int offset = 1 + rand() % (max_offset-1);
        cds.diagonal_offsets[2*n] = -offset;
        cds.diagonal_offsets[2*n+1] = offset;
    }
    // for odd number of diagonals add main diagonal
    if (num_diagonals%2 == 1)
        cds.diagonal_offsets[num_diagonals-1]=0;

    std::sort(cds.diagonal_offsets.begin(),cds.diagonal_offsets.end());

    size_t num_entries = 0;

    for(size_t n = 0; n < num_diagonals; n++)
    {
        const int offset = cds.diagonal_offsets[n];
        const size_t num_blocks = num_rows/block_size-std::abs(offset);
        num_entries += num_blocks*block_size*block_size;
        const size_t first = std::max(0, -offset*(int)block_size);
        for(size_t block = 0; block < num_blocks; block++)
        {
            for(size_t i = 0; i < block_size; i++)
            {
                for(size_t j = 0; j < block_size; j++)
                {
                    cds.values(first+block*block_size+i, n*block_size+j) = ValueType(1);
                }
            }
        }
    }

    cds.resize(num_rows, num_entries, block_size, num_diagonals);
    
    output = cds;
}
/*! \}
 */

} // end namespace gallery
} // end namespace cusp

