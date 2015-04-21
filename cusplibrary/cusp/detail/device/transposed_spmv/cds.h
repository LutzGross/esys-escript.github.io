/*
 *  Copyright 2014-2015 The University of Queensland
 *  http://www.uq.edu.au
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

#pragma once

#include <thrust/extrema.h>

#include <cusp/detail/device/arch.h>
#include <cusp/detail/device/common.h>
#include <cusp/detail/device/utils.h>
#include <cusp/detail/device/texture.h>

#include <thrust/device_ptr.h>

namespace cusp
{
namespace detail
{
namespace device
{

////////////////////////////////////////////////////////////////////////
// CDS transposed SpMV kernels
///////////////////////////////////////////////////////////////////////
//
// transposed_spmv_cds
//   Each thread computes y[i] += A[:,i] * x
//   (the dot product of the i-th column of A with the x vector, i.e. A^T*x)
//
// transposed_spmv_cds_tex
//   Same as transposed_spmv_cds, except x is accessed via texture cache.
//


template <typename IndexType, typename ValueType, unsigned int BLOCK_SIZE, bool UseCache>
__launch_bounds__(BLOCK_SIZE,1)
__global__ void
transposed_spmv_cds_kernel(const IndexType num_rows,
                           const IndexType num_cols,
                           const IndexType num_diagonals,
                           const IndexType block_size,
                           const IndexType pitch,
                           const IndexType* diagonal_offsets,
                           const ValueType* values,
                           const ValueType* x,
                                 ValueType* y)
{
    __shared__ IndexType offsets[BLOCK_SIZE];

    const IndexType thread_id = BLOCK_SIZE * blockIdx.x + threadIdx.x;
    const IndexType grid_size = BLOCK_SIZE * gridDim.x;
    const IndexType total_diags = num_diagonals * block_size;

    for (IndexType base = 0; base < total_diags; base += BLOCK_SIZE)
    {
        // read a chunk of the diagonal offsets into shared memory
        const IndexType chunk_size = thrust::min(IndexType(BLOCK_SIZE), total_diags - base);

        if (threadIdx.x < chunk_size) {
            offsets[threadIdx.x] = -diagonal_offsets[(base + threadIdx.x)/block_size] * block_size + threadIdx.x%block_size;
        }

        __syncthreads();

        // process chunk
        for (IndexType row = thread_id; row < num_cols; row += grid_size)
        {
            ValueType sum = (base == 0) ? ValueType(0) : y[row];
            const IndexType rowbase = block_size*(row/block_size);

            // loop over all diagonals
            for (IndexType n = 0; n < chunk_size; n++)
            {
                const IndexType col = rowbase + offsets[n];

                if (col >= 0 && col < num_rows)
                {
                    // index into values array
                    const IndexType diagbase = block_size*(n/block_size)
                                               + threadIdx.x%block_size;
                    IndexType idx = pitch * (base + diagbase) + col;

                    const ValueType A_ij = values[idx];
                    sum += A_ij * fetch_x<UseCache>(col, x);
                }
            }

            y[row] = sum;
        }

        // wait until all threads are done reading offsets
        __syncthreads();
    }
}


template <bool UseCache,
          typename Matrix,
          typename ValueType>
void __transposed_spmv_cds(const Matrix&    A,
                           const ValueType* x,
                                 ValueType* y)
{
    typedef typename Matrix::index_type IndexType;

    const size_t BLOCK_SIZE = 256;
    const size_t MAX_BLOCKS = cusp::detail::device::arch::max_active_blocks(
            transposed_spmv_cds_kernel<IndexType, ValueType, BLOCK_SIZE, UseCache>, BLOCK_SIZE, (size_t) sizeof(IndexType) * BLOCK_SIZE);
    const size_t NUM_BLOCKS = std::min<size_t>(MAX_BLOCKS, DIVIDE_INTO(A.num_cols, BLOCK_SIZE));

    const IndexType num_diagonals = A.values.num_cols / A.block_size;
    const IndexType pitch         = A.values.pitch;

    if (UseCache)
        bind_x(x);

    transposed_spmv_cds_kernel<IndexType, ValueType, BLOCK_SIZE, UseCache> <<<NUM_BLOCKS, BLOCK_SIZE>>>
        (A.num_rows, A.num_cols, num_diagonals, A.block_size, pitch,
         thrust::raw_pointer_cast(&A.diagonal_offsets[0]),
         thrust::raw_pointer_cast(&A.values.values[0]),
         x, y);

    if (UseCache)
        unbind_x(x);
}

template <typename Matrix,
          typename ValueType>
void transposed_spmv_cds(const Matrix&    A,
                         const ValueType* x,
                               ValueType* y)
{
    __transposed_spmv_cds<false>(A, x, y);
}

template <typename Matrix,
          typename ValueType>
void transposed_spmv_cds_tex(const Matrix&    A,
                             const ValueType* x,
                                   ValueType* y)
{
    __transposed_spmv_cds<true>(A, x, y);
}

} // end namespace device
} // end namespace detail
} // end namespace cusp

