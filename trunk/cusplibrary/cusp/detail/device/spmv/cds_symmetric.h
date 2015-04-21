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
// CDS symmetric SpMV kernels
///////////////////////////////////////////////////////////////////////

// generic version

template <typename IndexType, typename ValueType, unsigned int BLOCK_SIZE, bool UseCache>
__launch_bounds__(BLOCK_SIZE,1)
__global__ void
spmv_cds_symmetric_kernel(const IndexType  num_rows,
                          const IndexType  num_diagonals,
                          const IndexType  block_size,
                          const IndexType  pitch,
                          const IndexType* diagonal_offsets,
                          const ValueType* values,
                          const ValueType* x,
                                ValueType* y)
{
    __shared__ IndexType offsets[BLOCK_SIZE];
    __shared__ IndexType offsets2[BLOCK_SIZE];

    const IndexType thread_id = BLOCK_SIZE * blockIdx.x + threadIdx.x;
    const IndexType grid_size = BLOCK_SIZE * gridDim.x;
    //const float bsf = (float)block_size;
    const IndexType num_cols = num_diagonals * block_size;
    const IndexType t_mod_bs = threadIdx.x % block_size;

    // for BLOCK_SIZE=256 this is most likely only executed once
    for (IndexType base = 0; base < num_cols; base += BLOCK_SIZE)
    {
        // read a chunk of the diagonal offsets into shared memory
        const IndexType chunk_size = thrust::min(IndexType(BLOCK_SIZE), num_cols - base);

        if (threadIdx.x < chunk_size) {
            offsets[threadIdx.x] = diagonal_offsets[(base + threadIdx.x)/block_size] * block_size + t_mod_bs;
            offsets2[threadIdx.x] = -diagonal_offsets[(base + threadIdx.x)/block_size] * block_size + t_mod_bs;
        }

        __syncthreads();

        // process chunk
        for (IndexType row = thread_id; row < num_rows; row += grid_size)
        {
            ValueType sum = (base == 0) ? ValueType(0) : y[row];

            // for sm_10 this is faster than integer division
            // block_size*(row/block_size)
            //const IndexType colbase = block_size*(int)((row/bsf)+0.001f);
            const IndexType colbase = block_size*(row/block_size);

            // process subdiagonal blocks
            for (IndexType n = block_size; n < chunk_size; n++)
            {
                const IndexType col = colbase + offsets2[n];

                if (col >= 0 && col < num_rows)
                {
                    const IndexType diagbase = block_size*(n/block_size) + t_mod_bs;
                    const ValueType& A_ij = values[col + pitch*(base+diagbase)];
                    sum += A_ij * fetch_x<UseCache>(col, x);
                }
            }

            // index into values array
            IndexType idx = row + pitch * base;

            // process main and upper diagonal blocks
            for (IndexType n = 0; n < chunk_size; n++)
            {
                const IndexType col = colbase + offsets[n];

                if (col >= 0 && col < num_rows)
                {
                    const ValueType& A_ij = values[idx];
                    sum += A_ij * fetch_x<UseCache>(col, x);
                }
                idx += pitch;
            }

            y[row] = sum;
        }

        // wait until all threads are done reading offsets
        __syncthreads();
    }
}

//
// diagonal block size 2
//
template <typename IndexType, typename ValueType, unsigned int BLOCK_SIZE, bool UseCache>
__launch_bounds__(BLOCK_SIZE,1)
__global__ void
spmv_cds_symmetric_bs2_kernel(const IndexType  num_rows,
                              const IndexType  num_diagonals,
                              const IndexType  pitch,
                              const IndexType* diagonal_offsets,
                              const ValueType* values,
                              const ValueType* x,
                                    ValueType* y)
{
    __shared__ IndexType offsets[BLOCK_SIZE];
    __shared__ IndexType offsets2[BLOCK_SIZE];

    const IndexType thread_id = BLOCK_SIZE * blockIdx.x + threadIdx.x;
    const IndexType grid_size = BLOCK_SIZE * gridDim.x;
    const IndexType t_mod_bs = threadIdx.x % 2;

    // for BLOCK_SIZE=256 this is most likely only executed once
    for (IndexType base = 0; base < num_diagonals; base += BLOCK_SIZE)
    {
        // read a chunk of the diagonal offsets into shared memory
        const IndexType chunk_size = thrust::min(IndexType(BLOCK_SIZE), num_diagonals - base);

        if (threadIdx.x < chunk_size) {
            offsets[threadIdx.x] = 2*diagonal_offsets[base + threadIdx.x];
            offsets2[threadIdx.x] = -2*diagonal_offsets[num_diagonals - base - 1 - threadIdx.x];
        }

        __syncthreads();

        // process chunk
        for (IndexType row = thread_id; row < num_rows; row += grid_size)
        {
            const IndexType row_step = 2*(row/2);
            ValueType sum1 = (base == 0) ? ValueType(0) : y[row];

            // index into values array
            IndexType idx = 2 * pitch * (num_diagonals - base - 1) + t_mod_bs*pitch;

            // process subdiagonal blocks
            for (IndexType n = 0; n < chunk_size-1; n++)
            {
                const IndexType col = row_step + offsets2[n];

                if (col >= 0 && col <= num_rows-2)
                {
                    sum1 += values[col   + idx] * fetch_x<UseCache>(col, x);
                    sum1 += values[col+1 + idx] * fetch_x<UseCache>(col+1, x);
                }
                idx -= 2*pitch;
            }

            // index into values array
            idx = row + 2 * pitch * base;

            // process main and upper diagonal blocks
            for (IndexType n = 0; n < chunk_size; n++)
            {
                const IndexType col = row_step + offsets[n];

                if (col >= 0 && col <= num_rows-2)
                {
                    sum1 += values[idx] * fetch_x<UseCache>(col, x);
                    sum1 += values[idx+pitch] * fetch_x<UseCache>(col+1, x);
                }
                idx += 2*pitch;
            }

            y[row] = sum1;
        }

        // wait until all threads are done reading offsets
        __syncthreads();
    }
}

//
// diagonal block size 3
//
template <typename IndexType, typename ValueType, unsigned int BLOCK_SIZE, bool UseCache>
__launch_bounds__(BLOCK_SIZE,1)
__global__ void
spmv_cds_symmetric_bs3_kernel(const IndexType  num_rows,
                              const IndexType  num_diagonals,
                              const IndexType  pitch,
                              const IndexType* diagonal_offsets,
                              const ValueType* values,
                              const ValueType* x,
                                    ValueType* y)
{
    __shared__ IndexType offsets[BLOCK_SIZE];
    __shared__ IndexType offsets2[BLOCK_SIZE];

    const IndexType thread_id = BLOCK_SIZE * blockIdx.x + threadIdx.x;
    const IndexType grid_size = BLOCK_SIZE * gridDim.x;
    const IndexType t_mod_bs = threadIdx.x % 3;

    // for BLOCK_SIZE=256 this is most likely only executed once
    for (IndexType base = 0; base < num_diagonals; base += BLOCK_SIZE)
    {
        // read a chunk of the diagonal offsets into shared memory
        const IndexType chunk_size = thrust::min(IndexType(BLOCK_SIZE), num_diagonals - base);

        if (threadIdx.x < chunk_size) {
            offsets[threadIdx.x] = 3*diagonal_offsets[base + threadIdx.x];
            offsets2[threadIdx.x] = -3*diagonal_offsets[num_diagonals - base - 1 - threadIdx.x];
        }

        __syncthreads();

        // process chunk
        for (IndexType row = thread_id; row < num_rows; row += grid_size)
        {
            const IndexType row_step = 3*(row/3);
            ValueType sum1 = (base == 0) ? ValueType(0) : y[row];

            // index into values array
            IndexType idx = 3 * pitch * (num_diagonals - base - 1) + t_mod_bs*pitch;

            // process subdiagonal blocks
            for (IndexType n = 0; n < chunk_size-1; n++)
            {
                const IndexType col = row_step + offsets2[n];

                if (col >= 0 && col <= num_rows-3)
                {
                    sum1 += values[col   + idx] * fetch_x<UseCache>(col, x);
                    sum1 += values[col+1 + idx] * fetch_x<UseCache>(col+1, x);
                    sum1 += values[col+2 + idx] * fetch_x<UseCache>(col+2, x);
                }
                idx -= 3*pitch;
            }

            // index into values array
            idx = row + 3 * pitch * base;

            // process main and upper diagonal blocks
            for (IndexType n = 0; n < chunk_size; n++)
            {
                const IndexType col = row_step + offsets[n];

                if (col >= 0 && col <= num_rows-3)
                {
                    sum1 += values[idx] * fetch_x<UseCache>(col, x);
                    sum1 += values[idx+pitch] * fetch_x<UseCache>(col+1, x);
                    sum1 += values[idx+2*pitch] * fetch_x<UseCache>(col+2, x);
                }
                idx += 3*pitch;
            }

            y[row] = sum1;
        }

        // wait until all threads are done reading offsets
        __syncthreads();
    }
}

template <bool UseCache,
          typename Matrix,
          typename ValueType>
void __spmv_cds_symmetric(const Matrix&    A,
                          const ValueType* x,
                                ValueType* y)
{
    using cusp::detail::device::arch::max_active_blocks;
    typedef typename Matrix::index_type IndexType;

    const size_t BLOCK_SIZE = 192;
    const IndexType num_diagonals = A.values.num_cols / A.block_size;
    const IndexType pitch         = A.values.pitch;

    if (UseCache)
        bind_x(x);

    switch (A.block_size) {
        case 2: {
            const size_t MAX_BLOCKS = max_active_blocks(
                    spmv_cds_symmetric_bs2_kernel<IndexType, ValueType,
                    BLOCK_SIZE, UseCache>, BLOCK_SIZE,
                    (size_t) sizeof(IndexType) * BLOCK_SIZE);
            const size_t NUM_BLOCKS = std::min<size_t>(MAX_BLOCKS,
                    DIVIDE_INTO(A.num_rows, BLOCK_SIZE));
            spmv_cds_symmetric_bs2_kernel<IndexType, ValueType, BLOCK_SIZE, UseCache>
                <<<NUM_BLOCKS, BLOCK_SIZE>>>(A.num_rows, num_diagonals, pitch,
                    thrust::raw_pointer_cast(&A.diagonal_offsets[0]),
                    thrust::raw_pointer_cast(&A.values.values[0]), x, y);
        }
        break;

        case 3: {
            const size_t MAX_BLOCKS = max_active_blocks(
                    spmv_cds_symmetric_bs3_kernel<IndexType, ValueType,
                    BLOCK_SIZE, UseCache>, BLOCK_SIZE,
                    (size_t) sizeof(IndexType) * BLOCK_SIZE);
            const size_t NUM_BLOCKS = std::min<size_t>(MAX_BLOCKS,
                    DIVIDE_INTO(A.num_rows, BLOCK_SIZE));
            spmv_cds_symmetric_bs3_kernel<IndexType, ValueType, BLOCK_SIZE, UseCache>
                <<<NUM_BLOCKS, BLOCK_SIZE>>>(A.num_rows, num_diagonals, pitch,
                    thrust::raw_pointer_cast(&A.diagonal_offsets[0]),
                    thrust::raw_pointer_cast(&A.values.values[0]), x, y);
        }
        break;

        default: {
            const size_t MAX_BLOCKS = max_active_blocks(
                    spmv_cds_symmetric_kernel<IndexType, ValueType, BLOCK_SIZE,
                    UseCache>, BLOCK_SIZE,
                    (size_t) sizeof(IndexType) * BLOCK_SIZE);
            const size_t NUM_BLOCKS = std::min<size_t>(MAX_BLOCKS,
                    DIVIDE_INTO(A.num_rows, BLOCK_SIZE));
            spmv_cds_symmetric_kernel<IndexType, ValueType, BLOCK_SIZE, UseCache>
                <<<NUM_BLOCKS, BLOCK_SIZE>>>(A.num_rows, num_diagonals,
                    A.block_size, pitch,
                    thrust::raw_pointer_cast(&A.diagonal_offsets[0]),
                    thrust::raw_pointer_cast(&A.values.values[0]), x, y);
        }
    }

    if (UseCache)
        unbind_x(x);
}

} // end namespace device
} // end namespace detail
} // end namespace cusp

