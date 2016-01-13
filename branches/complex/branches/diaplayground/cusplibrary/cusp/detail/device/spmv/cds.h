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
// CDS SpMV kernels 
///////////////////////////////////////////////////////////////////////
//
// Diagonal matrices arise in grid-based discretizations using stencils.  
// For instance, the standard 5-point discretization of the two-dimensional 
// Laplacian operator has the stencil:
//      [  0  -1   0 ]
//      [ -1   4  -1 ]
//      [  0  -1   0 ]
// and the resulting CDS format has 5 diagonals.
//
// spmv_cds
//   Each thread computes y[i] += A[i,:] * x 
//   (the dot product of the i-th row of A with the x vector)
//
// spmv_cds_tex
//   Same as spmv_cds, except x is accessed via texture cache.
//


template <typename IndexType, typename ValueType, unsigned int BLOCK_SIZE, bool UseCache>
__launch_bounds__(BLOCK_SIZE,1)
__global__ void
spmv_cds_kernel(const IndexType num_rows, 
                const IndexType num_diagonals,
                const IndexType block_size,
                const IndexType pitch,
                const IndexType * diagonal_offsets,
                const ValueType * values,
                const ValueType * x, 
                      ValueType * y)
{
    __shared__ IndexType offsets[BLOCK_SIZE];

    const IndexType thread_id = BLOCK_SIZE * blockIdx.x + threadIdx.x;
    const IndexType grid_size = BLOCK_SIZE * gridDim.x;
    //const float bsf = (float)block_size;
    const IndexType num_cols = num_diagonals * block_size;

    // for BLOCK_SIZE=256 this is most likely only executed once
    for(IndexType base = 0; base < num_cols; base += BLOCK_SIZE)
    {
        // read a chunk of the diagonal offsets into shared memory
        const IndexType chunk_size = thrust::min(IndexType(BLOCK_SIZE), num_cols - base);

        if(threadIdx.x < chunk_size)
            offsets[threadIdx.x] = diagonal_offsets[(base + threadIdx.x)/block_size] * block_size + threadIdx.x%block_size;

        __syncthreads();

        // process chunk
        for(IndexType row = thread_id; row < num_rows; row += grid_size)
        {
            ValueType sum = (base == 0) ? ValueType(0) : y[row];

            // index into values array
            IndexType idx = row + pitch * base;
            // for sm_10 this is faster than integer division
            // block_size*(row/block_size)
            //const IndexType colbase = block_size*(int)((row/bsf)+0.001f);
            const IndexType colbase = block_size*(row/block_size);

            for(IndexType n = 0; n < chunk_size; n++)
            {
                const IndexType col = colbase + offsets[n];

                if(col >= 0 && col < num_rows)
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

    
template <bool UseCache,
          typename Matrix,
          typename ValueType>
void __spmv_cds(const Matrix&    A,
                const ValueType* x, 
                      ValueType* y)
{
    typedef typename Matrix::index_type IndexType;

    const size_t BLOCK_SIZE = 256;
    const size_t MAX_BLOCKS = cusp::detail::device::arch::max_active_blocks(spmv_cds_kernel<IndexType, ValueType, BLOCK_SIZE, UseCache>, BLOCK_SIZE, (size_t) sizeof(IndexType) * BLOCK_SIZE);
    const size_t NUM_BLOCKS = std::min<size_t>(MAX_BLOCKS, DIVIDE_INTO(A.num_rows, BLOCK_SIZE));
   
    const IndexType num_diagonals = A.values.num_cols / A.block_size;
    const IndexType pitch         = A.values.pitch;

    if (UseCache)
        bind_x(x);
  
    spmv_cds_kernel<IndexType, ValueType, BLOCK_SIZE, UseCache> <<<NUM_BLOCKS, BLOCK_SIZE>>>
        (A.num_rows, num_diagonals, A.block_size, pitch,
         thrust::raw_pointer_cast(&A.diagonal_offsets[0]),
         thrust::raw_pointer_cast(&A.values.values[0]),
         x, y);

    if (UseCache)
        unbind_x(x);
}

template <typename Matrix,
          typename ValueType>
void spmv_cds(const Matrix&    A, 
              const ValueType* x, 
                    ValueType* y)
{
    __spmv_cds<false>(A, x, y);
}

template <typename Matrix,
          typename ValueType>
void spmv_cds_tex(const Matrix&    A,
                  const ValueType* x, 
                        ValueType* y)
{
    __spmv_cds<true>(A, x, y);
}

} // end namespace device
} // end namespace detail
} // end namespace cusp

