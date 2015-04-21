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

#include <thrust/functional.h>
#include <cusp/detail/functional.h>

#ifndef DIA_CHUNKSIZE
#define DIA_CHUNKSIZE 1024
#endif

namespace cusp
{
namespace detail
{
namespace host
{

/////////////////////////
// DIA transposed SpMV //
/////////////////////////
template <typename Matrix,
          typename Vector1,
          typename Vector2,
          typename UnaryFunction,
          typename BinaryFunction1,
          typename BinaryFunction2>
void transposed_spmv_dia(const Matrix&  A,
                         const Vector1& x,
                               Vector2& y,
                         UnaryFunction   initialize,
                         BinaryFunction1 combine,
                         BinaryFunction2 reduce)
{
    typedef typename Matrix::index_type  IndexType;
    //typedef typename Vector2::value_type ValueType;

    const size_t num_diagonals = A.values.num_cols;

#pragma omp parallel for
    for (size_t ch = 0; ch < A.num_cols; ch += DIA_CHUNKSIZE) {
        // initialize chunk
        for (size_t row = ch; row < std::min(ch+DIA_CHUNKSIZE,A.num_cols); row++)
        {
            y[row] = initialize(y[row]);
        }
        // for each diagonal
        for (size_t d = 0; d < num_diagonals; d++)
        {
            for (IndexType row=ch; row<std::min(ch+DIA_CHUNKSIZE,A.num_cols); row++)
            {
                const IndexType col = row - A.diagonal_offsets[d];
                if (col >= 0 && col < A.num_rows)
                {
                    y[row] = reduce(y[row], combine(A.values(col, d), x[col]));
                }
            }
        }
    }
}

template <typename Matrix,
          typename Vector1,
          typename Vector2>
void transposed_spmv_dia(const Matrix&  A,
                         const Vector1& x,
                               Vector2& y)
{
    typedef typename Vector2::value_type ValueType;

    transposed_spmv_dia(A, x, y,
             cusp::detail::zero_function<ValueType>(),
             thrust::multiplies<ValueType>(),
             thrust::plus<ValueType>());
}

template <typename Matrix,
          typename Vector1,
          typename Vector2,
          typename UnaryFunction,
          typename BinaryFunction1,
          typename BinaryFunction2>
void transposed_spmv_cds(const Matrix&  A,
                         const Vector1& x,
                               Vector2& y,
                         UnaryFunction   initialize,
                         BinaryFunction1 combine,
                         BinaryFunction2 reduce)
{
    typedef typename Matrix::index_type  IndexType;
    typedef typename Vector2::value_type ValueType;

    const IndexType num_diagonals = A.diagonal_offsets.size();
    const IndexType block_size = (IndexType)A.block_size;
    const IndexType num_cols = (IndexType)A.num_cols;
    // make chunksize a multiple of block_size
    const IndexType chunksize = block_size*(DIA_CHUNKSIZE/block_size);

    // optimisation for special case
    if (block_size == 2) {
#pragma omp parallel for
        for (IndexType ch = 0; ch < num_cols; ch += chunksize) {
            for (IndexType row = ch; row < std::min(ch+chunksize,num_cols); row+=2)
            {
                ValueType sum1 = initialize(y[row]);
                ValueType sum2 = initialize(y[row+1]);
                // for each diagonal block
                for (IndexType d = 0; d < num_diagonals; d++)
                {
                    const IndexType col = row - A.diagonal_offsets[d]*2;
                    if (col >= 0 && col < A.num_rows)
                    {
                        sum1 = reduce(sum1, combine(A.values(col,   2*d),  x[col]));
                        sum2 = reduce(sum2, combine(A.values(col,   2*d+1),x[col]));
                        sum1 = reduce(sum1, combine(A.values(col+1, 2*d),  x[col+1]));
                        sum2 = reduce(sum2, combine(A.values(col+1, 2*d+1),x[col+1]));
                    }
                }
                y[row] = sum1;
                y[row+1] = sum2;
            }
        }
    } else { // block_size!=2
#pragma omp parallel for
        for (IndexType ch = 0; ch < num_cols; ch += chunksize) {
            for (IndexType row = ch; row < std::min(ch+chunksize,num_cols); row++)
            {
                y[row] = initialize(y[row]);

                // for each diagonal block
                for (IndexType d = 0; d < num_diagonals; d++)
                {
                    const IndexType k = A.diagonal_offsets[d]*block_size;
                    const IndexType col = block_size*(row/block_size) - k;
                    if (col >= 0 && col <= A.num_rows-block_size)
                    {
                        // for each column in block
                        for (IndexType i = 0; i < block_size; i++)
                        {
                            const ValueType& Aij = A.values(col+i, d*block_size+row%block_size);
                            const ValueType& xj = x[col + i];
                            y[row] = reduce(y[row], combine(Aij, xj));
                        }
                    }
                } // diagonals
            } // rows
        } // chunks
    } // block_size
}

template <typename Matrix,
          typename Vector1,
          typename Vector2>
void transposed_spmv_cds(const Matrix&  A,
                         const Vector1& x,
                               Vector2& y)
{
    typedef typename Vector2::value_type ValueType;

    if (A.block_size == 1) {
        transposed_spmv_dia(A, x, y,
                 cusp::detail::zero_function<ValueType>(),
                 thrust::multiplies<ValueType>(),
                 thrust::plus<ValueType>());
    } else {
        transposed_spmv_cds(A, x, y,
                  cusp::detail::zero_function<ValueType>(),
                  thrust::multiplies<ValueType>(),
                  thrust::plus<ValueType>());
    }
}

} // end namespace host
} // end namespace detail
} // end namespace cusp


