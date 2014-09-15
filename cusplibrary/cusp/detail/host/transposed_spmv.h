
/*****************************************************************************
*
* Copyright (c) 2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

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
    } else {
#pragma omp parallel for
        for (IndexType ch = 0; ch < num_cols; ch += chunksize) {
            for (IndexType row = ch; row < std::min(ch+chunksize,num_cols); row++)
            {
                y[row] = initialize(y[row]);
            }

            // for each diagonal block
            for (IndexType d = 0; d < num_diagonals; d++)
            {
                const IndexType k = A.diagonal_offsets[d]*block_size;
                for (IndexType row=ch; row<std::min(ch+chunksize,num_cols); row+=block_size)
                {
                    const IndexType col = row - k;
                    if (col >= 0 && col < A.num_rows-block_size)
                    {
                        //y[row] = reduce(y[row], combine(A.values(col, d), x[col]));
                    }
                }
            }
        }
    }
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


