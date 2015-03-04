#pragma once

#include <thrust/functional.h>
#include <cusp/detail/functional.h>

#ifndef DIA_CHUNKSIZE
#define DIA_CHUNKSIZE 1024
#endif

//MW: add some OpenMP pragmas
namespace cusp
{
namespace detail
{
namespace host
{

//////////////
// COO SpMV //
//////////////
template <typename Matrix,
          typename Vector1,
          typename Vector2,
          typename UnaryFunction,
          typename BinaryFunction1,
          typename BinaryFunction2>
void spmv_coo(const Matrix&  A,
              const Vector1& x,
                    Vector2& y,
              UnaryFunction   initialize,
              BinaryFunction1 combine,
              BinaryFunction2 reduce)
{
    typedef typename Matrix::index_type  IndexType;
    typedef typename Vector2::value_type ValueType;

    for(size_t i = 0; i < A.num_rows; i++)
        y[i] = initialize(y[i]);

    for(size_t n = 0; n < A.num_entries; n++)
    {
        const IndexType& i   = A.row_indices[n];
        const IndexType& j   = A.column_indices[n];
        const ValueType& Aij = A.values[n];
        const ValueType& xj  = x[j];

        y[i] = reduce(y[i], combine(Aij, xj));
    }
}

template <typename Matrix,
          typename Vector1,
          typename Vector2>
void spmv_coo(const Matrix&  A,
              const Vector1& x,
                    Vector2& y)
{
    typedef typename Vector2::value_type ValueType;

    spmv_coo(A, x, y,
             cusp::detail::zero_function<ValueType>(),
             thrust::multiplies<ValueType>(),
             thrust::plus<ValueType>());
}


//////////////
// CSR SpMV //
//////////////
template <typename Matrix,
          typename Vector1,
          typename Vector2,
          typename UnaryFunction,
          typename BinaryFunction1,
          typename BinaryFunction2>
void spmv_csr(const Matrix&  A,
              const Vector1& x,
                    Vector2& y,
              UnaryFunction   initialize,
              BinaryFunction1 combine,
              BinaryFunction2 reduce)
{
    typedef typename Matrix::index_type  IndexType;
    typedef typename Vector2::value_type ValueType;

#pragma omp parallel for 
    for(size_t i = 0; i < A.num_rows; i++)
    {
        const IndexType& row_start = A.row_offsets[i];
        const IndexType& row_end   = A.row_offsets[i+1];
 
        ValueType accumulator = initialize(y[i]);
 
        for (IndexType jj = row_start; jj < row_end; jj++)
        {
            const IndexType& j   = A.column_indices[jj];
            const ValueType& Aij = A.values[jj];
            const ValueType& xj  = x[j];
 
            accumulator = reduce(accumulator, combine(Aij, xj));
        }
 
        y[i] = accumulator;
    }
}


template <typename Matrix,
          typename Vector1,
          typename Vector2>
void spmv_csr(const Matrix&  A,
              const Vector1& x,
                    Vector2& y)
{
    typedef typename Vector2::value_type ValueType;

    spmv_csr(A, x, y,
             cusp::detail::zero_function<ValueType>(),
             thrust::multiplies<ValueType>(),
             thrust::plus<ValueType>());
}


//////////////
// DIA SpMV //
//////////////
template <typename Matrix,
          typename Vector1,
          typename Vector2,
          typename UnaryFunction,
          typename BinaryFunction1,
          typename BinaryFunction2>
void spmv_dia(const Matrix&  A,
              const Vector1& x,
                    Vector2& y,
              UnaryFunction   initialize,
              BinaryFunction1 combine,
              BinaryFunction2 reduce)
{
    typedef typename Matrix::index_type  IndexType;
    //typedef typename Vector2::value_type ValueType;

    const size_t num_diagonals = A.values.num_cols;

    if (A.symmetric) {
        // if matrix has a main diagonal it is the first in offsets and should
        // be skipped in the subdiagonal loop below. The main diagonal is
        // processed by the second loop
        const size_t d0 = (A.diagonal_offsets[0] == 0 ? 1 : 0);
#pragma omp parallel for
        for (size_t ch = 0; ch < A.num_rows; ch += DIA_CHUNKSIZE) {
            // initialize chunk
            for (size_t row = ch; row < std::min(ch+DIA_CHUNKSIZE,A.num_rows); row++)
            {
                y[row] = initialize(y[row]);
            }

            // process subdiagonals
            for (size_t d = 0; d < num_diagonals-d0; d++)
            {
                const size_t diag = num_diagonals-d-1;
                for (size_t row = ch; row < std::min(ch+DIA_CHUNKSIZE,A.num_rows); row++)
                {
                    const IndexType col = row - A.diagonal_offsets[diag];
                    if (col >= 0 && col < A.num_rows)
                    {
                        y[row] = reduce(y[row], combine(A.values(col, diag), x[col]));
                    }
                }
            }
            // process main and upper diagonals
            for (size_t d = 0; d < num_diagonals; d++)
            {
                for (size_t row = ch; row < std::min(ch+DIA_CHUNKSIZE,A.num_rows); row++)
                {
                    const IndexType col = row + A.diagonal_offsets[d];
                    if (col >= 0 && col < A.num_cols)
                    {
                        y[row] = reduce(y[row], combine(A.values(row, d), x[col]));
                    }
                }
            }
        }
    } else { // !A.symmetric
#pragma omp parallel for
        for (size_t ch = 0; ch < A.num_rows; ch += DIA_CHUNKSIZE) {
            // initialize chunk
            for (size_t row = ch; row < std::min(ch+DIA_CHUNKSIZE,A.num_rows); row++)
            {
                y[row] = initialize(y[row]);
            }
            // for each diagonal
            for (size_t d = 0; d < num_diagonals; d++)
            {
                for (IndexType row=ch; row<std::min(ch+DIA_CHUNKSIZE,A.num_rows); row++)
                {
                    const IndexType col = row + A.diagonal_offsets[d];
                    if (col >= 0 && col < A.num_cols)
                    {
                        y[row] = reduce(y[row], combine(A.values(row, d), x[col]));
                    }
                }
            }
        }
    }
}

template <typename Matrix,
          typename Vector1,
          typename Vector2>
void spmv_dia(const Matrix&  A,
              const Vector1& x,
                    Vector2& y)
{
    typedef typename Vector2::value_type ValueType;

    spmv_dia(A, x, y,
             cusp::detail::zero_function<ValueType>(),
             thrust::multiplies<ValueType>(),
             thrust::plus<ValueType>());
}


//////////////
// CDS SpMV //
//////////////
template <typename Matrix,
          typename Vector1,
          typename Vector2,
          typename UnaryFunction,
          typename BinaryFunction1,
          typename BinaryFunction2>
void spmv_cds(const Matrix&  A,
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
    const IndexType num_rows = (IndexType)A.num_rows;
    // make chunksize a multiple of block_size
    const IndexType chunksize = block_size*(DIA_CHUNKSIZE/block_size);

    // optimization for special case
    if (block_size == 2) {
        if (A.symmetric) {
            // if there is a main diagonal block, it is the first in offsets
            // and should be skipped in the first loop below since the main
            // diagonal is processed in the second loop
            const IndexType d0 = (A.diagonal_offsets[0] == 0 ? 1 : 0);
#pragma omp parallel for
            for (IndexType ch = 0; ch < num_rows; ch+=chunksize)
            {
                for (IndexType row = ch; row<std::min(ch+chunksize,num_rows); row++)
                {
                    y[row] = initialize(y[row]);
                }

                // process subdiagonal blocks
                for (IndexType d = 0; d < num_diagonals-d0; d++)
                {
                    const IndexType diag = num_diagonals-d-1;
                    const IndexType k = -2*A.diagonal_offsets[diag];
                    for (IndexType row = ch; row<std::min(ch+chunksize,num_rows); row+=2)
                    {
                        const IndexType col = row + k;
                        if (col >= 0 && col <= num_rows-2)
                        {
                            y[row]   = reduce(y[row],  combine(A.values(col,  2*diag),  x[col]));
                            y[row]   = reduce(y[row],  combine(A.values(col+1,2*diag),  x[col+1]));
                            y[row+1] = reduce(y[row+1],combine(A.values(col,  2*diag+1),x[col]));
                            y[row+1] = reduce(y[row+1],combine(A.values(col+1,2*diag+1),x[col+1]));
                        }
                    }
                }
                // process main and upper diagonal blocks
                for (IndexType d = 0; d < num_diagonals; d++)
                {
                    const IndexType k = 2*A.diagonal_offsets[d];
                    for (IndexType row = ch; row<std::min(ch+chunksize,num_rows); row+=2)
                    {
                        const IndexType col = row + k;
                        if (col >= 0 && col <= num_rows-2)
                        {
                            y[row]   = reduce(y[row],  combine(A.values(row,  2*d),  x[col]));
                            y[row+1] = reduce(y[row+1],combine(A.values(row+1,2*d),  x[col]));
                            y[row]   = reduce(y[row],  combine(A.values(row,  2*d+1),x[col+1]));
                            y[row+1] = reduce(y[row+1],combine(A.values(row+1,2*d+1),x[col+1]));
                        }
                    }
                }
            }
        } else { // !A.symmetric
#pragma omp parallel for
            for (IndexType ch = 0; ch < num_rows; ch+=chunksize)
            {
                for (IndexType row = ch; row<std::min(ch+chunksize,num_rows); row+=2)
                {
                    ValueType sum1 = initialize(y[row]);
                    ValueType sum2 = initialize(y[row+1]);
                    // for each diagonal block
                    for (IndexType d = 0; d < num_diagonals; d++)
                    {
                        const IndexType col = row + A.diagonal_offsets[d]*2;
                        if (col >= 0 && col <= num_rows-2)
                        {
                            sum1 = reduce(sum1,combine(A.values(row,  2*d),  x[col]));
                            sum2 = reduce(sum2,combine(A.values(row+1,2*d),  x[col]));
                            sum1 = reduce(sum1,combine(A.values(row,  2*d+1),x[col+1]));
                            sum2 = reduce(sum2,combine(A.values(row+1,2*d+1),x[col+1]));
                        }
                    }
                    y[row] = sum1;
                    y[row+1] = sum2;
                }
            }
        } // A.symmetric
    } else { // block size
        if (A.symmetric) {
            // if there is a main diagonal block, it is the first in offsets
            // and should be skipped in the first loop below since the main
            // diagonal is processed in the second loop
            const IndexType d0 = (A.diagonal_offsets[0] == 0 ? 1 : 0);
            const ValueType* values = thrust::raw_pointer_cast(&A.values.values[0]);
            const IndexType pitch = A.values.pitch;
#pragma omp parallel for
            for (IndexType ch = 0; ch < num_rows; ch+=chunksize)
            {
                for (IndexType row = ch; row<std::min(ch+chunksize,num_rows); row++)
                {
                    y[row] = initialize(y[row]);
                }

                IndexType idx = pitch*block_size*(num_diagonals-1);
                // process subdiagonal blocks
                for (IndexType d = 0; d < num_diagonals-d0; d++)
                {
                    const IndexType diag = num_diagonals-d-1;
                    const IndexType k = -block_size*A.diagonal_offsets[diag];
                    for (IndexType row = ch; row<std::min(ch+chunksize,num_rows); row+=block_size)
                    {
                        const IndexType col = row + k;
                        if (col >= 0 && col <= num_rows-block_size)
                        {
                            // for each row in block
                            for (IndexType j = 0; j < block_size; j++)
                            {
                                // for each column in block
                                for (IndexType i = 0; i < block_size; i++)
                                {
                                    const ValueType& Aij = values[idx+col+i+j*pitch];
                                    const ValueType& xj = x[col + i];

                                    y[row+j] = reduce(y[row+j], combine(Aij, xj));
                                }
                            }
                        }
                    }
                    idx -= block_size*pitch;
                }
                // process main and upper diagonal blocks
                for (IndexType d = 0; d < num_diagonals; d++)
                {
                    const IndexType k = A.diagonal_offsets[d]*block_size;
                    for (IndexType row = ch; row<std::min(ch+chunksize,num_rows); row+=block_size)
                    {
                        const IndexType col = row + k;
                        if (col >= 0 && col <= num_rows-block_size)
                        {
                            // for each column in block
                            for (IndexType i = 0; i < block_size; i++)
                            {
                                // for each row in block
                                for (IndexType j = 0; j < block_size; j++)
                                {
                                    const ValueType& Aij = values[row+j+(d*block_size+i)*pitch];
                                    const ValueType& xj = x[col + i];

                                    y[row+j] = reduce(y[row+j], combine(Aij, xj));
                                }
                            }
                        }
                    }
                } // diagonals
            }
        } else { // !A.symmetric
#pragma omp parallel for
            for (IndexType ch = 0; ch < num_rows; ch+=chunksize)
            {
                for (IndexType row = ch; row<std::min(ch+chunksize,num_rows); row++)
                {
                    y[row] = initialize(y[row]);
                }

                // for each diagonal block
                for (IndexType d = 0; d < num_diagonals; d++)
                {
                    const IndexType k = A.diagonal_offsets[d]*block_size;

                    for (IndexType row=ch; row<std::min(ch+chunksize,num_rows); row+=block_size)
                    {
                        const IndexType col = row + k;
                        if (col >= 0 && col <= num_rows-block_size)
                        {
                            // for each column in block
                            for (IndexType i = 0; i < block_size; i++)
                            {
                                // for each row in block
                                for (IndexType j = 0; j < block_size; j++)
                                {
                                    const ValueType& Aij = A.values(row+j, d*block_size+i);
                                    const ValueType& xj = x[col + i];

                                    y[row+j] = reduce(y[row+j], combine(Aij, xj));
                                }
                            }
                        }
                    }
                } // diagonals
            } // row chunks
        } // A.symmetric
    } // block size
}

template <typename Matrix,
          typename Vector1,
          typename Vector2>
void spmv_cds(const Matrix&  A,
              const Vector1& x,
                    Vector2& y)
{
    typedef typename Vector2::value_type ValueType;

    if (A.block_size == 1) {
        spmv_dia(A, x, y,
                 cusp::detail::zero_function<ValueType>(),
                 thrust::multiplies<ValueType>(),
                 thrust::plus<ValueType>());
    } else {
        spmv_cds(A, x, y,
                  cusp::detail::zero_function<ValueType>(),
                  thrust::multiplies<ValueType>(),
                  thrust::plus<ValueType>());
    }
}


//////////////
// ELL SpMV //
//////////////
template <typename Matrix,
          typename Vector1,
          typename Vector2,
          typename UnaryFunction,
          typename BinaryFunction1,
          typename BinaryFunction2>
void spmv_ell(const Matrix&  A,
              const Vector1& x,
                    Vector2& y,
              UnaryFunction   initialize,
              BinaryFunction1 combine,
              BinaryFunction2 reduce)
{
    typedef typename Matrix::index_type  IndexType;
    typedef typename Vector2::value_type ValueType;

    const size_t& num_entries_per_row = A.column_indices.num_cols;

    const IndexType invalid_index = Matrix::invalid_index;
    
    for(size_t i = 0; i < A.num_rows; i++)
        y[i] = initialize(y[i]);

    for(size_t n = 0; n < num_entries_per_row; n++)
    {
        for(size_t i = 0; i < A.num_rows; i++)
        {
            const IndexType& j   = A.column_indices(i, n);
            const ValueType& Aij = A.values(i,n);

            if (j != invalid_index)
            {
                const ValueType& xj = x[j];
                y[i] = reduce(y[i], combine(Aij, xj));
            }
        }
    }
}


template <typename Matrix,
          typename Vector1,
          typename Vector2>
void spmv_ell(const Matrix&  A,
              const Vector1& x,
                    Vector2& y)
{
    typedef typename Vector2::value_type ValueType;

    spmv_ell(A, x, y,
             cusp::detail::zero_function<ValueType>(),
             thrust::multiplies<ValueType>(),
             thrust::plus<ValueType>());
}

} // end namespace host
} // end namespace detail
} // end namespace cusp


