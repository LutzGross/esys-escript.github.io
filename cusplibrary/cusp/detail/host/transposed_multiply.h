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

#include <cusp/format.h>

#include <cusp/detail/functional.h>

#include <cusp/detail/host/spmv.h>
#include <cusp/detail/host/transposed_spmv.h>


namespace cusp
{
namespace detail
{
namespace host
{

/////////////////////////////////////////////
// Dense Transposed Matrix-Vector Multiply //
/////////////////////////////////////////////
template <typename Matrix,
         typename Vector1,
         typename Vector2>
void transposed_multiply(const Matrix&  A,
                         const Vector1& B,
                               Vector2& C,
                         cusp::array2d_format,
                         cusp::array1d_format,
                         cusp::array1d_format)
{
    typedef typename Vector2::value_type ValueType;

    for(size_t i = 0; i < A.num_cols; i++)
    {
        ValueType sum = 0;
        for(size_t j = 0; j < A.num_rows; j++)
        {
            sum += A(j,i) * B[j];
        }
        C[i] = sum;
    }
}

//////////////////////////////////////////////
// Sparse Transposed Matrix-Vector Multiply //
//////////////////////////////////////////////
template <typename Matrix,
         typename Vector1,
         typename Vector2>
void transposed_multiply(const Matrix&  A,
                         const Vector1& B,
                               Vector2& C,
                         cusp::cds_format,
                         cusp::array1d_format,
                         cusp::array1d_format)
{
    if (A.symmetric)
        cusp::detail::host::spmv_cds(A, B, C);
    else
        cusp::detail::host::transposed_spmv_cds(A, B, C);
}

template <typename Matrix,
         typename Vector1,
         typename Vector2>
void transposed_multiply(const Matrix&  A,
                         const Vector1& B,
                               Vector2& C,
                         cusp::dia_format,
                         cusp::array1d_format,
                         cusp::array1d_format)
{
    if (A.symmetric)
        cusp::detail::host::spmv_dia(A, B, C);
    else
        cusp::detail::host::transposed_spmv_dia(A, B, C);
}


/////////////////
// Entry Point //
/////////////////
template <typename Matrix,
         typename MatrixOrVector1,
         typename MatrixOrVector2>
void transposed_multiply(const Matrix& A,
                         const MatrixOrVector1& B,
                               MatrixOrVector2& C)
{
    cusp::detail::host::transposed_multiply(A, B, C,
                                 typename Matrix::format(),
                                 typename MatrixOrVector1::format(),
                                 typename MatrixOrVector2::format());
}

} // end namespace host
} // end namespace detail
} // end namespace cusp

