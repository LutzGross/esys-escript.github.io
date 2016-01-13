
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

