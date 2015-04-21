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

// SpMV
#include <cusp/detail/device/transposed_spmv/cds.h>
#include <cusp/detail/device/transposed_spmv/dia.h>

namespace cusp
{
namespace detail
{
namespace device
{

//////////////////////////////////
// Dense Matrix-Vector Multiply //
//////////////////////////////////
//// TODO implement this for both row and column-major ordering
//template <typename Matrix,
//          typename Vector1,
//          typename Vector2>
//void multiply(const Matrix&  A,
//              const Vector1& B,
//                    Vector2& C,
//              cusp::array2d_format,
//              cusp::array1d_format,
//              cusp::array1d_format)
//{
//}

///////////////////////////////////
// Sparse Matrix-Vector Multiply //
///////////////////////////////////
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
#ifdef CUSP_USE_TEXTURE_MEMORY
    cusp::detail::device::transposed_spmv_cds_tex(A, thrust::raw_pointer_cast(&B[0]), thrust::raw_pointer_cast(&C[0]));
#else
    cusp::detail::device::transposed_spmv_cds(A, thrust::raw_pointer_cast(&B[0]), thrust::raw_pointer_cast(&C[0]));
#endif
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
#ifdef CUSP_USE_TEXTURE_MEMORY
    cusp::detail::device::transposed_spmv_dia_tex(A, thrust::raw_pointer_cast(&B[0]), thrust::raw_pointer_cast(&C[0]));
#else
    cusp::detail::device::transposed_spmv_dia(A, thrust::raw_pointer_cast(&B[0]), thrust::raw_pointer_cast(&C[0]));
#endif
}


/////////////////
// Entry Point //
/////////////////
template <typename Matrix,
         typename MatrixOrVector1,
         typename MatrixOrVector2>
void transposed_multiply(const Matrix&  A,
                         const MatrixOrVector1& B,
                               MatrixOrVector2& C)
{
    cusp::detail::device::transposed_multiply(A, B, C,
                                   typename Matrix::format(),
                                   typename MatrixOrVector1::format(),
                                   typename MatrixOrVector2::format());
}

} // end namespace device
} // end namespace detail
} // end namespace cusp

