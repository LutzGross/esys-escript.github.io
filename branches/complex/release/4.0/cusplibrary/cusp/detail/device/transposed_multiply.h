
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

