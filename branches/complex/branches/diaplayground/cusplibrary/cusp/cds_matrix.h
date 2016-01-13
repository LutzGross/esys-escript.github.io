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

/*! \file cds_matrix.h
 *  \brief Column Diagonal Block matrix format.
 */

#pragma once

#include <cusp/detail/config.h>

#include <cusp/array1d.h>
#include <cusp/format.h>
#include <cusp/detail/matrix_base.h>
#include <cusp/detail/utils.h>

namespace cusp
{

/*! \addtogroup sparse_matrices Sparse Matrices
 */

/*! \addtogroup sparse_matrix_containers Sparse Matrix Containers
 *  \ingroup sparse_matrices
 *  \{
 */

    // Forward definitions
    struct column_major;
    template<typename ValueType, class MemorySpace, class Orientation> class array2d;
    template<typename Array, class Orientation>                        class array2d_view;
    template <typename Array1, typename Array2, typename IndexType, typename ValueType, typename MemorySpace> class cds_matrix_view;

/*! \p cds_matrix : Diagonal block matrix container (always square)
 *
 * \tparam IndexType Type used for matrix indices (e.g. \c int).
 * \tparam ValueType Type used for matrix values (e.g. \c float).
 * \tparam MemorySpace A memory space (e.g. \c cusp::host_memory or cusp::device_memory)
 *
 * \note The diagonal offsets should not contain duplicate entries.
 *
 *  The following code snippet demonstrates how to create a 8-by-8
 *  \p cds_matrix on the host with 2 diagonal 2x2 blocks (16 total nonzeros)
 *  and then copies the matrix to the device.
 *
 *  \code
 *  #include <cusp/cds_matrix.h>
 *  ...
 *
 *  // allocate storage for (8,8) matrix with 16 nonzeros in 2 diagonal blocks
 *  // of size 2x2
 *  cusp::cds_matrix<int,float,cusp::host_memory> A(8,16,2,2);
 *
 *  // initialize diagonal offsets
 *  A.diagonal_offsets[0] = -2;
 *  A.diagonal_offsets[1] =  2;
 *
 *  // initialize diagonal values
 *
 *  // first diagonal block
 *  A.values( 0,0) =  0;  // outside matrix
 *  A.values( 1,0) =  0;  // outside matrix
 *  A.values( 2,0) =  0;  // outside matrix
 *  A.values( 3,0) =  0;  // outside matrix
 *  A.values( 4,0) =  0;  // outside matrix
 *  A.values( 5,0) =  0;  // outside matrix
 *  A.values( 6,0) =  0;  // outside matrix
 *  A.values( 7,0) =  0;  // outside matrix
 *  A.values( 8,0) = -10;
 *  A.values( 9,0) = -20;
 *  A.values(10,0) = -30;
 *  A.values(11,0) = -40;
 *  A.values(12,0) = -50;
 *  A.values(13,0) = -60;
 *  A.values(14,0) = -70;
 *  A.values(15,0) = -80;
 *  
 *  // second diagonal block
 *  A.values( 0,1) =  10;
 *  A.values( 1,1) =  20;
 *  A.values( 2,1) =  30;
 *  A.values( 3,1) =  40;
 *  A.values( 4,1) =  50;
 *  A.values( 5,1) =  60;
 *  A.values( 6,1) =  70;
 *  A.values( 7,1) =  80;
 *  A.values( 8,1) =  0;  // outside matrix
 *  A.values( 9,1) =  0;  // outside matrix
 *  A.values(10,1) =  0;  // outside matrix
 *  A.values(11,1) =  0;  // outside matrix
 *  A.values(12,1) =  0;  // outside matrix
 *  A.values(13,1) =  0;  // outside matrix
 *  A.values(14,1) =  0;  // outside matrix
 *  A.values(15,1) =  0;  // outside matrix
 *
 *  // A now represents the following matrix
 *  //    [  0   0   0   0  10  30   0   0]
 *  //    [  0   0   0   0  20  40   0   0]
 *  //    [  0   0   0   0   0   0  50  70]
 *  //    [  0   0   0   0   0   0  60  80]
 *  //    [-10 -30   0   0   0   0   0   0]
 *  //    [-20 -40   0   0   0   0   0   0]
 *  //    [  0  0  -50 -70   0   0   0   0]
 *  //    [  0  0  -60 -80   0   0   0   0]
 *
 *  // copy to the device
 *  cusp::cds_matrix<int,float,cusp::device_memory> B = A;
 *  \endcode
 *
 */
template <typename IndexType, typename ValueType, class MemorySpace>
class cds_matrix : public detail::matrix_base<IndexType,ValueType,MemorySpace,cusp::cds_format>
{
  typedef cusp::detail::matrix_base<IndexType,ValueType,MemorySpace,cusp::cds_format> Parent;
  public:
    // TODO statically assert is_signed<IndexType>
    
    /*! rebind matrix to a different MemorySpace
     */
    template<typename MemorySpace2>
    struct rebind { typedef cusp::cds_matrix<IndexType, ValueType, MemorySpace2> type; };

    /*! type of diagonal offsets array
     */
    typedef typename cusp::array1d<IndexType, MemorySpace>                     diagonal_offsets_array_type;
    
    /*! type of values array
     */
    typedef typename cusp::array2d<ValueType, MemorySpace, cusp::column_major> values_array_type;

    /*! equivalent container type
     */
    typedef typename cusp::cds_matrix<IndexType, ValueType, MemorySpace> container;

    /*! equivalent view type
     */
    typedef typename cusp::cds_matrix_view<typename diagonal_offsets_array_type::view,
                                           typename values_array_type::view,
                                           IndexType, ValueType, MemorySpace> view;
    
    /*! equivalent const_view type
     */
    typedef typename cusp::cds_matrix_view<typename diagonal_offsets_array_type::const_view,
                                           typename values_array_type::const_view,
                                           IndexType, ValueType, MemorySpace> const_view;

    /*! Storage for the diagonal block offsets.
     */
    diagonal_offsets_array_type diagonal_offsets;
    
    /*! Storage for the nonzero entries of the CDS data structure.
     */
    values_array_type values;

    /*! Diagonal block size.
     */
    size_t block_size;

    /*! Construct an empty \p cds_matrix.
     */
    cds_matrix() {}

    /*! Construct a \p cds_matrix with a specific shape, number of nonzero entries,
     *  and number of occupied diagonals.
     *
     *  \param num_rows Number of rows & columns.
     *  \param num_entries Number of nonzero matrix entries.
     *  \param num_diagonals Number of occupied diagonals.
     *  \param alignment Amount of padding used to align the data structure (default 32).
     */
    cds_matrix(size_t num_rows, size_t num_entries,
               size_t num_diagonals, size_t blocksize, size_t alignment = 32)
      : Parent(num_rows, num_rows, num_entries),
        diagonal_offsets(num_diagonals),
        block_size(blocksize)
      {
        // TODO use array2d constructor when it can accept pitch
        values.resize(num_rows, num_diagonals*blocksize, detail::round_up(num_rows, alignment));
      }
    
    /*! Construct a \p cds_matrix from another matrix.
     *
     *  \param matrix Another sparse or dense matrix.
     */
    template <typename MatrixType>
    cds_matrix(const MatrixType& matrix);
    
    /*! Resize matrix dimensions and underlying storage
     */
    void resize(size_t num_rows, size_t num_entries, size_t num_diagonals,
                size_t blocksize)
    {
      Parent::resize(num_rows, num_rows, num_entries);
      diagonal_offsets.resize(num_diagonals);
      block_size = blocksize;
      values.resize(num_rows, num_diagonals*blocksize);
    }
               
    /*! Resize matrix dimensions and underlying storage
     */
    void resize(size_t num_rows, size_t num_entries, size_t num_diagonals,
                size_t blocksize, size_t alignment)
    {
      Parent::resize(num_rows, num_rows, num_entries);
      diagonal_offsets.resize(num_diagonals);
      block_size = blocksize;
      values.resize(num_rows, num_diagonals*blocksize, detail::round_up(num_rows, alignment));
    }
    
    /*! Swap the contents of two \p cds_matrix objects.
     *
     *  \param matrix Another \p cds_matrix with the same IndexType and ValueType.
     */
    void swap(cds_matrix& matrix)
    {
      Parent::swap(matrix);
      diagonal_offsets.swap(matrix.diagonal_offsets);
      thrust::swap(block_size, matrix.block_size);
      values.swap(matrix.values);
    }
    
    /*! Assignment from another matrix.
     *
     *  \param matrix Another sparse or dense matrix.
     */
    template <typename MatrixType>
    cds_matrix& operator=(const MatrixType& matrix);
}; // class cds_matrix
/*! \}
 */
    
/*! \addtogroup sparse_matrix_views Sparse Matrix Views
 *  \ingroup sparse_matrices
 *  \{
 */

/*! \p cds_matrix_view : Diagonal block matrix view
 *
 * \tparam Array1 Type of \c diagonal_offsets
 * \tparam Array2 Type of \c values array view
 * \tparam IndexType Type used for matrix indices (e.g. \c int).
 * \tparam ValueType Type used for matrix values (e.g. \c float).
 * \tparam MemorySpace A memory space (e.g. \c cusp::host_memory or cusp::device_memory)
 *
 */
template <typename Array1,
          typename Array2,
          typename IndexType   = typename Array1::value_type,
          typename ValueType   = typename Array2::value_type,
          typename MemorySpace = typename cusp::minimum_space<typename Array1::memory_space, typename Array2::memory_space>::type >
class cds_matrix_view : public detail::matrix_base<IndexType,ValueType,MemorySpace,cusp::cds_format>
{
  typedef cusp::detail::matrix_base<IndexType,ValueType,MemorySpace,cusp::cds_format> Parent;
  public:
    /*! type of \c diagonal_offsets array
     */
    typedef Array1 diagonal_offsets_array_type;
    
    /*! type of \c column_indices array
     */
    typedef Array2 values_array_type;

    /*! equivalent container type
     */
    typedef typename cusp::cds_matrix<IndexType, ValueType, MemorySpace> container;

    /*! equivalent view type
     */
    typedef typename cusp::cds_matrix_view<Array1, Array2, IndexType, ValueType, MemorySpace> view;

    /*! Storage for the diagonal offsets.
     */
    diagonal_offsets_array_type diagonal_offsets;

    /*! Storage for the nonzero entries of the DIA data structure.
     */
    values_array_type values;

    /*! Size of the diagonal blocks.
     */
    size_t block_size;

    /*! Construct an empty \p cds_matrix_view.
     */
    cds_matrix_view() {}

    template <typename OtherArray1, typename OtherArray2>
    cds_matrix_view(size_t num_rows, size_t num_entries, size_t blocksize,
                    OtherArray1& diagonal_offsets, OtherArray2& values)
    : Parent(num_rows, num_rows, num_entries), block_size(blocksize),
             diagonal_offsets(diagonal_offsets), values(values) {}

    template <typename OtherArray1, typename OtherArray2>
    cds_matrix_view(size_t num_rows, size_t num_entries, size_t blocksize,
                    const OtherArray1& diagonal_offsets, const OtherArray2& values)
    : Parent(num_rows, num_rows, num_entries), block_size(blocksize),
             diagonal_offsets(diagonal_offsets), values(values) {}
    
    template <typename Matrix>
    cds_matrix_view(Matrix& A)
    : Parent(A), block_size(A.block_size), diagonal_offsets(A.diagonal_offsets), values(A.values) {}
    
    template <typename Matrix>
    cds_matrix_view(const Matrix& A)
    : Parent(A), block_size(A.block_size), diagonal_offsets(A.diagonal_offsets), values(A.values) {}
    
    /*! Resize matrix dimensions and underlying storage
     */
    void resize(size_t num_rows, size_t num_entries, size_t num_diagonals,
                size_t blocksize)
    {
      Parent::resize(num_rows, num_rows, num_entries);
      diagonal_offsets.resize(num_diagonals);
      block_size = blocksize;
      values.resize(num_rows, num_diagonals*blocksize);
    }
               
    /*! Resize matrix dimensions and underlying storage
     */
    void resize(size_t num_rows, size_t num_entries, size_t num_diagonals,
                size_t blocksize, size_t alignment)
    {
      Parent::resize(num_rows, num_rows, num_entries);
      diagonal_offsets.resize(num_diagonals);
      block_size = blocksize;
      values.resize(num_rows, num_diagonals*blocksize, detail::round_up(num_rows, alignment));
    }
}; // class cds_matrix_view


template <typename Array1,
          typename Array2>
cds_matrix_view<Array1,Array2>
make_cds_matrix_view(size_t num_rows,
                     size_t num_entries,
                     size_t block_size,
                     Array1 diagonal_offsets,
                     Array2 values);

template <typename Array1,
          typename Array2,
          typename IndexType,
          typename ValueType,
          typename MemorySpace>
cds_matrix_view<Array1,Array2,IndexType,ValueType,MemorySpace>
make_cds_matrix_view(const cds_matrix_view<Array1,Array2,IndexType,ValueType,MemorySpace>& m);
    
template <typename IndexType, typename ValueType, class MemorySpace>
typename cds_matrix<IndexType,ValueType,MemorySpace>::view
make_cds_matrix_view(cds_matrix<IndexType,ValueType,MemorySpace>& m);

template <typename IndexType, typename ValueType, class MemorySpace>
typename cds_matrix<IndexType,ValueType,MemorySpace>::const_view
make_cds_matrix_view(const cds_matrix<IndexType,ValueType,MemorySpace>& m);
/*! \} // end Views
 */
    
} // end namespace cusp

#include <cusp/array2d.h>
#include <cusp/detail/cds_matrix.inl>

