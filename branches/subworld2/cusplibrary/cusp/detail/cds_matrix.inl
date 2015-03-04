/*****************************************************************************
 *
 * Copyright (c) 2014-2015 by University of Queensland
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

#include <cusp/convert.h>
#include <cusp/detail/utils.h>

namespace cusp
{

//////////////////
// Constructors //
//////////////////
        
// construct from a different matrix
template <typename IndexType, typename ValueType, class MemorySpace>
template <typename MatrixType>
cds_matrix<IndexType,ValueType,MemorySpace>
    ::cds_matrix(const MatrixType& matrix)
    {
        cusp::convert(matrix, *this);
    }

//////////////////////
// Member Functions //
//////////////////////

// copy a matrix in a different format
template <typename IndexType, typename ValueType, class MemorySpace>
template <typename MatrixType>
    cds_matrix<IndexType,ValueType,MemorySpace>&
    cds_matrix<IndexType,ValueType,MemorySpace>
    ::operator=(const MatrixType& matrix)
    {
        cusp::convert(matrix, *this);
        
        return *this;
    }

///////////////////////////
// Convenience Functions //
///////////////////////////

template <typename Array1,
          typename Array2>
cds_matrix_view<Array1,Array2>
make_cds_matrix_view(size_t num_rows,
                     size_t num_entries,
                     size_t block_size,
                     Array1 diagonal_offsets,
                     Array2 values)
{
  return cds_matrix_view<Array1,Array2>
    (num_rows, num_entries, block_size,
     diagonal_offsets, values);
}

template <typename Array1,
          typename Array2,
          typename IndexType,
          typename ValueType,
          typename MemorySpace>
cds_matrix_view<Array1,Array2,IndexType,ValueType,MemorySpace>
make_cds_matrix_view(const cds_matrix_view<Array1,Array2,IndexType,ValueType,MemorySpace>& m)
{
  return cds_matrix_view<Array1,Array2,IndexType,ValueType,MemorySpace>(m);
}
    
template <typename IndexType, typename ValueType, class MemorySpace>
typename cds_matrix<IndexType,ValueType,MemorySpace>::view
make_cds_matrix_view(cds_matrix<IndexType,ValueType,MemorySpace>& m)
{
  return make_cds_matrix_view
    (m.num_rows, m.num_entries, m.block_size,
     cusp::make_array1d_view(m.diagonal_offsets),
     cusp::make_array2d_view(m.values));
}

template <typename IndexType, typename ValueType, class MemorySpace>
typename cds_matrix<IndexType,ValueType,MemorySpace>::const_view
make_cds_matrix_view(const cds_matrix<IndexType,ValueType,MemorySpace>& m)
{
  return make_cds_matrix_view
    (m.num_rows, m.num_entries, m.block_size,
     cusp::make_array1d_view(m.diagonal_offsets),
     cusp::make_array2d_view(m.values));
}

} // end namespace cusp

