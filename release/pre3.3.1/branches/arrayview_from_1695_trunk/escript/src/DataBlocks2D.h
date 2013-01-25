
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#if !defined escript_DataBlocks2D_20040405_H
#define escript_DataBlocks2D_20040405_H
#include "system_dep.h"

#include "DataVector.h"

#include <sstream>
#include <iostream>

namespace escript {

/**
   \brief
   DataBlocks2D manages a 2D array of multi-dimensional data points.

   Description:
   This class is used to manage the data held by instances of
   the DataExpanded class.
*/

class DataBlocks2D {

 public:

  //
  // The type of the underlying data array under management.
  // The multi-dimensional data points are flattened and stored
  // serially as a vector of doubles.
  typedef DataVector ValueType;

  /**
     \brief
     Default constructor for DataBlocks2D.

     Description:
     Default constructor for DataBlocks2D.
     Creates an empty DataBlocks2D object.
  */
  ESCRIPT_DLL_API
  DataBlocks2D();

  /**
     \brief
     Copy constructor for DataBlocks2D.

     Description:
     Copy constructor for DataBlocks2D.
  */
  ESCRIPT_DLL_API
  DataBlocks2D(const DataBlocks2D& other);

  /**
     \brief
     Constructor for DataBlocks2D.

     Description:
     Constructor for DataBlocks2D.

     \param numRows - Input - Number of rows(samples).
     \param numCols - Input - Number of columns(data-points per sample).
     \param blockSize - Input - Number of elements per block(per data-point).

     All parameters must be >0, else an exception will be thrown.
  */
  ESCRIPT_DLL_API
  DataBlocks2D(int numRows, int numCols, int blockSize);

  /**
     \brief
     Default destructor for DataBlocks2D.

     Description:
     Default destructor for DataBlocks2D.
  */
  ESCRIPT_DLL_API
  ~DataBlocks2D();

  /**
     \brief
     Return the size of the underlying data array.
     ie: Number of rows * Number of columns * Number of elements per data point.
  */
  ESCRIPT_DLL_API
  inline
  ValueType::size_type
  size() const;

  /**
     \brief
     Return the number of rows in this DataBlocks2D array.
  */
  ESCRIPT_DLL_API
  inline
  ValueType::size_type
  getNumRows() const;

  /**
     \brief
     Return the number of columns in this DataBlocks2D array.
  */
  ESCRIPT_DLL_API
  inline
  ValueType::size_type
  getNumCols() const;

  /**
     \brief
     Return the data point size for this DataBlocks2D array.
  */
  ESCRIPT_DLL_API
  inline
  ValueType::size_type
  getBlockSize() const;

  /**
     \brief
     Resize the underlying data array. All current data is lost.
     The new data elements are initialised to 0.

     \param numRows - Input - Number of rows.
     \param numCols - Input - Number of columns.
     \param blockSize - Input - Number of elements per block.

     All parameters must be >0, else an exception will be thrown.
  */
  ESCRIPT_DLL_API
  void
  resize(int numRows, int numCols, int blockSize);

  /**
     \brief
     DataBlocks2D assignment operator =
     Assign the given DataBlocks2D object to this one.
  */
  ESCRIPT_DLL_API
  DataBlocks2D&
  operator=(const DataBlocks2D& other);

  /**
     \brief
     Swap all the values managed by the given DataBlocks2D objects.
  */
  ESCRIPT_DLL_API
  void
  Swap(DataBlocks2D& other);

  /**
    \brief
    Return the 1 dimensional index of the first element for data-point (i,j)
    within the underlying data array.
    Provides an index for accessing this data value via the [] operator.
    Subsequent elements of this data point can be accessed by manually
    incrementing the returned index value.
  */
  ESCRIPT_DLL_API
  inline
  ValueType::size_type
  index(int row, int col) const;

  /**
    \brief
    Return a reference to the first element for the data-point with index i
    within the underlying data array as determined by the index(i,j) method.
  */
  ESCRIPT_DLL_API
  inline
  ValueType::reference
  operator[](ValueType::size_type i);

  ESCRIPT_DLL_API
  inline
  ValueType::const_reference
  operator[](ValueType::size_type i) const;

  /**
    \brief
    Return a reference to the first element for the data-point (i,j).
  */
  ESCRIPT_DLL_API
  inline
  ValueType::reference
  operator()(int row, int col);

  ESCRIPT_DLL_API
  inline
  ValueType::const_reference
  operator()(int row, int col) const;

  /**
     \brief
     Return a reference to the underlying data array.
     Data returned is an array type object that can be indexed via indexes generated
     by DataBlocks2D::index.
  */
  ESCRIPT_DLL_API
  inline
  ValueType&
  getData();

  ESCRIPT_DLL_API
  inline
  const ValueType&
  getData() const;


 protected:

 private:

  //
  // The underlying array of data values.
  // The two dimensional array of multi-dimensional data points is flattened
  // and serialised within this one dimensional array of doubles.
  ValueType m_data;

  //
  // The dimensions of the 2D array of data points.
  ValueType::size_type m_numRows;
  ValueType::size_type m_numCols; 

  //
  // The number of values per data point.
  ValueType::size_type m_blockSize;

};

inline
DataBlocks2D::ValueType::size_type
DataBlocks2D::size() const
{
    EsysAssert(((m_numRows >= 0) && (m_numCols >= 0) && (m_blockSize >= 0)), "(DataBlocks2D) Invalid object.");
    return m_data.size();
}

inline
DataBlocks2D::ValueType::size_type
DataBlocks2D::getNumRows() const
{
    EsysAssert(((m_numRows >= 0) && (m_numCols >= 0) && (m_blockSize >= 0)), "(DataBlocks2D) Invalid object.");
    return m_numRows;
}

inline
DataBlocks2D::ValueType::size_type
DataBlocks2D::getNumCols() const
{
    EsysAssert(((m_numRows >= 0) && (m_numCols >= 0) && (m_blockSize >= 0)), "(DataBlocks2D) Invalid object.");
    return m_numCols;
}

inline
DataBlocks2D::ValueType::size_type
DataBlocks2D::getBlockSize() const
{
    EsysAssert(((m_numRows >= 0) && (m_numCols >= 0) && (m_blockSize >= 0)), "(DataBlocks2D) Invalid object.");
    return m_blockSize;
}

inline
DataBlocks2D::ValueType::size_type
DataBlocks2D::index(int row, int col) const
{
    EsysAssert(((m_numRows >= 0) && (m_numCols >= 0) && (m_blockSize >= 0)), "(DataBlocks2D) Invalid object.");
    EsysAssert(((row >= 0) && (col >= 0) && (m_data.size() > 0)), "(DataBlocks2D) Index value out of range.");
    ValueType::size_type temp=(row*m_numCols+col)*m_blockSize;
    EsysAssert((temp <= (m_data.size()-m_blockSize)), "(DataBlocks2D) Index value out of range.");
    return (temp);
}

inline
DataBlocks2D::ValueType::reference
DataBlocks2D::operator[](DataBlocks2D::ValueType::size_type i)
{
    EsysAssert(((m_numRows >= 0) && (m_numCols >= 0) && (m_blockSize >= 0)), "(DataBlocks2D) Invalid object.");
    return m_data[i];
}

inline
DataBlocks2D::ValueType::const_reference
DataBlocks2D::operator[](DataBlocks2D::ValueType::size_type i) const
{
    EsysAssert(((m_numRows >= 0) && (m_numCols >= 0) && (m_blockSize >= 0)), "(DataBlocks2D) Invalid object.");
    return m_data[i];
}

inline
DataBlocks2D::ValueType::reference
DataBlocks2D::operator()(int row, int col)
{
    EsysAssert(((m_numRows >= 0) && (m_numCols >= 0) && (m_blockSize >= 0)), "(DataBlocks2D) Invalid object.");
    return m_data[index(row,col)];
}

inline
DataBlocks2D::ValueType::const_reference
DataBlocks2D::operator()(int row, int col) const
{
    EsysAssert(((m_numRows >= 0) && (m_numCols >= 0) && (m_blockSize >= 0)), "(DataBlocks2D) Invalid object.");
    return m_data[index(row,col)];
}

inline
DataBlocks2D::ValueType&
DataBlocks2D::getData()
{
    EsysAssert(((m_numRows >= 0) && (m_numCols >= 0) && (m_blockSize >= 0)), "(DataBlocks2D) Invalid object.");
    return m_data;
}

inline
const DataBlocks2D::ValueType&
DataBlocks2D::getData() const
{
    EsysAssert(((m_numRows >= 0) && (m_numCols >= 0) && (m_blockSize >= 0)), "(DataBlocks2D) Invalid object.");
    return m_data;
}

} // end of namespace

#endif
