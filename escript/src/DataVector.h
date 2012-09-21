
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


#if !defined escript_DataVector_20050324_H
#define escript_DataVector_20050324_H
#include "system_dep.h"

#include "esysUtils/EsysAssert.h"

#include <vector>
#include <iostream>
#include <fstream>

namespace escript {

class WrappedArray;

/**
   \brief
   DataVector implements an arbitrarily long vector of data values.
   DataVector is the underlying data container for Data objects.

   Description:
   DataVector provides an implementation of a vector of data values for use
   by DataBlocks2D and DataArrayView. Hiding the vector in this container
   allows different implementations to be swapped in without disrupting the
   client classes.
*/

class ESCRIPT_DLL_API DataVector {

 public:

  //
  // The type of the elements stored in the vector.
  typedef double ElementType;

  //
  // The underlying type used to implement the vector.
  typedef ElementType *  ValueType;
  typedef const ElementType * ConstValueType;

  //
  // Various types exported to clients of this class.
  typedef ElementType          value_type;
  typedef long                 size_type;
  typedef ElementType &        reference;
  typedef const ElementType &  const_reference;

  /**
     \brief
     Default constructor for DataVector.

     Description:
     Constructs an empty DataVector object.
  */
  DataVector();

  /**
     \brief
     Copy constructor for DataVector.

     Description:
     Constructs a DataVector object which is a copy of the
     given DataVector object.
  */
  DataVector(const DataVector& other);

  /**
     \brief
     Constructor for DataVector.

     Description:
     Constructs a DataVector object of length "size" with all elements
     initilised to "val".

     \param size - Input - Number of elements in the vector.
     \param val - Input - Initial value for all elements in the vector. Default is 0.0.
     \param blockSize - Input - size of blocks within the vector, overall vector
                size must be a precise multiple of the block size. Default is 1.

     In escript::Data, blocksize corresponds to the number of elements required to hold all
     the data-points for a sample, ie: the product of the dimensions of a data-point and the
     number of data-points per sample. Size is the total number of elements required to hold
     all elements for all data-points in the given object, ie: number of samples * blocksize.
  */
  DataVector(const size_type size,
             const value_type val=0.0,
             const size_type blockSize=1);

  /**
     \brief
     Default destructor for DataVector.

     Description:
     Destroys the current DataVector object.
  */
  ~DataVector();

  /**
     \brief
     Resize the DataVector to the given length "newSize".
     All current data is lost. All elements in the new DataVector are
     initialised to "newVal".

     \param newSize - Input - New size for the vector.
     \param newVal - Input - New initial value for all elements in the vector.
     \param newBlockSize - Input - New block size for the vector.
  */
  void
  resize(const size_type newSize,
         const value_type newVal=0.0,
         const size_type newBlockSize=1);

  /**
    \brief 
    Populates the vector with the data from value.
    This method currently throws an exception if the specified number of copies won't fit.
    \warning This function does not attempt to perform shape checking.
  */
  void
  copyFromArray(const escript::WrappedArray& value, size_type copies);

  void 
  copyFromArrayToOffset(const WrappedArray& value, size_type offset, size_type copies);


  /**
     \brief
     Return the number of elements in this DataVector.
  */
  inline
  size_type
  size() const;

  /**
     \brief
     DataVector assignment operator "=".
     Assign the given DataVector object to this.
  */
  DataVector&
  operator=(const DataVector& other);

  /**
     \brief
     DataVector equality comparison operator "==".
     Return true if the given DataVector is equal to this.
  */
  bool
  operator==(const DataVector& other) const;

  /**
     \brief
     DataVector inequality comparison operator "!=".
     Return true if the given DataVector is not equal to this.
  */
  bool
  operator!=(const DataVector& other) const;

  /**
    \brief
    Return a reference to the element at position i in this DataVector.
    Will throw an exception if an invalid index "i" is given.

    NB: access to the element one past the end of the vector is permitted
    in order to provide a facility equivalent to an end() pointer.
  */
  inline
  reference
  operator[](const size_type i);

  inline
  const_reference
  operator[](const size_type i) const;


 protected:

 private:

  size_type m_size;
  size_type m_dim;
  size_type m_N;

  //
  // The container for the elements contained in this DataVector.
  ValueType m_array_data;
};

/**
  \brief
  releases unused memory in the memory manager.
*/
                                                                                                                                                                                                     
ESCRIPT_DLL_API void releaseUnusedMemory();
                                                                                                                                                                                                     


inline
DataVector::size_type
DataVector::size() const
{
  return m_size;
}

inline
DataVector::reference
DataVector::operator[](const DataVector::size_type i)
{
  EsysAssert(i<size(),"DataVector: invalid index specified. " << i << " of " << size());
  return m_array_data[i];
}

inline
DataVector::const_reference
DataVector::operator[](const DataVector::size_type i) const
{
  EsysAssert(i<size(),"DataVector: invalid index specified. " << i << " of " << size());
  return m_array_data[i];
}

} // end of namespace

#endif
