
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#if !defined escript_DataVectorTaipan_H
#define escript_DataVectorTaipan_H
#include "system_dep.h"

#include "Assert.h"
#include "DataTypes.h"
#include "WrappedArray.h"

namespace escript
{

namespace DataTypes
{
  
  /**
   \brief
   DataVectorTaipan implements an arbitrarily long vector of data values.
   DataVectorTaipan is the underlying data container for Data objects.

   Description:
   DataVectorTaipan provides an implementation of a vector of data values for use
   by DataBlocks2D and DataArrayView. Hiding the vector in this container
   allows different implementations to be swapped in without disrupting the
   client classes.
*/

class ESCRIPT_DLL_API DataVectorTaipan {

 public:

  //
  // The type of the elements stored in the vector.
  typedef double ElementType;

  //
  // The underlying type used to implement the vector.
  typedef ElementType *  VectorStorageType;


  //
  // Various types exported to clients of this class.
  typedef const ElementType *  const_pointer;  
  typedef ElementType          value_type;
  typedef long                 size_type;
  typedef ElementType &        reference;
  typedef const ElementType &  const_reference;

  /**
     \brief
     Default constructor for DataVectorTaipan.

     Description:
     Constructs an empty DataVectorTaipan object.
  */
  DataVectorTaipan();

  /**
     \brief
     Copy constructor for DataVectorTaipan.

     Description:
     Constructs a DataVectorTaipan object which is a copy of the
     given DataVectorTaipan object.
  */
  DataVectorTaipan(const DataVectorTaipan& other);

  /**
     \brief
     Constructor for DataVectorTaipan.

     Description:
     Constructs a DataVectorTaipan object of length "size" with all elements
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
  DataVectorTaipan(const size_type size,
             const value_type val=0.0,
             const size_type blockSize=1);

  /**
     \brief
     Default destructor for DataVectorTaipan.

     Description:
     Destroys the current DataVectorTaipan object.
  */
  ~DataVectorTaipan();

  /**
     \brief
     Resize the DataVectorTaipan to the given length "newSize".
     All current data is lost. All elements in the new DataVectorTaipan are
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
  copyFromArrayToOffset(const escript::WrappedArray& value, size_type offset, size_type copies);


  /**
     \brief
     Return the number of elements in this DataVectorTaipan.
  */
  inline
  size_type
  size() const;

  /**
     \brief
     DataVectorTaipan assignment operator "=".
     Assign the given DataVectorTaipan object to this.
  */
  DataVectorTaipan&
  operator=(const DataVectorTaipan& other);

  /**
     \brief
     DataVectorTaipan equality comparison operator "==".
     Return true if the given DataVectorTaipan is equal to this.
  */
  bool
  operator==(const DataVectorTaipan& other) const;

  /**
     \brief
     DataVectorTaipan inequality comparison operator "!=".
     Return true if the given DataVectorTaipan is not equal to this.
  */
  bool
  operator!=(const DataVectorTaipan& other) const;

  /**
    \brief
    Return a reference to the element at position i in this DataVectorTaipan.
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
  // The container for the elements contained in this DataVectorTaipan.
  VectorStorageType m_array_data;
};


/**
  \brief
  releases unused memory in the memory manager.
*/
                                                                                                                                                                                                     
ESCRIPT_DLL_API void releaseUnusedMemory();
                                                                                                                                                                                                     


inline
DataVectorTaipan::size_type
DataVectorTaipan::size() const
{
  return m_size;
}

inline
DataVectorTaipan::reference
DataVectorTaipan::operator[](const DataVectorTaipan::size_type i)
{
  ESYS_ASSERT(i<size(), "DataVectorTaipan: invalid index specified. " << i << " of " << size());
  return m_array_data[i];
}

inline
DataVectorTaipan::const_reference
DataVectorTaipan::operator[](const DataVectorTaipan::size_type i) const
{
  ESYS_ASSERT(i<size(),"DataVectorTaipan: invalid index specified. " << i << " of " << size());
  return m_array_data[i];
}


} // end of namespace 
} // end of namespace



#endif
