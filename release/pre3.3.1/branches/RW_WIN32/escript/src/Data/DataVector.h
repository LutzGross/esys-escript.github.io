/* 
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

#if !defined escript_DataVector_20050324_H
#define escript_DataVector_20050324_H
#ifdef MSVC
#ifdef ESCRIPT_EXPORTS
#define ESCRIPT_DLL __declspec(dllexport)
#else
#define ESCRIPT_DLL __declspec(dllimport)
#endif
#else
#define ESCRIPT_DLL
#endif

#include <vector>

#include "esysUtils/EsysAssert.h"

namespace escript {

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

	class ESCRIPT_DLL DataVector {

 public:

  //
  // The type of the elements stored in the vector.
  typedef double ElementType;

  //
  // The underlying type used to implement the vector.
  typedef ElementType *  ValueType;

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

  /**
    \brief
    Archive the data managed by this DataVector to the file referenced
    by ofstream. A count of the number of values expected to be written
    is provided as a cross-check.

    The return value indicates success (0) or otherwise (1).
  */
  int
  archiveData(std::ofstream& archiveFile,
              const size_type noValues) const;

  /**
    \brief
    Extract the number of values specified by noValues from the file
    referenced by ifstream to this DataVector.

    The return value indicates success (0) or otherwise (1).
  */
  int
  extractData(std::ifstream& archiveFile,
              const size_type noValues);

 protected:

 private:

  size_type m_size;
  size_type m_dim;
  size_type m_N;

  //
  // The container for the elements contained in this DataVector.
  ValueType m_array_data;
};

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
  EsysAssert(i<size(),"DataVector: invalid index specified.");
  return m_array_data[i];
}

inline
DataVector::const_reference
DataVector::operator[](const DataVector::size_type i) const
{
  EsysAssert(i<size(),"DataVector: invalid index specified.");
  return m_array_data[i];
}

} // end of namespace

#endif
