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

#include <vector>

#include "esysUtils/EsysAssert.h"

namespace escript {

/**
   \brief
   DataVector implements an arbitrarily long vector of data values.

   Description:
   DataVector provides an implementation of a vector of data values for use
   by DataBlocks2D and DataArrayView. Hiding the vector in this container
   allows different implementations to be swapped in without disrupting the
   client classes. This is the underlying data container for Data objects.
*/

class DataVector {

 public:

  //
  // The type of the elements stored in the vector.
  typedef double ElementType;

  //
  // The underlying type used to implement the vector.
  typedef std::vector<ElementType> ValueType;

  //
  // Various types needed by clients of this class.
  typedef ValueType::value_type       value_type;
  typedef ValueType::size_type        size_type;
  typedef ValueType::reference        reference;
  typedef ValueType::const_reference  const_reference;
  typedef ValueType::iterator         iterator;

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
     initilised to "val". Default for "val" is zero.

     \param size - Input - Number of elements in the vector.
  */
  DataVector(ValueType::size_type size, ValueType::value_type val=0.0);

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
     initialised to "newVal", which defaults to zero.

     \param newSize - Input - New size for the vector.
  */
  inline
  void
  resize(ValueType::size_type newSize, ValueType::value_type newVal=0.0);

  /**
     \brief
     Return the number of elements in this DataVector.
  */
  inline
  ValueType::size_type
  size() const;

  /**
     \brief
     DataVector assignment operator "=".
     Assign the given DataVector object to this.
  */
  inline
  DataVector&
  operator=(const DataVector& other);

  /**
     \brief
     DataVector equality comparison operator "==".
     Return true if the given DataVector is equal to this.
  */
  inline
  bool
  operator==(const DataVector& other);

  /**
     \brief
     DataVector inequality comparison operator "!=".
     Return true if the given DataVector is not equal to this.
  */
  inline
  bool
  operator!=(const DataVector& other);

  /**
    \brief
    Return a reference to the element at position i in this DataVector.
    Will throw an exception if an invalid index "i" is given.

    NB: access to the element one past the end of the vector is permitted
    in order to provide a facility equivalent to an end() pointer.
  */
  inline
  ValueType::reference
  operator[](ValueType::size_type i);

  inline
  ValueType::const_reference
  operator[](ValueType::size_type i) const;

  /**
    \brief
    Add the element x at the end of the vector.
  */
  inline
  void
  push_back(const ValueType::value_type& x);

  /**
    \brief
    Insert.
  */
  inline
  void
  insert(iterator pos, const ValueType::value_type* first, const ValueType::value_type* last);

  /**
    \brief
    End.
  */
  inline
  iterator
  end();

 protected:

 private:

  //
  // The container for the elements contained in this DataVector.
  ValueType m_data;

};

inline
DataVector::ValueType::size_type
DataVector::size() const
{
  return m_data.size();
}

inline
DataVector::ValueType::reference
DataVector::operator[](DataVector::ValueType::size_type i)
{
  // Allow access to element one beyond end of vector to simulate end().
  EsysAssert(i<=size(),"DataVector: invalid index specified.");
  return m_data[i];
}

inline
DataVector::ValueType::const_reference
DataVector::operator[](DataVector::ValueType::size_type i) const
{
  // Allow access to element one beyond end of vector to simulate end().
  EsysAssert(i<=size(),"DataVector: invalid index specified.");
  return m_data[i];
}

inline
DataVector&
DataVector::operator=(const DataVector& other)
{
  DataVector temp(other);
  swap(m_data,temp.m_data);
  return *this;
}

inline
bool
DataVector::operator==(const DataVector& other)
{
  return m_data==other.m_data;
}

inline
bool
DataVector::operator!=(const DataVector& other)
{
  return m_data!=other.m_data;
}

inline
void
DataVector::resize(DataVector::ValueType::size_type newSize, DataVector::ValueType::value_type newValue)
{
  m_data.resize(newSize,newValue);
}

inline
void
DataVector::push_back(const ValueType::value_type& x)
{
  m_data.push_back(x);
}

inline
void
DataVector::insert(DataVector::iterator pos, const ValueType::value_type* first, const ValueType::value_type* last)
{
  m_data.insert(pos, first, last);
}

inline
DataVector::iterator
DataVector::end()
{
  return m_data.end();
}

} // end of namespace

#endif
