
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


#ifndef __ESCRIPT_DATAVECTORALT_H__
#define __ESCRIPT_DATAVECTORALT_H__

#include "DataTypes.h"
#include "system_dep.h"
#include "Assert.h"
#include "DataException.h"
#include "WrappedArray.h"

#include <sstream>

namespace escript
{

namespace DataTypes
{

template <class T>
class /* ESCRIPT_DLL_API */ DataVectorAlt {

 public:

  //
  // The type of the elements stored in the vector.
  typedef T ElementType;

  //
  // Various types exported to clients of this class.
  
  typedef const ElementType * const_pointer;
  typedef ElementType          value_type;
  typedef DataTypes::vec_size_type size_type;
  typedef ElementType &        reference;
  typedef const ElementType &  const_reference;

  /**
     \brief
     Default constructor for DataVectorAlt.

     Description:
     Constructs an empty DataVectorAlt object.
  */
  DataVectorAlt();

  /**
     \brief
     Copy constructor for DataVectorAlt.

     Description:
     Constructs a DataVectorAlt object which is a copy of the
     given DataVectorAlt object.
  */
  DataVectorAlt(const DataVectorAlt<T>& other);

  /**
     \brief
     Constructor for DataVectorAlt.

     Description:
     Constructs a DataVectorAlt object of length "size" with all elements
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
  DataVectorAlt(const size_type size,
             const value_type val=0.0,
             const size_type blockSize=1);

  /**
     \brief
     Default destructor for DataVectorAlt.

     Description:
     Destroys the current DataVectorAlt object.
  */
  ESCRIPT_INLINE_DLL_API ~DataVectorAlt();

  /**
     \brief
     Resize the DataVectorAlt to the given length "newSize".
     All current data is lost. All elements in the new DataVectorAlt are
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
  copyFromArray(const WrappedArray& value, size_type copies);

  
  // Please make sure that any implementation changes here are reflected in the specialised 
  // version in the .cpp file
  void 
  copyFromArrayToOffset(const WrappedArray& value, size_type offset, size_type copies);


  /**
     \brief
     Return the number of elements in this DataVectorAlt.
  */
  inline
  ESCRIPT_INLINE_DLL_API size_type
  size() const;

  /**
     \brief
     DataVectorAlt assignment operator "=".
     Assign the given DataVectorAlt object to this.
  */
  DataVectorAlt&
  operator=(const DataVectorAlt<T>& other);

  /**
     \brief
     DataVectorAlt equality comparison operator "==".
     Return true if the given DataVectorAlt is equal to this.
  */
  bool
  operator==(const DataVectorAlt<T>& other) const;

  /**
     \brief
     DataVectorAlt inequality comparison operator "!=".
     Return true if the given DataVectorAlt is not equal to this.
  */
  bool
  operator!=(const DataVectorAlt<T>& other) const;

  /**
    \brief
    Return a reference to the element at position i in this DataVectorAlt.
    Will throw an exception if an invalid index "i" is given.

    NB: access to the element one past the end of the vector is permitted
    in order to provide a facility equivalent to an end() pointer.
  */
  inline
  ESCRIPT_INLINE_DLL_API reference
  operator[](const size_type i);

  inline
  ESCRIPT_INLINE_DLL_API const_reference
  operator[](const size_type i) const;

    // for compatibility with std::vector
  ESCRIPT_INLINE_DLL_API ElementType* data();
  ESCRIPT_INLINE_DLL_API const ElementType* data() const; 
  
 protected:

 private:

  size_type m_size;
  size_type m_dim;
  size_type m_N;

  ElementType* m_array_data;
};

template <class T>
inline
typename DataVectorAlt<T>::ElementType* DataVectorAlt<T>::data()
{
    return m_array_data;
}

template <class T>
inline
const typename DataVectorAlt<T>::ElementType* DataVectorAlt<T>::data() const
{
    return m_array_data;  
}

template <class T>
inline
typename DataVectorAlt<T>::size_type
DataVectorAlt<T>::size() const
{
  return m_size;
}

template <class T>
inline
typename DataVectorAlt<T>::reference
DataVectorAlt<T>::operator[](const DataVectorAlt::size_type i)
{
  ESYS_ASSERT(i<size(), "DataVectorAlt: invalid index specified, " << i << " of " << size());
  return m_array_data[i];
}

template <class T>
inline
typename DataVectorAlt<T>::const_reference
DataVectorAlt<T>::operator[](const DataVectorAlt::size_type i) const
{
  ESYS_ASSERT(i<size(), "DataVectorAlt: invalid index specified. " << i << " of " << size());
  return m_array_data[i];
}


  
template <class T>
DataTypes::DataVectorAlt<T>::DataVectorAlt() :
  m_size(0),
  m_dim(0),
  m_N(0),
  m_array_data(0)
{
}

template <class T>
DataTypes::DataVectorAlt<T>::DataVectorAlt(const DataVectorAlt& other) :
  m_size(other.m_size),
  m_dim(other.m_dim),
  m_N(other.m_N),
  m_array_data(0)
{
  m_array_data=reinterpret_cast<T*>(malloc(sizeof(T)*m_size));  
  int i;
  #pragma omp parallel for private(i) schedule(static)
  for (i=0; i<m_size; i++) {
    m_array_data[i] = other.m_array_data[i];
  }
}

template <class T>
DataTypes::DataVectorAlt<T>::DataVectorAlt(const DataVectorAlt<T>::size_type size,
                       const DataVectorAlt<T>::value_type val,
                       const DataVectorAlt<T>::size_type blockSize) :
  m_size(size),
  m_dim(blockSize),
  m_array_data(0)
{
  resize(size, val, blockSize);
}

template <class T>
DataTypes::DataVectorAlt<T>::~DataVectorAlt()
{
  // clear data members
  m_size = -1;
  m_dim = -1;
  m_N = -1;
  if (m_array_data!=0)
  {
      free(m_array_data);
  }
  m_array_data=0;
}

template <class T>
void
DataVectorAlt<T>::resize(const DataVectorAlt<T>::size_type newSize,
                   const DataVectorAlt<T>::value_type newValue,
                   const DataVectorAlt<T>::size_type newBlockSize)
{
        // The < 1 is to catch both ==0 and negatives
  if ( newBlockSize < 1) {
    std::ostringstream oss;
    oss << "DataVectorAlt: invalid blockSize specified (" << newBlockSize << ')';    
    throw DataException(oss.str());
  }

  if ( newSize < 0 ) {
    std::ostringstream oss;
    oss << "DataVectorAlt: invalid new size specified (" << newSize << ')';
    throw DataException(oss.str());
  }
  if ( (newSize % newBlockSize) != 0) {
    std::ostringstream oss;
    oss << "DataVectorAlt: newSize is not a multiple of blockSize: (" << newSize << ", " << newBlockSize<< ')';
    throw DataException(oss.str());
  }

  m_size = newSize;
  m_dim = newBlockSize;
  m_N = newSize / newBlockSize;

  if (m_array_data!=0)
  {
     free(m_array_data);
  } 
  m_array_data=reinterpret_cast<T*>(malloc(sizeof(T)*m_size));  
  long i;
  #pragma omp parallel for private(i) schedule(static)
  for (i=0; i<m_size; i++) {
    m_array_data[i] = newValue;
  }
}

template <class T>
DataVectorAlt<T>&
DataVectorAlt<T>::operator=(const DataVectorAlt& other)
{
  assert(m_size >= 0);


  m_size = other.m_size;
  m_dim = other.m_dim;
  m_N = other.m_N;

  if (m_array_data!=0)
  {
      free(m_array_data);
  }
  m_array_data=reinterpret_cast<T*>(malloc(sizeof(T)*m_size));
  int i;
  #pragma omp parallel for private(i) schedule(static)
  for (i=0; i<m_size; i++) {
    m_array_data[i] = other.m_array_data[i];
  }

  return *this;
}

template <class T>
bool
DataVectorAlt<T>::operator==(const DataVectorAlt& other) const
{
  assert(m_size >= 0);

  if (m_size!=other.m_size) {
    return false;
  }
  if (m_dim!=other.m_dim) {
    return false;
  }
  if (m_N!=other.m_N) {
    return false;
  }
  for (int i=0; i<m_size; i++) {
    if (m_array_data[i] != other.m_array_data[i]) {
      return false;
    }
  }
  return true;
}

template <class T>
bool
DataVectorAlt<T>::operator!=(const DataVectorAlt& other) const
{
  return !(*this==other);
}

template <class T>
void 
DataVectorAlt<T>::copyFromArrayToOffset(const WrappedArray& value, size_type offset, size_type copies)
{
  const DataTypes::ShapeType& tempShape=value.getShape();
  size_type len=DataTypes::noValues(tempShape);
  if (offset+len*copies>size())
  {
     std::ostringstream ss;
     ss << "Error - not enough room for that DataPoint at that offset. (";
     ss << "offset=" << offset << " + " << " len=" << len << " >= " << size();
     throw DataException(ss.str());
  }
  size_type si=0,sj=0,sk=0,sl=0;
  switch (value.getRank())
  {
  case 0:       
        for (size_type z=0;z<copies;++z)
        {
           m_array_data[offset+z]=value.getElt();
        }
        break;
  case 1:
        for (size_type z=0;z<copies;++z)
        {
           for (size_t i=0;i<tempShape[0];++i)
           {
              m_array_data[offset+i]=value.getElt(i);
           }
           offset+=len;
        }
        break;
  case 2:
        si=tempShape[0];
        sj=tempShape[1];
        for (size_type z=0;z<copies;++z)
        {
           for (size_type i=0;i<si;i++)
           {
              for (size_type j=0;j<sj;j++)
              {
                 m_array_data[offset+DataTypes::getRelIndex(tempShape,i,j)]=value.getElt(i,j);
              }
           }
           offset+=len;
        }
        break;
  case 3:
        si=tempShape[0];
        sj=tempShape[1];
        sk=tempShape[2];
        for (size_type z=0;z<copies;++z) 
        {
          for (size_type i=0;i<si;i++)
          {
            for (size_type j=0;j<sj;j++)
            {
              for (size_type k=0;k<sk;k++)
              {
                 m_array_data[offset+DataTypes::getRelIndex(tempShape,i,j,k)]=value.getElt(i,j,k);
              }
            }
          }
          offset+=len;
        }
        break;
  case 4:
        si=tempShape[0];
        sj=tempShape[1];
        sk=tempShape[2];
        sl=tempShape[3];
        for (size_type z=0;z<copies;++z)
        {
          for (size_type i=0;i<si;i++)
          {
            for (size_type j=0;j<sj;j++)
            {
              for (size_type k=0;k<sk;k++)
              {
                 for (size_type l=0;l<sl;l++)
                 {
                    m_array_data[offset+DataTypes::getRelIndex(tempShape,i,j,k,l)]=value.getElt(i,j,k,l);
                 }
              }
            }
          }
          offset+=len;
        }
        break;
  default:
        std::ostringstream oss;
        oss << "Error - unknown rank. Rank=" << value.getRank();
        throw DataException(oss.str());
  }
}

template <class T>
void
DataVectorAlt<T>::copyFromArray(const WrappedArray& value, size_type copies)
{
  DataTypes::ShapeType tempShape=value.getShape();
  DataVectorAlt<T>::size_type nelements=DataTypes::noValues(tempShape)*copies;
  if (m_array_data!=0)
  {
    free(m_array_data);
  }
  m_array_data=reinterpret_cast<T*>(malloc(sizeof(T)*nelements));
  m_size=nelements;     // total amount of elements
  m_dim=m_size;         // elements per sample
  m_N=1;                        // number of samples
  copyFromArrayToOffset(value,0,copies);
}



} // end of namespace
} // end of namespace

#endif // __ESCRIPT_DATAVECTORALT_H__
