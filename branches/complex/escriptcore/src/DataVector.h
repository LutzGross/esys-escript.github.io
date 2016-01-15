
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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


#if !defined escript_DataVectorTaipan_20050324_H
#define escript_DataVectorTaipan_20050324_H
#include "system_dep.h"

#include "esysUtils/EsysAssert.h"

#include <vector>
#include <iostream>
#include <fstream>

#include "DataTypes.h"

namespace escript {

class WrappedArray;

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
  copyFromArrayToOffset(const WrappedArray& value, size_type offset, size_type copies);


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



class ESCRIPT_DLL_API DataVectorAlt {

 public:

  //
  // The type of the elements stored in the vector.
  typedef double ElementType;

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
  DataVectorAlt(const DataVectorAlt& other);

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
  ~DataVectorAlt();

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
  copyFromArray(const escript::WrappedArray& value, size_type copies);

  void 
  copyFromArrayToOffset(const WrappedArray& value, size_type offset, size_type copies);


  /**
     \brief
     Return the number of elements in this DataVectorAlt.
  */
  inline
  size_type
  size() const;

  /**
     \brief
     DataVectorAlt assignment operator "=".
     Assign the given DataVectorAlt object to this.
  */
  DataVectorAlt&
  operator=(const DataVectorAlt& other);

  /**
     \brief
     DataVectorAlt equality comparison operator "==".
     Return true if the given DataVectorAlt is equal to this.
  */
  bool
  operator==(const DataVectorAlt& other) const;

  /**
     \brief
     DataVectorAlt inequality comparison operator "!=".
     Return true if the given DataVectorAlt is not equal to this.
  */
  bool
  operator!=(const DataVectorAlt& other) const;

  /**
    \brief
    Return a reference to the element at position i in this DataVectorAlt.
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
  std::vector<ElementType> m_array_data;
};



// This is the main version we had
//typedef DataVectorTaipan DataVector;

typedef DataVectorAlt DataVector;


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
  EsysAssert(i<size(),"DataVectorTaipan: invalid index specified. " << i << " of " << size());
  return m_array_data[i];
}

inline
DataVectorTaipan::const_reference
DataVectorTaipan::operator[](const DataVectorTaipan::size_type i) const
{
  EsysAssert(i<size(),"DataVectorTaipan: invalid index specified. " << i << " of " << size());
  return m_array_data[i];
}




inline
DataVectorAlt::size_type
DataVectorAlt::size() const
{
  return m_size;
}

inline
DataVectorAlt::reference
DataVectorAlt::operator[](const DataVectorAlt::size_type i)
{
  EsysAssert(i<size(),"DataVectorAlt: invalid index specified. " << i << " of " << size());
  return m_array_data[i];
}

inline
DataVectorAlt::const_reference
DataVectorAlt::operator[](const DataVectorAlt::size_type i) const
{
  EsysAssert(i<size(),"DataVectorAlt: invalid index specified. " << i << " of " << size());
  return m_array_data[i];
}


namespace DataTypes
{
  typedef escript::DataVector               FloatVectorType;//!< Vector to store underlying data.
  
  
  
  

  /**
     \brief
     Copy a data slice specified by the given region and offset from the
     "other" view into the "left" view at the given offset.
     
     \param left - vector to copy into
     \param leftShape - shape of datapoints for the left vector
     \param leftOffset - location within left to start copying to
     \param other - vector to copy from
     \param otherShape - shape of datapoints for the other vector
     \param otherOffset - location within other vector to start copying from
     \param region - Input -
                      Region in other view to copy data from.
  */
   ESCRIPT_DLL_API
   void
   copySlice(FloatVectorType& left,
			    const ShapeType& leftShape,
			    vec_size_type leftOffset,
                            const FloatVectorType& other,
			    const ShapeType& otherShape,
                            vec_size_type otherOffset,
                            const RegionLoopRangeType& region);

  /**
     \brief
     Copy data into a slice specified by the given region and offset in
     the left vector from the other vector at the given offset.

     \param left - vector to copy into
     \param leftShape - shape of datapoints for the left vector
     \param leftOffset - location within left to start copying to
     \param other - vector to copy from
     \param otherShape - shape of datapoints for the other vector
     \param otherOffset - location within other vector to start copying from
     \param region - Input -
                      Region in the left vector to copy data to.
  */
   ESCRIPT_DLL_API
   void
   copySliceFrom(FloatVectorType& left,
				const ShapeType& leftShape,
				vec_size_type leftOffset,
                                const FloatVectorType& other,
				const ShapeType& otherShape,
                                vec_size_type otherOffset,
                                const RegionLoopRangeType& region);


   /**
      \brief Display a single value (with the specified shape) from the data.

     Despite its similar name this function behaves differently to pointToString.
     There are no prefixes or (i,j,k) identifiers on each field. each datapoint is printed without
     new lines.
     It also works with double* rather than vectors so be careful what you pass it.

     \param os - stream to write to
     \param data - vector containing the datapoint
     \param shape - shape of the datapoint
     \param offset - start of the datapoint within data
     \param needsep - Does this output need to start with a separator
     \param sep - separator string to print between components
   */
   void
   pointToStream(std::ostream& os, const FloatVectorType::ElementType* data,const ShapeType& shape, int offset, bool needsep=true, const std::string& sep=",");

   /**
      \brief Display a single value (with the specified shape) from the data.

     \param data - vector containing the datapoint
     \param shape - shape of the datapoint
     \param offset - start of the datapoint within data
     \param prefix - string to prepend to the output
   */
   std::string
   pointToString(const FloatVectorType& data,const ShapeType& shape, int offset, const std::string& prefix);


   /**
      \brief  Copy a point from one vector to another. Note: This version does not check to see if shapes are the same.

   \param dest - vector to copy to
   \param doffset - beginning of the target datapoint in dest
   \param nvals - the number of values comprising the datapoint
   \param src - vector to copy from
   \param soffset - beginning of the datapoint in src
   */
   void copyPoint(FloatVectorType& dest, vec_size_type doffset, vec_size_type nvals, const FloatVectorType& src, vec_size_type soffset);  
  
  
  
  
}








  
} // end of namespace

#endif
