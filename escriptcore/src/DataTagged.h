
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#ifndef __ESCRIPT_DATATAGGED_H__
#define __ESCRIPT_DATATAGGED_H__

#include "system_dep.h"

#include "DataReady.h"
#include "DataTypes.h"

#include <map>
#include <vector>

namespace escript {

class DataConstant;

/**
   \brief
   Simulates a full dataset accessible via sampleNo and dataPointNo. 

   Description:
   Each data-point has an associated tag number, and a given tag can represent a
   range of dataPointNo and sampleNo. Each tag indexes only a single data-point.
   Thus only a single data-point needs to be stored for a range of sampleNo and
   dataPointNo values.
*/

class ESCRIPT_DLL_API DataTagged : public DataReady
{
    typedef DataReady parent;
 public:

  //
  // Types for the lists of tags and values.
  typedef std::vector<int> TagListType;
  typedef std::vector<DataTypes::RealVectorType::ElementType> FloatBatchType;
  typedef std::vector<DataTypes::CplxVectorType::ElementType> CplxBatchType;

  //
  // Map from a tag to an offset into the data array. 
  typedef std::map<int, int> DataMapType;

  /**
     \brief
     Alternative Constructor for DataTagged.

     Description:
     Alternative Constructor for DataTagged.
     \param what - Input - A description of what this data object represents.
     \param shape - Input - The shape of each data-point.
     \param tags - Input - An array of tags, one for each sample number (starts at tag[1]).
     \param data - The data values for each tag.
    NB: no unit testing yet
  */
  explicit DataTagged(const FunctionSpace& what,
             const DataTypes::ShapeType &shape,
             const int tags[],
             const DataTypes::RealVectorType& data);
  
  
  explicit DataTagged(const FunctionSpace& what,
             const DataTypes::ShapeType &shape,
             const int tags[],
             const DataTypes::CplxVectorType& data);  
  

 /**
     \brief
     Alternative Constructor for DataTagged.

     Description:
     Alternative Constructor for DataTagged.
     \param what - Input - A description of what this data object represents.
     \param shape - Input - The shape of each data-point.
     \param tags - Input - An vector of tags, one for each sample number.
     \param data - The data values for each tag.
TODO Make sure to document the relationship between tags and data, ie: data also contains the default value
 */
  explicit DataTagged(const FunctionSpace& what,
             const DataTypes::ShapeType &shape,
             const TagListType& tags,
             const DataTypes::RealVectorType& data);
  
  explicit DataTagged(const FunctionSpace& what,
             const DataTypes::ShapeType &shape,
             const TagListType& tags,
             const DataTypes::CplxVectorType& data);  
  

  /**
     \brief
     Copy Constructor for DataTagged.
     Performs a deep copy from the given DataTagged object.
    T
  */
  DataTagged(const DataTagged& other);

  /**
     \brief
     Copy Constructor for DataTagged.
     Construct a DataTagged object from a DataConstant object.
     The default value will be the value of the DataConstant object.
    T
  */
  explicit DataTagged(const DataConstant& other);

  /**
     \brief
     Copies the tags from a DataTagged into a new Data Tagged and assigns them the default value. ** Not unit tested **

     This is different from a deep copy because we are not copying shape or other information, just tags.
     \param what - Input - FunctionSpace for the new DataTagged
     \param shape - Input - Shape for points in the new DataTagged
     \param defaultvalue - Input - Default value for new DataTagged
     \param tagsource - Input - A DataTagged object which supplies the tags. 
  */
  explicit DataTagged(const FunctionSpace& what,
             const DataTypes::ShapeType& shape,
             const DataTypes::RealVectorType& defaultvalue,
             const DataTagged* tagsource=0);

  explicit DataTagged(const FunctionSpace& what,
             const DataTypes::ShapeType& shape,
             const DataTypes::CplxVectorType& defaultvalue,
             const DataTagged* tagsource=0);  
  
  
  /**
     \brief
     Destructor
  */
  inline virtual
  ~DataTagged() {};

  bool
  isTagged() const 
  {
    return true;
  };

  /**
  \brief Return true if any one of the datapoints contains a NaN.
  */
  bool
  hasNaN() const;

  /**
  \brief replaces all NaN values with value 
  */
  void
  replaceNaN(DataTypes::real_t value);
  
  /**
  \brief replaces all NaN values with value 
  */
  void
  replaceNaN(DataTypes::cplx_t value);

  /**
   \brief Return true if data contains Inf or -Inf 
  */
  // ESCRIPT_DLL_API class already exported
  virtual bool
  hasInf() const;

  /**
  \brief replaces all (+/-)Inf values with value 
  */
  // ESCRIPT_DLL_API class already exported
  virtual void
  replaceInf(DataTypes::real_t value);
  
  /**
  \brief replaces all (+/-)Inf values with value 
  */
  // ESCRIPT_DLL_API class already exported
  virtual void
  replaceInf(DataTypes::cplx_t value);    
  
  
  /**
     \brief Return a deep copy of the current object.
  */
  virtual
  DataAbstract*
  deepCopy() const;

  
  /**
     \brief Return an object with the same type, domain (and tags if appropriate)
     as this, but all values are zeroed.
  */  
  // ESCRIPT_DLL_API class already exported
  virtual
  DataAbstract*
  zeroedCopy() const;  
  
  
  /**
     \brief
     getSampleDataByTag

     Description:
     Return the data-point for the given tag. All of the data for the
     sample will be visible via the returned pointer.

     ** This provides an interface into the data suitable for legacy C code.
     ** NB: need to do array bounds checking when accessing returned value!
    T
  */
  virtual
  DataTypes::real_t*
  getSampleDataByTag(int tag, DataTypes::real_t dummy=0);  

  virtual
  DataTypes::cplx_t*
  getSampleDataByTag(int tag, DataTypes::cplx_t dummy);  
  
  
  /**
     \brief
     Write the data as a string.
     Writes out each tag, including the default, and the data-point which is
     associated with each tag.
    T
  */
  virtual
  std::string
  toString() const;
 /**
     \brief
     dumps the object into am HDF5 file
  */
  #ifdef ESYS_HAVE_HDF5
  virtual
  void
  dump_hdf5(const H5::Group h5_grp) const;
  #endif


  /**
    \brief invert square matricies
    \param out - Where to store the results
    \return errorcode (0 indicates success)
  */
  virtual int
  matrixInverse(DataAbstract* out) const;

 /**
     \brief
    sets all values to zero
  */
  virtual
  void
  setToZero();

  /**
     \brief
     Return the tag number associated with the given data-point number
     according to the associated function space.
    T
  */
  virtual
  int
  getTagNumber(int dpno);

  /**
     \brief
     getPointOffset

     Description:
     Return the offset to the given data-point value in the underlying
     data vector.

     \param sampleNo - Input - sample number.
     \param dataPointNo - Input - data-point number.
    T
  */
  virtual
  DataTypes::RealVectorType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

 /**
     \brief
     addTaggedValues

     Description:
     Add the given tags and values to this DataTagged object.
     \param tagKeys - Input - A vector of integer tags.
     \param values - Input - A vector of doubles. If this is empty, the default value for
                     this DataTagged will be used for all tags.
                     If it contains one value all tag values will be assigned
                     this value. Otherwise consecutive tags will be assigned 
                     consecutive values.  If there is a mismatch between  the
                     number of keys and the number of values an exception 
                     will be generated.
     \param vShape - shape of the datapoints in "values"
    T
 */
  void
  addTaggedValues(const TagListType& tagKeys,
                            const FloatBatchType& values,
                            const ShapeType& vShape);


  /**
   Description: 
   Add the given tags and values to this DataTagged object.
   \param tagKeys - Input - A vector of integer tags.
   \param values - Input - A DataVector containing the datapoints.
                     If this is empty, the default value for
                     this DataTagged will be used for all tags.
                     If it contains one value all tag values will be assigned
                     this value. Otherwise consecutive tags will be assigned 
                     consecutive values.  If there is a mismatch between  the
                     number of keys and the number of values an exception 
                     will be generated.
    \param vShape - shape of the datapoints in "values"

   TODO Makesure this is properly unit tested
  */
  void
  addTaggedValues(const TagListType& tagKeys,
                            const DataTypes::RealVectorType& values,
                            const ShapeType& vShape);




  /**
     \brief
     addTaggedValue

     Description:
     Add a single tag and value to this DataTagged object. If this tag already has
     a value associated with it, setTaggedValue will be used to update this value.
     \param tagKey - Input - Integer tag.
     \param pointshape - Shape of the value parameter
     \param value - Input - Single DataArrayView value to be assigned to the tag. 
     \param dataOffset - Input - Offset of the beginning of the point in the value parameter 
  */
  void
  addTaggedValue(int tagKey,
                 const DataTypes::ShapeType& pointshape,
                 const DataTypes::RealVectorType& value,
                 int dataOffset=0);
  
  void
  addTaggedValue(int tagKey,
                 const DataTypes::ShapeType& pointshape,
                 const DataTypes::CplxVectorType& value,
                 int dataOffset=0);  
  

  /**
     \brief
     addTag - does not modify the default value for this object. ** Not unit tested **

     Description:
     Add a single tag. The default value for this DataTagged will be associated with the tag.
     If this tag already has a value associated with it, then no change will be made.
     \param tagKey - Input - Integer tag.
    TODO: Make sure this is unit tested
  */
  void
  addTag(int tagKey);

  /**
     \brief
     setTaggedValue

     Description:
     Assign the given value to the given tag.
     \param tagKey - Input - Integer tag.
     \param pointshape - the shape of the value parameter
     \param value - Input - Vector storing the datapoint to be assigned to the tag.
     \param dataOffset - beginning of the datapoint within "value".
    T
  */
  void
  setTaggedValue(int tagKey,
                 const DataTypes::ShapeType& pointshape,
                 const DataTypes::RealVectorType& value,
                 int dataOffset=0);
  
  void
  setTaggedValue(int tagKey,
                 const DataTypes::ShapeType& pointshape,
                 const DataTypes::CplxVectorType& value,
                 int dataOffset=0);  
  

  /**
     \brief
     getDataByTag

     Return a pointer to the beginning of the datapoint with the specified tag.
     TODO Eventually these should be inlined.
     \param tag - Input - Integer key.
     \param i - position in the underlying datastructure
  */

  DataTypes::RealVectorType::reference
  getDataByTagRW(int tag, DataTypes::RealVectorType::size_type i, DataTypes::real_t dummy=0);

  DataTypes::RealVectorType::const_reference
  getDataByTagRO(int tag, DataTypes::RealVectorType::size_type i, DataTypes::real_t dummy=0) const;


  DataTypes::CplxVectorType::reference
  getDataByTagRW(int tag, DataTypes::CplxVectorType::size_type i, DataTypes::cplx_t dummy);

  DataTypes::CplxVectorType::const_reference
  getDataByTagRO(int tag, DataTypes::CplxVectorType::size_type i, DataTypes::cplx_t dummy) const;

  /**
      \brief 
      getOffsetForTag

      \param tag
      \return the offset of the beginning of the datapoint corresponding to tag.

      Note: If the tag is not valid, the offset of the default value is returned instead.
  */
  DataTypes::RealVectorType::size_type
  getOffsetForTag(int tag) const;


  /**
     \brief
     Return a reference to the underlying DataVector.
  */

  DataTypes::RealVectorType&
  getVectorRW();

  const DataTypes::RealVectorType&
  getVectorRO() const;


  DataTypes::CplxVectorType&
  getVectorRWC();

  const DataTypes::CplxVectorType&
  getVectorROC() const;
  

  virtual DataTypes::RealVectorType&
  getTypedVectorRW(DataTypes::real_t dummy);  
  
  virtual const DataTypes::RealVectorType&
  getTypedVectorRO(DataTypes::real_t dummy) const;

  virtual DataTypes::CplxVectorType&
  getTypedVectorRW(DataTypes::cplx_t dummy);
  
  virtual const DataTypes::CplxVectorType&
  getTypedVectorRO(DataTypes::cplx_t dummy) const;  

  
  


  /**
     \brief 
     getTagLookup

     Description:
     Return a reference to the tag offset lookup table.
    T
  */
  const DataMapType&
  getTagLookup() const;

  /**
     \brief
     isCurrentTag

     Description:
     Return true if the given tag exists within the DataTagged tag map.

     *** NB: The DataTagged tag map does not necessarily coincide with the tag
     keys in the associated function space.
    T
  */
  bool
  isCurrentTag(int tag) const;

  /**
     \brief
     getDefaultValue

     Description:
     Return the default value. This value is associated with any tag which
     is not explicitly recorded in this DataTagged object's tag map.
     \param i - position in the underlying datastructure
  */
  DataTypes::RealVectorType::reference
  getDefaultValueRW(DataTypes::RealVectorType::size_type i, DataTypes::real_t dummy=0);

  DataTypes::RealVectorType::const_reference
  getDefaultValueRO(DataTypes::RealVectorType::size_type i, DataTypes::real_t dummy=0) const;

  DataTypes::CplxVectorType::reference
  getDefaultValueRW(DataTypes::CplxVectorType::size_type i, DataTypes::cplx_t dummy);

  DataTypes::CplxVectorType::const_reference
  getDefaultValueRO(DataTypes::CplxVectorType::size_type i, DataTypes::cplx_t dummy) const;




  /**
     \brief
     getLength

     Description:
     Return the total number of doubles stored for this DataTagged object.
    T
  */
  virtual
  DataTypes::RealVectorType::size_type
  getLength() const;

  /**
     \brief
     getSlice

     Description:
     Factory method that returns a newly created DataTagged object generated
     by taking the specified slice from this DataTagged object.
     The caller is reponsible for managing the returned object.
    T
  */
  virtual
  DataAbstract*
  getSlice(const DataTypes::RegionType& region) const;

  /**
     \brief
     Slice Constructor for DataTagged.

     Description:
     Creates a DataTagged object which is the specified slice
     from the given DataTagged object.
     \param other - Input - DataTagged object to slice from.
     \param region - Input - Region to slice.
    T
  */
  DataTagged(const DataTagged& other, 
             const DataTypes::RegionType& region);

  /**
     \brief
     setSlice

     Description:
     Copy the given Data object into the specified region in this object.
     \param other - Input - Data object to copy from.
     \param region - Input - Region to copy into (NB: must have same shape as other!).
    T
  */
  virtual
  void
  setSlice(const DataAbstract* other,
           const DataTypes::RegionType& region);


  /**
     \brief
     Computes a symmetric matrix (A + AT) / 2

     \param ev - Output - symmetric matrix

  */
  virtual void
  symmetric(DataAbstract* ev);

  /**
     \brief
     Computes a antisymmetric matrix (A - AT) / 2

     \param ev - Output - antisymmetric matrix

  */
  virtual void
  antisymmetric(DataAbstract* ev);

  /**
     \brief
     Computes an hermitian matrix (A + A*) / 2

     \param ev - Output - hermitian matrix

  */
  virtual void
  hermitian(DataAbstract* ev);

  /**
     \brief
     Computes an antihermitian matrix (A - A*) / 2

     \param ev - Output - anti-hermitian matrix

  */
  virtual void
  antihermitian(DataAbstract* ev);

  /**
     \brief
     Computes the trace of a matrix

     \param ev - Output - the trace of a matrix
     \param axis_offset
  */
  virtual void
  trace(DataAbstract* ev, int axis_offset);

  /**
     \brief
     swaps components axis0 and axis1

     \param ev - Output - swapped components
     \param axis0
     \param axis1
  */
  virtual void
  swapaxes(DataAbstract* ev, int axis0, int axis1);

  /**
     \brief
     Transpose each data point of this Data object around the given axis.

     \param ev - Output - the transpose of a matrix
     \param axis_offset
  */
  virtual void
  transpose(DataAbstract* ev, int axis_offset);

  /**
     \brief
     solves the eigenvalue problem this*V=ev*V for the eigenvalues ev

     \param ev - Output - eigenvalues in increasing order at each data point

  */
  virtual void
  eigenvalues(DataAbstract* ev);

  /**
     \brief
     solves the eigenvalue problem this*V=ev*V for the eigenvalues ev and eigenvectors V

     \param ev - Output - eigenvalues in increasing order at each data point
     \param V - Output - corresponding eigenvectors. They are normalized such that their length is one
                         and the first nonzero component is positive.
     \param tol - Input - eigenvalue with relative distance tol are treated as equal.

  */

  virtual void
  eigenvalues_and_eigenvectors(DataAbstract* ev,DataAbstract* V,const double tol=1.e-13);


  /**
     \brief  Returns the offset in the structure which stores the default value
  */
  DataTypes::RealVectorType::size_type
  getDefaultOffset() const;
  
  /**
   \brief Return the number of tags which have been given values (+the default) 
  */ 
  size_t
  getTagCount() const;
  
  void
  complicate();
  
 protected:

 private:

  //
  // The offset lookup table
  DataMapType m_offsetLookup;

  //
  // the offset to the default value
  static const int m_defaultValueOffset = 0;
  
  // the actual data
  DataTypes::RealVectorType m_data_r;
  DataTypes::CplxVectorType m_data_c;  
  

};

inline
bool
DataTagged::isCurrentTag(int tag) const
{
  DataMapType::const_iterator pos(m_offsetLookup.find(tag));
  return (pos!=m_offsetLookup.end());
}

inline 
DataTypes::RealVectorType::size_type
DataTagged::getDefaultOffset() const
{
  return m_defaultValueOffset;  
}

inline
DataTypes::RealVectorType::reference
DataTagged::getDefaultValueRW(DataTypes::RealVectorType::size_type i, DataTypes::real_t dummy)
{       
        return getVectorRW()[i];                // getVectorRW has exclusive write checks
}

inline
DataTypes::RealVectorType::const_reference
DataTagged::getDefaultValueRO(DataTypes::RealVectorType::size_type i, DataTypes::real_t dummy) const
{
        return getVectorRO()[i];
}

inline
DataTypes::CplxVectorType::reference
DataTagged::getDefaultValueRW(DataTypes::RealVectorType::size_type i, DataTypes::cplx_t dummy)
{       
        return getVectorRWC()[i];                // getVectorRW has exclusive write checks
}

inline
DataTypes::CplxVectorType::const_reference
DataTagged::getDefaultValueRO(DataTypes::CplxVectorType::size_type i, DataTypes::cplx_t dummy) const
{
        return getVectorROC()[i];
}




inline
const DataTagged::DataMapType&
DataTagged::getTagLookup() const
{
  return m_offsetLookup;
}

inline
DataTypes::RealVectorType::size_type
DataTagged::getLength() const
{
  return std::max(m_data_c.size(), m_data_r.size());
}

} // end of namespace

#endif // __ESCRIPT_DATATAGGED_H__

