
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined escript_DataTagged_20040615_H
#define escript_DataTagged_20040615_H
#include "system_dep.h"

#include "DataReady.h"
#include "DataTypes.h"

#include <vector>
#include <map>

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

class DataTagged : public DataReady {
typedef DataReady parent;
 public:

  //
  // Types for the lists of tags and values.
  typedef std::vector<int>           TagListType;
//  typedef std::vector<DataArrayView> ValueListType;
  typedef DataTypes::ValueType   ValueType;
  typedef std::vector<ValueType::ElementType> ValueBatchType;

  //
  // Map from a tag to an offset into the data array. 
  typedef std::map<int, int> DataMapType;

  /**
     \brief
     Default constructor for DataTagged.

     Description:
     Default constructor for DataTagged. Creates a DataTagged object for which
     the default data-point is a scalar data-point with value 0.0, and no other
     tag values are stored.
    T
  */
  ESCRIPT_DLL_API
  DataTagged();

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
  ESCRIPT_DLL_API
  DataTagged(const FunctionSpace& what,
             const DataTypes::ShapeType &shape,
             const int tags[],
             const ValueType& data);

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
  ESCRIPT_DLL_API
  DataTagged(const FunctionSpace& what,
             const DataTypes::ShapeType &shape,
             const TagListType& tags,
             const ValueType& data);

  /**
     \brief
     Copy Constructor for DataTagged.
     Performs a deep copy from the given DataTagged object.
    T
  */
  ESCRIPT_DLL_API
  DataTagged(const DataTagged& other);

  /**
     \brief
     Copy Constructor for DataTagged.
     Construct a DataTagged object from a DataConstant object.
     The default value will be the value of the DataConstant object.
    T
  */
  ESCRIPT_DLL_API
  DataTagged(const DataConstant& other);

  /**
     \brief
     Copies the tags from a DataTagged into a new Data Tagged and assigns them the default value. ** Not unit tested **

     This is different from a deep copy because we are not copying shape or other information, just tags.
     \param what - Input - FunctionSpace for the new DataTagged
     \param shape - Input - Shape for points in the new DataTagged
     \param defaultvalue - Input - Default value for new DataTagged
     \param tagsource - Input - A DataTagged object which supplies the tags. 
  */
  ESCRIPT_DLL_API
  DataTagged(const FunctionSpace& what,
             const DataTypes::ShapeType& shape,
	     const DataTypes::ValueType& defaultvalue,
             const DataTagged* tagsource=0);

  /**
     \brief
     Destructor
  */
  ESCRIPT_DLL_API
  inline virtual
  ~DataTagged() {};

  ESCRIPT_DLL_API
  bool
  isTagged() const 
  {
    return true;
  };


  /**
     \brief Return a deep copy of the current object.
  */
  ESCRIPT_DLL_API
  virtual
  DataAbstract*
  deepCopy();


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
  ESCRIPT_DLL_API
  virtual
  double*
  getSampleDataByTag(int tag);

  /**
     \brief
     Write the data as a string.
     Writes out each tag, including the default, and the data-point which is
     associated with each tag.
    T
  */
  ESCRIPT_DLL_API
  virtual
  std::string
  toString() const;
 /**
     \brief
     dumps the object into a netCDF file
  */
  ESCRIPT_DLL_API
  virtual
  void
  dump(const std::string fileName) const;

 /**
     \brief
    sets all values to zero
  */
  ESCRIPT_DLL_API
  virtual
  void
  setToZero();

  /**
     \brief
     Return the tag number associated with the given data-point number
     according to the associated function space.
    T
  */
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
  virtual
  ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

  ESCRIPT_DLL_API
  virtual
  ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo);

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
  ESCRIPT_DLL_API
  void
  addTaggedValues(const TagListType& tagKeys,
                            const ValueBatchType& values,
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
  ESCRIPT_DLL_API
  void
  addTaggedValues(const TagListType& tagKeys,
                            const ValueType& values,
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
  ESCRIPT_DLL_API
  void
  addTaggedValue(int tagKey,
		 const DataTypes::ShapeType& pointshape,
                 const ValueType& value,
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
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
  void
  setTaggedValue(int tagKey,
		 const DataTypes::ShapeType& pointshape,
                 const ValueType& value,
		 int dataOffset=0);

  /**
     \brief
     getDataByTag

     Return a pointer to the beginning of the datapoint with the specified tag.
     TODO Eventually these should be inlined.
     \param tag - Input - Integer key.
     \param i - position in the underlying datastructure
  */

  ESCRIPT_DLL_API
  DataTypes::ValueType::reference
  getDataByTagRW(int tag, DataTypes::ValueType::size_type i);

  ESCRIPT_DLL_API
  DataTypes::ValueType::const_reference
  getDataByTagRO(int tag, DataTypes::ValueType::size_type i) const;



  /**
      \brief 
      getOffsetForTag

      \param tag
      \return the offset of the beginning of the datapoint corresponding to tag.

      Note: If the tag is not valid, the offset of the default value is returned instead.
  */
  ESCRIPT_DLL_API
  DataTypes::ValueType::size_type
  getOffsetForTag(int tag) const;


  /**
     \brief
     Return a reference to the underlying DataVector.
  */

  ESCRIPT_DLL_API
  DataTypes::ValueType&
  getVectorRW();

  ESCRIPT_DLL_API
  const DataTypes::ValueType&
  getVectorRO() const;



  /**
     \brief 
     getTagLookup

     Description:
     Return a reference to the tag offset lookup table.
    T
  */
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
  DataTypes::ValueType::reference
  getDefaultValueRW(DataTypes::ValueType::size_type i);

  ESCRIPT_DLL_API
  DataTypes::ValueType::const_reference
  getDefaultValueRO(DataTypes::ValueType::size_type i) const;





  /**
     \brief
     getLength

     Description:
     Return the total number of doubles stored for this DataTagged object.
    T
  */
  ESCRIPT_DLL_API
  virtual
  ValueType::size_type
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
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
  virtual
  void
  setSlice(const DataAbstract* other,
           const DataTypes::RegionType& region);


  /**
     \brief
     Computes a symmetric matrix (A + AT) / 2

     \param ev - Output - symmetric matrix

  */
  ESCRIPT_DLL_API
  virtual void
  symmetric(DataAbstract* ev);

  /**
     \brief
     Computes a nonsymmetric matrix (A - AT) / 2

     \param ev - Output - nonsymmetric matrix

  */
  ESCRIPT_DLL_API
  virtual void
  nonsymmetric(DataAbstract* ev);

  /**
     \brief
     Computes the trace of a matrix

     \param ev - Output - the trace of a matrix
     \param axis_offset
  */
  ESCRIPT_DLL_API
  virtual void
  trace(DataAbstract* ev, int axis_offset);

  /**
     \brief
     swaps components axis0 and axis1

     \param ev - Output - swapped components
     \param axis0
     \param axis1
  */
  ESCRIPT_DLL_API
  virtual void
  swapaxes(DataAbstract* ev, int axis0, int axis1);

  /**
     \brief
     Transpose each data point of this Data object around the given axis.

     \param ev - Output - the transpose of a matrix
     \param axis_offset
  */
  ESCRIPT_DLL_API
  virtual void
  transpose(DataAbstract* ev, int axis_offset);

  /**
     \brief
     solves the eigenvalue problem this*V=ev*V for the eigenvalues ev

     \param ev - Output - eigenvalues in increasing order at each data point

  */
  ESCRIPT_DLL_API
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

  ESCRIPT_DLL_API
  virtual void
  eigenvalues_and_eigenvectors(DataAbstract* ev,DataAbstract* V,const double tol=1.e-13);


  /**
     \brief  Returns the offset in the structure which stores the default value
  */
  ESCRIPT_DLL_API
  DataTypes::ValueType::size_type
  getDefaultOffset() const;

 protected:

 private:

  //
  // The offset lookup table
  DataMapType m_offsetLookup;

  //
  // the offset to the default value
  static const int m_defaultValueOffset = 0;

  //
  // The actual data
  ValueType m_data;

};

inline
bool
DataTagged::isCurrentTag(int tag) const
{
  DataMapType::const_iterator pos(m_offsetLookup.find(tag));
  return (pos!=m_offsetLookup.end());
}

inline 
DataTypes::ValueType::size_type
DataTagged::getDefaultOffset() const
{
  return m_defaultValueOffset;	
}

inline
DataTypes::ValueType::reference
DataTagged::getDefaultValueRW(DataTypes::ValueType::size_type i)
{	
	return getVectorRW()[i];		// getVectorRW has exclusive write checks
}

inline
DataTypes::ValueType::const_reference
DataTagged::getDefaultValueRO(DataTypes::ValueType::size_type i) const
{
	return getVectorRO()[i];
}

inline
const DataTagged::DataMapType&
DataTagged::getTagLookup() const
{
  return m_offsetLookup;
}

inline
DataTypes::ValueType::size_type
DataTagged::getLength() const
{
  return m_data.size();
}

} // end of namespace

#endif
