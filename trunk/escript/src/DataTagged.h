// $Id$

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

#if !defined escript_DataTagged_20040615_H
#define escript_DataTagged_20040615_H

#include "DataAbstract.h"
#include "DataArrayView.h"

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

class DataTagged : public DataAbstract {

 public:

  //
  // Types for the lists of tags and values.
  typedef std::vector<int>           TagListType;
  typedef std::vector<DataArrayView> ValueListType;
  typedef DataArrayView::ValueType   ValueType;

  //
  // Map from a tag to an offset into the data array. 
  typedef std::map<int, int> DataMapType;

  /**
     \brief
     Default constructor for DataTagged.

     Description:
     Default constructor for DataTagged. Creates a DataTagged object for which
     the only data-point is a scalar data-point with value 0.0. All tags
     will map to this single data-point.
     T
  */
  DataTagged();

  /**
     \brief
     Constructor for DataTagged.

     Description:
     Constructor for DataTagged.
     \param tagKeys - Input - A vector of integer keys.
     \param values - Input - A vector of DataArrayViews. If this is empty
                     all tag values will be assigned a value of zero. If 
                     it contains one value all tag values will be assigned the 
		     same value. Otherwise if there is a mismatch between 
		     the number of keys and the number of values an exception 
		     will be generated.
     \param defaultValue - Input - Value returned if a requested tag doesn't exist.
     \param what - Input - A description of what this data represents.
     T
  */
  DataTagged(const TagListType& tagKeys,
             const ValueListType& values,
	     const DataArrayView& defaultValue,
	     const FunctionSpace& what);

  /**
     \brief
     Alternative Constructor for DataTagged.

     Description:
     Alternative Constructor for DataTagged.
     \param what - Input - A description of what this data object represents.
     \param shape - Input - The shape of each data-point.
     \param tags - Input - An array of tags, one for each sample number.
     \param data - The data values for each tag.
  */
  DataTagged(const FunctionSpace& what,
             const DataArrayView::ShapeType &shape,
             const int tags[],
             const ValueType& data);

  /**
     \brief
     Slice Constructor for DataTagged.

     Description:
     Slice Constructor for DataTagged.
     Copies a slice from another DataTagged object.
     \param other - Input - DataTagged object to copy from.
     \param region - Input - Region to copy.
  */
  DataTagged(const DataTagged& other, 
	     const DataArrayView::RegionType& region);

  /**
     \brief
     Copy Constructor for DataTagged.
     Performs a deep copy from the given DataTagged object.
  */
  DataTagged(const DataTagged& other);

  /**
     \brief
     Copy Constructor for DataTagged.
     Construct a tagged data from a DataConstant object.
     The default value will be that held by the DataConstant object.
  */
  DataTagged(const DataConstant& other);

  /**
     \brief
     getSampleDataByTag

     Description:
     Return the data-point for the given tag. All of the data for the
     sample will be visible via the returned pointer.

     ** This provides an interface into the data suitable for legacy C code.
  */
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
  virtual
  std::string
  toString() const;

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
  ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

  /**
     \brief
     addTaggedValues

     Description:
     Add the given tags and values to this DataTagged object, by repeatedly
     using addTaggedValue for each given tag/value pair.
     \param tagKeys - Input - A vector of integer keys.
     \param values - Input - A vector of DataArrayViews. If this is empty
                      then all given tags will be assigned a value of zero. If 
                      it contains one value all tags will be assigned the same value.
                      Otherwise if there is a mismatch between the number of tags and 
                      the number of values an exception will be generated.
  */
  void
  addTaggedValues(const TagListType& tagKeys,
                  const ValueListType& values);  

  /**
     \brief
     addTaggedValue

     Description:
     Add a single tag and value to this DataTagged object. If this tag already has
     a value associated with it, setTaggedValue will be used to update this value.
     \param tagKey - Input - Integer key.
     \param value - Input - Single DataArrayView value to be assigned to the tag.
  */
  void
  addTaggedValue(int tagKey,
                 const DataArrayView& value);

  /**
     \brief
     setTaggedValues

     Description:
     Set the given tags to the given values in this DataTagged object, by repeatedly
     using setTaggedValue for each given tag/value pair.
     \param tagKeys - Input - A vector of integer keys.
     \param values - Input - A vector of DataArrayViews. If this is empty
                      then all given tags will be assigned a value of zero. If 
                      it contains one value all tag values will be assigned the same value.
                      Otherwise if there is a mismatch between the number of keys and 
                      the number of values an exception will be generated.
  */
  void
  setTaggedValues(const TagListType& tagKeys,
                  const ValueListType& values); 

  /**
     \brief
     setTaggedValue

     Description:
     Assign the given value to the given tag. If this tag does not already have a value
     associated with it, addTaggedValue will be used to add this tag/value pair.
     \param tagKey - Input - Integer key.
     \param value - Input - Single DataArrayView value to be assigned to the tag.
  */
  virtual
  void
  setTaggedValue(int tagKey,
                 const DataArrayView& value);

  /**
     \brief
     getDataPointByTag

     Description:
     Return a view into the data-point associated with the given tag.
     \param tag - Input - Integer key.
     T
  */
  DataArrayView
  getDataPointByTag(int tag) const;

  /**
     \brief
     getDataPoint

     Description:
     Return a view into the data-point specified by the given sample
     and data-point numbers.
     \param sampleNo - Input.
     \param dataPointNo - Input.
     T
  */
  virtual
  DataArrayView
  getDataPoint(int sampleNo,
               int dataPointNo);

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

     NOTE: The DataTagged tag map does not necessarily coincide with the tag
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
     T
  */
  DataArrayView&
  getDefaultValue();

  const DataArrayView&
  getDefaultValue() const;

  /**
     \brief
     getLength

     Description:
     Return the number of doubles stored for the Data.
     T
  */
  virtual
  ValueType::size_type
  getLength() const;

  /**
     \brief
     getSlice

     Description:
     Factory method that returns a newly created DataTagged object.
     The caller is reponsible for managing the object created.
  */
  virtual
  DataAbstract*
  getSlice(const DataArrayView::RegionType& region) const;

  /**
     \brief
     setSlice

     Description:
     Copy the specified region from the given value into this object.
     \param value - Input - Data to copy from.
     \param region - Input - Region to copy.
  */
  virtual
  void
  setSlice(const DataAbstract* value,
           const DataArrayView::RegionType& region);

  /**
     \brief
     reshapeDataPoint

     Description:
     Reshape the data point only if the data-point is currently rank 0.
     An exception is thrown if the data-point has rank other than 0.
     The original data point value is used for all values of the new
     data point.
  */
  void
  reshapeDataPoint(const DataArrayView::ShapeType& shape);

  /**
    \brief
    Archive the underlying data values to the file referenced
    by ofstream. A count of the number of values expected to be written
    is provided as a cross-check.

    The return value indicates success (0) or otherwise (1).
  */
  int
  archiveData(std::ofstream& archiveFile,
              const DataArrayView::ValueType::size_type noValues) const;

  /**
    \brief
    Extract the number of values specified by noValues from the file
    referenced by ifstream to the underlying data structure.

    The return value indicates success (0) or otherwise (1).
  */
  int
  extractData(std::ifstream& archiveFile,
              const DataArrayView::ValueType::size_type noValues);

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
DataArrayView&
DataTagged::getDefaultValue()
{
  // The default value is always the first value.
  return getPointDataView();
}

inline
const DataArrayView&
DataTagged::getDefaultValue() const
{
  // The default value is always the first value.
  return getPointDataView();
}

inline
const DataTagged::DataMapType&
DataTagged::getTagLookup() const
{
  return m_offsetLookup;
}

inline
DataArrayView::ValueType::size_type
DataTagged::getLength() const
{
  return m_data.size();
}

} // end of namespace

#endif
