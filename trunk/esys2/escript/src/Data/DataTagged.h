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

#include "escript/Data/DataAbstract.h"
#include "escript/Data/DataArray.h"
#include "escript/Data/DataArrayView.h"

#include <vector>
#include <map>

namespace escript {

class DataConstant;

/**
   \brief
   Creates the illusion of a full dataset accessible via sampleNo and dataPointNo. 

   Description:
   Creates the illusion of a full dataset accessible via sampleNo and
   dataPointNo. In reality a much smaller number of data-points is stored.
   Each data-point has an associated key, thus a given key represents a specific
   range of dataPointNo and sampleNo. Each key indexes a single data-point.
   Thus only a single data-point needs to be stored for a range of sampleNo and
   dataPointNo values.
*/

class DataTagged : public DataAbstract {

 public:

  typedef std::vector<int> TagListType;
  typedef std::vector<DataArrayView> ValueListType;

  //
  // Map type, maps the key to an array offset. 
  typedef std::map<int, int> DataMapType;

  /**
     \brief
     Default constructor for DataTagged.

     Description:
     Default constructor for DataTagged. Creates a DataTagged object for which
     the default data-point is a scalar data-point with value 0.0. Any given tag
     will map to this data-point.
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
  */
  DataTagged(const TagListType& tagKeys,
             const ValueListType& values,
	     const DataArrayView& defaultValue,
	     const FunctionSpace& what);

  /**
     \brief
     Slice constructor for DataTagged.

     Description:
     Slice constructor for DataTagged.
     Copies a slice from another DataTagged object.
     \param other - Input - DataTagged object to copy from.
     \param region - Input - region to copy.
  */
  DataTagged(const DataTagged& other, 
	     const DataArrayView::RegionType& region);

  /**
     \brief
     Copy constructorfor DataTagged.
     Performs a deep copy from the given DataTagged object.
  */
  DataTagged(const DataTagged& other);

  /**
     \brief
     Construct a tagged data from a DataConstant object.
     The default data-point will be that held by the DataConstant object.
  */
  DataTagged(const DataConstant& other);

  /**
     \brief
     getSampleDataByTag

     Description:
     Return the data-point for the given tag. All of the data for the entire
     sample should be visable via the returned pointer. This provides an 
     interface into the data suitable for legacy C code.
  */
  virtual
  double*
  getSampleDataByTag(int tag);

  /**
     \brief
     Write the data as a string.
     Writes out each tag, including the default, and the data-point which is
     associated with each tag.
  */
  virtual
  std::string
  toString() const;

  /**
     \brief
     getPointOffset

     Description:
     Return the offset to the given data-point. This is somewhat artificial,
     but returns the offset for the given point in the DataTagged object.
     Only really necessary to avoid many DataArrayView objects.

     \param sampleNo - Input - sample number.
     \param dataPointNo - Input - data-point number.
   */
  virtual
  DataArrayView::ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

  /**
     \brief
     addTaggedValues

     Description:
     Add the given tags and values to this DataTagged object. 
     \param tagKeys - Input - A vector of integer keys.
     \param values - Input - A vector of DataArrayViews. If this is empty
                   then all given tags will be assigned a value of zero. If 
                   it contains one value all tags will be assigned the same value.
                   Otherwise if there is a mismatch between the number of tags and 
                   the number of values an exception will be generated.
     NB: If a tag given here already exists in this object, this attempt to add the given
     value will have no effect. setTaggedValues is more versatile and should be
     used in most cases in preference to addTaggedValues.
  */
  void
  addTaggedValues(const TagListType& tagKeys,
                  const ValueListType& values);  

  /**
     \brief
     addTaggedValue

     Description:
     Add a single tag and value to this DataTagged object.
     \param tagKey - Input - Integer key.
     \param value - Input - DataArrayView.
     NB: If this tag already exists in this object, this attempt to add the given
     value will have no effect. setTaggedValue is more versatile and should be
     used in most cases in preference to addTaggedValue.
  */
  void
  addTaggedValue(int tagKey,
                 const DataArrayView& value);

  /**
     \brief
     setTaggedValues

     Description:
     Assign the given values to the tags listed.
     \param tagKeys - Input - A vector of integer keys.
     \param values - Input - A vector of DataArrayViews. If this is empty
                      then all given tags will be assigned a value of zero. If 
                      it contains one value all tag values will be assigned the same value.
                      Otherwise if there is a mismatch between the number of keys and 
                      the number of values an exception will be generated.
     NB: If a given tag does not yet exist in this DataTagged object, it will be added.
  */
  void
  setTaggedValues(const TagListType& tagKeys,
                  const ValueListType& values); 

  /**
     \brief
     setTaggedValue

     Description:
     Assign the given value to the given tag.
     \param tagKey - Input - Integer key.
     \param value - Input - Single DataArrayView value to be assigned to the tag.
     NB: If the given tag does not yet exist in this DataTagged object, it will be added.
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
  */
  DataArrayView
  getDataPointByTag(int tag) const;

  /**
     \brief
     getDataPoint

     Description:
     Return a view into the data-point specified by the given sample
     and data-point numbers.
     NOTE: Construction of the DataArrayView is a relatively expensive 
     operation.
     \param sampleNo - Input.
     \param dataPointNo - Input.
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
  */
  const DataMapType&
  getTagLookup() const;

  /**
     \brief
     isCurrentTag

     Description:
     Return true if the given tag exists within the DataTagged tag keys.
     NOTE: The DataTagged keys do not necessarily coincide with the tag
     keys for the function space.
   */
  bool
  isCurrentTag(int tag) const;

  /**
     \brief
     getDefaultValue

     Description:
     Return the default value. This value is associated with any tag which
     is not explicitly recorded in this DataTagged object.
  */
  DataArrayView&
  getDefaultValue();

  /**
     \brief
     getDefaultValue

     Description:
     Return the default value, const version.
  */
  const DataArrayView&
  getDefaultValue() const;

  /**
     \brief
     getLength

     Description:
     Return the number of doubles stored for the Data.
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

 protected:

 private:

  //
  // The offset lookup table
  DataMapType m_offsetLookup;

  //
  // The actual data
  DataArrayView::ValueType m_data;

  //
  // the default value offset
  static const int m_defaultValueOffset = 0;

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
  //
  // The default value is always the first value.
  return getPointDataView();
}

inline
const DataArrayView&
DataTagged::getDefaultValue() const
{
  //
  // The default value is always the first value.
  return getPointDataView();
}

} // end of namespace

#endif
