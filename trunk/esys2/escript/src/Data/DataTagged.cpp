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

#include "escript/Data/DataTagged.h"
#include "escript/Data/DataConstant.h"
#include "escript/Data/DataExpanded.h"
#include "escript/Data/DataException.h"

#include <sstream>

using namespace std;

namespace escript {

DataTagged::DataTagged():
  DataAbstract(FunctionSpace())
{
  //
  // create a scalar default value
  m_data.push_back(0.0);
  DataArrayView temp(m_data,DataArrayView::ShapeType());
  setPointDataView(temp);
}

DataTagged::DataTagged(const TagListType& tagKeys, 
		       const ValueListType& values,
		       const DataArrayView& defaultValue,
		       const FunctionSpace& what)
  : DataAbstract(what)
{
  //
  // Initialise the array of data values
  // The default value is always the first item in the values list
  m_data.insert(m_data.end(), &defaultValue.getData(0), &defaultValue.getData(defaultValue.noValues()) );
  // create the data view
  DataArrayView temp(m_data,defaultValue.getShape());
  setPointDataView(temp);
  // add remaining tags and values
  addTaggedValues(tagKeys,values);
}

DataTagged::DataTagged(const FunctionSpace& what,
                       const DataArrayView::ShapeType &shape,
                       const int tags[],
                       const DataArrayView::ValueType &data)
  : DataAbstract(what)
{
  //
  // copy the data in the correct format
  m_data=data;
  //
  // create the view of the data
  DataArrayView tempView(m_data,shape);
  setPointDataView(tempView);
  //
  // create the tag lookup map
  for (int sampleNo=0; sampleNo<getNumSamples(); sampleNo++) {
    m_offsetLookup.insert(DataMapType::value_type(sampleNo,tags[sampleNo]));
  }
}

DataTagged::DataTagged(const DataTagged& other)
  : DataAbstract(other.getFunctionSpace()),
  m_data(other.m_data),
  m_offsetLookup(other.m_offsetLookup)
{
  // create the data view
  DataArrayView temp(m_data,other.getPointDataView().getShape());
  setPointDataView(temp);
}

DataTagged::DataTagged(const DataConstant& other)
  : DataAbstract(other.getFunctionSpace())
{
  //
  // Fill the default value with the constant value item from other
  const DataArrayView& value=other.getPointDataView();
  m_data.insert(m_data.end(), &value.getData(0), &value.getData(value.noValues()) );
  // create the data view
  DataArrayView temp(m_data,value.getShape());
  setPointDataView(temp);
}

DataTagged::DataTagged(const DataTagged& other, 
		       const DataArrayView::RegionType& region)
  : DataAbstract(other.getFunctionSpace())
{
  //
  // get the shape of the slice to copy from other
  DataArrayView::ShapeType shape(DataArrayView::getResultSliceShape(region));
  DataArrayView::RegionLoopRangeType region_loop_range=getSliceRegionLoopRange(region);
  // allocate enough space for all values
  m_data.resize(DataArrayView::noValues(shape)*(other.m_offsetLookup.size()+1));
  // create the data view
  DataArrayView temp(m_data,shape);
  setPointDataView(temp);
  // copy the default value
  getDefaultValue().copySlice(other.getDefaultValue(),region_loop_range);
  //
  // Loop through the tag values copying these
  DataMapType::const_iterator pos;
  DataArrayView::ValueType::size_type tagOffset=getPointDataView().noValues();
  for (pos=other.m_offsetLookup.begin();pos!=other.m_offsetLookup.end();++pos){
    getPointDataView().copySlice(tagOffset,other.getPointDataView(), pos->second,region_loop_range);
    m_offsetLookup.insert(DataMapType::value_type(pos->first,tagOffset));
    tagOffset+=getPointDataView().noValues();
  } 
}

void
DataTagged::reshapeDataPoint(const DataArrayView::ShapeType& shape) 
{
  //
  // can only reshape a rank zero data point
  if (getPointDataView().getRank()!=0) {
    stringstream temp;
    temp << "Error - Can only reshape Data with data points of rank 0. "
	 << "This Data has data points with rank: " 
	 << getPointDataView().getRank();
    throw DataException(temp.str());
  }
  //
  // allocate enough space for all values
  DataArrayView::ValueType newData(DataArrayView::noValues(shape)*(m_offsetLookup.size()+1));
  DataArrayView newView(newData,shape);
  newView.copy(0,getDefaultValue()());
  //
  // Loop through the tag values
  DataMapType::iterator pos;
  DataArrayView::ValueType::size_type tagOffset=DataArrayView::noValues(shape);
  for (pos=m_offsetLookup.begin();pos!=m_offsetLookup.end();++pos){
    newView.copy(tagOffset,m_data[pos->second]);
    pos->second=tagOffset;
    tagOffset+=DataArrayView::noValues(shape);
  }
  m_data=newData;
  DataArrayView temp(m_data,shape);
  setPointDataView(temp);
}

DataAbstract*
DataTagged::getSlice(const DataArrayView::RegionType& region) const 
{
  return new DataTagged(*this,region);
}

void
DataTagged::setSlice(const DataAbstract* value,
                     const DataArrayView::RegionType& region) 
{
  const DataTagged* tempDataTag=dynamic_cast<const DataTagged*>(value);
  if (tempDataTag==0) {
    throw DataException("Programming error - casting to DataTagged.");
  }
  //
  DataArrayView::ShapeType shape(DataArrayView::getResultSliceShape(region));
  DataArrayView::RegionLoopRangeType region_loop_range=getSliceRegionLoopRange(region);
  //
  if (getPointDataView().getRank()!=region.size()) {
    throw DataException("Error - Invalid slice region.");
  }
  if (tempDataTag->getPointDataView().getRank()>0 && !value->getPointDataView().checkShape(shape)) {
    throw DataException (value->getPointDataView().createShapeErrorMessage(
                "Error - Couldn't copy slice due to shape mismatch.",shape));
  }
  //
  getDefaultValue().copySliceFrom(tempDataTag->getDefaultValue(),region_loop_range);
  //
  // Loop through the tag values
  DataMapType::const_iterator pos;
  for (pos=m_offsetLookup.begin();pos!=m_offsetLookup.end();++pos){
    getDataPointByTag(pos->first).copySliceFrom(tempDataTag->getDataPointByTag(pos->first),region_loop_range);
  } 
}

void
DataTagged::setTaggedValue(int tagKey,
                           const DataArrayView& value)
{
  DataMapType::iterator pos(m_offsetLookup.find(tagKey));
  if (pos==m_offsetLookup.end()) {
    //
    // tag couldn't be found so add as a new tag
    addTaggedValue(tagKey,value);
  } else {
    if (!getPointDataView().checkShape(value.getShape())) {
      throw DataException(getPointDataView().createShapeErrorMessage(
		 "Error - Cannot setTaggedValue due to shape mismatch.", value.getShape()));
    }
    //
    // copy the values into tagged data storage
    copy(&value.getData(0), &value.getData(getPointDataView().noValues()), &m_data[pos->second]);
  }
}

void
DataTagged::setTaggedValues(const TagListType& tagKeys,
                            const ValueListType& values)
{
  for (int i=0;i<tagKeys.size();++i) {
    setTaggedValue(tagKeys[i],values[i]);
  }
}

void
DataTagged::addTaggedValue(int tagKey,
                           const DataArrayView& value)
{
  if (!getPointDataView().checkShape(value.getShape())) {
    throw DataException(getPointDataView().createShapeErrorMessage(
		 "Error - Cannot addTaggedValue due to shape mismatch.", value.getShape()));
  }
  //
  // save the key and the location of its data
  m_offsetLookup.insert( DataMapType::value_type(tagKey,m_data.size()) );
  //
  // insert the data given in value at the end of m_data
  m_data.insert( m_data.end(), &(value.getData(0)), &(value.getData(value.noValues())) );
}

void
DataTagged::addTaggedValues(const TagListType& tagKeys,
                            const ValueListType& values)
{
  if (values.size()==0) {
    //
    // Copy the default value for each of the tags
    TagListType::const_iterator iT;
    for (iT=tagKeys.begin();iT!=tagKeys.end();++iT) {
      //
      // the point data view for DataTagged points at the default value
      addTaggedValue(*iT,getPointDataView());
    }
  } else if (values.size()==1 && tagKeys.size()>1) {
    //
    // assume the one value will be used for all tag values
    // Copy the input data
    TagListType::const_iterator iT;
    for (iT=tagKeys.begin();iT!=tagKeys.end();++iT) {
      addTaggedValue(*iT,values[0]);
    }
  } else {
    if (tagKeys.size()!=values.size()) {
      stringstream temp;
      temp << "Error - (addTaggedValue) Number of tags: " << tagKeys.size()
	   << " doesn't match the number of values: " << values.size();
      throw DataException(temp.str());
    } else {
      for (int i=0;i<tagKeys.size();++i) {
        addTaggedValue(tagKeys[i],values[i]);
      }
    }
  }
}

double*
DataTagged::getSampleDataByTag(int tag)
{
  DataMapType::iterator pos(m_offsetLookup.find(tag));
  if (pos==m_offsetLookup.end()) {
    //
    // tag couldn't be found so return the default value
    return &(m_data[0]);
  } else {
    //
    // return the data-point corresponding to the given tag
    return &(m_data[pos->second]);
  }
}

string
DataTagged::toString() const
{
  stringstream temp;
  DataMapType::const_iterator i;
  temp << "Tag(Default)" << endl;
  temp << getDefaultValue().toString() << endl;
  //
  // create a temporary view as the offset will be changed
  DataArrayView tempView(getPointDataView().getData(), getPointDataView().getShape());
  for (i=m_offsetLookup.begin();i!=m_offsetLookup.end();++i) {
    temp << "Tag(" << i->first << ")" << endl;
    tempView.setOffset(i->second);
    temp << tempView.toString() << endl;
  }
  return temp.str();
}

DataArrayView
DataTagged::getDataPointByTag(int tag) const
{
  DataMapType::const_iterator pos(m_offsetLookup.find(tag));
  DataArrayView::ValueType::size_type offset=m_defaultValueOffset;
  if (pos!=m_offsetLookup.end()) {
    offset=pos->second;
  }
  DataArrayView temp(getPointDataView());
  temp.setOffset(offset);
  return temp;
}

DataArrayView::ValueType::size_type 
DataTagged::getPointOffset(int sampleNo,
                           int dataPointNo) const
{
  int tagKey=getFunctionSpace().getTagFromSampleNo(sampleNo);
  DataMapType::const_iterator pos(m_offsetLookup.find(tagKey));
  DataArrayView::ValueType::size_type offset=m_defaultValueOffset;
  if (pos!=m_offsetLookup.end()) {
    offset=pos->second;
  }
  return offset;
}

DataArrayView
DataTagged::getDataPoint(int sampleNo,
                         int dataPointNo)
{
  EsysAssert(validSampleNo(sampleNo),"(getDataPoint) Invalid sampleNo: " << sampleNo);
  int tagKey=getFunctionSpace().getTagFromSampleNo(sampleNo);
  return getDataPointByTag(tagKey);
}

const DataTagged::DataMapType&
DataTagged::getTagLookup() const
{
  return m_offsetLookup;
}

DataArrayView::ValueType::size_type
DataTagged::getLength() const
{
  return m_data.size();
}

}  // end of namespace
