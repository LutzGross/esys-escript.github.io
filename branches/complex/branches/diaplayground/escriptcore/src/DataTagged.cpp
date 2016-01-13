
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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


#include "Data.h"
#include "DataTagged.h"
#include "DataConstant.h"
#include "DataException.h"
#include "esysUtils/Esys_MPI.h"

#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

#include "DataMaths.h"


#define CHECK_FOR_EX_WRITE if (!checkNoSharing()) {throw DataException("Attempt to modify shared object");}

// #define CHECK_FOR_EX_WRITE if (!checkNoSharing()) {std::ostringstream ss; ss << " Attempt to modify shared object. line " << __LINE__ << " of " << __FILE__; throw DataException(ss.str());}

using namespace std;

namespace escript {

DataTagged::DataTagged()
  : parent(FunctionSpace(),DataTypes::scalarShape)
{
  // default constructor

  // create a scalar default value
  m_data.resize(1,0.,1);
}

DataTagged::DataTagged(const FunctionSpace& what,
                       const DataTypes::ShapeType &shape,
                       const int tags[],
                       const ValueType& data)
  : parent(what,shape)
{
  // alternative constructor
  // not unit_tested tested yet
  // It is not explicitly unit tested yet, but it is called from DataFactory

  if (!what.canTag())
  {
    throw DataException("Programming error - DataTag created with a non-taggable FunctionSpace.");
  }
  // copy the data
  m_data=data;

  // we can't rely on the tag array to give us the number of tags so 
  // use the data we have been passed
  int valsize=DataTypes::noValues(shape);
  int ntags=data.size()/valsize;

  // create the tag lookup map
  // we assume that the first value and first tag are the default value so we skip
  for (int i=1;i<ntags;++i)
  {
    m_offsetLookup.insert(DataMapType::value_type(tags[i],i*valsize));
  }
}

DataTagged::DataTagged(const FunctionSpace& what,
                       const DataTypes::ShapeType &shape,
                       const TagListType& tags,
                       const ValueType& data)
  : parent(what,shape)
{
  // alternative constructor

  if (!what.canTag())
  {
    throw DataException("Programming error - DataTag created with a non-taggable FunctionSpace.");
  }

  // copy the data
  m_data=data;

  // create the view of the data
//   DataArrayView tempView(m_data,shape);
//   setPointDataView(tempView);

  // create the tag lookup map

//   for (int sampleNo=0; sampleNo<getNumSamples(); sampleNo++) {
//     m_offsetLookup.insert(DataMapType::value_type(sampleNo,tags[sampleNo]));
//   }

  // The above code looks like it will create a map the wrong way around

  int valsize=DataTypes::noValues(shape);
  int npoints=(data.size()/valsize)-1;
  int ntags=tags.size();
  if (ntags>npoints)
  {		// This throw is not unit tested yet
	throw DataException("Programming error - Too many tags for the supplied values.");
  }

  // create the tag lookup map
  // we assume that the first value is the default value so we skip it (hence the i+1 below)
  for (int i=0;i<ntags;++i)
  {
    m_offsetLookup.insert(DataMapType::value_type(tags[i],(i+1)*valsize));
  }
}


DataTagged::DataTagged(const DataTagged& other)
  : parent(other.getFunctionSpace(),other.getShape()),
  m_offsetLookup(other.m_offsetLookup),
  m_data(other.m_data)
{
  // copy constructor

  // create the data view
//   DataArrayView temp(m_data,other.getPointDataView().getShape());
//   setPointDataView(temp);
}

DataTagged::DataTagged(const DataConstant& other)
  : parent(other.getFunctionSpace(),other.getShape())
{
  // copy constructor

  if (!other.getFunctionSpace().canTag())
  {
    throw DataException("Programming error - DataTag created with a non-taggable FunctionSpace.");
  }

  // fill the default value with the constant value item from "other"
  int len = other.getNoValues();
  m_data.resize(len,0.,len);
  for (int i=0; i<len; i++) {
    m_data[i]=other.getVectorRO()[i];
  }
}


// Create a new object by copying tags
DataTagged::DataTagged(const FunctionSpace& what,
             const DataTypes::ShapeType& shape,
	     const DataTypes::ValueType& defaultvalue,
             const DataTagged* tagsource)
 : parent(what,shape)
{
// This constructor has not been unit tested yet

  if (defaultvalue.size()!=DataTypes::noValues(shape)) {
    throw DataException("Programming error - defaultvalue does not match supplied shape.");
  }


  if (!what.canTag())
  {
    throw DataException("Programming error - DataTag created with a non-taggable FunctionSpace.");
  }

  if (tagsource!=0)
  {
       m_data.resize(defaultvalue.size(),0.);	// since this is tagged data, we should have blocksize=1

       DataTagged::DataMapType::const_iterator i;
       for (i=tagsource->getTagLookup().begin();i!=tagsource->getTagLookup().end();i++) {
	  addTag(i->first);
       }
  }
  else
  {
	m_data.resize(defaultvalue.size());
  }

  // need to set the default value ....
  for (int i=0; i<defaultvalue.size(); i++) {
     m_data[i]=defaultvalue[i];
  }
}

DataAbstract*
DataTagged::deepCopy()
{
  return new DataTagged(*this);
}

DataAbstract*
DataTagged::getSlice(const DataTypes::RegionType& region) const 
{
  return new DataTagged(*this, region);
}

DataTagged::DataTagged(const DataTagged& other, 
		       const DataTypes::RegionType& region)
  : parent(other.getFunctionSpace(),DataTypes::getResultSliceShape(region))
{
  // slice constructor

  // get the shape of the slice to copy from other
  DataTypes::ShapeType regionShape(DataTypes::getResultSliceShape(region));
  DataTypes::RegionLoopRangeType regionLoopRange=DataTypes::getSliceRegionLoopRange(region);

  // allocate enough space in this for all values
  // (need to add one to allow for the default value)
  int len = DataTypes::noValues(regionShape)*(other.m_offsetLookup.size()+1);
  m_data.resize(len,0.0,len);

  // copy the default value from other to this
  const DataTypes::ShapeType& otherShape=other.getShape();
  const DataTypes::ValueType& otherData=other.getVectorRO();
  DataTypes::copySlice(getVectorRW(),getShape(),getDefaultOffset(),otherData,otherShape,other.getDefaultOffset(), regionLoopRange);

  // loop through the tag values copying these
  DataMapType::const_iterator pos;
  DataTypes::ValueType::size_type tagOffset=getNoValues();
  for (pos=other.m_offsetLookup.begin();pos!=other.m_offsetLookup.end();pos++){
    DataTypes::copySlice(m_data,getShape(),tagOffset,otherData, otherShape, pos->second, regionLoopRange);
    m_offsetLookup.insert(DataMapType::value_type(pos->first,tagOffset));
    tagOffset+=getNoValues();
  }
}

void
DataTagged::setSlice(const DataAbstract* other,
                     const DataTypes::RegionType& region)
{

  // other must be another DataTagged object
  // Data:setSlice implementation should ensure this
  const DataTagged* otherTemp=dynamic_cast<const DataTagged*>(other);
  if (otherTemp==0) {
    throw DataException("Programming error - casting to DataTagged.");
  }

  CHECK_FOR_EX_WRITE

  // determine shape of the specified region
  DataTypes::ShapeType regionShape(DataTypes::getResultSliceShape(region));

  // modify region specification as needed to match rank of this object
  DataTypes::RegionLoopRangeType regionLoopRange=DataTypes::getSliceRegionLoopRange(region);

  // ensure rank/shape of this object is compatible with specified region
  if (getRank()!=region.size()) {
    throw DataException("Error - Invalid slice region.");
  }
  if (otherTemp->getRank()>0 && !DataTypes::checkShape(other->getShape(),regionShape)) {
    throw DataException (DataTypes::createShapeErrorMessage(
                         "Error - Couldn't copy slice due to shape mismatch.",regionShape,other->getShape()));
  }

  const DataTypes::ValueType& otherData=otherTemp->getVectorRO();
  const DataTypes::ShapeType& otherShape=otherTemp->getShape();
  // copy slice from other default value to this default value
  DataTypes::copySliceFrom(m_data,getShape(),getDefaultOffset(),otherData,otherShape,otherTemp->getDefaultOffset(),regionLoopRange);

  // loop through tag values in other, adding any which aren't in this, using default value
  DataMapType::const_iterator pos;
  for (pos=otherTemp->m_offsetLookup.begin();pos!=otherTemp->m_offsetLookup.end();pos++) {
    if (!isCurrentTag(pos->first)) {
      addTag(pos->first);
    }
  }

  // loop through the tag values copying slices from other to this
  for (pos=m_offsetLookup.begin();pos!=m_offsetLookup.end();pos++) {
    DataTypes::copySliceFrom(m_data,getShape(),getOffsetForTag(pos->first),otherData, otherShape, otherTemp->getOffsetForTag(pos->first), regionLoopRange);

  }

}

int
DataTagged::getTagNumber(int dpno)
{
  //
  // Get the number of samples and data-points per sample
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  int numDataPoints = numSamples * numDataPointsPerSample;

  if (numDataPointsPerSample==0) {
    throw DataException("DataTagged::getTagNumber error: no data-points associated with this object.");
  }

  if (dpno<0 || dpno>numDataPoints-1) {
    throw DataException("DataTagged::getTagNumber error: invalid data-point number supplied.");
  }

  //
  // Determine the sample number which corresponds to this data-point number
  int sampleNo = dpno / numDataPointsPerSample;

  //
  // Determine the tag number which corresponds to this sample number
  int tagNo = getFunctionSpace().getTagFromSampleNo(sampleNo);

  //
  // return the tag number
  return(tagNo);
}

void
DataTagged::setTaggedValue(int tagKey,
			   const DataTypes::ShapeType& pointshape,
                           const ValueType& value,
			   int dataOffset)
{
  if (!DataTypes::checkShape(getShape(), pointshape)) {
      throw DataException(DataTypes::createShapeErrorMessage(
                          "Error - Cannot setTaggedValue due to shape mismatch.", pointshape,getShape()));
  }
  CHECK_FOR_EX_WRITE
  DataMapType::iterator pos(m_offsetLookup.find(tagKey));
  if (pos==m_offsetLookup.end()) {
    // tag couldn't be found so use addTaggedValue
    addTaggedValue(tagKey,pointshape, value, dataOffset);
  } else {
    // copy the values into the data array at the offset determined by m_offsetLookup
    int offset=pos->second;
    for (unsigned int i=0; i<getNoValues(); i++) {
      m_data[offset+i]=value[i+dataOffset];
    }
  }
}


void
DataTagged::addTaggedValues(const TagListType& tagKeys,
                            const ValueBatchType& values,
                            const ShapeType& vShape)
{
  DataTypes::ValueType t(values.size(),0);
  for (size_t i=0;i<values.size();++i)
  {
	t[i]=values[i];
  }
  addTaggedValues(tagKeys,t,vShape);
}


// Note: The check to see if vShape==our shape is done in the addTaggedValue method
void
DataTagged::addTaggedValues(const TagListType& tagKeys,
                            const ValueType& values,
                            const ShapeType& vShape)
{
  unsigned int n=getNoValues();
  unsigned int numVals=values.size()/n;
  if (values.size()==0) {
    // copy the current default value for each of the tags
    TagListType::const_iterator iT;
    for (iT=tagKeys.begin();iT!=tagKeys.end();iT++) {
      // the point data view for DataTagged points at the default value
      addTag(*iT);
    }
  } else if (numVals==1 && tagKeys.size()>1) {
    // assume the one given value will be used for all tag values
    TagListType::const_iterator iT;
    for (iT=tagKeys.begin();iT!=tagKeys.end();iT++) {
      addTaggedValue(*iT, vShape, values,0);
    }
  } else {
    if (tagKeys.size()!=numVals) {
      stringstream temp;
      temp << "Error - (addTaggedValue) Number of tags: " << tagKeys.size()
	   << " doesn't match number of values: " << values.size();
      throw DataException(temp.str());
    } else {
      unsigned int i;
      int offset=0;
      for (i=0;i<tagKeys.size();i++ ,offset+=n) {
        addTaggedValue(tagKeys[i],vShape,values,offset);
      }
    }
  }
}




void
DataTagged::addTaggedValue(int tagKey,
			   const DataTypes::ShapeType& pointshape,
                           const ValueType& value,
			   int dataOffset)
{
  if (!DataTypes::checkShape(getShape(), pointshape)) {
    throw DataException(DataTypes::createShapeErrorMessage(
                        "Error - Cannot addTaggedValue due to shape mismatch.", pointshape,getShape()));
  }
  CHECK_FOR_EX_WRITE
  DataMapType::iterator pos(m_offsetLookup.find(tagKey));
  if (pos!=m_offsetLookup.end()) {
    // tag already exists so use setTaggedValue
    setTaggedValue(tagKey,pointshape, value, dataOffset);
  } else {
    // save the key and the location of its data in the lookup tab
    m_offsetLookup.insert(DataMapType::value_type(tagKey,m_data.size()));
    // add the data given in "value" at the end of m_data
    // need to make a temp copy of m_data, resize m_data, then copy
    // all the old values plus the value to be added back into m_data
    ValueType m_data_temp(m_data);
    int oldSize=m_data.size();
    int newSize=m_data.size()+getNoValues();
    m_data.resize(newSize,0.,newSize);
    for (int i=0;i<oldSize;i++) {
      m_data[i]=m_data_temp[i];
    }
    for (unsigned int i=0;i<getNoValues();i++) {
      m_data[oldSize+i]=value[i+dataOffset];
    }
  }
}

void
DataTagged::addTag(int tagKey)
{
  CHECK_FOR_EX_WRITE
  DataMapType::iterator pos(m_offsetLookup.find(tagKey));
  if (pos!=m_offsetLookup.end()) {
    // tag already exists so use setTaggedValue
//    setTaggedValue(tagKey,value);
  } else {
    // save the key and the location of its data in the lookup tab
    m_offsetLookup.insert(DataMapType::value_type(tagKey,m_data.size()));
    // add the data given in "value" at the end of m_data
    // need to make a temp copy of m_data, resize m_data, then copy
    // all the old values plus the value to be added back into m_data
    ValueType m_data_temp(m_data);
    int oldSize=m_data.size();
    int newSize=m_data.size()+getNoValues();
    m_data.resize(newSize,0.,newSize);
    for (int i=0;i<oldSize;i++) {
      m_data[i]=m_data_temp[i];
    }
    for (unsigned int i=0;i<getNoValues();i++) {
      m_data[oldSize+i]=m_data[m_defaultValueOffset+i];
    }
  }
}


double*
DataTagged::getSampleDataByTag(int tag)
{
  CHECK_FOR_EX_WRITE
  DataMapType::iterator pos(m_offsetLookup.find(tag));
  if (pos==m_offsetLookup.end()) {
    // tag couldn't be found so return the default value
    return &(m_data[0]);
  } else {
    // return the data-point corresponding to the given tag
    return &(m_data[pos->second]);
  }
}


bool
DataTagged::hasNaN() const
{
  bool haveNaN=false;
  #pragma omp parallel for
	for (ValueType::size_type i=0;i<m_data.size();++i)
	{
		if (nancheck(m_data[i]))	// can't assume we have new standard NaN checking
		{
        #pragma omp critical 
        {
            haveNaN=true;
        }
		}
	}
	return haveNaN;
}

void
DataTagged::replaceNaN(double value) {
  #pragma omp parallel for
  for (ValueType::size_type i=0;i<m_data.size();++i)
  {
    if (nancheck(m_data[i]))  
    {
      m_data[i] = value;
    }
  }
}


string
DataTagged::toString() const
{
  using namespace escript::DataTypes;
  string empty="";
  stringstream temp;
  DataMapType::const_iterator i;
  temp << "Tag(Default)" << endl;
  temp << pointToString(m_data,getShape(),getDefaultOffset(),empty) << endl;
  // create a temporary view as the offset will be changed
//   DataArrayView tempView(getPointDataView().getData(), getPointDataView().getShape());
  for (i=m_offsetLookup.begin();i!=m_offsetLookup.end();++i) {
    temp << "Tag(" << i->first << ")" << endl;
    temp << pointToString(m_data,getShape(),i->second,empty) << endl;
//     tempView.setOffset(i->second);
//     temp << tempView.toString() << endl;
  }
  return temp.str();
}

DataTypes::ValueType::size_type 
DataTagged::getPointOffset(int sampleNo,
                           int dataPointNo) const
{
  int tagKey=getFunctionSpace().getTagFromSampleNo(sampleNo);
  DataMapType::const_iterator pos(m_offsetLookup.find(tagKey));
  DataTypes::ValueType::size_type offset=m_defaultValueOffset;
  if (pos!=m_offsetLookup.end()) {
    offset=pos->second;
  }
  return offset;
}

DataTypes::ValueType::size_type 
DataTagged::getPointOffset(int sampleNo,
                           int dataPointNo)
{
  int tagKey=getFunctionSpace().getTagFromSampleNo(sampleNo);
  DataMapType::const_iterator pos(m_offsetLookup.find(tagKey));
  DataTypes::ValueType::size_type offset=m_defaultValueOffset;
  if (pos!=m_offsetLookup.end()) {
    offset=pos->second;
  }
  return offset;
}

DataTypes::ValueType::size_type
DataTagged::getOffsetForTag(int tag) const
{
  DataMapType::const_iterator pos(m_offsetLookup.find(tag));
  DataTypes::ValueType::size_type offset=m_defaultValueOffset;
  if (pos!=m_offsetLookup.end()) {
    offset=pos->second;
  }
  return offset;
}

DataTypes::ValueType::const_reference
DataTagged::getDataByTagRO(int tag, DataTypes::ValueType::size_type i) const
{
  DataMapType::const_iterator pos(m_offsetLookup.find(tag));
  DataTypes::ValueType::size_type offset=m_defaultValueOffset;
  if (pos!=m_offsetLookup.end()) {
    offset=pos->second;
  }
  return m_data[offset+i];
}

DataTypes::ValueType::reference
DataTagged::getDataByTagRW(int tag, DataTypes::ValueType::size_type i)
{
  CHECK_FOR_EX_WRITE
  DataMapType::const_iterator pos(m_offsetLookup.find(tag));
  DataTypes::ValueType::size_type offset=m_defaultValueOffset;
  if (pos!=m_offsetLookup.end()) {
    offset=pos->second;
  }
  return m_data[offset+i];
}

void
DataTagged::symmetric(DataAbstract* ev)
{
  DataTagged* temp_ev=dynamic_cast<DataTagged*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataTagged::symmetric casting to DataTagged failed (probably a programming error).");
  }
  const DataTagged::DataMapType& thisLookup=getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
  ValueType& evVec=temp_ev->getVectorRW();
  const ShapeType& evShape=temp_ev->getShape();
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTag(i->first);
      DataTypes::ValueType::size_type offset=getOffsetForTag(i->first);
      DataTypes::ValueType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
      DataMaths::symmetric(m_data,getShape(),offset,evVec, evShape, evoffset);
  }
  DataMaths::symmetric(m_data,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset());
}


void
DataTagged::nonsymmetric(DataAbstract* ev)
{
  DataTagged* temp_ev=dynamic_cast<DataTagged*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataTagged::nonsymmetric casting to DataTagged failed (probably a programming error).");
  }
  const DataTagged::DataMapType& thisLookup=getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
  ValueType& evVec=temp_ev->getVectorRW();
  const ShapeType& evShape=temp_ev->getShape();
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTag(i->first);
      DataTypes::ValueType::size_type offset=getOffsetForTag(i->first);
      DataTypes::ValueType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
      DataMaths::nonsymmetric(m_data,getShape(),offset,evVec, evShape, evoffset);
  }
  DataMaths::nonsymmetric(m_data,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset());
}


void
DataTagged::trace(DataAbstract* ev, int axis_offset)
{
  DataTagged* temp_ev=dynamic_cast<DataTagged*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataTagged::trace casting to DataTagged failed (probably a programming error).");
  }
  const DataTagged::DataMapType& thisLookup=getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
  ValueType& evVec=temp_ev->getVectorRW();
  const ShapeType& evShape=temp_ev->getShape();
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTag(i->first);
      DataTypes::ValueType::size_type offset=getOffsetForTag(i->first);
      DataTypes::ValueType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
      DataMaths::trace(m_data,getShape(),offset,evVec, evShape, evoffset, axis_offset);
  }
  DataMaths::trace(m_data,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset(),axis_offset);
}

void
DataTagged::transpose(DataAbstract* ev, int axis_offset)
{
  DataTagged* temp_ev=dynamic_cast<DataTagged*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataTagged::transpose casting to DataTagged failed (probably a programming error).");
  }
  const DataTagged::DataMapType& thisLookup=getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
  ValueType& evVec=temp_ev->getVectorRW();
  const ShapeType& evShape=temp_ev->getShape();
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTag(i->first);
      DataTypes::ValueType::size_type offset=getOffsetForTag(i->first);
      DataTypes::ValueType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
      DataMaths::transpose(m_data,getShape(),offset,evVec, evShape, evoffset, axis_offset);
  }
  DataMaths::transpose(m_data,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset(),axis_offset);
}

void
DataTagged::swapaxes(DataAbstract* ev, int axis0, int axis1)
{
  DataTagged* temp_ev=dynamic_cast<DataTagged*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataTagged::swapaxes casting to DataTagged failed (probably a programming error).");
  }
  const DataTagged::DataMapType& thisLookup=getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
  ValueType& evVec=temp_ev->getVectorRW();
  const ShapeType& evShape=temp_ev->getShape();
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTag(i->first);
      DataTypes::ValueType::size_type offset=getOffsetForTag(i->first);
      DataTypes::ValueType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
      DataMaths::swapaxes(m_data,getShape(),offset,evVec, evShape, evoffset,axis0,axis1);
  }
  DataMaths::swapaxes(m_data,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset(),axis0,axis1);
}

void
DataTagged::eigenvalues(DataAbstract* ev)
{
  DataTagged* temp_ev=dynamic_cast<DataTagged*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataTagged::eigenvalues casting to DataTagged failed (probably a programming error).");
  }
  const DataTagged::DataMapType& thisLookup=getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
  ValueType& evVec=temp_ev->getVectorRW();
  const ShapeType& evShape=temp_ev->getShape();
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTag(i->first);
//       DataArrayView thisView=getDataPointByTag(i->first);
//       DataArrayView evView=temp_ev->getDataPointByTag(i->first);
      DataTypes::ValueType::size_type offset=getOffsetForTag(i->first);
      DataTypes::ValueType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
      DataMaths::eigenvalues(m_data,getShape(),offset,evVec, evShape, evoffset);
  }
  DataMaths::eigenvalues(m_data,getShape(),getDefaultOffset(),evVec, evShape, temp_ev->getDefaultOffset());
}
void
DataTagged::eigenvalues_and_eigenvectors(DataAbstract* ev,DataAbstract* V,const double tol)
{
  DataTagged* temp_ev=dynamic_cast<DataTagged*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataTagged::eigenvalues_and_eigenvectors casting to DataTagged failed (probably a programming error).");
  }
  DataTagged* temp_V=dynamic_cast<DataTagged*>(V);
  if (temp_V==0) {
    throw DataException("Error - DataTagged::eigenvalues_and_eigenvectors casting to DataTagged failed (probably a programming error).");
  }
  const DataTagged::DataMapType& thisLookup=getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
  ValueType& evVec=temp_ev->getVectorRW();
  const ShapeType& evShape=temp_ev->getShape();
  ValueType& VVec=temp_V->getVectorRW();
  const ShapeType& VShape=temp_V->getShape();
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTag(i->first);
      temp_V->addTag(i->first);
/*      DataArrayView thisView=getDataPointByTag(i->first);
      DataArrayView evView=temp_ev->getDataPointByTag(i->first);
      DataArrayView VView=temp_V->getDataPointByTag(i->first);*/
      DataTypes::ValueType::size_type offset=getOffsetForTag(i->first);
      DataTypes::ValueType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
      DataTypes::ValueType::size_type Voffset=temp_V->getOffsetForTag(i->first);
/*      DataArrayView::eigenvalues_and_eigenvectors(thisView,0,evView,0,VView,0,tol);*/
      DataMaths::eigenvalues_and_eigenvectors(m_data,getShape(),offset,evVec, evShape, evoffset,VVec,VShape,Voffset,tol);

  }
  DataMaths::eigenvalues_and_eigenvectors(m_data,getShape(),getDefaultOffset(),evVec, evShape,
					  temp_ev->getDefaultOffset(),VVec,VShape,
					  temp_V->getDefaultOffset(), tol);


}

int
DataTagged::matrixInverse(DataAbstract* out) const
{
  DataTagged* temp=dynamic_cast<DataTagged*>(out);
  if (temp==0)
  {
	throw DataException("Error - DataTagged::matrixInverse: casting to DataTagged failed (probably a programming error).");
  }
  if (getRank()!=2)
  {
	throw DataException("Error - DataExpanded::matrixInverse: input must be rank 2.");
  }
  const DataTagged::DataMapType& thisLookup=getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
  ValueType& outVec=temp->getVectorRW();
  const ShapeType& outShape=temp->getShape();
  LapackInverseHelper h(getShape()[0]);
  int err=0;
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp->addTag(i->first);
      DataTypes::ValueType::size_type inoffset=getOffsetForTag(i->first);
      DataTypes::ValueType::size_type outoffset=temp->getOffsetForTag(i->first);

      err=DataMaths::matrix_inverse(m_data, getShape(), inoffset, outVec, outShape, outoffset, 1, h);
      if (!err) break;
  }
  if (!err)
  {
      DataMaths::matrix_inverse(m_data, getShape(), getDefaultOffset(), outVec, outShape, temp->getDefaultOffset(), 1, h);
  }
  return err;
}

void
DataTagged::setToZero(){
    CHECK_FOR_EX_WRITE
    DataTypes::ValueType::size_type n=m_data.size();
    for (int i=0; i<n ;++i) m_data[i]=0.;
}

void
DataTagged::dump(const std::string fileName) const
{
   #ifdef USE_NETCDF
   const int ldims=DataTypes::maxRank+1;
   const NcDim* ncdims[ldims];
   NcVar *var, *tags_var;
   int rank = getRank();
   int type=  getFunctionSpace().getTypeCode();
   int ndims =0;
   long dims[ldims];
   const double* d_ptr=&(m_data[0]);
   DataTypes::ShapeType shape = getShape();
   int mpi_iam=getFunctionSpace().getDomain()->getMPIRank();
   int mpi_num=getFunctionSpace().getDomain()->getMPISize();
#ifdef ESYS_MPI
   MPI_Status status;
#endif

#ifdef ESYS_MPI
   /* Serialize NetCDF I/O */
   if (mpi_iam>0) MPI_Recv(&ndims, 0, MPI_INT, mpi_iam-1, 81803, MPI_COMM_WORLD, &status);
#endif

   // netCDF error handler
   NcError err(NcError::verbose_nonfatal);
   // Create the file.
   const std::string newFileName(esysUtils::appendRankToFileName(fileName,
                                                            mpi_num, mpi_iam));
   NcFile dataFile(newFileName.c_str(), NcFile::Replace);
   // check if writing was successful
   if (!dataFile.is_valid())
        throw DataException("Error - DataTagged:: opening of netCDF file for output failed.");
   if (!dataFile.add_att("type_id",1) )
        throw DataException("Error - DataTagged:: appending data type to netCDF file failed.");
   if (!dataFile.add_att("rank",rank) )
        throw DataException("Error - DataTagged:: appending rank attribute to netCDF file failed.");
   if (!dataFile.add_att("function_space_type",type))
        throw DataException("Error - DataTagged:: appending function space attribute to netCDF file failed.");
   ndims=rank+1;
   if ( rank >0 ) {
       dims[0]=shape[0];
       if (! (ncdims[0] = dataFile.add_dim("d0",shape[0])) )
            throw DataException("Error - DataTagged:: appending ncdimension 0 to netCDF file failed.");
   }
   if ( rank >1 ) {
       dims[1]=shape[1];
       if (! (ncdims[1] = dataFile.add_dim("d1",shape[1])) )
            throw DataException("Error - DataTagged:: appending ncdimension 1 to netCDF file failed.");
   }
   if ( rank >2 ) {
       dims[2]=shape[2];
       if (! (ncdims[2] = dataFile.add_dim("d2", shape[2])) )
            throw DataException("Error - DataTagged:: appending ncdimension 2 to netCDF file failed.");
   }
   if ( rank >3 ) {
       dims[3]=shape[3];
       if (! (ncdims[3] = dataFile.add_dim("d3", shape[3])) )
            throw DataException("Error - DataTagged:: appending ncdimension 3 to netCDF file failed.");
   }
   const DataTagged::DataMapType& thisLookup=getTagLookup();
   DataTagged::DataMapType::const_iterator i;
   DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
   std::vector<int> tags;
   tags.push_back(-1);
   for (i=thisLookup.begin();i!=thisLookupEnd;i++)
       tags.push_back(i->first);
   dims[rank]=tags.size();
   if (! (ncdims[rank] = dataFile.add_dim("num_tags", dims[rank])) )
   {
           throw DataException("Error - DataTagged:: appending num_tags to netCDF file failed.");
   }
   if (! ( tags_var = dataFile.add_var("tags", ncInt, ncdims[rank])) )
   {
        throw DataException("Error - DataTagged:: appending tags to netCDF file failed.");
   }
   if (! (tags_var->put(&tags[0], dims[rank])) )
   {
        throw DataException("Error - DataTagged:: copy tags to netCDF buffer failed.");
   }
   if (! ( var = dataFile.add_var("data", ncDouble, ndims, ncdims)) )
   {
        throw DataException("Error - DataTagged:: appending variable to netCDF file failed.");
   }
   if (! (var->put(d_ptr,dims)) )
   {
        throw DataException("Error - DataTagged:: copy data to netCDF buffer failed.");
   }
#ifdef ESYS_MPI
   if (mpi_iam<mpi_num-1) MPI_Send(&ndims, 0, MPI_INT, mpi_iam+1, 81803, MPI_COMM_WORLD);
#endif
   #else
   throw DataException("Error - DataTagged:: dump is not configured with netCDF. Please contact your installation manager.");
   #endif
}

DataTypes::ValueType&
DataTagged::getVectorRW()
{
    CHECK_FOR_EX_WRITE
    return m_data;
}

const DataTypes::ValueType&
DataTagged::getVectorRO() const
{
	return m_data;
}

}  // end of namespace
