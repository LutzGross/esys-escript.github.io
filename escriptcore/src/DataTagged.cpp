
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "Data.h"
#include "DataConstant.h"
#include "DataException.h"
#include "DataVectorOps.h"
#include "DataTagged.h"
#include "Utils.h"

#include <complex>


#ifdef ESYS_HAVE_HDF5
  #include <H5Cpp.h>
#endif

#ifdef SLOWSHARECHECK
  #define CHECK_FOR_EX_WRITE if (isShared()) {throw DataException("Attempt to modify shared object");}
#else
  #define CHECK_FOR_EX_WRITE
#endif

using namespace std;

namespace escript {

DataTagged::DataTagged(const FunctionSpace& what,
                       const DataTypes::ShapeType &shape,
                       const int tags[],
                       const DataTypes::RealVectorType& data)
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
  m_data_r=data;

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
                       const int tags[],
                       const DataTypes::CplxVectorType& data)
  : parent(what,shape)
{
  // alternative constructor
  // not unit_tested tested yet
  // It is not explicitly unit tested yet, but it is called from DataFactory

  m_iscompl=true;
  if (!what.canTag())
  {
    throw DataException("Programming error - DataTag created with a non-taggable FunctionSpace.");
  }
  // copy the data
  m_data_c=data;

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
                       const DataTypes::RealVectorType& data)
  : parent(what,shape)
{
  // alternative constructor

  if (!what.canTag())
  {
    throw DataException("Programming error - DataTag created with a non-taggable FunctionSpace.");
  }

  // copy the data
  m_data_r=data;

  // The above code looks like it will create a map the wrong way around

  int valsize=DataTypes::noValues(shape);
  int npoints=(data.size()/valsize)-1;
  int ntags=tags.size();
  if (ntags>npoints)
  {     // This throw is not unit tested yet
        throw DataException("Programming error - Too many tags for the supplied values.");
  }

  // create the tag lookup map
  // we assume that the first value is the default value so we skip it (hence the i+1 below)
  for (int i=0;i<ntags;++i)
  {
    m_offsetLookup.insert(DataMapType::value_type(tags[i],(i+1)*valsize));
  }
}


DataTagged::DataTagged(const FunctionSpace& what,
                       const DataTypes::ShapeType &shape,
                       const TagListType& tags,
                       const DataTypes::CplxVectorType& data)
  : parent(what,shape)
{
  // alternative constructor
  m_iscompl=true;
  if (!what.canTag())
  {
    throw DataException("Programming error - DataTag created with a non-taggable FunctionSpace.");
  }

  // copy the data
  m_data_c=data;

  // The above code looks like it will create a map the wrong way around

  int valsize=DataTypes::noValues(shape);
  int npoints=(data.size()/valsize)-1;
  int ntags=tags.size();
  if (ntags>npoints)
  {             // This throw is not unit tested yet
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
  m_data_r(other.m_data_r), m_data_c(other.m_data_c)
{
  // copy constructor
    m_iscompl=other.m_iscompl;
}

DataTagged::DataTagged(const DataConstant& other)
  : parent(other.getFunctionSpace(),other.getShape())
{
  // copy constructor
  m_iscompl=other.isComplex();
  if (!other.getFunctionSpace().canTag())
  {
    throw DataException("Programming error - DataTag created with a non-taggable FunctionSpace.");
  }

  // fill the default value with the constant value item from "other"
  int len = other.getNoValues();
  if (m_iscompl)
  {
      DataTypes::cplx_t dummy=0;
      m_data_c.resize(len,0.,len);
      for (int i=0; i<len; i++) {
        m_data_c[i]=other.getTypedVectorRO(dummy)[i];
      }
  }
  else
  {
      DataTypes::real_t dummy=0;
      m_data_r.resize(len,0.,len);
      for (int i=0; i<len; i++) {
        m_data_r[i]=other.getTypedVectorRO(dummy)[i];
      }
  }
}


// Create a new object by copying tags
DataTagged::DataTagged(const FunctionSpace& what,
             const DataTypes::ShapeType& shape,
             const DataTypes::RealVectorType& defaultvalue,
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
       m_data_r.resize(defaultvalue.size(),0.); // since this is tagged data, we should have blocksize=1

       DataTagged::DataMapType::const_iterator i;
       for (i=tagsource->getTagLookup().begin();i!=tagsource->getTagLookup().end();i++) {
          addTag(i->first);
       }
  }
  else
  {
        m_data_r.resize(defaultvalue.size());
  }

  // need to set the default value ....
  for (int i=0; i<defaultvalue.size(); i++) {
     m_data_r[i]=defaultvalue[i];
  }
}

// Create a new object by copying tags
DataTagged::DataTagged(const FunctionSpace& what,
             const DataTypes::ShapeType& shape,
             const DataTypes::CplxVectorType& defaultvalue,
             const DataTagged* tagsource)
 : parent(what,shape)
{
// This constructor has not been unit tested yet
  m_iscompl=true;
  
  if (defaultvalue.size()!=DataTypes::noValues(shape)) {
    throw DataException("Programming error - defaultvalue does not match supplied shape.");
  }

  if (!what.canTag())
  {
    throw DataException("Programming error - DataTag created with a non-taggable FunctionSpace.");
  }

  if (tagsource!=0)
  {
       m_data_r.resize(defaultvalue.size(),0.); // since this is tagged data, we should have blocksize=1

       DataTagged::DataMapType::const_iterator i;
       for (i=tagsource->getTagLookup().begin();i!=tagsource->getTagLookup().end();i++) {
          addTag(i->first);
       }
  }
  else
  {
        m_data_r.resize(defaultvalue.size());
  }

  // need to set the default value ....
  for (int i=0; i<defaultvalue.size(); i++) {
     m_data_c[i]=defaultvalue[i];
  }
}


DataAbstract*
DataTagged::deepCopy() const
{
  return new DataTagged(*this);
}


DataAbstract*
DataTagged::zeroedCopy() const
{
    DataTagged* p=0;
    if (isComplex())
    {
        DataTypes::CplxVectorType v(1);
	v[0]=0;
        p=new DataTagged(this->getFunctionSpace(), this->getShape(), v, this);
    }
    else
    {
        DataTypes::RealVectorType v(1);
	v[0]=0;      
        p=new DataTagged(this->getFunctionSpace(), this->getShape(), v, this);
    }   
  return p;
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
  m_iscompl=other.isComplex();
  
  
  // get the shape of the slice to copy from other
  DataTypes::ShapeType regionShape(DataTypes::getResultSliceShape(region));
  DataTypes::RegionLoopRangeType regionLoopRange=DataTypes::getSliceRegionLoopRange(region);

  // allocate enough space in this for all values
  // (need to add one to allow for the default value)
  int len = DataTypes::noValues(regionShape)*(other.m_offsetLookup.size()+1);
  if (m_iscompl)
  {
      m_data_c.resize(len,0.0,len);
      // copy the default value from other to this
      const DataTypes::ShapeType& otherShape=other.getShape();
      const DataTypes::CplxVectorType& otherData=other.getTypedVectorRO((DataTypes::cplx_t)0);
      DataTypes::copySlice(getTypedVectorRW((DataTypes::cplx_t)0),getShape(),getDefaultOffset(),otherData,otherShape,other.getDefaultOffset(), regionLoopRange);

      // loop through the tag values copying these
      DataMapType::const_iterator pos;
      DataTypes::CplxVectorType::size_type tagOffset=getNoValues();
      for (pos=other.m_offsetLookup.begin();pos!=other.m_offsetLookup.end();pos++){
        DataTypes::copySlice(m_data_c,getShape(),tagOffset,otherData, otherShape, pos->second, regionLoopRange);
        m_offsetLookup.insert(DataMapType::value_type(pos->first,tagOffset));
        tagOffset+=getNoValues();
      }      
      
      
  }
  else
  {
      m_data_r.resize(len,0.0,len);    
      // copy the default value from other to this
      const DataTypes::ShapeType& otherShape=other.getShape();
      const DataTypes::RealVectorType& otherData=other.getTypedVectorRO((DataTypes::real_t)0);
      DataTypes::copySlice(getTypedVectorRW((DataTypes::real_t)0),getShape(),getDefaultOffset(),otherData,otherShape,other.getDefaultOffset(), regionLoopRange);

      // loop through the tag values copying these
      DataMapType::const_iterator pos;
      DataTypes::RealVectorType::size_type tagOffset=getNoValues();
      for (pos=other.m_offsetLookup.begin();pos!=other.m_offsetLookup.end();pos++){
        DataTypes::copySlice(m_data_r,getShape(),tagOffset,otherData, otherShape, pos->second, regionLoopRange);
        m_offsetLookup.insert(DataMapType::value_type(pos->first,tagOffset));
        tagOffset+=getNoValues();
      }      
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
  if (isComplex()!=other->isComplex())
  {
    throw DataException("Error - cannot copy between slices of different complexity.");
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


  const DataTypes::ShapeType& otherShape=otherTemp->getShape();
  if (isComplex())      // from check earlier, other will have the same complexity
  {
      // copy slice from other default value to this default value
      DataTypes::copySliceFrom(m_data_c,getShape(),getDefaultOffset(),otherTemp->getTypedVectorRO((DataTypes::cplx_t)0),
                               otherShape,otherTemp->getDefaultOffset(),regionLoopRange);
  } 
  else
  {
      // copy slice from other default value to this default value
      DataTypes::copySliceFrom(m_data_r,getShape(),getDefaultOffset(),otherTemp->getTypedVectorRO((DataTypes::real_t)0),
                               otherShape,otherTemp->getDefaultOffset(),regionLoopRange);
  }

  // loop through tag values in other, adding any which aren't in this, using default value
  DataMapType::const_iterator pos;
  for (pos=otherTemp->m_offsetLookup.begin();pos!=otherTemp->m_offsetLookup.end();pos++) {
    if (!isCurrentTag(pos->first)) {
      addTag(pos->first);
    }
  }
  if (isComplex())
  {
    // loop through the tag values copying slices from other to this
    for (pos=m_offsetLookup.begin();pos!=m_offsetLookup.end();pos++) {
      DataTypes::copySliceFrom(m_data_c,getShape(),getOffsetForTag(pos->first),otherTemp->getTypedVectorRO((DataTypes::cplx_t)0),
                               otherShape, otherTemp->getOffsetForTag(pos->first), regionLoopRange);
    }
  }
  else
  {
    // loop through the tag values copying slices from other to this
    for (pos=m_offsetLookup.begin();pos!=m_offsetLookup.end();pos++) {
      DataTypes::copySliceFrom(m_data_r,getShape(),getOffsetForTag(pos->first),otherTemp->getTypedVectorRO((DataTypes::real_t)0),
                               otherShape, otherTemp->getOffsetForTag(pos->first), regionLoopRange);
    }
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
                           const DataTypes::RealVectorType& value,
                           int dataOffset)
{
  if (!DataTypes::checkShape(getShape(), pointshape)) {
      throw DataException(DataTypes::createShapeErrorMessage(
                          "Error - Cannot setTaggedValue due to shape mismatch.", pointshape,getShape()));
  }
  if (isComplex())
  {
      throw DataException("Programming Error - attempt to set real value on complex data.");
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
      m_data_r[offset+i]=value[i+dataOffset];
    }
  }
}

void
DataTagged::setTaggedValue(int tagKey,
                           const DataTypes::ShapeType& pointshape,
                           const DataTypes::CplxVectorType& value,
                           int dataOffset)
{
  if (!DataTypes::checkShape(getShape(), pointshape)) {
      throw DataException(DataTypes::createShapeErrorMessage(
                          "Error - Cannot setTaggedValue due to shape mismatch.", pointshape,getShape()));
  }
  if (!isComplex())
  {
      throw DataException("Programming Error - attempt to set a complex value on real data");
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
      m_data_c[offset+i]=value[i+dataOffset];
    }
  }
}


void
DataTagged::addTaggedValues(const TagListType& tagKeys,
                            const FloatBatchType& values,
                            const ShapeType& vShape)
{
  DataTypes::RealVectorType t(values.size(),0);
  for (size_t i=0;i<values.size();++i)
  {
        t[i]=values[i];
  }
  addTaggedValues(tagKeys,t,vShape);
}


// Note: The check to see if vShape==our shape is done in the addTaggedValue method
void
DataTagged::addTaggedValues(const TagListType& tagKeys,
                            const DataTypes::RealVectorType& values,
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
                           const DataTypes::RealVectorType& value,
                           int dataOffset)
{
  if (!DataTypes::checkShape(getShape(), pointshape)) {
    throw DataException(DataTypes::createShapeErrorMessage(
                        "Error - Cannot addTaggedValue due to shape mismatch.", pointshape,getShape()));
  }
  if (isComplex())
  {
      throw DataException("Programming Error - attempt to set a real value on complex data");
  }
  CHECK_FOR_EX_WRITE
  DataMapType::iterator pos(m_offsetLookup.find(tagKey));
  if (pos!=m_offsetLookup.end()) {
    // tag already exists so use setTaggedValue
    setTaggedValue(tagKey,pointshape, value, dataOffset);
  } else {
    // save the key and the location of its data in the lookup tab
    m_offsetLookup.insert(DataMapType::value_type(tagKey,m_data_r.size()));
    // add the data given in "value" at the end of m_data_r
    // need to make a temp copy of m_data_r, resize m_data_r, then copy
    // all the old values plus the value to be added back into m_data_r
    DataTypes::RealVectorType m_data_r_temp(m_data_r);
    int oldSize=m_data_r.size();
    int newSize=m_data_r.size()+getNoValues();
    m_data_r.resize(newSize,0.,newSize);
    for (int i=0;i<oldSize;i++) {
      m_data_r[i]=m_data_r_temp[i];
    }
    for (unsigned int i=0;i<getNoValues();i++) {
      m_data_r[oldSize+i]=value[i+dataOffset];
    }
  }
}


void
DataTagged::addTaggedValue(int tagKey,
                           const DataTypes::ShapeType& pointshape,
                           const DataTypes::CplxVectorType& value,
                           int dataOffset)
{
  if (!DataTypes::checkShape(getShape(), pointshape)) {
    throw DataException(DataTypes::createShapeErrorMessage(
                        "Error - Cannot addTaggedValue due to shape mismatch.", pointshape,getShape()));
  }
  if (!isComplex())
  {
      throw DataException("Programming error - attempt to set a complex value on real data.");
  }
  CHECK_FOR_EX_WRITE
  DataMapType::iterator pos(m_offsetLookup.find(tagKey));
  if (pos!=m_offsetLookup.end()) {
    // tag already exists so use setTaggedValue
    setTaggedValue(tagKey,pointshape, value, dataOffset);
  } else {
    // save the key and the location of its data in the lookup tab
    m_offsetLookup.insert(DataMapType::value_type(tagKey,m_data_c.size()));
    // add the data given in "value" at the end of m_data_c
    // need to make a temp copy of m_data_c, resize m_data_c, then copy
    // all the old values plus the value to be added back into m_data_c
    DataTypes::CplxVectorType m_data_c_temp(m_data_c);
    int oldSize=m_data_c.size();
    int newSize=m_data_c.size()+getNoValues();
    m_data_c.resize(newSize,0.,newSize);
    for (int i=0;i<oldSize;i++) {
      m_data_c[i]=m_data_c_temp[i];
    }
    for (unsigned int i=0;i<getNoValues();i++) {
      m_data_c[oldSize+i]=value[i+dataOffset];
    }
  }
}

void
DataTagged::addTag(int tagKey)
{
  CHECK_FOR_EX_WRITE
  DataMapType::iterator pos(m_offsetLookup.find(tagKey));
  if (pos==m_offsetLookup.end()) {
    if (isComplex())
    {
	// save the key and the location of its data in the lookup tab
	m_offsetLookup.insert(DataMapType::value_type(tagKey,m_data_c.size()));
	// add the data given in "value" at the end of m_data_c
	// need to make a temp copy of m_data_c, resize m_data_c, then copy
	// all the old values plus the value to be added back into m_data_c
	DataTypes::CplxVectorType m_data_c_temp(m_data_c);
	int oldSize=m_data_c.size();
	int newSize=m_data_c.size()+getNoValues();
	m_data_c.resize(newSize,0.,newSize);
	for (int i=0;i<oldSize;i++) {
	  m_data_c[i]=m_data_c_temp[i];
	}
	for (unsigned int i=0;i<getNoValues();i++) {
	  m_data_c[oldSize+i]=m_data_c[m_defaultValueOffset+i];
	}
    }
    else
    {
	// save the key and the location of its data in the lookup tab
	m_offsetLookup.insert(DataMapType::value_type(tagKey,m_data_r.size()));
	// add the data given in "value" at the end of m_data_r
	// need to make a temp copy of m_data_r, resize m_data_r, then copy
	// all the old values plus the value to be added back into m_data_r
	DataTypes::RealVectorType m_data_r_temp(m_data_r);
	int oldSize=m_data_r.size();
	int newSize=m_data_r.size()+getNoValues();
	m_data_r.resize(newSize,0.,newSize);
	for (int i=0;i<oldSize;i++) {
	  m_data_r[i]=m_data_r_temp[i];
	}
	for (unsigned int i=0;i<getNoValues();i++) {
	  m_data_r[oldSize+i]=m_data_r[m_defaultValueOffset+i];
	}
    }
  }
}


DataTypes::real_t*
DataTagged::getSampleDataByTag(int tag, DataTypes::real_t dummy)
{
  CHECK_FOR_EX_WRITE
  DataMapType::iterator pos(m_offsetLookup.find(tag));
  if (pos==m_offsetLookup.end()) {
    // tag couldn't be found so return the default value
    return &(m_data_r[0]);
  } else {
    // return the data-point corresponding to the given tag
    return &(m_data_r[pos->second]);
  }
}

DataTypes::cplx_t*
DataTagged::getSampleDataByTag(int tag, DataTypes::cplx_t dummy)
{
  CHECK_FOR_EX_WRITE
  DataMapType::iterator pos(m_offsetLookup.find(tag));
  if (pos==m_offsetLookup.end()) {
    // tag couldn't be found so return the default value
    return &(m_data_c[0]);
  } else {
    // return the data-point corresponding to the given tag
    return &(m_data_c[pos->second]);
  }
}


bool
DataTagged::hasNaN() const
{
  bool haveNaN=false;
  if (isComplex())
  {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
          if (std::isnan(m_data_c[i].real()) || std::isnan(m_data_c[i].imag()))
          {
              #pragma omp critical 
              {
                  haveNaN=true;
              }
          }
      }
  }
  else
  {
      #pragma omp parallel for
      for (DataTypes::RealVectorType::size_type i=0;i<m_data_r.size();++i)
      {
          if (std::isnan(m_data_r[i]))
          {
              #pragma omp critical 
              {
                  haveNaN=true;
              }
          }
      }
  }
  return haveNaN;
}

void
DataTagged::replaceNaN(double value) {
  CHECK_FOR_EX_WRITE  
  if (isComplex())
  {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
        if (std::isnan(m_data_c[i].real()) || std::isnan(m_data_c[i].imag()))  
        {
          m_data_c[i] = value;
        }
      }
  }
  else
  {
      #pragma omp parallel for
      for (DataTypes::RealVectorType::size_type i=0;i<m_data_r.size();++i)
      {
        if (std::isnan(m_data_r[i]))  
        {
          m_data_r[i] = value;
        }
      }    
  }
}

void
DataTagged::replaceNaN(DataTypes::cplx_t value) {
  CHECK_FOR_EX_WRITE  
  if (isComplex())
  {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
        if (std::isnan(m_data_c[i].real()) || std::isnan(m_data_c[i].imag())) 
        {
          m_data_c[i] = value;
        }
      }
  }
  else
  {
      complicate();
      replaceNaN(value);
  }
}

bool
DataTagged::hasInf() const
{
  bool haveNaN=false;
  if (isComplex())
  {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
          if (std::isinf(m_data_c[i].real()) || std::isinf(m_data_c[i].imag()))
          {
              #pragma omp critical 
              {
                  haveNaN=true;
              }
          }
      }
  }
  else
  {
      #pragma omp parallel for
      for (DataTypes::RealVectorType::size_type i=0;i<m_data_r.size();++i)
      {
          if (std::isinf(m_data_r[i]))
          {
              #pragma omp critical 
              {
                  haveNaN=true;
              }
          }
      }
  }
  return haveNaN;
}

void
DataTagged::replaceInf(double value) {
  CHECK_FOR_EX_WRITE  
  if (isComplex())
  {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
        if (std::isinf(m_data_c[i].real()) || std::isinf(m_data_c[i].imag()))  
        {
          m_data_c[i] = value;
        }
      }
  }
  else
  {
      #pragma omp parallel for
      for (DataTypes::RealVectorType::size_type i=0;i<m_data_r.size();++i)
      {
        if (std::isinf(m_data_r[i]))  
        {
          m_data_r[i] = value;
        }
      }    
  }
}

void
DataTagged::replaceInf(DataTypes::cplx_t value) {
  CHECK_FOR_EX_WRITE  
  if (isComplex())
  {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
        if (std::isinf(m_data_c[i].real()) || std::isinf(m_data_c[i].imag())) 
        {
          m_data_c[i] = value;
        }
      }
  }
  else
  {
      complicate();
      replaceInf(value);
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
  
  const FunctionSpace& fs=getFunctionSpace();
  int usedCount=fs.getNumberOfTagsInUse();
  const int* usedTags=fs.borrowListOfTagsInUse();
  
  if (isComplex())
  {
  
      temp << pointToString(m_data_c,getShape(),getDefaultOffset(),empty) << endl;
      for (i=m_offsetLookup.begin();i!=m_offsetLookup.end();++i) {
        temp << "Tag(" << i->first << ")";
        bool found=false;
        for (int j=0;j<usedCount;++j) {
            if (i->first==usedTags[j]) {
                found=true;
            }
        }
        if (!found) {
            temp << " - Unused";
        }
        temp << endl;
        temp << pointToString(m_data_c,getShape(),i->second,empty) << endl;
      }
  }
  else
  {
  
      temp << pointToString(m_data_r,getShape(),getDefaultOffset(),empty) << endl;
      for (i=m_offsetLookup.begin();i!=m_offsetLookup.end();++i) {
        temp << "Tag(" << i->first << ")";
        bool found=false;
        for (int j=0;j<usedCount;++j) {
            if (i->first==usedTags[j]) {
                found=true;
            }
        }
        if (!found) {
            temp << " - Unused";
        }
        temp << endl;
        temp << pointToString(m_data_r,getShape(),i->second,empty) << endl;
      }    
  }
  return temp.str();
}

DataTypes::RealVectorType::size_type 
DataTagged::getPointOffset(int sampleNo,
                           int dataPointNo) const
{
  int tagKey=getFunctionSpace().getTagFromSampleNo(sampleNo);
  DataMapType::const_iterator pos(m_offsetLookup.find(tagKey));
  DataTypes::RealVectorType::size_type offset=m_defaultValueOffset;
  if (pos!=m_offsetLookup.end()) {
    offset=pos->second;
  }
  return offset;
}

DataTypes::RealVectorType::size_type
DataTagged::getOffsetForTag(int tag) const
{
  DataMapType::const_iterator pos(m_offsetLookup.find(tag));
  DataTypes::RealVectorType::size_type offset=m_defaultValueOffset;
  if (pos!=m_offsetLookup.end()) {
    offset=pos->second;
  }
  return offset;
}

DataTypes::RealVectorType::const_reference
DataTagged::getDataByTagRO(int tag, DataTypes::RealVectorType::size_type i, DataTypes::real_t dummy) const
{
  DataMapType::const_iterator pos(m_offsetLookup.find(tag));
  DataTypes::RealVectorType::size_type offset=m_defaultValueOffset;
  if (pos!=m_offsetLookup.end()) {
    offset=pos->second;
  }
  return m_data_r[offset+i];
}

DataTypes::RealVectorType::reference
DataTagged::getDataByTagRW(int tag, DataTypes::RealVectorType::size_type i, DataTypes::real_t dummy)
{
  CHECK_FOR_EX_WRITE
  DataMapType::const_iterator pos(m_offsetLookup.find(tag));
  DataTypes::RealVectorType::size_type offset=m_defaultValueOffset;
  if (pos!=m_offsetLookup.end()) {
    offset=pos->second;
  }
  return m_data_r[offset+i];
}

DataTypes::CplxVectorType::const_reference
DataTagged::getDataByTagRO(int tag, DataTypes::RealVectorType::size_type i, DataTypes::cplx_t dummy) const
{
  DataMapType::const_iterator pos(m_offsetLookup.find(tag));
  DataTypes::CplxVectorType::size_type offset=m_defaultValueOffset;
  if (pos!=m_offsetLookup.end()) {
    offset=pos->second;
  }
  return m_data_c[offset+i];
}

DataTypes::CplxVectorType::reference
DataTagged::getDataByTagRW(int tag, DataTypes::RealVectorType::size_type i, DataTypes::cplx_t dummy)
{
  CHECK_FOR_EX_WRITE
  DataMapType::const_iterator pos(m_offsetLookup.find(tag));
  DataTypes::CplxVectorType::size_type offset=m_defaultValueOffset;
  if (pos!=m_offsetLookup.end()) {
    offset=pos->second;
  }
  return m_data_c[offset+i];
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
  const ShapeType& evShape=temp_ev->getShape();

  if (isComplex())
  {
      DataTypes::CplxVectorType& evVec=temp_ev->getTypedVectorRW((DataTypes::cplx_t)0);
      for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
	  temp_ev->addTag(i->first);
	  DataTypes::CplxVectorType::size_type offset=getOffsetForTag(i->first);
	  DataTypes::CplxVectorType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
	  escript::symmetric(m_data_c,getShape(),offset,evVec, evShape, evoffset);
      }
      escript::symmetric(m_data_c,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset());      
  }
  else
  {
      DataTypes::RealVectorType& evVec=temp_ev->getTypedVectorRW(0.0);
      for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
	  temp_ev->addTag(i->first);
	  DataTypes::RealVectorType::size_type offset=getOffsetForTag(i->first);
	  DataTypes::RealVectorType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
	  escript::symmetric(m_data_r,getShape(),offset,evVec, evShape, evoffset);
      }
      escript::symmetric(m_data_r,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset());      
  }
}


void
DataTagged::antisymmetric(DataAbstract* ev)
{
  DataTagged* temp_ev=dynamic_cast<DataTagged*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataTagged::antisymmetric casting to DataTagged failed (probably a programming error).");
  }
  const DataTagged::DataMapType& thisLookup=getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
  const ShapeType& evShape=temp_ev->getShape();
  if (isComplex())
  {
      DataTypes::CplxVectorType& evVec=temp_ev->getTypedVectorRW((DataTypes::cplx_t)0);
      for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
	  temp_ev->addTag(i->first);
	  DataTypes::CplxVectorType::size_type offset=getOffsetForTag(i->first);
	  DataTypes::CplxVectorType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
	  escript::antisymmetric(m_data_c,getShape(),offset,evVec, evShape, evoffset);
      }
      escript::antisymmetric(m_data_c,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset());      
  }
  else
  {
      DataTypes::RealVectorType& evVec=temp_ev->getTypedVectorRW(0.0);
      for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
	  temp_ev->addTag(i->first);
	  DataTypes::RealVectorType::size_type offset=getOffsetForTag(i->first);
	  DataTypes::RealVectorType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
	  escript::antisymmetric(m_data_r,getShape(),offset,evVec, evShape, evoffset);
      }
      escript::antisymmetric(m_data_r,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset());      
  }  
}

void
DataTagged::hermitian(DataAbstract* ev)
{
  DataTagged* temp_ev=dynamic_cast<DataTagged*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataTagged::hermitian casting to DataTagged failed (probably a programming error).");
  }
  if (!isComplex() || !temp_ev->isComplex())
  {
      throw DataException("DataTagged::hermitian: do not call this method with real data");
  }  
  
  const DataTagged::DataMapType& thisLookup=getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
  const ShapeType& evShape=temp_ev->getShape();

  DataTypes::CplxVectorType& evVec=temp_ev->getTypedVectorRW((DataTypes::cplx_t)0);
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTag(i->first);
      DataTypes::CplxVectorType::size_type offset=getOffsetForTag(i->first);
      DataTypes::CplxVectorType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
      escript::hermitian(m_data_c,getShape(),offset,evVec, evShape, evoffset);
  }
  escript::hermitian(m_data_c,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset());      
}


void
DataTagged::antihermitian(DataAbstract* ev)
{
  DataTagged* temp_ev=dynamic_cast<DataTagged*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataTagged::antihermitian casting to DataTagged failed (probably a programming error).");
  }
  if (!isComplex() || !temp_ev->isComplex())
  {
      throw DataException("DataTagged::antihermitian: do not call this method with real data");
  }  
  const DataTagged::DataMapType& thisLookup=getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
  const ShapeType& evShape=temp_ev->getShape();
  DataTypes::CplxVectorType& evVec=temp_ev->getTypedVectorRW((DataTypes::cplx_t)0);
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTag(i->first);
      DataTypes::CplxVectorType::size_type offset=getOffsetForTag(i->first);
      DataTypes::CplxVectorType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
      escript::antihermitian(m_data_c,getShape(),offset,evVec, evShape, evoffset);
  }
  escript::antihermitian(m_data_c,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset());      
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
  const ShapeType& evShape=temp_ev->getShape();
  if (isComplex())
  {
      DataTypes::CplxVectorType& evVec=temp_ev->getVectorRWC();  
      for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
	  temp_ev->addTag(i->first);
	  DataTypes::CplxVectorType::size_type offset=getOffsetForTag(i->first);
	  DataTypes::CplxVectorType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
	  escript::trace(m_data_c,getShape(),offset,evVec, evShape, evoffset, axis_offset);
      }
      escript::trace(m_data_c,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset(),axis_offset);
  }
  else
  {
      DataTypes::RealVectorType& evVec=temp_ev->getVectorRW();  
      for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
	  temp_ev->addTag(i->first);
	  DataTypes::RealVectorType::size_type offset=getOffsetForTag(i->first);
	  DataTypes::RealVectorType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
	  escript::trace(m_data_r,getShape(),offset,evVec, evShape, evoffset, axis_offset);
      }
      escript::trace(m_data_r,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset(),axis_offset);
  }
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
  const ShapeType& evShape=temp_ev->getShape();
  if (isComplex())
  {
      DataTypes::CplxVectorType& evVec=temp_ev->getVectorRWC();  
      for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
	  temp_ev->addTag(i->first);
	  DataTypes::CplxVectorType::size_type offset=getOffsetForTag(i->first);
	  DataTypes::CplxVectorType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
	  escript::transpose(m_data_c,getShape(),offset,evVec, evShape, evoffset, axis_offset);
      }
      escript::transpose(m_data_c,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset(),axis_offset);
  }
  else
  {
      DataTypes::RealVectorType& evVec=temp_ev->getVectorRW();  
      for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
	  temp_ev->addTag(i->first);
	  DataTypes::RealVectorType::size_type offset=getOffsetForTag(i->first);
	  DataTypes::RealVectorType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
	  escript::transpose(m_data_r,getShape(),offset,evVec, evShape, evoffset, axis_offset);
      }
      escript::transpose(m_data_r,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset(),axis_offset);
  }
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
  const ShapeType& evShape=temp_ev->getShape();
  if (isComplex())
  {
      DataTypes::CplxVectorType& evVec=temp_ev->getVectorRWC();  
      for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
	  temp_ev->addTag(i->first);
	  DataTypes::CplxVectorType::size_type offset=getOffsetForTag(i->first);
	  DataTypes::CplxVectorType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
	  escript::swapaxes(m_data_c,getShape(),offset,evVec, evShape, evoffset,axis0,axis1);
      }
      escript::swapaxes(m_data_c,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset(),axis0,axis1);    
  }
  else
  {
      DataTypes::RealVectorType& evVec=temp_ev->getVectorRW();  
      for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
	  temp_ev->addTag(i->first);
	  DataTypes::RealVectorType::size_type offset=getOffsetForTag(i->first);
	  DataTypes::RealVectorType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
	  escript::swapaxes(m_data_r,getShape(),offset,evVec, evShape, evoffset,axis0,axis1);
      }
      escript::swapaxes(m_data_r,getShape(),getDefaultOffset(),evVec,evShape,temp_ev->getDefaultOffset(),axis0,axis1);
  }
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
  const ShapeType& evShape=temp_ev->getShape();
  if (isComplex())
  {
      DataTypes::CplxVectorType& evVec=temp_ev->getVectorRWC();

      for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
	  temp_ev->addTag(i->first);
	  DataTypes::CplxVectorType::size_type offset=getOffsetForTag(i->first);
	  DataTypes::CplxVectorType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
	  escript::eigenvalues(m_data_c,getShape(),offset,evVec, evShape, evoffset);
      }
      escript::eigenvalues(m_data_c,getShape(),getDefaultOffset(),evVec, evShape, temp_ev->getDefaultOffset());
  }
  else
  {
      DataTypes::RealVectorType& evVec=temp_ev->getVectorRW();

      for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
	  temp_ev->addTag(i->first);
	  DataTypes::RealVectorType::size_type offset=getOffsetForTag(i->first);
	  DataTypes::RealVectorType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
	  escript::eigenvalues(m_data_r,getShape(),offset,evVec, evShape, evoffset);
      }
      escript::eigenvalues(m_data_r,getShape(),getDefaultOffset(),evVec, evShape, temp_ev->getDefaultOffset());
  }
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
  DataTypes::RealVectorType& evVec=temp_ev->getVectorRW();
  const ShapeType& evShape=temp_ev->getShape();
  DataTypes::RealVectorType& VVec=temp_V->getVectorRW();
  const ShapeType& VShape=temp_V->getShape();
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTag(i->first);
      temp_V->addTag(i->first);
/*      DataArrayView thisView=getDataPointByTag(i->first);
      DataArrayView evView=temp_ev->getDataPointByTag(i->first);
      DataArrayView VView=temp_V->getDataPointByTag(i->first);*/
      DataTypes::RealVectorType::size_type offset=getOffsetForTag(i->first);
      DataTypes::RealVectorType::size_type evoffset=temp_ev->getOffsetForTag(i->first);
      DataTypes::RealVectorType::size_type Voffset=temp_V->getOffsetForTag(i->first);
/*      DataArrayView::eigenvalues_and_eigenvectors(thisView,0,evView,0,VView,0,tol);*/
      escript::eigenvalues_and_eigenvectors(m_data_r,getShape(),offset,evVec, evShape, evoffset,VVec,VShape,Voffset,tol);

  }
  escript::eigenvalues_and_eigenvectors(m_data_r,getShape(),getDefaultOffset(),evVec, evShape,
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
  DataTypes::RealVectorType& outVec=temp->getVectorRW();
  const ShapeType& outShape=temp->getShape();
  LapackInverseHelper h(getShape()[0]);
  int err=0;
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp->addTag(i->first);
      DataTypes::RealVectorType::size_type inoffset=getOffsetForTag(i->first);
      DataTypes::RealVectorType::size_type outoffset=temp->getOffsetForTag(i->first);

      err=escript::matrix_inverse(m_data_r, getShape(), inoffset, outVec, outShape, outoffset, 1, h);
      if (!err) break;
  }
  if (!err)
  {
      escript::matrix_inverse(m_data_r, getShape(), getDefaultOffset(), outVec, outShape, temp->getDefaultOffset(), 1, h);
  }
  return err;
}

void
DataTagged::setToZero(){
    CHECK_FOR_EX_WRITE
    if (isComplex())
    {
        DataTypes::CplxVectorType::size_type n=m_data_c.size();
        for (int i=0; i<n ;++i) m_data_c[i]=0.;
    }
    else 
    {    
        DataTypes::RealVectorType::size_type n=m_data_r.size();
        for (int i=0; i<n ;++i) m_data_r[i]=0.;
    }
}

#ifdef ESYS_HAVE_HDF5
void DataTagged::dump_hdf5(const H5::Group h5_grp) const
{
    uint rank = getRank();
    int fs_type=  getFunctionSpace().getTypeCode();
    const DataTypes::ShapeType& shape = getShape();
    if (isComplex())
    {
        throw DataException("Error - DataTagged::dump_hdf5: complex data are not supported. Split into real and imaginary part.");
    }


#ifdef ESYS_MPI
   /* Serialize I/O */
   JMPI mpiInfo = getFunctionSpace().getDomain()->getMPI();
   const int mpi_iam = mpiInfo->rank;
   const int mpi_num = mpiInfo->size;
   MPI_Comm comm = mpiInfo->comm;
   MPI_Status status;
   int dummy = 0;
   if (mpi_iam > 0)
       MPI_Recv(&dummy, 1, MPI_INT, mpi_iam-1, 81802, comm, &status);
#endif
    try
    {
        // .... add meta data ............
        uint h5_shape[DataTypes::maxRank]; // dataset dimensions
        for (uint i = 0; i < rank; i++) {
            h5_shape[i]= shape[i];
        }
        hsize_t h5_shape_dims[1] = {rank};
        H5::DataSet h5_dmeta = h5_grp.createDataSet("Meta", H5::PredType::NATIVE_UINT, H5::DataSpace(1, h5_shape_dims ) );
        h5_dmeta.write( h5_shape, H5::PredType::NATIVE_UINT);
        // data type
        hsize_t h5_typeid_dims[1] = { 1 };
        int h5_type_id[1] = { 1 };
        H5::Attribute h5_typeid_attr = h5_dmeta.createAttribute("type_id", H5::PredType::NATIVE_INT, H5::DataSpace(1,h5_typeid_dims ) );
        h5_typeid_attr.write( H5::PredType::NATIVE_INT , h5_type_id );

        uint h5_rank[1] = { rank };
        H5::Attribute h5_rank_attr = h5_dmeta.createAttribute("rank", H5::PredType::NATIVE_UINT, H5::DataSpace(1,h5_typeid_dims ) );
        h5_rank_attr.write( H5::PredType::NATIVE_UINT , h5_rank );

        int dfs_type[1] = { fs_type };
        H5::Attribute h5_fs_type_attr = h5_dmeta.createAttribute("function_space_type", H5::PredType::NATIVE_INT, H5::DataSpace(1,h5_typeid_dims) );
        h5_fs_type_attr.write( H5::PredType::NATIVE_INT , dfs_type );
        // .... end meta data ............
        // collect tags:
        const DataTagged::DataMapType& thisLookup=getTagLookup();
        DataTagged::DataMapType::const_iterator i;
        DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
        std::vector<int> tags;
        tags.push_back(-1);
        for (i=thisLookup.begin();i!=thisLookupEnd;i++)
            tags.push_back(i->first);
        // ... add list of tags ...
        hsize_t h5_tag_length[1] = { static_cast<hsize_t>(tags.size()) };
        H5::DataSet h5_dataset_tags = h5_grp.createDataSet("tags", H5::PredType::NATIVE_INT, H5::DataSpace(1 , h5_tag_length ) );
        h5_dataset_tags.write(&(tags[0]), H5::PredType::NATIVE_INT);
        // ... add data ....
        hsize_t h5_data_length[1] = { static_cast<hsize_t>(m_data_r.size()) };
        const double* d_ptr=&(m_data_r[0]);
        H5::DataSet h5_dataset_data = h5_grp.createDataSet("data", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1 , h5_data_length ) );
        h5_dataset_data.write(d_ptr, H5::PredType::NATIVE_DOUBLE);
    }
    // catch failure caused by the H5File operations
    catch (H5::Exception& error)
    {
        #ifdef ESYS_MPI
            if ( mpi_iam < mpi_num-1 ) MPI_Send(&dummy, 1, MPI_INT, mpi_iam+1, 81802, comm);
        #endif
        error.printErrorStack();
        throw DataException("Error - DataConstant:: creating HDF5 file failed.");
    }
    #ifdef ESYS_MPI
       if ( mpi_iam < mpi_num-1 ) MPI_Send(&dummy, 1, MPI_INT, mpi_iam+1, 81802, comm);
    #endif

}
#endif

DataTypes::RealVectorType&
DataTagged::getVectorRW()
{
    CHECK_FOR_EX_WRITE
    return m_data_r;
}

const DataTypes::RealVectorType&
DataTagged::getVectorRO() const
{
        return m_data_r;
}

DataTypes::CplxVectorType&
DataTagged::getVectorRWC()
{
    CHECK_FOR_EX_WRITE
    return m_data_c;
}

const DataTypes::CplxVectorType&
DataTagged::getVectorROC() const
{
        return m_data_c;
}

DataTypes::RealVectorType&
DataTagged::getTypedVectorRW(DataTypes::real_t dummy)
{
  CHECK_FOR_EX_WRITE
  return m_data_r;
}

const DataTypes::RealVectorType&
DataTagged::getTypedVectorRO(DataTypes::real_t dummy) const
{
  return m_data_r;
}

DataTypes::CplxVectorType&
DataTagged::getTypedVectorRW(DataTypes::cplx_t dummy)
{
  CHECK_FOR_EX_WRITE
  return m_data_c;
}

const DataTypes::CplxVectorType&
DataTagged::getTypedVectorRO(DataTypes::cplx_t dummy) const
{
  return m_data_c;
}

size_t
DataTagged::getTagCount() const
{
    return m_offsetLookup.size();
}


void DataTagged::complicate()
{
    if (!isComplex())
    {
        fillComplexFromReal(m_data_r, m_data_c);
        this->m_iscompl=true;
        m_data_r.resize(0,0,1);
    }
}

}  // end of namespace

