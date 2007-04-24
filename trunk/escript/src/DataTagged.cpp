// $Id$

/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

#include "DataTagged.h"

#include "DataConstant.h"
#include "DataException.h"
#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

using namespace std;

namespace escript {

DataTagged::DataTagged()
  : DataAbstract(FunctionSpace())
{
  // default constructor

  // create a scalar default value
  m_data.resize(1,0.,1);
  DataArrayView temp(m_data,DataArrayView::ShapeType());
  setPointDataView(temp);
}

DataTagged::DataTagged(const TagListType& tagKeys, 
		       const ValueListType& values,
		       const DataArrayView& defaultValue,
		       const FunctionSpace& what)
  : DataAbstract(what)
{
  // constructor

  // initialise the array of data values
  // the default value is always the first item in the values list
  int len = defaultValue.noValues();
  m_data.resize(len,0.,len);
  for (int i=0; i<defaultValue.noValues(); i++) {
    m_data[i]=defaultValue.getData(i);
  }

  // create the data view
  DataArrayView temp(m_data,defaultValue.getShape());
  setPointDataView(temp);

  // add remaining tags and values
  addTaggedValues(tagKeys,values);
}

DataTagged::DataTagged(const FunctionSpace& what,
                       const DataArrayView::ShapeType &shape,
                       const int tags[],
                       const ValueType& data)
  : DataAbstract(what)
{
  // alternative constructor
  // not unit_tested tested yet

  // copy the data
  m_data=data;

  // create the view of the data
  DataArrayView tempView(m_data,shape);
  setPointDataView(tempView);

  // create the tag lookup map
  for (int sampleNo=0; sampleNo<getNumSamples(); sampleNo++) {
    m_offsetLookup.insert(DataMapType::value_type(sampleNo,tags[sampleNo]));
  }
}

DataTagged::DataTagged(const FunctionSpace& what,
                       const DataArrayView::ShapeType &shape,
                       const TagListType& tags,
                       const ValueType& data)
  : DataAbstract(what)
{
  // alternative constructor
  // not unit_tested tested yet

  // copy the data
  m_data=data;

  // create the view of the data
  DataArrayView tempView(m_data,shape);
  setPointDataView(tempView);

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
  // copy constructor

  // create the data view
  DataArrayView temp(m_data,other.getPointDataView().getShape());
  setPointDataView(temp);
}

DataTagged::DataTagged(const DataConstant& other)
  : DataAbstract(other.getFunctionSpace())
{
  // copy constructor

  // fill the default value with the constant value item from "other"
  const DataArrayView& value=other.getPointDataView();
  int len = value.noValues();
  m_data.resize(len,0.,len);
  for (int i=0; i<value.noValues(); i++) {
    m_data[i]=value.getData(i);
  }

  // create the data view
  DataArrayView temp(m_data,value.getShape());
  setPointDataView(temp);
}

DataAbstract*
DataTagged::getSlice(const DataArrayView::RegionType& region) const 
{
  return new DataTagged(*this, region);
}

DataTagged::DataTagged(const DataTagged& other, 
		       const DataArrayView::RegionType& region)
  : DataAbstract(other.getFunctionSpace())
{
  // slice constructor

  // get the shape of the slice to copy from other
  DataArrayView::ShapeType regionShape(DataArrayView::getResultSliceShape(region));
  DataArrayView::RegionLoopRangeType regionLoopRange=getSliceRegionLoopRange(region);

  // allocate enough space in this for all values
  // (need to add one to allow for the default value)
  int len = DataArrayView::noValues(regionShape)*(other.m_offsetLookup.size()+1);
  m_data.resize(len,0.0,len);

  // create the data view
  DataArrayView temp(m_data,regionShape);
  setPointDataView(temp);

  // copy the default value from other to this
  getDefaultValue().copySlice(other.getDefaultValue(), regionLoopRange);

  // loop through the tag values copying these
  DataMapType::const_iterator pos;
  DataArrayView::ValueType::size_type tagOffset=getPointDataView().noValues();
  for (pos=other.m_offsetLookup.begin();pos!=other.m_offsetLookup.end();pos++){
    getPointDataView().copySlice(tagOffset,other.getPointDataView(),pos->second,regionLoopRange);
    m_offsetLookup.insert(DataMapType::value_type(pos->first,tagOffset));
    tagOffset+=getPointDataView().noValues();
  }
}

void
DataTagged::setSlice(const DataAbstract* other,
                     const DataArrayView::RegionType& region)
{

  // other must be another DataTagged object
  // Data:setSlice implementation should ensure this
  const DataTagged* otherTemp=dynamic_cast<const DataTagged*>(other);
  if (otherTemp==0) {
    throw DataException("Programming error - casting to DataTagged.");
  }

  // determine shape of the specified region
  DataArrayView::ShapeType regionShape(DataArrayView::getResultSliceShape(region));

  // modify region specification as needed to match rank of this object
  DataArrayView::RegionLoopRangeType regionLoopRange=getSliceRegionLoopRange(region);

  // ensure rank/shape of this object is compatible with specified region
  if (getPointDataView().getRank()!=region.size()) {
    throw DataException("Error - Invalid slice region.");
  }
  if (otherTemp->getPointDataView().getRank()>0 && !other->getPointDataView().checkShape(regionShape)) {
    throw DataException (other->getPointDataView().createShapeErrorMessage(
                         "Error - Couldn't copy slice due to shape mismatch.",regionShape));
  }

  // copy slice from other default value to this default value
  getDefaultValue().copySliceFrom(otherTemp->getDefaultValue(), regionLoopRange);

  // loop through tag values in other, adding any which aren't in this, using default value
  DataMapType::const_iterator pos;
  for (pos=otherTemp->m_offsetLookup.begin();pos!=otherTemp->m_offsetLookup.end();pos++) {
    if (!isCurrentTag(pos->first)) {
      addTaggedValue(pos->first,getDefaultValue());
    }
  }

  // loop through the tag values copying slices from other to this
  for (pos=m_offsetLookup.begin();pos!=m_offsetLookup.end();pos++) {
    getDataPointByTag(pos->first).copySliceFrom(otherTemp->getDataPointByTag(pos->first), regionLoopRange);
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
DataTagged::setTaggedValues(const TagListType& tagKeys,
                            const ValueListType& values)
{
  addTaggedValues(tagKeys,values);
}

void
DataTagged::setTaggedValue(int tagKey,
                           const DataArrayView& value)
{
  if (!getPointDataView().checkShape(value.getShape())) {
      throw DataException(getPointDataView().createShapeErrorMessage(
                          "Error - Cannot setTaggedValue due to shape mismatch.", value.getShape()));
  }
  DataMapType::iterator pos(m_offsetLookup.find(tagKey));
  if (pos==m_offsetLookup.end()) {
    // tag couldn't be found so use addTaggedValue
    addTaggedValue(tagKey,value);
  } else {
    // copy the values into the data array at the offset determined by m_offsetLookup
    int offset=pos->second;
    for (int i=0; i<getPointDataView().noValues(); i++) {
      m_data[offset+i]=value.getData(i);
    }
  }
}

void
DataTagged::addTaggedValues(const TagListType& tagKeys,
                            const ValueListType& values)
{
  if (values.size()==0) {
    // copy the current default value for each of the tags
    TagListType::const_iterator iT;
    for (iT=tagKeys.begin();iT!=tagKeys.end();iT++) {
      // the point data view for DataTagged points at the default value
      addTaggedValue(*iT,getPointDataView());
    }
  } else if (values.size()==1 && tagKeys.size()>1) {
    // assume the one given value will be used for all tag values
    TagListType::const_iterator iT;
    for (iT=tagKeys.begin();iT!=tagKeys.end();iT++) {
      addTaggedValue(*iT,values[0]);
    }
  } else {
    if (tagKeys.size()!=values.size()) {
      stringstream temp;
      temp << "Error - (addTaggedValue) Number of tags: " << tagKeys.size()
	   << " doesn't match number of values: " << values.size();
      throw DataException(temp.str());
    } else {
      for (int i=0;i<tagKeys.size();i++) {
        addTaggedValue(tagKeys[i],values[i]);
      }
    }
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
  DataMapType::iterator pos(m_offsetLookup.find(tagKey));
  if (pos!=m_offsetLookup.end()) {
    // tag already exists so use setTaggedValue
    setTaggedValue(tagKey,value);
  } else {
    // save the key and the location of its data in the lookup tab
    m_offsetLookup.insert(DataMapType::value_type(tagKey,m_data.size()));
    // add the data given in "value" at the end of m_data
    // need to make a temp copy of m_data, resize m_data, then copy
    // all the old values plus the value to be added back into m_data
    ValueType m_data_temp(m_data);
    int oldSize=m_data.size();
    int newSize=m_data.size()+value.noValues();
    m_data.resize(newSize,0.,newSize);
    for (int i=0;i<oldSize;i++) {
      m_data[i]=m_data_temp[i];
    }
    for (int i=0;i<value.noValues();i++) {
      m_data[oldSize+i]=value.getData(i);
    }
  }
}

double*
DataTagged::getSampleDataByTag(int tag)
{
  DataMapType::iterator pos(m_offsetLookup.find(tag));
  if (pos==m_offsetLookup.end()) {
    // tag couldn't be found so return the default value
    return &(m_data[0]);
  } else {
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
  // create a temporary view as the offset will be changed
  DataArrayView tempView(getPointDataView().getData(), getPointDataView().getShape());
  for (i=m_offsetLookup.begin();i!=m_offsetLookup.end();++i) {
    temp << "Tag(" << i->first << ")" << endl;
    tempView.setOffset(i->second);
    temp << tempView.toString() << endl;
  }
  return temp.str();
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

DataArrayView
DataTagged::getDataPoint(int sampleNo,
                         int dataPointNo)
{
  EsysAssert(validSampleNo(sampleNo),"(getDataPoint) Invalid sampleNo: " << sampleNo);
  int tagKey=getFunctionSpace().getTagFromSampleNo(sampleNo);
  return getDataPointByTag(tagKey);
}

int
DataTagged::archiveData(ofstream& archiveFile,
                        const DataArrayView::ValueType::size_type noValues) const
{
  return(m_data.archiveData(archiveFile, noValues));
}

int
DataTagged::extractData(ifstream& archiveFile,
                        const DataArrayView::ValueType::size_type noValues)
{
  return(m_data.extractData(archiveFile, noValues));
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
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTaggedValue(i->first,temp_ev->getDefaultValue());
      DataArrayView thisView=getDataPointByTag(i->first);
      DataArrayView evView=temp_ev->getDataPointByTag(i->first);
      DataArrayView::symmetric(thisView,0,evView,0);
  }
  DataArrayView::symmetric(getDefaultValue(),0,temp_ev->getDefaultValue(),0);
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
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTaggedValue(i->first,temp_ev->getDefaultValue());
      DataArrayView thisView=getDataPointByTag(i->first);
      DataArrayView evView=temp_ev->getDataPointByTag(i->first);
      DataArrayView::nonsymmetric(thisView,0,evView,0);
  }
  DataArrayView::nonsymmetric(getDefaultValue(),0,temp_ev->getDefaultValue(),0);
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
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTaggedValue(i->first,temp_ev->getDefaultValue());
      DataArrayView thisView=getDataPointByTag(i->first);
      DataArrayView evView=temp_ev->getDataPointByTag(i->first);
      DataArrayView::trace(thisView,0,evView,0, axis_offset);
  }
  DataArrayView::trace(getDefaultValue(),0,temp_ev->getDefaultValue(),0,axis_offset);
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
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTaggedValue(i->first,temp_ev->getDefaultValue());
      DataArrayView thisView=getDataPointByTag(i->first);
      DataArrayView evView=temp_ev->getDataPointByTag(i->first);
      DataArrayView::transpose(thisView,0,evView,0, axis_offset);
  }
  DataArrayView::transpose(getDefaultValue(),0,temp_ev->getDefaultValue(),0,axis_offset);
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
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTaggedValue(i->first,temp_ev->getDefaultValue());
      DataArrayView thisView=getDataPointByTag(i->first);
      DataArrayView evView=temp_ev->getDataPointByTag(i->first);
      DataArrayView::swapaxes(thisView,0,evView,0,axis0,axis1);
  }
  DataArrayView::swapaxes(getDefaultValue(),0,temp_ev->getDefaultValue(),0,axis0,axis1);
}

void
DataTagged::eigenvalues(DataAbstract* ev)
{
  DataTagged* temp_ev=dynamic_cast<DataTagged*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataTagged::eigenvalues casting to DataTagged failed (propably a programming error).");
  }
  const DataTagged::DataMapType& thisLookup=getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTaggedValue(i->first,temp_ev->getDefaultValue());
      DataArrayView thisView=getDataPointByTag(i->first);
      DataArrayView evView=temp_ev->getDataPointByTag(i->first);
      DataArrayView::eigenvalues(thisView,0,evView,0);
  }
  DataArrayView::eigenvalues(getDefaultValue(),0,temp_ev->getDefaultValue(),0);
}
void
DataTagged::eigenvalues_and_eigenvectors(DataAbstract* ev,DataAbstract* V,const double tol)
{
  DataTagged* temp_ev=dynamic_cast<DataTagged*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataTagged::eigenvalues_and_eigenvectors casting to DataTagged failed (propably a programming error).");
  }
  DataTagged* temp_V=dynamic_cast<DataTagged*>(V);
  if (temp_V==0) {
    throw DataException("Error - DataTagged::eigenvalues_and_eigenvectors casting to DataTagged failed (propably a programming error).");
  }
  const DataTagged::DataMapType& thisLookup=getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
  for (i=thisLookup.begin();i!=thisLookupEnd;i++) {
      temp_ev->addTaggedValue(i->first,temp_ev->getDefaultValue());
      temp_V->addTaggedValue(i->first,temp_V->getDefaultValue());
      DataArrayView thisView=getDataPointByTag(i->first);
      DataArrayView evView=temp_ev->getDataPointByTag(i->first);
      DataArrayView VView=temp_V->getDataPointByTag(i->first);
      DataArrayView::eigenvalues_and_eigenvectors(thisView,0,evView,0,VView,0,tol);
  }
  DataArrayView::eigenvalues_and_eigenvectors(getDefaultValue(),0,
                                              temp_ev->getDefaultValue(),0,
                                              temp_V->getDefaultValue(),0,
                                              tol);


}

void
DataTagged::setToZero(){
    DataArrayView::ValueType::size_type n=m_data.size();
    for (int i=0; i<n ;++i) m_data[i]=0.;
}

void
DataTagged::dump(const std::string fileName) const
{
   #ifdef PASO_MPI
   throw DataException("Error - DataTagged:: dump is not implemented for MPI yet.")
   #endif
   #ifdef USE_NETCDF
   const int ldims=DataArrayView::maxRank+1;
   const NcDim* ncdims[ldims];
   NcVar *var, *tags_var;
   int rank = getPointDataView().getRank();
   int type=  getFunctionSpace().getTypeCode();
   int ndims =0;
   long dims[ldims];
   DataArrayView::ShapeType shape = getPointDataView().getShape();

   // netCDF error handler
   NcError err(NcError::verbose_nonfatal);
   // Create the file.
   NcFile dataFile(fileName.c_str(), NcFile::Replace);
   // check if writing was successful
   if (!dataFile.is_valid())
        throw DataException("Error - DataTagged:: opening of netCDF file for output failed.");
   if (!dataFile.add_att("type","tagged") )
        throw DataException("Error - DataTagged:: appending data type to netCDF file failed.");
   if (!dataFile.add_att("rank",rank) )
        throw DataException("Error - DataTagged:: appending rank attribute to netCDF file failed.");
   if (!dataFile.add_att("function_space_type",type))
        throw DataException("Error - DataTagged:: appending function space attribute to netCDF file failed.");
   ndims=rank+1;
   if ( rank >0 ) {
       dims[0]=shape[0];
       if (! (ncdims[0] = dataFile.add_dim("d0",shape[0])) )
            throw DataException("Error - DataTagged:: appending ncdimsion 0 to netCDF file failed.");
   }
   if ( rank >1 ) {
       dims[1]=shape[1];
       if (! (ncdims[1] = dataFile.add_dim("d1",shape[1])) )
            throw DataException("Error - DataTagged:: appending ncdimsion 1 to netCDF file failed.");
   }
   if ( rank >2 ) {
       dims[2]=shape[2];
       if (! (ncdims[2] = dataFile.add_dim("d2", shape[2])) )
            throw DataException("Error - DataTagged:: appending ncdimsion 2 to netCDF file failed.");
   }
   if ( rank >3 ) {
       dims[3]=shape[3];
       if (! (ncdims[3] = dataFile.add_dim("d3", shape[3])) )
            throw DataException("Error - DataTagged:: appending ncdimsion 3 to netCDF file failed.");
   }
   const DataTagged::DataMapType& thisLookup=getTagLookup();
   DataTagged::DataMapType::const_iterator i;
   DataTagged::DataMapType::const_iterator thisLookupEnd=thisLookup.end();
   int ntags=1;
   for (i=thisLookup.begin();i!=thisLookupEnd;i++) ntags++;
   int* tags =(int*) malloc(ntags*sizeof(int));
   int c=1;
   tags[0]=-1;
   for (i=thisLookup.begin();i!=thisLookupEnd;i++) tags[c++]=i->first;
   dims[rank]=ntags;
   if (! (ncdims[rank] = dataFile.add_dim("num_tags", dims[rank])) )
   {
	   free(tags);
           throw DataException("Error - DataTagged:: appending num_tags to netCDF file failed.");
   }
   if (! ( tags_var = dataFile.add_var("tags", ncInt, ncdims[rank])) )
   {
	free(tags);
        throw DataException("Error - DataTagged:: appending tags to netCDF file failed.");
   }
   if (! (tags_var->put(tags,dims[rank])) )
   {
	free(tags);
        throw DataException("Error - DataTagged:: copy tags to netCDF buffer failed.");
   }
   if (! ( var = dataFile.add_var("data", ncDouble, ndims, ncdims)) )
   {
	free(tags);
        throw DataException("Error - DataTagged:: appending variable to netCDF file failed.");
   }
   if (! (var->put(&m_data[0],dims)) )
   {
	free(tags);
        throw DataException("Error - DataTagged:: copy data to netCDF buffer failed.");
   }
   #else
   throw DataException("Error - DataTagged:: dump is not configured with netCDF. Please contact your installation manager.");
   #endif
}
}  // end of namespace
