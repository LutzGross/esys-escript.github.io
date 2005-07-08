//$Id$
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

#include "escript/Data/DataConstant.h"
#include "escript/Data/DataException.h"
#include "esysUtils/EsysAssert.h"

#include <iostream>
#include <boost/python/extract.hpp>

using namespace std;

namespace escript {

DataConstant::DataConstant(const boost::python::numeric::array& value,
                           const FunctionSpace& what)
  : DataAbstract(what)
{
  DataArray temp(value);
  //
  // copy the data in the correct format
  m_data=temp.getData();
  //
  // create the view of the data
  DataArrayView tempView(m_data,temp.getView().getShape());
  setPointDataView(tempView);
}

DataConstant::DataConstant(const DataArrayView& value,
                           const FunctionSpace& what)
  : DataAbstract(what)
{
  //
  // copy the data in the correct format
  m_data=value.getData();
  //
  // create the view of the data
  DataArrayView tempView(m_data,value.getShape());
  setPointDataView(tempView);
}

DataConstant::DataConstant(const DataConstant& other)
  : DataAbstract(other.getFunctionSpace())
{
  // 
  // copy the data in the correct format
  m_data=other.m_data;
  //
  // create the view of the data
  DataArrayView tempView(m_data,other.getPointDataView().getShape());
  setPointDataView(tempView);
}

DataConstant::DataConstant(const DataConstant& other,
                           const DataArrayView::RegionType& region)
  : DataAbstract(other.getFunctionSpace())
{
  //
  // get the shape of the slice to copy from
  DataArrayView::ShapeType shape(DataArrayView::getResultSliceShape(region));
  //
  // allocate space for this new DataConstant's data
  m_data.resize(DataArrayView::noValues(shape));
  //
  // create a view of the data with the correct shape
  DataArrayView tempView(m_data,shape);
  DataArrayView::RegionLoopRangeType region_loop_range=getSliceRegionLoopRange(region);
  //
  // load the view with the data from the slice
  tempView.copySlice(other.getPointDataView(),region_loop_range);
  setPointDataView(tempView);
}

DataConstant::DataConstant(const FunctionSpace& what,
                           const DataArrayView::ShapeType &shape,
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
}

string
DataConstant::toString() const
{
  return getPointDataView().toString("");
}

DataArrayView::ValueType::size_type
DataConstant::getPointOffset(int sampleNo,
                             int dataPointNo) const
{
  EsysAssert((validSamplePointNo(dataPointNo) && validSampleNo(sampleNo)),
              "Invalid index, sampleNo: " << sampleNo << " dataPointNo: " << dataPointNo);
  //
  // Whatever the coord's always return the same value as this is constant data.
  return 0;
}

DataArrayView::ValueType::size_type
DataConstant::getLength() const
{
  return m_data.size();
}

DataArrayView
DataConstant::getDataPoint(int sampleNo,
                           int dataPointNo)
{
  EsysAssert((validSamplePointNo(dataPointNo) && validSampleNo(sampleNo)),
             "Invalid index, sampleNo: " << sampleNo << " dataPointNo: " << dataPointNo);
  //
  // Whatever the coord's always return the same value as this is constant data.
  return getPointDataView();
}
  
DataAbstract*
DataConstant::getSlice(const DataArrayView::RegionType& region) const
{
  return new DataConstant(*this,region);
}

void
DataConstant::setSlice(const DataAbstract* value,
                       const DataArrayView::RegionType& region) 
{
  const DataConstant* tempDataConst=dynamic_cast<const DataConstant*>(value);
  if (tempDataConst==0) {
    throw DataException("Programming error - casting to DataConstant.");
  }
  // 
  DataArrayView::ShapeType shape(DataArrayView::getResultSliceShape(region));
  DataArrayView::RegionLoopRangeType region_loop_range=getSliceRegionLoopRange(region);
  //
  // check shape:
  if (getPointDataView().getRank()!=region.size()) {
    throw DataException("Error - Invalid slice region.");
  }
  if (tempDataConst->getPointDataView().getRank()>0 and !value->getPointDataView().checkShape(shape)) {
    throw DataException (value->getPointDataView().createShapeErrorMessage(
                "Error - Couldn't copy slice due to shape mismatch.",shape));
  }
  //
  getPointDataView().copySliceFrom(tempDataConst->getPointDataView(),region_loop_range);
}

void
DataConstant::reshapeDataPoint(const DataArrayView::ShapeType& shape) 
{
  if (getPointDataView().getRank()!=0) {
    stringstream temp;
    temp << "Error - Can only reshape Data with data points of rank 0. "
         << "This Data has data points with rank: " << getPointDataView().getRank();
    throw DataException(temp.str());
  }
  m_data.resize(DataArrayView::noValues(shape),getPointDataView()());
  DataArrayView newView(m_data,shape);
  setPointDataView(newView);
}

int
DataConstant::archiveData(ofstream& archiveFile,
                          const DataArrayView::ValueType::size_type noValues) const
{
  return(m_data.archiveData(archiveFile, noValues));
}

int
DataConstant::extractData(ifstream& archiveFile,
                          const DataArrayView::ValueType::size_type noValues)
{
  return(m_data.extractData(archiveFile, noValues));
}

}  // end of namespace
