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

  DataConstant::DataConstant(const boost::python::numeric::array& value, const FunctionSpace& what): DataAbstract(what)
  {
    //cout << "Calling DataConstant constructor 1." << endl;
    DataArray temp(value);
    //
    // copy the data in the correct format
    m_data=temp.getData();
    DataArrayView tempView(m_data,temp.getView().getShape());
    //
    // copy the view of the data
    setPointDataView(tempView);
  }

  DataConstant::DataConstant(const DataArrayView& value, const FunctionSpace& what): DataAbstract(what)
  {
    //cout << "Calling DataConstant constructor 2." << endl;
    //
    // copy the data in the correct format
    m_data=value.getData();
    DataArrayView tempView(m_data,value.getShape());
    //
    // copy the view of the data
    setPointDataView(tempView);
  }

  DataConstant::DataConstant(const DataConstant& other): DataAbstract(other.getFunctionSpace()), m_data(other.m_data)
  {
    //cout << "Calling DataConstant copy constructor." << endl;
    // 
    DataArrayView tempView(m_data,other.getPointDataView().getShape());
    //
    // copy the view of the data
    setPointDataView(tempView);
  }

  DataConstant::DataConstant(const DataConstant& other, const DataArrayView::RegionType& region): DataAbstract(other.getFunctionSpace())
  {
    //cout << "Calling DataConstant slice constructor." << endl;
    // 
    //
    // get the shape of the slice
    DataArrayView::ShapeType shape(DataArrayView::getResultSliceShape(region));
    //
    // allocate space for this DataConstant
    m_data.resize(DataArrayView::noValues(shape));
    DataArrayView tempView(m_data,shape);
    tempView.copySlice(other.getPointDataView(),region);
    //
    // copy the view of the data
    setPointDataView(tempView);
  }

  string DataConstant::toString() const
  {
    return getPointDataView().toString("");
  }

  DataArrayView::ValueType::size_type DataConstant::getPointOffset(int sampleNo, int dataPointNo) const
  {
    EsysAssert((validSamplePointNo(dataPointNo) && validSampleNo(sampleNo)),
	       "Invalid index, sampleNo: " << sampleNo 
	       << " dataPointNo: " << dataPointNo);
    return 0;
  }

  DataArrayView::ValueType::size_type DataConstant::getLength() const
  {
    return m_data.size();
  }

  DataArrayView DataConstant::getDataPoint(int sampleNo, int dataPointNo)
  {
    EsysAssert((validSamplePointNo(dataPointNo) && validSampleNo(sampleNo)),
	       "Invalid index, sampleNo: " << sampleNo 
	       << " dataPointNo: " << dataPointNo);
    //
    // Whatever the coord's always return the same value
    return getPointDataView();
  }
  
  DataAbstract* DataConstant::getSlice(const DataArrayView::RegionType& region) const
  {
    return new DataConstant(*this,region);
  }

  void DataConstant::setSlice(const DataAbstract* value, const DataArrayView::RegionType& region) 
  {
    const DataConstant* tempDataConst=dynamic_cast<const DataConstant*>(value);
    if (tempDataConst==0)
    {
      throw DataException("Programming error - casting to DataConstant.");
    }
    getPointDataView().copySliceFrom(tempDataConst->getPointDataView(),region);
  }

  void DataConstant::reshapeDataPoint(const DataArrayView::ShapeType& shape) 
  {
    //
    // reshape a rank zero data point
    if (getPointDataView().getRank()!=0) {
      stringstream temp;
      temp << "Error - Can only reshape Data with data points of rank 0. "
	   << "This Data has data points with rank: " 
	   << getPointDataView().getRank();
      throw DataException(temp.str());
    }
    m_data.resize(DataArrayView::noValues(shape),getPointDataView()());
    DataArrayView newView(m_data,shape);
    setPointDataView(newView);
  }

}  // end of namespace




