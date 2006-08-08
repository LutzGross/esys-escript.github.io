//$Id$
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

#include "DataConstant.h"
#include "DataException.h"
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
  int len = DataArrayView::noValues(shape);
  m_data.resize(len,0.,len);
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
  if (tempDataConst->getPointDataView().getRank()>0 && !value->getPointDataView().checkShape(shape)) {
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
  int len = DataArrayView::noValues(shape);
  m_data.resize(len,getPointDataView()(),len);
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

void
DataConstant::symmetric(DataAbstract* ev)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::symmetric: casting to DataConstant failed (propably a programming error).");
  }
  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();
  DataArrayView::symmetric(thisView,0,evView,0);
}

void
DataConstant::nonsymmetric(DataAbstract* ev)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::nonsymmetric: casting to DataConstant failed (propably a programming error).");
  }
  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();
  DataArrayView::nonsymmetric(thisView,0,evView,0);
}

void
DataConstant::trace(DataAbstract* ev, int axis_offset)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::trace: casting to DataConstant failed (propably a programming error).");
  }
  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();
  DataArrayView::trace(thisView,0,evView,0,axis_offset);
}

void
DataConstant::swap(DataAbstract* ev, int axis_offset)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::swap: casting to DataConstant failed (propably a programming error).");
  }
  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();
  DataArrayView::swap(thisView,0,evView,0,axis_offset);
}

void
DataConstant::transpose(DataAbstract* ev, int axis_offset)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::transpose: casting to DataConstant failed (propably a programming error).");
  }
  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();
  DataArrayView::transpose(thisView,0,evView,0,axis_offset);
}

void
DataConstant::eigenvalues(DataAbstract* ev)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::eigenvalues: casting to DataConstant failed (propably a programming error).");
  }
  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();
  DataArrayView::eigenvalues(thisView,0,evView,0);
}
void
DataConstant::eigenvalues_and_eigenvectors(DataAbstract* ev,DataAbstract* V,const double tol)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::eigenvalues_and_eigenvectors: casting to DataConstant failed (propably a programming error).");
  }
  DataConstant* temp_V=dynamic_cast<DataConstant*>(V);
  if (temp_V==0) {
    throw DataException("Error - DataConstant::eigenvalues_and_eigenvectors: casting to DataConstant failed (propably a programming error).");
  }
  DataArrayView thisView=getPointDataView();
  DataArrayView evView=ev->getPointDataView();
  DataArrayView VView=V->getPointDataView();

  DataArrayView::eigenvalues_and_eigenvectors(thisView,0,evView,0,VView,tol);
}


}  // end of namespace
