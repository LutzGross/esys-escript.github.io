// $Id$
/*=============================================================================

 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT ACcESS 2004 -  All Rights Reserved                         *
 *                                                                            *
 * This software is the property of ACcESS.  No part of this code             *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that                          *
 * person has a software license agreement with ACcESS.                       *
 *                                                                            *
 ******************************************************************************

******************************************************************************/

#include "escript/Data/Data.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <functional>
#include <math.h>

#include <boost/python/str.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/long.hpp>

#include "escript/Data/DataException.h"

#include "escript/Data/DataExpanded.h"
#include "escript/Data/DataConstant.h"
#include "escript/Data/DataTagged.h"
#include "escript/Data/DataEmpty.h"
#include "escript/Data/DataArray.h"
#include "escript/Data/DataAlgorithm.h"
#include "escript/Data/FunctionSpaceFactory.h"
#include "escript/Data/AbstractContinuousDomain.h"
#include "escript/Data/UnaryFuncs.h"

using namespace std;
using namespace boost::python;
using namespace boost;
using namespace escript;

Data::Data()
{
  //
  // Default data is type DataEmpty
  DataAbstract* temp=new DataEmpty();
  shared_ptr<DataAbstract> temp_data(temp);
  m_data=temp_data;
}

Data::Data(double value,
           const tuple& shape,
           const FunctionSpace& what,
           bool expanded)
{
  DataArrayView::ShapeType dataPointShape;
  for (int i = 0; i < shape.attr("__len__")(); ++i) {
    dataPointShape.push_back(extract<const int>(shape[i]));
  }
  DataArray temp(dataPointShape,value);
  initialise(temp.getView(),what,expanded);
}

Data::Data(double value,
	   const DataArrayView::ShapeType& dataPointShape,
	   const FunctionSpace& what,
           bool expanded)
{
  DataArray temp(dataPointShape,value);
  pair<int,int> dataShape=what.getDataShape();
  initialise(temp.getView(),what,expanded);
}

Data::Data(const Data& inData)
{
  m_data=inData.m_data;
}

Data::Data(const Data& inData,
           const DataArrayView::RegionType& region)
{
  //
  // Create Data which is a slice of another Data
  DataAbstract* tmp = inData.m_data->getSlice(region);
  shared_ptr<DataAbstract> temp_data(tmp);
  m_data=temp_data;
}

Data::Data(const Data& inData,
           const FunctionSpace& functionspace)
{
  if (inData.getFunctionSpace()==functionspace) {
    m_data=inData.m_data;
  } else {
    Data tmp(0,inData.getPointDataView().getShape(),functionspace,true);
    // Note for Lutz, Must use a reference or pointer to a derived object
    // in order to get polymorphic behaviour. Shouldn't really
    // be able to create an instance of AbstractDomain but that was done
    // as a boost python work around which may no longer be required.
    const AbstractDomain& inDataDomain=inData.getDomain();
    if  (inDataDomain==functionspace.getDomain()) {
      inDataDomain.interpolateOnDomain(tmp,inData);
    } else {
      inDataDomain.interpolateACross(tmp,inData);
    }
    m_data=tmp.m_data;
  }
}

Data::Data(const DataTagged::TagListType& tagKeys,
           const DataTagged::ValueListType & values,
           const DataArrayView& defaultValue,
           const FunctionSpace& what,
           bool expanded)
{
  DataAbstract* temp=new DataTagged(tagKeys,values,defaultValue,what);
  shared_ptr<DataAbstract> temp_data(temp);
  m_data=temp_data;
  if (expanded) {
    expand();
  }
}

Data::Data(const numeric::array& value,
	   const FunctionSpace& what,
           bool expanded)
{
  initialise(value,what,expanded);
}

Data::Data(const DataArrayView& value,
	   const FunctionSpace& what,
           bool expanded)
{
  initialise(value,what,expanded);
}

Data::Data(const object& value,
	   const FunctionSpace& what,
           bool expanded)
{
  numeric::array asNumArray(value);
  initialise(asNumArray,what,expanded);
}

Data::Data(const object& value,
           const Data& other)
{
  //
  // Create DataConstant using the given value and all other parameters
  // copied from other. If value is a rank 0 object this Data
  // will assume the point data shape of other.
  DataArray temp(value);
  if (temp.getView().getRank()==0) {
    //
    // Create a DataArray with the scalar value for all elements
    DataArray temp2(other.getPointDataView().getShape(),temp.getView()());
    initialise(temp2.getView(),other.getFunctionSpace(),false);
  } else {
    //
    // Create a DataConstant with the same sample shape as other
    initialise(temp.getView(),other.getFunctionSpace(),false);
  }
}

escriptDataC
Data::getDataC()
{
  escriptDataC temp;
  temp.m_dataPtr=(void*)this;
  return temp;
}

escriptDataC
Data::getDataC() const
{
  escriptDataC temp;
  temp.m_dataPtr=(void*)this;
  return temp;
}

tuple
Data::getShapeTuple() const
{
  const DataArrayView::ShapeType& shape=getDataPointShape();
  switch(getDataPointRank()) {
     case 0:
        return make_tuple();
     case 1:
        return make_tuple(long_(shape[0]));
     case 2:
        return make_tuple(long_(shape[0]),long_(shape[1]));
     case 3:
        return make_tuple(long_(shape[0]),long_(shape[1]),long_(shape[2]));
     case 4:
        return make_tuple(long_(shape[0]),long_(shape[1]),long_(shape[2]),long_(shape[3]));
     default:
        throw DataException("Error - illegal Data rank.");
  }
}

void
Data::copy(const Data& other)
{
  //
  // Perform a deep copy
  {
    DataExpanded* temp=dynamic_cast<DataExpanded*>(other.m_data.get());
    if (temp!=0) {
      //
      // Construct a DataExpanded copy
      DataAbstract* newData=new DataExpanded(*temp);
      shared_ptr<DataAbstract> temp_data(newData);
      m_data=temp_data;
      return;
    }
  }
  {
    DataTagged* temp=dynamic_cast<DataTagged*>(other.m_data.get());
    if (temp!=0) {
      //
      // Construct a DataTagged copy
      DataAbstract* newData=new DataTagged(*temp);
      shared_ptr<DataAbstract> temp_data(newData);
      m_data=temp_data;
      return;
    }
  }
  {
    DataConstant* temp=dynamic_cast<DataConstant*>(other.m_data.get());
    if (temp!=0) {
      //
      // Construct a DataConstant copy
      DataAbstract* newData=new DataConstant(*temp);
      shared_ptr<DataAbstract> temp_data(newData);
      m_data=temp_data;
      return;
    }
  }
  {
    DataEmpty* temp=dynamic_cast<DataEmpty*>(other.m_data.get());
    if (temp!=0) {
      //
      // Construct a DataEmpty copy
      DataAbstract* newData=new DataEmpty();
      shared_ptr<DataAbstract> temp_data(newData);
      m_data=temp_data;
      return;
    }
  }
  throw DataException("Error - Copy not implemented for this Data type.");
}

void
Data::copyWithMask(const Data& other,
                   const Data& mask)
{
  Data mask1;
  Data mask2;

  mask1 = mask.wherePositive();
  mask2.copy(mask1);

  mask1 *= other;
  mask2 *= *this;
  mask2 = *this - mask2;

  *this = mask1 + mask2;
}

bool
Data::isExpanded() const
{
  DataExpanded* temp=dynamic_cast<DataExpanded*>(m_data.get());
  return (temp!=0);
}

bool
Data::isTagged() const
{
  DataTagged* temp=dynamic_cast<DataTagged*>(m_data.get());
  return (temp!=0);
}

bool
Data::isEmpty() const
{
  DataEmpty* temp=dynamic_cast<DataEmpty*>(m_data.get());
  return (temp!=0);
}

bool
Data::isConstant() const
{
  DataConstant* temp=dynamic_cast<DataConstant*>(m_data.get());
  return (temp!=0);
}

void
Data::expand()
{
  if (isConstant()) {
    DataConstant* tempDataConst=dynamic_cast<DataConstant*>(m_data.get());
    DataAbstract* temp=new DataExpanded(*tempDataConst);
    shared_ptr<DataAbstract> temp_data(temp);
    m_data=temp_data;
  } else if (isTagged()) {
    DataTagged* tempDataTag=dynamic_cast<DataTagged*>(m_data.get());
    DataAbstract* temp=new DataExpanded(*tempDataTag);
    shared_ptr<DataAbstract> temp_data(temp);
    m_data=temp_data;
  } else if (isExpanded()) {
    //
    // do nothing
  } else if (isEmpty()) {
    throw DataException("Error - Expansion of DataEmpty not possible.");
  } else {
    throw DataException("Error - Expansion not implemented for this Data type.");
  }
}

void
Data::tag()
{
  if (isConstant()) {
    DataConstant* tempDataConst=dynamic_cast<DataConstant*>(m_data.get());
    DataAbstract* temp=new DataTagged(*tempDataConst);
    shared_ptr<DataAbstract> temp_data(temp);
    m_data=temp_data;
  } else if (isTagged()) {
    // do nothing
  } else if (isExpanded()) {
    throw DataException("Error - Creating tag data from DataExpanded not possible.");
  } else if (isEmpty()) {
    throw DataException("Error - Creating tag data from DataEmpty not possible.");
  } else {
    throw DataException("Error - Tagging not implemented for this Data type.");
  }
}

void
Data::reshapeDataPoint(const DataArrayView::ShapeType& shape) 
{
  m_data->reshapeDataPoint(shape);
}

Data
Data::wherePositive() const
{
  return escript::unaryOp(*this,bind2nd(greater<double>(),0.0));
}

Data
Data::whereNegative() const
{
  return escript::unaryOp(*this,bind2nd(less<double>(),0.0));
}

Data
Data::whereNonNegative() const
{
  return escript::unaryOp(*this,bind2nd(greater_equal<double>(),0.0));
}

Data
Data::whereNonPositive() const
{
  return escript::unaryOp(*this,bind2nd(less_equal<double>(),0.0));
}

Data
Data::whereZero() const
{
  return escript::unaryOp(*this,bind2nd(equal_to<double>(),0.0));
}

Data
Data::whereNonZero() const
{
  return escript::unaryOp(*this,bind2nd(not_equal_to<double>(),0.0));
}

Data
Data::interpolate(const FunctionSpace& functionspace) const
{
  return Data(*this,functionspace);
}

bool
Data::probeInterpolation(const FunctionSpace& functionspace) const
{
  if (getFunctionSpace()==functionspace) {
    return true;
  } else {
    const AbstractDomain& domain=getDomain();
    if  (domain==functionspace.getDomain()) {
      return domain.probeInterpolationOnDomain(getFunctionSpace().getTypeCode(),functionspace.getTypeCode());
    } else {
      return domain.probeInterpolationACross(getFunctionSpace().getTypeCode(),functionspace.getDomain(),functionspace.getTypeCode());
    }
  }
}

Data
Data::gradOn(const FunctionSpace& functionspace) const
{
  if (functionspace.getDomain()!=getDomain())
    throw DataException("Error - gradient cannot be calculated on different domains.");
  DataArrayView::ShapeType grad_shape=getPointDataView().getShape();
  grad_shape.push_back(functionspace.getDim());
  Data out(0.0,grad_shape,functionspace,true);
  getDomain().setToGradient(out,*this);
  return out;
}

Data
Data::grad() const
{
  return gradOn(escript::function(getDomain()));
}

int
Data::getDataPointSize() const
{
  return getPointDataView().noValues();
}

DataArrayView::ValueType::size_type
Data::getLength() const
{
  return m_data->getLength();
}

const DataArrayView::ShapeType&
Data::getDataPointShape() const
{
  return getPointDataView().getShape();
}

boost::python::numeric::array
Data::integrate() const
{
  int index;
  int rank = getDataPointRank();
  DataArrayView::ShapeType shape = getDataPointShape();

  //
  // calculate the integral values
  vector<double> integrals(getDataPointSize());
  AbstractContinuousDomain::asAbstractContinuousDomain(getDomain()).setToIntegrals(integrals,*this);

  //
  // create the numeric array to be returned
  // and load the array with the integral values
  boost::python::numeric::array bp_array(1.0);
  if (rank==0) {
    index = 0;
    bp_array[0] = integrals[index];
  }
  if (rank==1) {
    bp_array.resize(shape[0]);
    for (int i=0; i<shape[0]; i++) {
      index = i;
      bp_array[i] = integrals[index];
    }
  }
  if (rank==2) {
    bp_array.resize(shape[0],shape[1]);
    for (int i=0; i<shape[0]; i++) {
      for (int j=0; j<shape[1]; j++) {
        index = i + shape[0] * j;
        bp_array[i,j] = integrals[index];
      }
    }
  }
  if (rank==3) {
    bp_array.resize(shape[0],shape[1],shape[2]);
    for (int i=0; i<shape[0]; i++) {
      for (int j=0; j<shape[1]; j++) {
        for (int k=0; k<shape[2]; k++) {
          index = i + shape[0] * ( j + shape[1] * k );
          bp_array[i,j,k] = integrals[index];
        }
      }
    }
  }
  if (rank==4) {
    bp_array.resize(shape[0],shape[1],shape[2],shape[3]);
    for (int i=0; i<shape[0]; i++) {
      for (int j=0; j<shape[1]; j++) {
        for (int k=0; k<shape[2]; k++) {
          for (int l=0; l<shape[3]; l++) {
            index = i + shape[0] * ( j + shape[1] * ( k + shape[2] * l ) );
            bp_array[i,j,k,l] = integrals[index];
          }
        }
      }
    }
  }

  //
  // return the loaded array
  return bp_array;
}

Data
Data::sin() const
{
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::sin);
}

Data
Data::cos() const
{
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::cos);
}

Data
Data::tan() const
{
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::tan);
}

Data
Data::log() const
{
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::log10);
}

Data
Data::ln() const
{
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::log);
}

double
Data::Lsup() const
{
  //
  // set the initial absolute maximum value to zero
  return algorithm(DataAlgorithmAdapter<AbsMax>(0));
}

double
Data::sup() const
{
  //
  // set the initial maximum value to min possible double
  return algorithm(DataAlgorithmAdapter<FMax>(numeric_limits<double>::min()));
}

double
Data::inf() const
{
  //
  // set the initial minimum value to max possible double
  return algorithm(DataAlgorithmAdapter<FMin>(numeric_limits<double>::max()));
}

Data
Data::maxval() const
{
  // not implemented - will use dp_algorithm
  return (*this);
}

Data
Data::minval() const
{
  // not implemented - will use dp_algorithm
  return (*this);
}

Data
Data::length() const
{
  // not implemented - will use dp_algorithm
  return (*this);
}

Data
Data::trace() const
{
  // not implemented - will use dp_algorithm
  return (*this);
}

Data
Data::transpose(int axis) const
{
  // not implemented
  return (*this);
}

Data
Data::sign() const
{
  return escript::unaryOp(*this,escript::fsign);
}

Data
Data::abs() const
{
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::fabs);
}

Data
Data::neg() const
{
  return escript::unaryOp(*this,negate<double>());
}

Data
Data::pos() const
{
  return (*this);
}

Data
Data::exp() const
{
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::exp);
}

Data
Data::sqrt() const
{
  return escript::unaryOp(*this,(Data::UnaryDFunPtr)::sqrt);
}

Data&
Data::operator+=(const Data& right)
{
  binaryOp(right,plus<double>());
  return (*this);
}

Data&
Data::operator+=(const boost::python::object& right)
{
  binaryOp(right,plus<double>());
  return (*this);
}

Data&
Data::operator-=(const Data& right)
{
  binaryOp(right,minus<double>());
  return (*this);
}

Data&
Data::operator-=(const boost::python::object& right)
{
  binaryOp(right,minus<double>());
  return (*this);
}

Data&
Data::operator*=(const Data& right)
{
  binaryOp(right,multiplies<double>());
  return (*this);
}

Data&
Data::operator*=(const boost::python::object& right)
{
  binaryOp(right,multiplies<double>());
  return (*this);
}

Data&
Data::operator/=(const Data& right)
{
  binaryOp(right,divides<double>());
  return (*this);
}

Data&
Data::operator/=(const boost::python::object& right)
{
  binaryOp(right,divides<double>());
  return (*this);
}

Data
Data::powO(const boost::python::object& right) const
{
  Data result;
  result.copy(*this);
  result.binaryOp(right,(Data::BinaryDFunPtr)::pow);
  return result;
}

Data
Data::powD(const Data& right) const
{
  Data result;
  result.copy(*this);
  result.binaryOp(right,(Data::BinaryDFunPtr)::pow);
  return result;
}

//
// NOTE: It is essential to specify the namepsace this operator belongs to
Data
escript::operator+(const Data& left, const Data& right)
{
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result+=right;
  return result;
}

//
// NOTE: It is essential to specify the namepsace this operator belongs to
Data
escript::operator-(const Data& left, const Data& right)
{
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result-=right;
  return result;
}

//
// NOTE: It is essential to specify the namepsace this operator belongs to
Data
escript::operator*(const Data& left, const Data& right)
{
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result*=right;
  return result;
}

//
// NOTE: It is essential to specify the namepsace this operator belongs to
Data
escript::operator/(const Data& left, const Data& right)
{
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result/=right;
  return result;
}

//
// NOTE: It is essential to specify the namepsace this operator belongs to
Data
escript::operator+(const Data& left, const boost::python::object& right)
{
  //
  // Convert to DataArray format if possible
  DataArray temp(right);
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result+=right;
  return result;
}

//
// NOTE: It is essential to specify the namepsace this operator belongs to
Data
escript::operator-(const Data& left, const boost::python::object& right)
{
  //
  // Convert to DataArray format if possible
  DataArray temp(right);
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result-=right;
  return result;
}

//
// NOTE: It is essential to specify the namepsace this operator belongs to
Data
escript::operator*(const Data& left, const boost::python::object& right)
{
  //
  // Convert to DataArray format if possible
  DataArray temp(right);
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result*=right;
  return result;
}

//
// NOTE: It is essential to specify the namepsace this operator belongs to
Data
escript::operator/(const Data& left, const boost::python::object& right)
{
  //
  // Convert to DataArray format if possible
  DataArray temp(right);
  Data result;
  //
  // perform a deep copy
  result.copy(left);
  result/=right;
  return result;
}

//
// NOTE: It is essential to specify the namepsace this operator belongs to
Data
escript::operator+(const boost::python::object& left, const Data& right)
{
  //
  // Construct the result using the given value and the other parameters
  // from right
  Data result(left,right);
  result+=right;
  return result;
}

//
// NOTE: It is essential to specify the namepsace this operator belongs to
Data
escript::operator-(const boost::python::object& left, const Data& right)
{
  //
  // Construct the result using the given value and the other parameters
  // from right
  Data result(left,right);
  result-=right;
  return result;
}

//
// NOTE: It is essential to specify the namepsace this operator belongs to
Data
escript::operator*(const boost::python::object& left, const Data& right)
{
  //
  // Construct the result using the given value and the other parameters
  // from right
  Data result(left,right);
  result*=right;
  return result;
}

//
// NOTE: It is essential to specify the namepsace this operator belongs to
Data
escript::operator/(const boost::python::object& left, const Data& right)
{
  //
  // Construct the result using the given value and the other parameters
  // from right
  Data result(left,right);
  result/=right;
  return result;
}

//
// NOTE: It is essential to specify the namepsace this operator belongs to
//bool escript::operator==(const Data& left, const Data& right)
//{
//  /*
//  NB: this operator does very little at this point, and isn't to
//  be relied on. Requires further implementation.
//  */
//
//  bool ret;
//
//  if (left.isEmpty()) {
//    if(!right.isEmpty()) {
//      ret = false;
//    } else {
//      ret = true;
//    }
//  }
//
//  if (left.isConstant()) {
//    if(!right.isConstant()) {
//      ret = false;
//    } else {
//      ret = true;
//    }
// }
//
//  if (left.isTagged()) {
//   if(!right.isTagged()) {
//      ret = false;
//    } else {
//      ret = true;
//    }
//  }
//
//  if (left.isExpanded()) {
//    if(!right.isExpanded()) {
//      ret = false;
//    } else {
//      ret = true;
//    }
//  }
//
//  return ret;
//}

Data
Data::getItem(const boost::python::object& key) const 
{
  const DataArrayView& view=getPointDataView();

  DataArrayView::RegionType slice_region=view.getSliceRegion(key);

  if (slice_region.size()!=view.getRank()) {
    throw DataException("Error - slice size does not match Data rank.");
  }

  return getSlice(slice_region);
}

Data
Data::getSlice(const DataArrayView::RegionType& region) const
{
  return Data(*this,region);
}

void
Data::setItemO(const boost::python::object& key,
               const boost::python::object& value)
{
  Data tempData(value,getFunctionSpace());
  setItemD(key,tempData);
}

void
Data::setItemD(const boost::python::object& key,
               const Data& value)
{
  const DataArrayView& view=getPointDataView();
  DataArrayView::RegionType slice_region=view.getSliceRegion(key);
  if (slice_region.size()!=view.getRank()) {
    throw DataException("Error - slice size does not match Data rank.");
  }
  setSlice(value,slice_region);
}

void
Data::setSlice(const Data& value,
               const DataArrayView::RegionType& region)
{
  Data tempValue(value);
  typeMatchLeft(tempValue);
  typeMatchRight(tempValue);
  m_data->setSlice(tempValue.m_data.get(),region);
}

void
Data::typeMatchLeft(Data& right) const
{
  if (isExpanded()){
    right.expand();
  } else if (isTagged()) {
    if (right.isConstant()) {
      right.tag();
    }
  }
}

void
Data::typeMatchRight(const Data& right)
{
  if (isTagged()) {
    if (right.isExpanded()) {
      expand();
    }
  } else if (isConstant()) {
    if (right.isExpanded()) {
      expand();
    } else if (right.isTagged()) {
      tag();
    }
  }
}

void
Data::setTaggedValue(int tagKey,
                     const boost::python::object& value)
{
  //
  // Ensure underlying data object is of type DataTagged
  tag();

  if (!isTagged()) {
    throw DataException("Error - DataTagged conversion failed!!");
  }

  //
  // Construct DataArray from boost::python::object input value
  DataArray valueDataArray(value);

  //
  // Call DataAbstract::setTaggedValue
  m_data->setTaggedValue(tagKey,valueDataArray.getView());
}

/*
Note: this version removed for now. Not needed, and breaks escript.cpp
void
Data::setTaggedValue(int tagKey,
                     const DataArrayView& value)
{
  //
  // Ensure underlying data object is of type DataTagged
  tag();

  if (!isTagged()) {
    throw DataException("Error - DataTagged conversion failed!!");
  }
                                                                                                               
  //
  // Call DataAbstract::setTaggedValue
  m_data->setTaggedValue(tagKey,value);
}
*/

ostream& escript::operator<<(ostream& o, const Data& data)
{
  o << data.toString();
  return o;
}
