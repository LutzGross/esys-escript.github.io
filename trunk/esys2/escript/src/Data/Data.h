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

#ifndef DATA_H
#define DATA_H

#include "escript/Data/DataAbstract.h"
#include "escript/Data/DataTagged.h"
#include "escript/Data/FunctionSpace.h"
#include "escript/Data/BinaryOp.h"
#include "escript/Data/UnaryOp.h"

extern "C" {
#include "escript/Data/DataC.h"
}

#include <iostream>
#include <string>
#include <memory>
#include <algorithm>

#include <boost/shared_ptr.hpp>
#include <boost/python/object.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/numeric.hpp>

/**
   \brief
   Data is essentially a factory class which creates the appropriate Data
   object for the given construction arguments. It retains control over
   the object created for the lifetime of the object.
   The type of Data object referred to may change during the lifetime of
   the Data object.

   Description:
   Data is essentially a factory class which creates the appropriate Data
   object for the given construction arguments. It retains control over
   the object created for the lifetime of the object.
   The type of Data object referred to may change during the lifetime of
   the Data object.
*/

namespace escript {

//
// Forward declaration for various implimentations of Data.
class DataEmpty;
class DataConstant;
class DataTagged;
class DataExpanded;

class Data {

  public:

  typedef double (*UnaryDFunPtr)(double);
  typedef double (*BinaryDFunPtr)(double,double);

  /**
     \brief
     Default constructor.
     Creates a DataEmpty object.
  */
  Data();

  /**
     \brief
     Copy constructor.
     WARNING: Only performs a shallow copy.
  */
  Data(const Data& inData);

  /**
     \brief
     Constructor from another Data object. If "what" is different from the
     function space of inData the inData are tried to be interpolated to what
     otherwise a shallow copy of inData is returned.
  */
  Data(const Data& inData,
       const FunctionSpace& what);

  /**
     \brief
     Constructor which copies data from a DataArrayView.

     \param value - Input - Data value for a single point.
     \param what - Input - A description of what this data represents.
     \param expanded - Input - Flag, if true fill the entire container with
                       the value. Otherwise a more efficient storage
                       mechanism will be used.
  */
  Data(const DataArrayView& value,
       const FunctionSpace& what=FunctionSpace(),
       bool expanded=false);

  /**
     \brief
     Constructor which creates a Data from a DataArrayView shape.

     \param value - Input - Single value applied to all Data.
     \param dataPointShape - Input - The shape of each data point.
     \param what - Input - A description of what this data represents.
     \param expanded - Input - Flag, if true fill the entire container with
                       the given value. Otherwise a more efficient storage
                       mechanism will be used.
  */
  Data(double value,
       const DataArrayView::ShapeType& dataPointShape=DataArrayView::ShapeType(),
       const FunctionSpace& what=FunctionSpace(),
       bool expanded=false);

  /**
     \brief
     Constructor which performs a deep copy of a region from another Data object.

     \param inData - Input - Input Data object.
     \param region - Input - Region to copy.
  */
  Data(const Data& inData,
       const DataArrayView::RegionType& region);

  /**
     \brief
     Constructor which will create Tagged data if expanded is false.
     No attempt is made to ensure the tag keys match the tag keys
     within the function space.

     \param tagKeys - Input - List of tag values.
     \param values - Input - List of values, one for each tag.
     \param defaultValue - Input - A default value, used if tag doesn't exist.
     \param what - Input - A description of what this data represents.
     \param expanded - Input - Flag, if true fill the entire container with
                       the appropriate values.
  */
  Data(const DataTagged::TagListType& tagKeys,
       const DataTagged::ValueListType& values,
       const DataArrayView& defaultValue,
       const FunctionSpace& what=FunctionSpace(),
       bool expanded=false);

  /**
     \brief
     Constructor which copies data from a python numarray.

     \param value - Input - Data value for a single point.
     \param what - Input - A description of what this data represents.
     \param expanded - Input - Flag, if true fill the entire container with
                       the value. Otherwise a more efficient storage
                       mechanism will be used.
  */
  Data(const boost::python::numeric::array& value,
       const FunctionSpace& what=FunctionSpace(),
       bool expanded=false);

  /**
     \brief
     Constructor which copies data from any object that can be converted into
     a python numarray.

     \param value - Input - Input data.
     \param what - Input - A description of what this data represents.
     \param expanded - Input - Flag, if true fill the entire container with
                       the value. Otherwise a more efficient storage
                       mechanism will be used.
  */
  Data(const boost::python::object& value,
       const FunctionSpace& what=FunctionSpace(),
       bool expanded=false);

  /**
     \brief
     Constructor which creates a DataConstant.
     Copies data from any object that can be converted
     into a numarray. All other parameters are copied from other.

     \param value - Input - Input data.
     \param other - Input - contains all other parameters.
  */
  Data(const boost::python::object& value,
       const Data& other);

  /**
     \brief
     Constructor which creates a DataConstant of "shape" with constant value.
  */
  Data(double value, 
       const boost::python::tuple& shape=boost::python::make_tuple(), 
       const FunctionSpace& what=FunctionSpace(),
       bool expanded=false);

  /**
     \brief
     Perform the specified algorithm on the Data and return result.
  */
  template <class UnaryFunction>
  inline double
  algorithm(UnaryFunction operation) const;

  /**
     \brief
     Perform the given unary operation on all of the data's elements.
  */
  template <class UnaryFunction>
  void
  unaryOp(UnaryFunction operation);

  /**
     \brief
     Perform the given binary operation on all of the data's elements.
     The underlying type of the right hand side (right) determines the final
     type of *this after the operation. For example if the right hand side
     is expanded *this will be expanded if necessary.
  */
  template <class BinaryFunction>
  void
  binaryOp(const Data& right,
           BinaryFunction operation);

  /**
     \brief
     Perform the given binary operation on all of the data's elements.
  */
  template <class BinaryFunction>
  void
  binaryOp(const boost::python::object& right,
           BinaryFunction operation);

  /**
     \brief
     Overloaded operator +=
     \param right - Input - The right hand side.
  */
  Data& operator+=(const Data& right);
  Data& operator+=(const boost::python::object& right);

  /**
     \brief
     Overloaded operator -=
     \param right - Input - The right hand side.
  */
  Data& operator-=(const Data& right);
  Data& operator-=(const boost::python::object& right);

 /**
     \brief
     Overloaded operator *=
     \param right - Input - The right hand side.
  */
  Data& operator*=(const Data& right);
  Data& operator*=(const boost::python::object& right);

 /**
     \brief
     Overloaded operator /=
     \param right - Input - The right hand side.
  */
  Data& operator/=(const Data& right);
  Data& operator/=(const boost::python::object& right);

  /**
     \brief
     Return the power of Data.
  */
  Data powD(const Data& right) const;
  Data powO(const boost::python::object& right) const;

  /**
     \brief
     Return the C wrapper for the Data object.
  */
  escriptDataC getDataC();

  /**
     \brief
     Return the C wrapper for the Data object - const version.
  */
  escriptDataC getDataC() const;

  /**
     \brief
     Write the data as a string.
  */
  std::string toString() const;

  /**
     \brief
     Return the DataArrayView of the point data. This essentially contains
     the shape information for each data point although it also may be used
     to manipulate the point data.
  */
  inline
  const DataArrayView&
  getPointDataView() const
  {
     return m_data->getPointDataView();
  }

  /**
     \brief
     Perform a deep copy.
  */
  void
  copy(const Data& other);

  /**
    \brief
    Copy other Data object into this Data object where mask is positive.
  */
  void
  copyWithMask(const Data& other,
               const Data& mask);

  /**
     \brief
     Whatever the current Data type make it expanded.
  */
  void
  expand();

  /**
     \brief
     If possible convert the Data type to tagged. This will only allow
     Constant data to be converted to tagged. An attempt to convert
     Expanded data to tagged will throw an exception.
  */
  void
  tag();

  /**
     \brief
     Return true if this Data is expanded.
  */
  bool
  isExpanded() const;

  /**
     \brief
     Return true if this Data is tagged.
  */
  bool
  isTagged() const;

  /**
     \brief
     Return true if this Data is empty.
  */
  bool
  isEmpty() const;

  /**
     \brief
     Return true if this Data is constant.
  */
  bool
  isConstant() const;

  /**
     \brief
     Return the function space.
  */
  inline
  const FunctionSpace&
  getFunctionSpace() const
  {
    return m_data->getFunctionSpace();
  }

  /**
     \brief
     Return a copy of the function space.
  */
  const FunctionSpace
  getCopyOfFunctionSpace() const;

  /**
     \brief
     Return the domain.
  */
  inline
  const AbstractDomain&
  getDomain() const
  {
     return getFunctionSpace().getDomain();
  }

  /**
     \brief
     Return a copy of the domain.
  */
  const AbstractDomain
  getCopyOfDomain() const;

  /**
     \brief
     Return the rank of the point data.
  */
  inline
  int
  getDataPointRank() const
  {
    return m_data->getPointDataView().getRank();
  }

  /**
     \brief
     Return the number of data points per sample.
  */
  inline
  int
  getNumDataPointsPerSample() const
  {
    return m_data->getNumDPPSample();
  }

  /**
     \brief
     Return the number of samples.
  */
  int
  getNumSamples() const;

  /**
     \brief
     Check *this and the right operand are compatible. Throws
     an exception if they aren't.
     \param right - Input - The right hand side.
  */
  void
  operandCheck(const Data& right) const;

  /**
     \brief
     Return the sample data for the given sample no. This is not the
     preferred interface but is provided for use by C code.
     \param sampleNo - Input - the given sample no.
  */
  DataAbstract::ValueType::value_type*
  getSampleData(DataAbstract::ShapeType::size_type sampleNo);

  /**
     \brief
     Return the sample data for the given tag. If an attempt is made to
     access data that isn't tagged an exception will be thrown.
     \param tag - Input - the tag key.
  */
  DataAbstract::ValueType::value_type*
  getSampleDataByTag(int tag);

  /**
     \brief
     Return a view into the data for the data point specified.
     NOTE: Construction of the DataArrayView is a relatively expensive
     operation.
     \param sampleNo - Input -
     \param dataPointNo - Input -
  */
  inline
  DataArrayView
  getDataPoint(int sampleNo,
               int dataPointNo)
  {
    return m_data->getDataPoint(sampleNo,dataPointNo);
  }

  /**
     \brief
     Return a reference to the data point shape.
  */
  const DataArrayView::ShapeType&
  getDataPointShape() const;

  /**
     \brief
     Return data point shape as a tuple of integers:
  */
  boost::python::tuple
  getShapeTuple() const;

  /**
     \brief
     Return the size of the data point. It is the product of the
     data point shape dimensions.
  */
  int
  getDataPointSize() const;

  /**
     \brief
     Return the number of doubles stored for Data.
  */
  DataArrayView::ValueType::size_type
  getLength() const;

  /**
     \brief
     Interpolates this onto the given functionspace and returns the result as a Data object.
  */
  Data
  interpolate(const FunctionSpace& functionspace) const;

  /**
     \brief
     Calculates the gradient of the data at the data points of functionspace.
     If functionspace is not present the function space of Function(getDomain()) is used.
  */
  Data
  gradOn(const FunctionSpace& functionspace) const;

  Data
  grad() const;

  /**
     \brief
     Calculate the integral over the function space domain.
  */
  boost::python::numeric::array
  integrate() const;

  /**
     \brief
     Return a Data with a 1 for +ive values and a 0 for 0 or -ive values.
  */
  Data
  wherePositive() const;

  /**
     \brief
     Return a Data with a 1 for +ive values and a 0 for -ive values.
  */
  Data
  whereNonNegative() const;

  /**
     \brief
     Return a Data with a 1 for -ive values and a 0 for +ive or 0 values.
  */
  Data
  whereNegative() const;

  /**
     \brief
     Return a Data with a 1 for 0 values a 0 for +ive or -ive.
  */
  Data
  whereZero() const;

  /**
     \brief
     Return the sin of Data.
  */
  Data
  sin() const;

  /**
     \brief
     Return the cos of Data.
  */
  Data
  cos() const;

  /**
     \brief
     Return the tan of Data.
  */
  Data
  tan() const;

  /**
     \brief
     Return the log to base 10 of Data.
  */
  Data
  log() const;

  /**
     \brief
     Return the natural log of Data.
  */
  Data
  ln() const;

  /**
     \brief
     Return a Data containing a slice of this Data.
  */
  Data
  getSlice(const DataArrayView::RegionType& region) const;

  /**
     \brief
     Copy the specified region from the given value.
     \param value - Input - Data to copy from.
     \param region - Input - Region to copy.
  */
  void
  setSlice(const Data& value,
           const DataArrayView::RegionType& region);

  /**
     \brief
     Return the maximum absolute value.
  */
  double
  Lsup() const;

  /**
     \brief
     Return the maximum value.
  */
  double
  sup() const;

  /**
     \brief
     Return the minimum value.
  */
  double
  inf() const;

  /**
     \brief
     Returns a slice from this.
  */
  Data
  getItem(const boost::python::object& key) const;

  /**
     \brief
     Copies slice from value into this.
  */
  void
  setItem(const boost::python::object& key,
          const Data& value);

  /**
     \brief
     Convert the underlying data type to match the RHS.
     \param right - Input - data type to match.
  */
  void
  typeMatch(const Data& right);

  /**
     \brief
     Returns true if this can be interpolated to functionspace.
  */
  bool
  probeInterpolation(const FunctionSpace& functionspace) const;

  /**
     \brief
     Assign the given value to the tag. Implicitly converts this
     object to type DataTagged. Throws an exception if this object
     cannot be converted to a DataTagged object.
     \param tagKey - Input - Integer key.
     \param value - Input - Value to associate with given key.
  */
  void
  setTaggedValue(int tagKey,
                 const boost::python::object& value);

  /**
     \brief
     Assign the given value to the tag. Implicitly converts this
     object to type DataTagged. Throws an exception if this object
     cannot be converted to a DataTagged object.
     \param tagKey - Input - Integer key.
     \param value - Input - Value to associate with given key.
     Note: removed for now - this version not needed, and breaks escript.cpp
  */
  /*
  void
  setTaggedValue(int tagKey,
                 const DataArrayView& value);
  */

 private:

  /**
     \brief
     Construct a Data object of the appropriate type.
  */
  template <class IValueType>
  void
  initialise(const IValueType& value,
             const FunctionSpace& what,
             bool expanded);

  /**
     \brief
     Reshape the data point if the data point is currently rank 0.
     Will throw an exception if the data points are not rank 0.
     The original data point value is used for all values of the new
     data point.
  */
  void
  reshapeDataPoint(const DataArrayView::ShapeType& shape);

  //
  // pointer to the actual data
  boost::shared_ptr<DataAbstract> m_data;

};

template <class IValueType>
void
Data::initialise(const IValueType& value,
                 const FunctionSpace& what,
                 bool expanded)
{
  //
  // Construct a Data object of the appropriate type.
  // Construct the object first as there seems to be a bug which causes
  // undefined behaviour if an exception is thrown during construction
  // within the shared_ptr constructor.
  if (expanded) {
    DataAbstract* temp=new DataExpanded(value,what);
    m_data=boost::shared_ptr<DataAbstract>(temp);
  } else {
    DataAbstract* temp=new DataConstant(value,what);
    m_data=boost::shared_ptr<DataAbstract>(temp);
  }
}

inline
DataAbstract::ValueType::value_type*
Data::getSampleData(DataAbstract::ValueType::size_type sampleNo)
{
  return m_data->getSampleData(sampleNo);
}

inline
DataAbstract::ValueType::value_type*
Data::getSampleDataByTag(int tag)
{
  return m_data->getSampleDataByTag(tag);
}

inline
void
Data::operandCheck(const Data& right) const
{
  return m_data->operandCheck(*(right.m_data.get()));
}

inline
int
Data::getNumSamples() const
{
  return m_data->getNumSamples();
}

inline
std::string
Data::toString() const
{
  return m_data->toString();
}

/**
  \brief
  Operator+
  Takes two Data objects.
*/
Data operator+(const Data& left, const Data& right);

/**
  \brief
  Operator-
  Takes two Data objects.
*/
Data operator-(const Data& left, const Data& right);

/**
  \brief
  Operator*
  Takes two Data objects.
*/
Data operator*(const Data& left, const Data& right);

/**
  \brief
  Operator/
  Takes two Data objects.
*/
Data operator/(const Data& left, const Data& right);

/**
  \brief
  Operator+
  Takes LHS Data object and RHS python::object.
  python::object must be convertable to Data type.
*/
Data operator+(const Data& left, const boost::python::object& right);

/**
  \brief
  Operator-
  Takes LHS Data object and RHS python::object.
  python::object must be convertable to Data type.
*/
Data operator-(const Data& left, const boost::python::object& right);

/**
  \brief
  Operator*
  Takes LHS Data object and RHS python::object.
  python::object must be convertable to Data type.
*/
Data operator*(const Data& left, const boost::python::object& right);

/**
  \brief
  Operator/
  Takes LHS Data object and RHS python::object.
  python::object must be convertable to Data type.
*/
Data operator/(const Data& left, const boost::python::object& right);

/**
  \brief
  Operator+
  Takes LHS python::object and RHS Data object.
  python::object must be convertable to Data type.
*/
Data operator+(const boost::python::object& left, const Data& right);

/**
  \brief
  Operator-
  Takes LHS python::object and RHS Data object.
  python::object must be convertable to Data type.
*/
Data operator-(const boost::python::object& left, const Data& right);

/**
  \brief
  Operator*
  Takes LHS python::object and RHS Data object.
  python::object must be convertable to Data type.
*/
Data operator*(const boost::python::object& left, const Data& right);

/**
  \brief
  Operator/
  Takes LHS python::object and RHS Data object.
  python::object must be convertable to Data type.
*/
Data operator/(const boost::python::object& left, const Data& right);

/**
  \brief
  Output operator
*/
std::ostream& operator<<(std::ostream& o, const Data& data);

/**
  \brief
  Return true if operands are equivalent, else return false.
  NB: this operator does very little at this point, and isn't to 
  be relied on. Requires further implementation.
*/
bool operator==(const Data& left, const Data& right);

template <class BinaryFunction>
inline
void
Data::binaryOp(const Data& right,
               BinaryFunction operation)
{
   //
   // if this has a rank of zero promote it to the rank of the RHS
   if (getPointDataView().getRank()==0 && right.getPointDataView().getRank()!=0) {
     reshapeDataPoint(right.getPointDataView().getShape());
   }
   //
   // initially make the temporary a shallow copy
   Data tempRight(right);
   if (getFunctionSpace()!=right.getFunctionSpace()) {
     if (right.probeInterpolation(getFunctionSpace())) {
       //
       // an interpolation is required so create a new Data 
       tempRight=Data(right,this->getFunctionSpace());
     } else if (probeInterpolation(right.getFunctionSpace())) {
       //
       // interpolate onto the RHS function space
       Data tempLeft(*this,right.getFunctionSpace());
       m_data=tempLeft.m_data;
     }
   }
   operandCheck(tempRight);
   //
   // ensure this has the right type for the RHS
   typeMatch(tempRight);
   //
   // Need to cast to the concrete types so that the correct binaryOp
   // is called.
   if (isExpanded()) {
     //
     // Expanded data will be done in parallel, the right hand side can be
     // of any data type
     DataExpanded* leftC=dynamic_cast<DataExpanded*>(m_data.get());
     EsysAssert((leftC!=0), "Programming error - casting to DataExpanded.");
     escript::binaryOp(*leftC,*(tempRight.m_data.get()),operation);
   } else if (isTagged()) {
     //
     // Tagged data is operated on serially, the right hand side can be
     // either DataConstant or DataTagged
     DataTagged* leftC=dynamic_cast<DataTagged*>(m_data.get());
     EsysAssert((leftC!=0), "Programming error - casting to DataTagged.");
     if (right.isTagged()) {
       DataTagged* rightC=dynamic_cast<DataTagged*>(tempRight.m_data.get());
       EsysAssert((rightC!=0), "Programming error - casting to DataTagged.");
       escript::binaryOp(*leftC,*rightC,operation);
     } else {
       DataConstant* rightC=dynamic_cast<DataConstant*>(tempRight.m_data.get());
       EsysAssert((rightC!=0), "Programming error - casting to DataConstant.");
       escript::binaryOp(*leftC,*rightC,operation);
     }
   } else {
     //
     // can only be DataConstant
     DataConstant* leftC=dynamic_cast<DataConstant*>(m_data.get());
     DataConstant* rightC=dynamic_cast<DataConstant*>(tempRight.m_data.get());
     EsysAssert((leftC!=0 && rightC!=0), "Programming error - casting to DataConstant.");
     escript::binaryOp(*leftC,*rightC,operation);
   }
}

template <class BinaryFunction>
inline
void
Data::binaryOp(const boost::python::object& right,
               BinaryFunction operation)
{
   DataArray temp(right);
   //
   // if this has a rank of zero promote it to the rank of the RHS.
   if (getPointDataView().getRank()==0 && temp.getView().getRank()!=0) {
      reshapeDataPoint(temp.getView().getShape());
   }
   //
   // Always allow scalar values for the RHS but check other shapes
   if (temp.getView().getRank()!=0) {
     if (!getPointDataView().checkShape(temp.getView().getShape())) {
       throw DataException(getPointDataView().createShapeErrorMessage(
                  "Error - RHS shape doesn't match LHS shape.",temp.getView().getShape()));
     }
   }

   if (isExpanded()) {
     //
     // Expanded data will be done in parallel
     DataExpanded* leftC=dynamic_cast<DataExpanded*>(m_data.get());
     EsysAssert((leftC!=0),"Programming error - casting to DataExpanded.");
     escript::binaryOp(*leftC,temp.getView(),operation);
   } else if (isTagged()) {
     DataTagged* leftC=dynamic_cast<DataTagged*>(m_data.get());
     EsysAssert((leftC!=0), "Programming error - casting to DataTagged.");
     escript::binaryOp(*leftC,temp.getView(),operation);
   } else {
     //
     // can only be DataConstant
     DataConstant* leftC=dynamic_cast<DataConstant*>(m_data.get());
     EsysAssert((leftC!=0),"Programming error - casting to DataConstant.");
     escript::binaryOp(*leftC,temp.getView(),operation);
   }
}

template <class UnaryFunction>
inline
void
Data::unaryOp(UnaryFunction operation)
{
  if (isExpanded()) {
    //
    // Expanded data will be done in parallel
    DataExpanded* leftC=dynamic_cast<DataExpanded*>(m_data.get());
    EsysAssert((leftC!=0), "Programming error - casting to DataExpanded.");
    escript::unaryOp(*leftC,operation);
  } else if (isTagged()) {
    DataTagged* leftC=dynamic_cast<DataTagged*>(m_data.get());
    EsysAssert((leftC!=0), "Programming error - casting to DataTagged.");
    escript::unaryOp(*leftC,operation);
  } else {
    DataConstant* leftC=dynamic_cast<DataConstant*>(m_data.get());
    EsysAssert((leftC!=0), "Programming error - casting to DataConstant.");
    escript::unaryOp(*leftC,operation);
  }
}

template <class UnaryFunction>
inline
double
Data::algorithm(UnaryFunction operation) const
{
  if (isExpanded()) {
    //
    // Expanded data will be done in parallel
    DataExpanded* leftC=dynamic_cast<DataExpanded*>(m_data.get());
    EsysAssert((leftC!=0), "Programming error - casting to DataExpanded.");
    return escript::algorithm(*leftC,operation);
  }
  else if (isTagged()) {
    DataTagged* leftC=dynamic_cast<DataTagged*>(m_data.get());
    EsysAssert((leftC!=0), "Programming error - casting to DataTagged.");
    return escript::algorithm(*leftC,operation);
  } else {
    DataConstant* leftC=dynamic_cast<DataConstant*>(m_data.get());
    EsysAssert((leftC!=0), "Programming error - casting to DataConstant.");
    return escript::algorithm(*leftC,operation);
  }
}

template <class UnaryFunction>
inline
Data
unaryOp(const Data& other,
        UnaryFunction operation)
{
  //
  // Perform the given operation on a copy of the input data and return the result
  Data result;
  //
  // perform a deep copy
  result.copy(other);
  result.unaryOp(operation);
  return result;
}

}

#endif
