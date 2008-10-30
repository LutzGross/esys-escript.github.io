
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/** \file Data.h */

#ifndef DATA_H
#define DATA_H
#include "system_dep.h"

#include "DataAbstract.h"
#include "DataAlgorithm.h"
#include "FunctionSpace.h"
#include "BinaryOp.h"
#include "UnaryOp.h"
#include "DataException.h"
#include "DataTypes.h"

extern "C" {
#include "DataC.h"
/* #include "paso/Paso.h" doesn't belong in this file...causes trouble for BruceFactory.cpp */
}

#include "esysmpi.h"
#include <string>
#include <algorithm>
#include <sstream>


#include <boost/shared_ptr.hpp>
#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/numeric.hpp>

namespace escript {

//
// Forward declaration for various implementations of Data.
class DataConstant;
class DataTagged;
class DataExpanded;

/**
   \brief
   Data represents a collection of datapoints.

   Description:
   Internally, the datapoints are actually stored by a DataAbstract object.
   The specific instance of DataAbstract used may vary over the lifetime
   of the Data object.
   Some methods on this class return references (eg getShape()).
   These references should not be used after an operation which changes the underlying DataAbstract object.
   Doing so will lead to invalid memory access.
   This should not affect any methods exposed via boost::python.
*/
class Data {

  public:

  // These typedefs allow function names to be cast to pointers
  // to functions of the appropriate type when calling unaryOp etc.
  typedef double (*UnaryDFunPtr)(double);
  typedef double (*BinaryDFunPtr)(double,double);


  /**
     Constructors.
  */

  /**
     \brief
     Default constructor.
     Creates a DataEmpty object.
  */
  ESCRIPT_DLL_API
  Data();

  /**
     \brief
     Copy constructor.
     WARNING: Only performs a shallow copy.
  */
  ESCRIPT_DLL_API
  Data(const Data& inData);

  /**
     \brief
     Constructor from another Data object. If "what" is different from the
     function space of inData the inData are tried to be interpolated to what,
     otherwise a shallow copy of inData is returned.
  */
  ESCRIPT_DLL_API
  Data(const Data& inData,
       const FunctionSpace& what);

  /**
	\brief Copy Data from an existing vector
  */ 

  ESCRIPT_DLL_API
  Data(const DataTypes::ValueType& value,
		 const DataTypes::ShapeType& shape,
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
  ESCRIPT_DLL_API
  Data(double value,
       const DataTypes::ShapeType& dataPointShape=DataTypes::ShapeType(),
       const FunctionSpace& what=FunctionSpace(),
       bool expanded=false);

  /**
     \brief
     Constructor which performs a deep copy of a region from another Data object.

     \param inData - Input - Input Data object.
     \param region - Input - Region to copy.
  */
  ESCRIPT_DLL_API
  Data(const Data& inData,
       const DataTypes::RegionType& region);

  /**
     \brief
     Constructor which copies data from a python numarray.

     \param value - Input - Data value for a single point.
     \param what - Input - A description of what this data represents.
     \param expanded - Input - Flag, if true fill the entire container with
                       the value. Otherwise a more efficient storage
                       mechanism will be used.
  */
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
  Data(const boost::python::object& value,
       const Data& other);

  /**
     \brief
     Constructor which creates a DataConstant of "shape" with constant value.
  */
  ESCRIPT_DLL_API
  Data(double value,
       const boost::python::tuple& shape=boost::python::make_tuple(),
       const FunctionSpace& what=FunctionSpace(),
       bool expanded=false);



  /**
	\brief Create a Data using an existing DataAbstract. Warning: The new object assumes ownership of the pointer!
	Once you have passed the pointer, do not delete it.
  */
  ESCRIPT_DLL_API
  explicit Data(DataAbstract* underlyingdata);

  /**
	\brief Create a Data based on the supplied DataAbstract
  */
  ESCRIPT_DLL_API
  explicit Data(DataAbstract_ptr underlyingdata);

  /**
     \brief
     Destructor
  */
  ESCRIPT_DLL_API
  ~Data();

  /**
     \brief Make this object a deep copy of "other".
  */
  ESCRIPT_DLL_API
  void
  copy(const Data& other);

  /**
     \brief Return a pointer to a deep copy of this object.
  */
  ESCRIPT_DLL_API
  Data*
  copySelf();


  /**
     \brief produce a delayed evaluation version of this Data.
  */
  ESCRIPT_DLL_API
  Data
  delay();

  /**
     \brief convert the current data into lazy data.
  */
  ESCRIPT_DLL_API
  void 
  delaySelf();


  /**
     Member access methods.
  */

  /**
     \brief
     switches on update protection

  */
  ESCRIPT_DLL_API
  void
  setProtection();

  /**
     \brief
     Returns trueif the data object is protected against update

  */
  ESCRIPT_DLL_API
  bool
  isProtected() const;

  /**
     \brief
     Return the values of a data point on this process
  */
  ESCRIPT_DLL_API
  const boost::python::numeric::array
  getValueOfDataPoint(int dataPointNo);

  /**
     \brief
     sets the values of a data-point from a python object on this process
  */
  ESCRIPT_DLL_API
  void
  setValueOfDataPointToPyObject(int dataPointNo, const boost::python::object& py_object);

  /**
     \brief
     sets the values of a data-point from a numarray object on this process
  */
  ESCRIPT_DLL_API
  void
  setValueOfDataPointToArray(int dataPointNo, const boost::python::numeric::array&);

  /**
     \brief
     sets the values of a data-point on this process
  */
  ESCRIPT_DLL_API
  void
  setValueOfDataPoint(int dataPointNo, const double);

  /**
     \brief
     Return the value of the specified data-point across all processors
  */
  ESCRIPT_DLL_API
  const boost::python::numeric::array
  getValueOfGlobalDataPoint(int procNo, int dataPointNo);

  /**
     \brief
     Return the tag number associated with the given data-point.

  */
  ESCRIPT_DLL_API
  int
  getTagNumber(int dpno);

  /**
     \brief
     Return the C wrapper for the Data object.
  */
  ESCRIPT_DLL_API
  escriptDataC
  getDataC();






// REMOVE ME
// ESCRIPT_DLL_API
// void
// CompareDebug(const Data& rd);


  /**
     \brief
     Return the C wrapper for the Data object - const version.
  */
  ESCRIPT_DLL_API
  escriptDataC
  getDataC() const;

  /**
     \brief
     Write the data as a string. For large amounts of data, a summary is printed.
  */
  ESCRIPT_DLL_API
  std::string
  toString() const;

  /**
     \brief
     Whatever the current Data type make this into a DataExpanded.
  */
  ESCRIPT_DLL_API
  void
  expand();

  /**
     \brief
     If possible convert this Data to DataTagged. This will only allow
     Constant data to be converted to tagged. An attempt to convert
     Expanded data to tagged will throw an exception.
  */
  ESCRIPT_DLL_API
  void
  tag();

  /**
    \brief If this data is lazy, then convert it to ready data.
    What type of ready data depends on the expression. For example, Constant+Tagged==Tagged.
  */
  ESCRIPT_DLL_API
  void
  resolve();


  /**
     \brief
     Return true if this Data is expanded.
  */
  ESCRIPT_DLL_API
  bool
  isExpanded() const;

  /**
     \brief
     Return true if this Data is tagged.
  */
  ESCRIPT_DLL_API
  bool
  isTagged() const;

  /**
     \brief
     Return true if this Data is constant.
  */
  ESCRIPT_DLL_API
  bool
  isConstant() const;

  /**
     \brief Return true if this Data is lazy.
  */
  ESCRIPT_DLL_API
  bool
  isLazy() const;

  /**
     \brief Return true if this data is ready.
  */
  ESCRIPT_DLL_API
  bool
  isReady() const;

  /**
     \brief
     Return true if this Data holds an instance of DataEmpty. This is _not_ the same as asking if the object 
contains datapoints.
  */
  ESCRIPT_DLL_API
  bool
  isEmpty() const;

  /**
     \brief
     Return the function space.
  */
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
  const FunctionSpace
  getCopyOfFunctionSpace() const;

  /**
     \brief
     Return the domain.
  */
  ESCRIPT_DLL_API
  inline
//   const AbstractDomain&
  const_Domain_ptr
  getDomain() const
  {
     return getFunctionSpace().getDomain();
  }


  /**
     \brief
     Return the domain.
     TODO: For internal use only.   This should be removed.
  */
  ESCRIPT_DLL_API
  inline
//   const AbstractDomain&
  Domain_ptr
  getDomainPython() const
  {
     return getFunctionSpace().getDomainPython();
  }

  /**
     \brief
     Return a copy of the domain.
  */
  ESCRIPT_DLL_API
  const AbstractDomain
  getCopyOfDomain() const;

  /**
     \brief
     Return the rank of the point data.
  */
  ESCRIPT_DLL_API
  inline
  int
  getDataPointRank() const
  {
//    return m_data->getPointDataView().getRank();
    return m_data->getRank();
  }

  /**
     \brief
     Return the number of data points
  */
  ESCRIPT_DLL_API
  inline
  int
  getNumDataPoints() const
  {
    return getNumSamples() * getNumDataPointsPerSample();
  }
  /**
     \brief
     Return the number of samples.
  */
  ESCRIPT_DLL_API
  inline
  int
  getNumSamples() const
  {
    return m_data->getNumSamples();
  }

  /**
     \brief
     Return the number of data points per sample.
  */
  ESCRIPT_DLL_API
  inline
  int
  getNumDataPointsPerSample() const
  {
    return m_data->getNumDPPSample();
  }


  /**
	\brief
	Return the number of values in the shape for this object.
  */
  ESCRIPT_DLL_API
  int
  getNoValues() const
  {
    return m_data->getNoValues();
  }


  /**
     \brief
     dumps the object into a netCDF file
  */
  ESCRIPT_DLL_API
  void
  dump(const std::string fileName) const;

  /**
     \brief
     Return the sample data for the given sample no. This is not the
     preferred interface but is provided for use by C code.
     \param sampleNo - Input - the given sample no.
  */
  ESCRIPT_DLL_API
  inline
  DataAbstract::ValueType::value_type*
  getSampleData(DataAbstract::ValueType::size_type sampleNo);

  /**
     \brief
     Return the sample data for the given tag. If an attempt is made to
     access data that isn't tagged an exception will be thrown.
     \param tag - Input - the tag key.
  */
  ESCRIPT_DLL_API
  inline
  DataAbstract::ValueType::value_type*
  getSampleDataByTag(int tag)
  {
    return m_data->getSampleDataByTag(tag);
  }

  /**
     \brief
     Return a view into the data for the data point specified.
     NOTE: Construction of the DataArrayView is a relatively expensive
     operation.
     \param sampleNo - Input -
     \param dataPointNo - Input -
  */
  ESCRIPT_DLL_API
  DataTypes::ValueType::const_reference
  getDataPoint(int sampleNo, int dataPointNo) const;


  ESCRIPT_DLL_API
  DataTypes::ValueType::reference
  getDataPoint(int sampleNo, int dataPointNo);



  /**
     \brief 
     Return the offset for the given sample and point within the sample
  */
  ESCRIPT_DLL_API
  inline
  DataTypes::ValueType::size_type
  getDataOffset(int sampleNo,
               int dataPointNo)
  {
      return m_data->getPointOffset(sampleNo,dataPointNo);
  }

  /**
     \brief
     Return a reference to the data point shape.
  */
  ESCRIPT_DLL_API
  inline
  const DataTypes::ShapeType&
  getDataPointShape() const
  {
	return m_data->getShape();
  }

  /**
     \brief
     Return the data point shape as a tuple of integers.
  */
  ESCRIPT_DLL_API
  const boost::python::tuple
  getShapeTuple() const;

  /**
     \brief
     Return the size of the data point. It is the product of the
     data point shape dimensions.
  */
  ESCRIPT_DLL_API
  int
  getDataPointSize() const;

  /**
     \brief
     Return the number of doubles stored for this Data.
  */
  ESCRIPT_DLL_API
  DataTypes::ValueType::size_type
  getLength() const;



  /**
     \brief
     Assign the given value to the tag assocciated with name. Implicitly converts this
     object to type DataTagged. Throws an exception if this object
     cannot be converted to a DataTagged object or name cannot be mapped onto a tag key.
     \param tagKey - Input - Integer key.
     \param value - Input - Value to associate with given key.
    ==>*
  */
  ESCRIPT_DLL_API
  void
  setTaggedValueByName(std::string name,
                       const boost::python::object& value);

  /**
     \brief
     Assign the given value to the tag. Implicitly converts this
     object to type DataTagged if it is constant.

     \param tagKey - Input - Integer key.
     \param value - Input - Value to associate with given key.
    ==>*
  */
  ESCRIPT_DLL_API
  void
  setTaggedValue(int tagKey,
                 const boost::python::object& value);

  /**
     \brief
     Assign the given value to the tag. Implicitly converts this
     object to type DataTagged if it is constant.

     \param tagKey - Input - Integer key.
     \param pointshape - Input - The shape of the value parameter
     \param value - Input - Value to associate with given key.
     \param dataOffset - Input - Offset of the begining of the point within the value parameter
  */
  ESCRIPT_DLL_API
  void
  setTaggedValueFromCPP(int tagKey,
			const DataTypes::ShapeType& pointshape,
                        const DataTypes::ValueType& value,
			int dataOffset=0);



  /**
    \brief
    Copy other Data object into this Data object where mask is positive.
  */
  ESCRIPT_DLL_API
  void
  copyWithMask(const Data& other,
               const Data& mask);

  /**
     Data object operation methods and operators.
  */

  /**
     \brief
     set all values to zero
     *
  */
  ESCRIPT_DLL_API
  void
  setToZero();

  /**
     \brief
     Interpolates this onto the given functionspace and returns
     the result as a Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  interpolate(const FunctionSpace& functionspace) const;
  /**
     \brief
     Calculates the gradient of the data at the data points of functionspace.
     If functionspace is not present the function space of Function(getDomain()) is used.
     *
  */
  ESCRIPT_DLL_API
  Data
  gradOn(const FunctionSpace& functionspace) const;

  ESCRIPT_DLL_API
  Data
  grad() const;

  /**
     \brief
     Calculate the integral over the function space domain.
     *
  */
  ESCRIPT_DLL_API
  boost::python::numeric::array
  integrate_const() const;

  ESCRIPT_DLL_API
  boost::python::numeric::array
  integrate();

  /**
     \brief
     Returns 1./ Data object
     *
  */
  ESCRIPT_DLL_API
  Data
  oneOver() const;
  /**
     \brief
     Return a Data with a 1 for +ive values and a 0 for 0 or -ive values.
     *
  */
  ESCRIPT_DLL_API
  Data
  wherePositive() const;

  /**
     \brief
     Return a Data with a 1 for -ive values and a 0 for +ive or 0 values.
     *
  */
  ESCRIPT_DLL_API
  Data
  whereNegative() const;

  /**
     \brief
     Return a Data with a 1 for +ive or 0 values and a 0 for -ive values.
     *
  */
  ESCRIPT_DLL_API
  Data
  whereNonNegative() const;

  /**
     \brief
     Return a Data with a 1 for -ive or 0 values and a 0 for +ive values.
     *
  */
  ESCRIPT_DLL_API
  Data
  whereNonPositive() const;

  /**
     \brief
     Return a Data with a 1 for 0 values and a 0 for +ive or -ive values.
     *
  */
  ESCRIPT_DLL_API
  Data
  whereZero(double tol=0.0) const;

  /**
     \brief
     Return a Data with a 0 for 0 values and a 1 for +ive or -ive values.
     *
  */
  ESCRIPT_DLL_API
  Data
  whereNonZero(double tol=0.0) const;

  /**
     \brief
     Return the maximum absolute value of this Data object.

     The method is not const because lazy data needs to be expanded before Lsup can be computed.
     The _const form can be used when the Data object is const, however this will only work for
     Data which is not Lazy.

     For Data which contain no samples (or tagged Data for which no tags in use have a value)
     zero is returned.
  */
  ESCRIPT_DLL_API
  double
  Lsup();

  ESCRIPT_DLL_API
  double
  Lsup_const() const;


  /**
     \brief
     Return the maximum value of this Data object.

     The method is not const because lazy data needs to be expanded before sup can be computed.
     The _const form can be used when the Data object is const, however this will only work for
     Data which is not Lazy.

     For Data which contain no samples (or tagged Data for which no tags in use have a value)
     a large negative value is returned.
  */
  ESCRIPT_DLL_API
  double
  sup();

  ESCRIPT_DLL_API
  double
  sup_const() const;


  /**
     \brief
     Return the minimum value of this Data object.

     The method is not const because lazy data needs to be expanded before inf can be computed.
     The _const form can be used when the Data object is const, however this will only work for
     Data which is not Lazy.

     For Data which contain no samples (or tagged Data for which no tags in use have a value)
     a large positive value is returned.
  */
  ESCRIPT_DLL_API
  double
  inf();

  ESCRIPT_DLL_API
  double
  inf_const() const;



  /**
     \brief
     Return the absolute value of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  abs() const;

  /**
     \brief
     Return the maximum value of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  maxval() const;

  /**
     \brief
     Return the minimum value of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  minval() const;

  /**
     \brief
     Return the (sample number, data-point number) of the data point with
     the minimum value in this Data object.
  */
  ESCRIPT_DLL_API
  const boost::python::tuple
  minGlobalDataPoint() const;

  ESCRIPT_DLL_API
  void
  calc_minGlobalDataPoint(int& ProcNo,  int& DataPointNo) const;
  /**
     \brief
     Return the sign of each data point of this Data object.
     -1 for negative values, zero for zero values, 1 for positive values.
     *
  */
  ESCRIPT_DLL_API
  Data
  sign() const;

  /**
     \brief
     Return the symmetric part of a matrix which is half the matrix plus its transpose.
     *
  */
  ESCRIPT_DLL_API
  Data
  symmetric() const;

  /**
     \brief
     Return the nonsymmetric part of a matrix which is half the matrix minus its transpose.
     *
  */
  ESCRIPT_DLL_API
  Data
  nonsymmetric() const;

  /**
     \brief
     Return the trace of a matrix
     *
  */
  ESCRIPT_DLL_API
  Data
  trace(int axis_offset) const;

  /**
     \brief
     Transpose each data point of this Data object around the given axis.
     *
  */
  ESCRIPT_DLL_API
  Data
  transpose(int axis_offset) const;

  /**
     \brief
     Return the eigenvalues of the symmetric part at each data point of this Data object in increasing values.
     Currently this function is restricted to rank 2, square shape, and dimension 3.
     *
  */
  ESCRIPT_DLL_API
  Data
  eigenvalues() const;

  /**
     \brief
     Return the eigenvalues and corresponding eigenvcetors of the symmetric part at each data point of this Data object.
     the eigenvalues are ordered in increasing size where eigenvalues with relative difference less than
     tol are treated as equal. The eigenvectors are orthogonal, normalized and the sclaed such that the
     first non-zero entry is positive.
     Currently this function is restricted to rank 2, square shape, and dimension 3
     *
  */
  ESCRIPT_DLL_API
  const boost::python::tuple
  eigenvalues_and_eigenvectors(const double tol=1.e-12) const;

  /**
     \brief
     swaps the components axis0 and axis1
     *
  */
  ESCRIPT_DLL_API
  Data
  swapaxes(const int axis0, const int axis1) const;

  /**
     \brief
     Return the error function erf of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  erf() const;

  /**
     \brief
     Return the sin of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  sin() const;

  /**
     \brief
     Return the cos of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  cos() const;

  /**
     \brief
     Return the tan of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  tan() const;

  /**
     \brief
     Return the asin of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  asin() const;

  /**
     \brief
     Return the acos of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  acos() const;

  /**
     \brief
     Return the atan of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  atan() const;

  /**
     \brief
     Return the sinh of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  sinh() const;

  /**
     \brief
     Return the cosh of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  cosh() const;

  /**
     \brief
     Return the tanh of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  tanh() const;

  /**
     \brief
     Return the asinh of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  asinh() const;

  /**
     \brief
     Return the acosh of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  acosh() const;

  /**
     \brief
     Return the atanh of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  atanh() const;

  /**
     \brief
     Return the log to base 10 of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  log10() const;

  /**
     \brief
     Return the natural log of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  log() const;

  /**
     \brief
     Return the exponential function of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  exp() const;

  /**
     \brief
     Return the square root of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  sqrt() const;

  /**
     \brief
     Return the negation of each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  neg() const;

  /**
     \brief
     Return the identity of each data point of this Data object.
     Simply returns this object unmodified.
     *
  */
  ESCRIPT_DLL_API
  Data
  pos() const;

  /**
     \brief
     Return the given power of each data point of this Data object.

     \param right Input - the power to raise the object to.
     *
  */
  ESCRIPT_DLL_API
  Data
  powD(const Data& right) const;

  /**
     \brief
     Return the given power of each data point of this boost python object.

     \param right Input - the power to raise the object to.
     *
   */
  ESCRIPT_DLL_API
  Data
  powO(const boost::python::object& right) const;

  /**
     \brief
     Return the given power of each data point of this boost python object.

     \param left Input - the bases
     *
   */

  ESCRIPT_DLL_API
  Data
  rpowO(const boost::python::object& left) const;

  /**
     \brief
     writes the object to a file in the DX file format
  */
  ESCRIPT_DLL_API
  void
  saveDX(std::string fileName) const;

  /**
     \brief
     writes the object to a file in the VTK file format
  */
  ESCRIPT_DLL_API
  void
  saveVTK(std::string fileName) const;

  /**
     \brief
     Overloaded operator +=
     \param right - Input - The right hand side.
     *
  */
  ESCRIPT_DLL_API
  Data& operator+=(const Data& right);
  ESCRIPT_DLL_API
  Data& operator+=(const boost::python::object& right);

  ESCRIPT_DLL_API
  Data& operator=(const Data& other);

  /**
     \brief
     Overloaded operator -=
     \param right - Input - The right hand side.
     *
  */
  ESCRIPT_DLL_API
  Data& operator-=(const Data& right);
  ESCRIPT_DLL_API
  Data& operator-=(const boost::python::object& right);

 /**
     \brief
     Overloaded operator *=
     \param right - Input - The right hand side.
     *
  */
  ESCRIPT_DLL_API
  Data& operator*=(const Data& right);
  ESCRIPT_DLL_API
  Data& operator*=(const boost::python::object& right);

 /**
     \brief
     Overloaded operator /=
     \param right - Input - The right hand side.
     *
  */
  ESCRIPT_DLL_API
  Data& operator/=(const Data& right);
  ESCRIPT_DLL_API
  Data& operator/=(const boost::python::object& right);

  /**
     \brief
     Returns true if this can be interpolated to functionspace.
  */
  ESCRIPT_DLL_API
  bool
  probeInterpolation(const FunctionSpace& functionspace) const;

  /**
     Data object slicing methods.
  */

  /**
     \brief
     Returns a slice from this Data object.

     /description
     Implements the [] get operator in python.
     Calls getSlice.

     \param key - Input - python slice tuple specifying
     slice to return.
  */
  ESCRIPT_DLL_API
  Data
  getItem(const boost::python::object& key) const;

  /**
     \brief
     Copies slice from value into this Data object.

     Implements the [] set operator in python.
     Calls setSlice.

     \param key - Input - python slice tuple specifying
     slice to copy from value.
     \param value - Input - Data object to copy from.
  */
  ESCRIPT_DLL_API
  void
  setItemD(const boost::python::object& key,
           const Data& value);

  ESCRIPT_DLL_API
  void
  setItemO(const boost::python::object& key,
           const boost::python::object& value);

  // These following public methods should be treated as private.

  /**
     \brief
     Perform the given unary operation on every element of every data point in
     this Data object.
  */
  template <class UnaryFunction>
  ESCRIPT_DLL_API
  inline
  void
  unaryOp2(UnaryFunction operation);

  /**
     \brief
     Return a Data object containing the specified slice of
     this Data object.
     \param region - Input - Region to copy.
     *
  */
  ESCRIPT_DLL_API
  Data
  getSlice(const DataTypes::RegionType& region) const;

  /**
     \brief
     Copy the specified slice from the given value into this
     Data object.
     \param value - Input - Data to copy from.
     \param region - Input - Region to copy.
     *
  */
  ESCRIPT_DLL_API
  void
  setSlice(const Data& value,
           const DataTypes::RegionType& region);

  /**
     \brief
     print the data values to stdout. Used for debugging
  */
  ESCRIPT_DLL_API
  void
        print(void);

  /**
     \brief
     return the MPI rank number of the local data
                 MPI_COMM_WORLD is assumed and the result of MPI_Comm_size()
                 is returned
  */
  ESCRIPT_DLL_API
        int
        get_MPIRank(void) const;

  /**
     \brief
     return the MPI rank number of the local data
                 MPI_COMM_WORLD is assumed and the result of MPI_Comm_rank()
                 is returned
  */
  ESCRIPT_DLL_API
        int
        get_MPISize(void) const;

  /**
     \brief
     return the MPI rank number of the local data
                 MPI_COMM_WORLD is assumed and returned.
  */
  ESCRIPT_DLL_API
        MPI_Comm
        get_MPIComm(void) const;

  /**
     \brief
     return the object produced by the factory, which is a DataConstant or DataExpanded
	TODO Ownership of this object should be explained in doco.
  */
  ESCRIPT_DLL_API
        DataAbstract*
        borrowData(void) const;

  ESCRIPT_DLL_API
        DataAbstract_ptr
        borrowDataPtr(void) const;

  ESCRIPT_DLL_API
        DataReady_ptr
        borrowReadyPtr(void) const;



  /**
     \brief
     Return a pointer to the beginning of the datapoint at the specified offset.
     TODO Eventually these should be inlined.
     \param i - position(offset) in the underlying datastructure
  */
  ESCRIPT_DLL_API
        DataTypes::ValueType::const_reference
        getDataAtOffset(DataTypes::ValueType::size_type i) const;


  ESCRIPT_DLL_API
        DataTypes::ValueType::reference
        getDataAtOffset(DataTypes::ValueType::size_type i);

 protected:

 private:

  double
  LsupWorker() const;

  double
  supWorker() const;

  double
  infWorker() const;

  boost::python::numeric::array
  integrateWorker() const;

  /**
     \brief
     Check *this and the right operand are compatible. Throws
     an exception if they aren't.
     \param right - Input - The right hand side.
  */
  inline
  void
  operandCheck(const Data& right) const
  {
    return m_data->operandCheck(*(right.m_data.get()));
  }

  /**
     \brief
     Perform the specified reduction algorithm on every element of every data point in
     this Data object according to the given function and return the single value result.
  */
  template <class BinaryFunction>
  inline
  double
  algorithm(BinaryFunction operation,
            double initial_value) const;

  /**
     \brief
     Reduce each data-point in this Data object using the given operation. Return a Data
     object with the same number of data-points, but with each data-point containing only
     one value - the result of the reduction operation on the corresponding data-point in
     this Data object
  */
  template <class BinaryFunction>
  inline
  Data
  dp_algorithm(BinaryFunction operation,
               double initial_value) const;

  /**
     \brief
     Perform the given binary operation on all of the data's elements.
     The underlying type of the right hand side (right) determines the final
     type of *this after the operation. For example if the right hand side
     is expanded *this will be expanded if necessary.
     RHS is a Data object.
  */
  template <class BinaryFunction>
  inline
  void
  binaryOp(const Data& right,
           BinaryFunction operation);

  /**
     \brief
     Convert the data type of the RHS to match this.
     \param right - Input - data type to match.
  */
  void
  typeMatchLeft(Data& right) const;

  /**
     \brief
     Convert the data type of this to match the RHS.
     \param right - Input - data type to match.
  */
  void
  typeMatchRight(const Data& right);

  /**
     \brief
     Construct a Data object of the appropriate type.
  */

  void
  initialise(const DataTypes::ValueType& value,
	     const DataTypes::ShapeType& shape,
             const FunctionSpace& what,
             bool expanded);

  void
  initialise(const boost::python::numeric::array& value,
                 const FunctionSpace& what,
                 bool expanded);

  //
  // flag to protect the data object against any update
  bool m_protected;

  //
  // pointer to the actual data object
//   boost::shared_ptr<DataAbstract> m_data;
  DataAbstract_ptr m_data;

// If possible please use getReadyPtr instead
  const DataReady*
  getReady() const;

  DataReady*
  getReady();

  DataReady_ptr
  getReadyPtr();

  const_DataReady_ptr
  getReadyPtr() const;


};

}   // end namespace escript


// No, this is not supposed to be at the top of the file
// DataAbstact needs to be declared first, then DataReady needs to be fully declared
// so that I can dynamic cast between them below.
#include "DataReady.h"

namespace escript
{

inline
const DataReady*
Data::getReady() const
{
   const DataReady* dr=dynamic_cast<const DataReady*>(m_data.get());
   EsysAssert((dr!=0), "Error - casting to DataReady.");
   return dr;
}

inline
DataReady*
Data::getReady()
{
   DataReady* dr=dynamic_cast<DataReady*>(m_data.get());
   EsysAssert((dr!=0), "Error - casting to DataReady.");
   return dr;
}

inline
DataReady_ptr
Data::getReadyPtr()
{
   DataReady_ptr dr=boost::dynamic_pointer_cast<DataReady>(m_data);
   EsysAssert((dr.get()!=0), "Error - casting to DataReady.");
   return dr;
}


inline
const_DataReady_ptr
Data::getReadyPtr() const
{
   const_DataReady_ptr dr=boost::dynamic_pointer_cast<const DataReady>(m_data);
   EsysAssert((dr.get()!=0), "Error - casting to DataReady.");
   return dr;
}

inline
DataAbstract::ValueType::value_type*
Data::getSampleData(DataAbstract::ValueType::size_type sampleNo)
{
   if (isLazy())
   {
	resolve();
   }
   return getReady()->getSampleData(sampleNo);
}


/**
   Modify a filename for MPI parallel output to multiple files
*/
char *Escript_MPI_appendRankToFileName(const char *, int, int);

/**
   Binary Data object operators.
*/
inline double rpow(double x,double y)
{
    return pow(y,x);
}

/**
  \brief
  Operator+
  Takes two Data objects.
*/
ESCRIPT_DLL_API Data operator+(const Data& left, const Data& right);

/**
  \brief
  Operator-
  Takes two Data objects.
*/
ESCRIPT_DLL_API Data operator-(const Data& left, const Data& right);

/**
  \brief
  Operator*
  Takes two Data objects.
*/
ESCRIPT_DLL_API Data operator*(const Data& left, const Data& right);

/**
  \brief
  Operator/
  Takes two Data objects.
*/
ESCRIPT_DLL_API Data operator/(const Data& left, const Data& right);

/**
  \brief
  Operator+
  Takes LHS Data object and RHS python::object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API Data operator+(const Data& left, const boost::python::object& right);

/**
  \brief
  Operator-
  Takes LHS Data object and RHS python::object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API Data operator-(const Data& left, const boost::python::object& right);

/**
  \brief
  Operator*
  Takes LHS Data object and RHS python::object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API Data operator*(const Data& left, const boost::python::object& right);

/**
  \brief
  Operator/
  Takes LHS Data object and RHS python::object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API Data operator/(const Data& left, const boost::python::object& right);

/**
  \brief
  Operator+
  Takes LHS python::object and RHS Data object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API Data operator+(const boost::python::object& left, const Data& right);

/**
  \brief
  Operator-
  Takes LHS python::object and RHS Data object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API Data operator-(const boost::python::object& left, const Data& right);

/**
  \brief
  Operator*
  Takes LHS python::object and RHS Data object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API Data operator*(const boost::python::object& left, const Data& right);

/**
  \brief
  Operator/
  Takes LHS python::object and RHS Data object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API Data operator/(const boost::python::object& left, const Data& right);



/**
  \brief
  Output operator
*/
ESCRIPT_DLL_API std::ostream& operator<<(std::ostream& o, const Data& data);

/**
  \brief
  Compute a tensor product of two Data objects
  \param arg0 - Input - Data object
  \param arg1 - Input - Data object
  \param axis_offset - Input - axis offset
  \param transpose - Input - 0: transpose neither, 1: transpose arg0, 2: transpose arg1
*/
ESCRIPT_DLL_API
Data
C_GeneralTensorProduct(Data& arg0,
                     Data& arg1,
                     int axis_offset=0,
                     int transpose=0);

/**
  \brief
  Perform the given binary operation with this and right as operands.
  Right is a Data object.
*/
template <class BinaryFunction>
inline
void
Data::binaryOp(const Data& right,
               BinaryFunction operation)
{
   //
   // if this has a rank of zero promote it to the rank of the RHS
   if (getDataPointRank()==0 && right.getDataPointRank()!=0) {
     throw DataException("Error - attempt to update rank zero object with object with rank bigger than zero.");
   }

   if (isLazy() || right.isLazy())
   {
     throw DataException("Programmer error - attempt to call binaryOp with Lazy Data.");
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
   typeMatchRight(tempRight);
   //
   // Need to cast to the concrete types so that the correct binaryOp
   // is called.
   if (isExpanded()) {
     //
     // Expanded data will be done in parallel, the right hand side can be
     // of any data type
     DataExpanded* leftC=dynamic_cast<DataExpanded*>(m_data.get());
     EsysAssert((leftC!=0), "Programming error - casting to DataExpanded.");
     escript::binaryOp(*leftC,*(tempRight.getReady()),operation);
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
   } else if (isConstant()) {
     DataConstant* leftC=dynamic_cast<DataConstant*>(m_data.get());
     DataConstant* rightC=dynamic_cast<DataConstant*>(tempRight.m_data.get());
     EsysAssert((leftC!=0 && rightC!=0), "Programming error - casting to DataConstant.");
     escript::binaryOp(*leftC,*rightC,operation);
   }
}

/**
  \brief
  Perform the given Data object reduction algorithm on this and return the result.
  Given operation combines each element of each data point, thus argument
  object (*this) is a rank n Data object, and returned object is a scalar.
  Calls escript::algorithm.
*/
template <class BinaryFunction>
inline
double
Data::algorithm(BinaryFunction operation, double initial_value) const
{
  if (isExpanded()) {
    DataExpanded* leftC=dynamic_cast<DataExpanded*>(m_data.get());
    EsysAssert((leftC!=0), "Programming error - casting to DataExpanded.");
    return escript::algorithm(*leftC,operation,initial_value);
  } else if (isTagged()) {
    DataTagged* leftC=dynamic_cast<DataTagged*>(m_data.get());
    EsysAssert((leftC!=0), "Programming error - casting to DataTagged.");
    return escript::algorithm(*leftC,operation,initial_value);
  } else if (isConstant()) {
    DataConstant* leftC=dynamic_cast<DataConstant*>(m_data.get());
    EsysAssert((leftC!=0), "Programming error - casting to DataConstant.");
    return escript::algorithm(*leftC,operation,initial_value);
  } else if (isEmpty()) {
    throw DataException("Error - Operations not permitted on instances of DataEmpty.");
  } else if (isLazy()) {
    throw DataException("Error - Operations not permitted on instances of DataLazy.");
  } else {
    throw DataException("Error - Data encapsulates an unknown type.");
  }
}

/**
  \brief
  Perform the given data point reduction algorithm on data and return the result.
  Given operation combines each element within each data point into a scalar,
  thus argument object is a rank n Data object, and returned object is a
  rank 0 Data object.
  Calls escript::dp_algorithm.
*/
template <class BinaryFunction>
inline
Data
Data::dp_algorithm(BinaryFunction operation, double initial_value) const
{
  if (isEmpty()) {
    throw DataException("Error - Operations not permitted on instances of DataEmpty.");
  } 
  else if (isExpanded()) {
    Data result(0,DataTypes::ShapeType(),getFunctionSpace(),isExpanded());
    DataExpanded* dataE=dynamic_cast<DataExpanded*>(m_data.get());
    DataExpanded* resultE=dynamic_cast<DataExpanded*>(result.m_data.get());
    EsysAssert((dataE!=0), "Programming error - casting data to DataExpanded.");
    EsysAssert((resultE!=0), "Programming error - casting result to DataExpanded.");
    escript::dp_algorithm(*dataE,*resultE,operation,initial_value);
    return result;
  }
  else if (isTagged()) {
    DataTagged* dataT=dynamic_cast<DataTagged*>(m_data.get());
    EsysAssert((dataT!=0), "Programming error - casting data to DataTagged.");
    DataTypes::ValueType defval(1);
    defval[0]=0;
    DataTagged* resultT=new DataTagged(getFunctionSpace(), DataTypes::scalarShape, defval, dataT);
    escript::dp_algorithm(*dataT,*resultT,operation,initial_value);
    return Data(resultT);   // note: the Data object now owns the resultT pointer
  } 
  else if (isConstant()) {
    Data result(0,DataTypes::ShapeType(),getFunctionSpace(),isExpanded());
    DataConstant* dataC=dynamic_cast<DataConstant*>(m_data.get());
    DataConstant* resultC=dynamic_cast<DataConstant*>(result.m_data.get());
    EsysAssert((dataC!=0), "Programming error - casting data to DataConstant.");
    EsysAssert((resultC!=0), "Programming error - casting result to DataConstant.");
    escript::dp_algorithm(*dataC,*resultC,operation,initial_value);
    return result;
  } else if (isLazy()) {
    throw DataException("Error - Operations not permitted on instances of DataLazy.");
  } else {
    throw DataException("Error - Data encapsulates an unknown type.");
  }
}

/**
  \brief
  Compute a tensor operation with two Data objects
  \param arg0 - Input - Data object
  \param arg1 - Input - Data object
  \param operation - Input - Binary op functor
*/
template <typename BinaryFunction>
inline
Data
C_TensorBinaryOperation(Data const &arg_0,
                        Data const &arg_1,
                        BinaryFunction operation)
{
  if (arg_0.isEmpty() || arg_1.isEmpty())
  {
     throw DataException("Error - Operations not permitted on instances of DataEmpty.");
  }
  if (arg_0.isLazy() || arg_1.isLazy())
  {
     throw DataException("Error - Operations not permitted on lazy data.");
  }
  // Interpolate if necessary and find an appropriate function space
  Data arg_0_Z, arg_1_Z;
  if (arg_0.getFunctionSpace()!=arg_1.getFunctionSpace()) {
    if (arg_0.probeInterpolation(arg_1.getFunctionSpace())) {
      arg_0_Z = arg_0.interpolate(arg_1.getFunctionSpace());
      arg_1_Z = Data(arg_1);
    }
    else if (arg_1.probeInterpolation(arg_0.getFunctionSpace())) {
      arg_1_Z=arg_1.interpolate(arg_0.getFunctionSpace());
      arg_0_Z =Data(arg_0);
    }
    else {
      throw DataException("Error - C_TensorBinaryOperation: arguments have incompatible function spaces.");
    }
  } else {
      arg_0_Z = Data(arg_0);
      arg_1_Z = Data(arg_1);
  }
  // Get rank and shape of inputs
  int rank0 = arg_0_Z.getDataPointRank();
  int rank1 = arg_1_Z.getDataPointRank();
  DataTypes::ShapeType shape0 = arg_0_Z.getDataPointShape();
  DataTypes::ShapeType shape1 = arg_1_Z.getDataPointShape();
  int size0 = arg_0_Z.getDataPointSize();
  int size1 = arg_1_Z.getDataPointSize();

  // Declare output Data object
  Data res;

  if (shape0 == shape1) {

    if (arg_0_Z.isConstant()   && arg_1_Z.isConstant()) {
      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace());      // DataConstant output
/*      double *ptr_0 = &((arg_0_Z.getPointDataView().getData())[0]);
      double *ptr_1 = &((arg_1_Z.getPointDataView().getData())[0]);
      double *ptr_2 = &((res.getPointDataView().getData())[0]);*/
      double *ptr_0 = &(arg_0_Z.getDataAtOffset(0));
      double *ptr_1 = &(arg_1_Z.getDataAtOffset(0));
      double *ptr_2 = &(res.getDataAtOffset(0));

      tensor_binary_operation(size0, ptr_0, ptr_1, ptr_2, operation);
    }
    else if (arg_0_Z.isConstant()   && arg_1_Z.isTagged()) {

      // Prepare the DataConstant input
      DataConstant* tmp_0=dynamic_cast<DataConstant*>(arg_0_Z.borrowData());

      // Borrow DataTagged input from Data object
      DataTagged* tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());

      // Prepare a DataTagged output 2
      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace());      // DataTagged output
      res.tag();
      DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());

      // Prepare offset into DataConstant
      int offset_0 = tmp_0->getPointOffset(0,0);
      double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
      // Get the views
//       DataArrayView view_1 = tmp_1->getDefaultValue();
//       DataArrayView view_2 = tmp_2->getDefaultValue();
//       // Get the pointers to the actual data
//       double *ptr_1 = &((view_1.getData())[0]);
//       double *ptr_2 = &((view_2.getData())[0]);

      // Get the pointers to the actual data
      double *ptr_1 = &(tmp_1->getDefaultValue(0));
      double *ptr_2 = &(tmp_2->getDefaultValue(0));

      // Compute a result for the default
      tensor_binary_operation(size0, ptr_0, ptr_1, ptr_2, operation);
      // Compute a result for each tag
      const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_1.begin();i!=lookup_1.end();i++) {
        tmp_2->addTag(i->first);
/*        DataArrayView view_1 = tmp_1->getDataPointByTag(i->first);
        DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
        double *ptr_1 = &view_1.getData(0);
        double *ptr_2 = &view_2.getData(0);*/
        double *ptr_1 = &(tmp_1->getDataByTag(i->first,0));
        double *ptr_2 = &(tmp_2->getDataByTag(i->first,0));

        tensor_binary_operation(size0, ptr_0, ptr_1, ptr_2, operation);
      }

    }
    else if (arg_0_Z.isConstant()   && arg_1_Z.isExpanded()) {

      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataConstant* tmp_0=dynamic_cast<DataConstant*>(arg_0_Z.borrowData());
      DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_1,dataPointNo_1;
      int numSamples_1 = arg_1_Z.getNumSamples();
      int numDataPointsPerSample_1 = arg_1_Z.getNumDataPointsPerSample();
      int offset_0 = tmp_0->getPointOffset(0,0);
      #pragma omp parallel for private(sampleNo_1,dataPointNo_1) schedule(static)
      for (sampleNo_1 = 0; sampleNo_1 < numSamples_1; sampleNo_1++) {
        for (dataPointNo_1 = 0; dataPointNo_1 < numDataPointsPerSample_1; dataPointNo_1++) {
          int offset_1 = tmp_1->getPointOffset(sampleNo_1,dataPointNo_1);
          int offset_2 = tmp_2->getPointOffset(sampleNo_1,dataPointNo_1);
//           double *ptr_0 = &((arg_0_Z.getPointDataView().getData())[offset_0]);
//           double *ptr_1 = &((arg_1_Z.getPointDataView().getData())[offset_1]);
//           double *ptr_2 = &((res.getPointDataView().getData())[offset_2]);

          double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
          double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1)); 
          double *ptr_2 = &(res.getDataAtOffset(offset_2)); 
          tensor_binary_operation(size0, ptr_0, ptr_1, ptr_2, operation);
        }
      }

    }
    else if (arg_0_Z.isTagged()     && arg_1_Z.isConstant()) {

      // Borrow DataTagged input from Data object
      DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());

      // Prepare the DataConstant input
      DataConstant* tmp_1=dynamic_cast<DataConstant*>(arg_1_Z.borrowData());

      // Prepare a DataTagged output 2
      res = Data(0.0, shape0, arg_0_Z.getFunctionSpace());      // DataTagged output
      res.tag();
      DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());

      // Prepare offset into DataConstant
      int offset_1 = tmp_1->getPointOffset(0,0);
//       double *ptr_1 = &((arg_1_Z.getPointDataView().getData())[offset_1]);
      double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
      // Get the views
//       DataArrayView view_0 = tmp_0->getDefaultValue();
//       DataArrayView view_2 = tmp_2->getDefaultValue();
//       // Get the pointers to the actual data
//       double *ptr_0 = &((view_0.getData())[0]);
//       double *ptr_2 = &((view_2.getData())[0]);
      // Get the pointers to the actual data
      double *ptr_0 = &(tmp_0->getDefaultValue(0));
      double *ptr_2 = &(tmp_2->getDefaultValue(0));
      // Compute a result for the default
      tensor_binary_operation(size0, ptr_0, ptr_1, ptr_2, operation);
      // Compute a result for each tag
      const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_0.begin();i!=lookup_0.end();i++) {
        tmp_2->addTag(i->first);
//         DataArrayView view_0 = tmp_0->getDataPointByTag(i->first);
//         DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
//         double *ptr_0 = &view_0.getData(0);
//         double *ptr_2 = &view_2.getData(0);
        double *ptr_0 = &(tmp_0->getDataByTag(i->first,0));
        double *ptr_2 = &(tmp_2->getDataByTag(i->first,0));
        tensor_binary_operation(size0, ptr_0, ptr_1, ptr_2, operation);
      }

    }
    else if (arg_0_Z.isTagged()     && arg_1_Z.isTagged()) {

      // Borrow DataTagged input from Data object
      DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());

      // Borrow DataTagged input from Data object
      DataTagged* tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());

      // Prepare a DataTagged output 2
      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace());
      res.tag();        // DataTagged output
      DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());

//       // Get the views
//       DataArrayView view_0 = tmp_0->getDefaultValue();
//       DataArrayView view_1 = tmp_1->getDefaultValue();
//       DataArrayView view_2 = tmp_2->getDefaultValue();
//       // Get the pointers to the actual data
//       double *ptr_0 = &((view_0.getData())[0]);
//       double *ptr_1 = &((view_1.getData())[0]);
//       double *ptr_2 = &((view_2.getData())[0]);

      // Get the pointers to the actual data
      double *ptr_0 = &(tmp_0->getDefaultValue(0));
      double *ptr_1 = &(tmp_1->getDefaultValue(0));
      double *ptr_2 = &(tmp_2->getDefaultValue(0));

      // Compute a result for the default
      tensor_binary_operation(size0, ptr_0, ptr_1, ptr_2, operation);
      // Merge the tags
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
      const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
      for (i=lookup_0.begin();i!=lookup_0.end();i++) {
        tmp_2->addTag(i->first); // use tmp_2 to get correct shape
      }
      for (i=lookup_1.begin();i!=lookup_1.end();i++) {
        tmp_2->addTag(i->first);
      }
      // Compute a result for each tag
      const DataTagged::DataMapType& lookup_2=tmp_2->getTagLookup();
      for (i=lookup_2.begin();i!=lookup_2.end();i++) {

//         DataArrayView view_0 = tmp_0->getDataPointByTag(i->first);
//         DataArrayView view_1 = tmp_1->getDataPointByTag(i->first);
//         DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
//         double *ptr_0 = &view_0.getData(0);
//         double *ptr_1 = &view_1.getData(0);
//         double *ptr_2 = &view_2.getData(0);

        double *ptr_0 = &(tmp_0->getDataByTag(i->first,0));
        double *ptr_1 = &(tmp_1->getDataByTag(i->first,0));
        double *ptr_2 = &(tmp_2->getDataByTag(i->first,0));

        tensor_binary_operation(size0, ptr_0, ptr_1, ptr_2, operation);
      }

    }
    else if (arg_0_Z.isTagged()     && arg_1_Z.isExpanded()) {

      // After finding a common function space above the two inputs have the same numSamples and num DPPS
      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataTagged*   tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());
      DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_0,dataPointNo_0;
      int numSamples_0 = arg_0_Z.getNumSamples();
      int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,0); // They're all the same, so just use #0
        double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);

//           double *ptr_1 = &((arg_1_Z.getPointDataView().getData())[offset_1]);
//           double *ptr_2 = &((res.getPointDataView().getData())[offset_2]);
          double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
          double *ptr_2 = &(res.getDataAtOffset(offset_2));


          tensor_binary_operation(size0, ptr_0, ptr_1, ptr_2, operation);
        }
      }

    }
    else if (arg_0_Z.isExpanded()   && arg_1_Z.isConstant()) {

      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
      DataConstant* tmp_1=dynamic_cast<DataConstant*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_0,dataPointNo_0;
      int numSamples_0 = arg_0_Z.getNumSamples();
      int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
      int offset_1 = tmp_1->getPointOffset(0,0);
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);

//           double *ptr_0 = &((arg_0_Z.getPointDataView().getData())[offset_0]);
//           double *ptr_1 = &((arg_1_Z.getPointDataView().getData())[offset_1]);
//           double *ptr_2 = &((res.getPointDataView().getData())[offset_2]);

          double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
          double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
          double *ptr_2 = &(res.getDataAtOffset(offset_2));


          tensor_binary_operation(size0, ptr_0, ptr_1, ptr_2, operation);
        }
      }

    }
    else if (arg_0_Z.isExpanded()   && arg_1_Z.isTagged()) {

      // After finding a common function space above the two inputs have the same numSamples and num DPPS
      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
      DataTagged*   tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_0,dataPointNo_0;
      int numSamples_0 = arg_0_Z.getNumSamples();
      int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        int offset_1 = tmp_1->getPointOffset(sampleNo_0,0);
        double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
          double *ptr_2 = &(res.getDataAtOffset(offset_2));
          tensor_binary_operation(size0, ptr_0, ptr_1, ptr_2, operation);
        }
      }

    }
    else if (arg_0_Z.isExpanded()   && arg_1_Z.isExpanded()) {

      // After finding a common function space above the two inputs have the same numSamples and num DPPS
      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
      DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_0,dataPointNo_0;
      int numSamples_0 = arg_0_Z.getNumSamples();
      int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
          double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
          double *ptr_2 = &(res.getDataAtOffset(offset_2));
          tensor_binary_operation(size0, ptr_0, ptr_1, ptr_2, operation);
        }
      }

    }
    else {
      throw DataException("Error - C_TensorBinaryOperation: unknown combination of inputs");
    }

  } else if (0 == rank0) {

    if (arg_0_Z.isConstant()   && arg_1_Z.isConstant()) {
      res = Data(0.0, shape1, arg_1_Z.getFunctionSpace());      // DataConstant output
      double *ptr_0 = &(arg_0_Z.getDataAtOffset(0));
      double *ptr_1 = &(arg_1_Z.getDataAtOffset(0));
      double *ptr_2 = &(res.getDataAtOffset(0));
      tensor_binary_operation(size1, ptr_0[0], ptr_1, ptr_2, operation);
    }
    else if (arg_0_Z.isConstant()   && arg_1_Z.isTagged()) {

      // Prepare the DataConstant input
      DataConstant* tmp_0=dynamic_cast<DataConstant*>(arg_0_Z.borrowData());

      // Borrow DataTagged input from Data object
      DataTagged* tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());

      // Prepare a DataTagged output 2
      res = Data(0.0, shape1, arg_1_Z.getFunctionSpace());      // DataTagged output
      res.tag();
      DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());

      // Prepare offset into DataConstant
      int offset_0 = tmp_0->getPointOffset(0,0);
      double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
      // Get the views
//       DataArrayView view_1 = tmp_1->getDefaultValue();
//       DataArrayView view_2 = tmp_2->getDefaultValue();
//       // Get the pointers to the actual data
//       double *ptr_1 = &((view_1.getData())[0]);
//       double *ptr_2 = &((view_2.getData())[0]);
       double *ptr_1 = &(tmp_1->getDefaultValue(0));
       double *ptr_2 = &(tmp_2->getDefaultValue(0));

      // Compute a result for the default
      tensor_binary_operation(size1, ptr_0[0], ptr_1, ptr_2, operation);
      // Compute a result for each tag
      const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_1.begin();i!=lookup_1.end();i++) {
        tmp_2->addTag(i->first);
//         DataArrayView view_1 = tmp_1->getDataPointByTag(i->first);
//         DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
//         double *ptr_1 = &view_1.getData(0);
//         double *ptr_2 = &view_2.getData(0);
        double *ptr_1 = &(tmp_1->getDataByTag(i->first,0));
        double *ptr_2 = &(tmp_2->getDataByTag(i->first,0));
        tensor_binary_operation(size1, ptr_0[0], ptr_1, ptr_2, operation);
      }

    }
    else if (arg_0_Z.isConstant()   && arg_1_Z.isExpanded()) {

      res = Data(0.0, shape1, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataConstant* tmp_0=dynamic_cast<DataConstant*>(arg_0_Z.borrowData());
      DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_1,dataPointNo_1;
      int numSamples_1 = arg_1_Z.getNumSamples();
      int numDataPointsPerSample_1 = arg_1_Z.getNumDataPointsPerSample();
      int offset_0 = tmp_0->getPointOffset(0,0);
      #pragma omp parallel for private(sampleNo_1,dataPointNo_1) schedule(static)
      for (sampleNo_1 = 0; sampleNo_1 < numSamples_1; sampleNo_1++) {
        for (dataPointNo_1 = 0; dataPointNo_1 < numDataPointsPerSample_1; dataPointNo_1++) {
          int offset_1 = tmp_1->getPointOffset(sampleNo_1,dataPointNo_1);
          int offset_2 = tmp_2->getPointOffset(sampleNo_1,dataPointNo_1);
          double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
          double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
          double *ptr_2 = &(res.getDataAtOffset(offset_2));
          tensor_binary_operation(size1, ptr_0[0], ptr_1, ptr_2, operation);

        }
      }

    }
    else if (arg_0_Z.isTagged()     && arg_1_Z.isConstant()) {

      // Borrow DataTagged input from Data object
      DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());

      // Prepare the DataConstant input
      DataConstant* tmp_1=dynamic_cast<DataConstant*>(arg_1_Z.borrowData());

      // Prepare a DataTagged output 2
      res = Data(0.0, shape1, arg_0_Z.getFunctionSpace());      // DataTagged output
      res.tag();
      DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());

      // Prepare offset into DataConstant
      int offset_1 = tmp_1->getPointOffset(0,0);
//       double *ptr_1 = &((arg_1_Z.getPointDataView().getData())[offset_1]);
      double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
      // Get the views
/*      DataArrayView view_0 = tmp_0->getDefaultValue();
      DataArrayView view_2 = tmp_2->getDefaultValue();
      // Get the pointers to the actual data
      double *ptr_0 = &((view_0.getData())[0]);
      double *ptr_2 = &((view_2.getData())[0]);*/

      // Get the pointers to the actual data
      double *ptr_0 = &(tmp_0->getDefaultValue(0));
      double *ptr_2 = &(tmp_2->getDefaultValue(0));


      // Compute a result for the default
      tensor_binary_operation(size1, ptr_0[0], ptr_1, ptr_2, operation);
      // Compute a result for each tag
      const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_0.begin();i!=lookup_0.end();i++) {
        tmp_2->addTag(i->first);
/*        DataArrayView view_0 = tmp_0->getDataPointByTag(i->first);
        DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
        double *ptr_0 = &view_0.getData(0);
        double *ptr_2 = &view_2.getData(0);*/
        double *ptr_0 = &(tmp_0->getDataByTag(i->first,0));
        double *ptr_2 = &(tmp_2->getDataByTag(i->first,0));

        tensor_binary_operation(size1, ptr_0[0], ptr_1, ptr_2, operation);
      }

    }
    else if (arg_0_Z.isTagged()     && arg_1_Z.isTagged()) {

      // Borrow DataTagged input from Data object
      DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());

      // Borrow DataTagged input from Data object
      DataTagged* tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());

      // Prepare a DataTagged output 2
      res = Data(0.0, shape1, arg_1_Z.getFunctionSpace());
      res.tag();        // DataTagged output
      DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());

      // Get the views
/*      DataArrayView view_0 = tmp_0->getDefaultValue();
      DataArrayView view_1 = tmp_1->getDefaultValue();
      DataArrayView view_2 = tmp_2->getDefaultValue();
      // Get the pointers to the actual data
      double *ptr_0 = &((view_0.getData())[0]);
      double *ptr_1 = &((view_1.getData())[0]);
      double *ptr_2 = &((view_2.getData())[0]);*/

      // Get the pointers to the actual data
      double *ptr_0 = &(tmp_0->getDefaultValue(0));
      double *ptr_1 = &(tmp_1->getDefaultValue(0));
      double *ptr_2 = &(tmp_2->getDefaultValue(0));


      // Compute a result for the default
      tensor_binary_operation(size1, ptr_0[0], ptr_1, ptr_2, operation);
      // Merge the tags
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
      const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
      for (i=lookup_0.begin();i!=lookup_0.end();i++) {
        tmp_2->addTag(i->first); // use tmp_2 to get correct shape
      }
      for (i=lookup_1.begin();i!=lookup_1.end();i++) {
        tmp_2->addTag(i->first);
      }
      // Compute a result for each tag
      const DataTagged::DataMapType& lookup_2=tmp_2->getTagLookup();
      for (i=lookup_2.begin();i!=lookup_2.end();i++) {

/*        DataArrayView view_0 = tmp_0->getDataPointByTag(i->first);
        DataArrayView view_1 = tmp_1->getDataPointByTag(i->first);
        DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
        double *ptr_0 = &view_0.getData(0);
        double *ptr_1 = &view_1.getData(0);
        double *ptr_2 = &view_2.getData(0);*/

        double *ptr_0 = &(tmp_0->getDataByTag(i->first,0));
        double *ptr_1 = &(tmp_1->getDataByTag(i->first,0));
        double *ptr_2 = &(tmp_2->getDataByTag(i->first,0));

        tensor_binary_operation(size1, ptr_0[0], ptr_1, ptr_2, operation);
      }

    }
    else if (arg_0_Z.isTagged()     && arg_1_Z.isExpanded()) {

      // After finding a common function space above the two inputs have the same numSamples and num DPPS
      res = Data(0.0, shape1, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataTagged*   tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());
      DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_0,dataPointNo_0;
      int numSamples_0 = arg_0_Z.getNumSamples();
      int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,0); // They're all the same, so just use #0
        double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
          double *ptr_2 = &(res.getDataAtOffset(offset_2));
          tensor_binary_operation(size1, ptr_0[0], ptr_1, ptr_2, operation);
        }
      }

    }
    else if (arg_0_Z.isExpanded()   && arg_1_Z.isConstant()) {

      res = Data(0.0, shape1, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
      DataConstant* tmp_1=dynamic_cast<DataConstant*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_0,dataPointNo_0;
      int numSamples_0 = arg_0_Z.getNumSamples();
      int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
      int offset_1 = tmp_1->getPointOffset(0,0);
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
          double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
          double *ptr_2 = &(res.getDataAtOffset(offset_2));
          tensor_binary_operation(size1, ptr_0[0], ptr_1, ptr_2, operation);
        }
      }


    }
    else if (arg_0_Z.isExpanded()   && arg_1_Z.isTagged()) {

      // After finding a common function space above the two inputs have the same numSamples and num DPPS
      res = Data(0.0, shape1, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
      DataTagged*   tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_0,dataPointNo_0;
      int numSamples_0 = arg_0_Z.getNumSamples();
      int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        int offset_1 = tmp_1->getPointOffset(sampleNo_0,0);
        double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
          double *ptr_2 = &(res.getDataAtOffset(offset_2));
          tensor_binary_operation(size1, ptr_0[0], ptr_1, ptr_2, operation);
        }
      }

    }
    else if (arg_0_Z.isExpanded()   && arg_1_Z.isExpanded()) {

      // After finding a common function space above the two inputs have the same numSamples and num DPPS
      res = Data(0.0, shape1, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
      DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_0,dataPointNo_0;
      int numSamples_0 = arg_0_Z.getNumSamples();
      int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
          double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
          double *ptr_2 = &(res.getDataAtOffset(offset_2));
          tensor_binary_operation(size1, ptr_0[0], ptr_1, ptr_2, operation);
        }
      }

    }
    else {
      throw DataException("Error - C_TensorBinaryOperation: unknown combination of inputs");
    }

  } else if (0 == rank1) {

    if (arg_0_Z.isConstant()   && arg_1_Z.isConstant()) {
      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace());      // DataConstant output
      double *ptr_0 = &(arg_0_Z.getDataAtOffset(0));
      double *ptr_1 = &(arg_1_Z.getDataAtOffset(0));
      double *ptr_2 = &(res.getDataAtOffset(0));
      tensor_binary_operation(size0, ptr_0, ptr_1[0], ptr_2, operation);
    }
    else if (arg_0_Z.isConstant()   && arg_1_Z.isTagged()) {

      // Prepare the DataConstant input
      DataConstant* tmp_0=dynamic_cast<DataConstant*>(arg_0_Z.borrowData());

      // Borrow DataTagged input from Data object
      DataTagged* tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());

      // Prepare a DataTagged output 2
      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace());      // DataTagged output
      res.tag();
      DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());

      // Prepare offset into DataConstant
      int offset_0 = tmp_0->getPointOffset(0,0);
      double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
      // Get the views
/*      DataArrayView view_1 = tmp_1->getDefaultValue();
      DataArrayView view_2 = tmp_2->getDefaultValue();
      // Get the pointers to the actual data
      double *ptr_1 = &((view_1.getData())[0]);
      double *ptr_2 = &((view_2.getData())[0]);*/
      //Get the pointers to the actual data
      double *ptr_1 = &(tmp_1->getDefaultValue(0));
      double *ptr_2 = &(tmp_2->getDefaultValue(0));

      // Compute a result for the default
      tensor_binary_operation(size0, ptr_0, ptr_1[0], ptr_2, operation);
      // Compute a result for each tag
      const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_1.begin();i!=lookup_1.end();i++) {
        tmp_2->addTag(i->first);
//         DataArrayView view_1 = tmp_1->getDataPointByTag(i->first);
//         DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
//         double *ptr_1 = &view_1.getData(0);
//         double *ptr_2 = &view_2.getData(0);
        double *ptr_1 = &(tmp_1->getDataByTag(i->first,0));
        double *ptr_2 = &(tmp_2->getDataByTag(i->first,0));


        tensor_binary_operation(size0, ptr_0, ptr_1[0], ptr_2, operation);
      }

    }
    else if (arg_0_Z.isConstant()   && arg_1_Z.isExpanded()) {

      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataConstant* tmp_0=dynamic_cast<DataConstant*>(arg_0_Z.borrowData());
      DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_1,dataPointNo_1;
      int numSamples_1 = arg_1_Z.getNumSamples();
      int numDataPointsPerSample_1 = arg_1_Z.getNumDataPointsPerSample();
      int offset_0 = tmp_0->getPointOffset(0,0);
      #pragma omp parallel for private(sampleNo_1,dataPointNo_1) schedule(static)
      for (sampleNo_1 = 0; sampleNo_1 < numSamples_1; sampleNo_1++) {
        for (dataPointNo_1 = 0; dataPointNo_1 < numDataPointsPerSample_1; dataPointNo_1++) {
          int offset_1 = tmp_1->getPointOffset(sampleNo_1,dataPointNo_1);
          int offset_2 = tmp_2->getPointOffset(sampleNo_1,dataPointNo_1);
          double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
          double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
          double *ptr_2 = &(res.getDataAtOffset(offset_2));
          tensor_binary_operation(size0, ptr_0, ptr_1[0], ptr_2, operation);
        }
      }

    }
    else if (arg_0_Z.isTagged()     && arg_1_Z.isConstant()) {

      // Borrow DataTagged input from Data object
      DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());

      // Prepare the DataConstant input
      DataConstant* tmp_1=dynamic_cast<DataConstant*>(arg_1_Z.borrowData());

      // Prepare a DataTagged output 2
      res = Data(0.0, shape0, arg_0_Z.getFunctionSpace());      // DataTagged output
      res.tag();
      DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());

      // Prepare offset into DataConstant
      int offset_1 = tmp_1->getPointOffset(0,0);
      double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
      // Get the views
//       DataArrayView view_0 = tmp_0->getDefaultValue();
//       DataArrayView view_2 = tmp_2->getDefaultValue();
//       // Get the pointers to the actual data
//       double *ptr_0 = &((view_0.getData())[0]);
//       double *ptr_2 = &((view_2.getData())[0]);
      // Get the pointers to the actual data
      double *ptr_0 = &(tmp_0->getDefaultValue(0));
      double *ptr_2 = &(tmp_2->getDefaultValue(0));
      // Compute a result for the default
      tensor_binary_operation(size0, ptr_0, ptr_1[0], ptr_2, operation);
      // Compute a result for each tag
      const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_0.begin();i!=lookup_0.end();i++) {
        tmp_2->addTag(i->first);
/*        DataArrayView view_0 = tmp_0->getDataPointByTag(i->first);
        DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
        double *ptr_0 = &view_0.getData(0);
        double *ptr_2 = &view_2.getData(0);*/
        double *ptr_0 = &(tmp_0->getDataByTag(i->first,0));
        double *ptr_2 = &(tmp_2->getDataByTag(i->first,0));
        tensor_binary_operation(size0, ptr_0, ptr_1[0], ptr_2, operation);
      }

    }
    else if (arg_0_Z.isTagged()     && arg_1_Z.isTagged()) {

      // Borrow DataTagged input from Data object
      DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());

      // Borrow DataTagged input from Data object
      DataTagged* tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());

      // Prepare a DataTagged output 2
      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace());
      res.tag();        // DataTagged output
      DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());

      // Get the views
//       DataArrayView view_0 = tmp_0->getDefaultValue();
//       DataArrayView view_1 = tmp_1->getDefaultValue();
//       DataArrayView view_2 = tmp_2->getDefaultValue();
//       // Get the pointers to the actual data
//       double *ptr_0 = &((view_0.getData())[0]);
//       double *ptr_1 = &((view_1.getData())[0]);
//       double *ptr_2 = &((view_2.getData())[0]);

      // Get the pointers to the actual data
      double *ptr_0 = &(tmp_0->getDefaultValue(0));
      double *ptr_1 = &(tmp_1->getDefaultValue(0));
      double *ptr_2 = &(tmp_2->getDefaultValue(0));

      // Compute a result for the default
      tensor_binary_operation(size0, ptr_0, ptr_1[0], ptr_2, operation);
      // Merge the tags
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
      const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
      for (i=lookup_0.begin();i!=lookup_0.end();i++) {
        tmp_2->addTag(i->first); // use tmp_2 to get correct shape
      }
      for (i=lookup_1.begin();i!=lookup_1.end();i++) {
        tmp_2->addTag(i->first);
      }
      // Compute a result for each tag
      const DataTagged::DataMapType& lookup_2=tmp_2->getTagLookup();
      for (i=lookup_2.begin();i!=lookup_2.end();i++) {
//         DataArrayView view_0 = tmp_0->getDataPointByTag(i->first);
//         DataArrayView view_1 = tmp_1->getDataPointByTag(i->first);
//         DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
//         double *ptr_0 = &view_0.getData(0);
//         double *ptr_1 = &view_1.getData(0);
//         double *ptr_2 = &view_2.getData(0);

        double *ptr_0 = &(tmp_0->getDataByTag(i->first,0));
        double *ptr_1 = &(tmp_1->getDataByTag(i->first,0));
        double *ptr_2 = &(tmp_2->getDataByTag(i->first,0));
        tensor_binary_operation(size0, ptr_0, ptr_1[0], ptr_2, operation);
      }

    }
    else if (arg_0_Z.isTagged()     && arg_1_Z.isExpanded()) {

      // After finding a common function space above the two inputs have the same numSamples and num DPPS
      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataTagged*   tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());
      DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_0,dataPointNo_0;
      int numSamples_0 = arg_0_Z.getNumSamples();
      int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,0); // They're all the same, so just use #0
        double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
          double *ptr_2 = &(res.getDataAtOffset(offset_2));
          tensor_binary_operation(size0, ptr_0, ptr_1[0], ptr_2, operation);
        }
      }

    }
    else if (arg_0_Z.isExpanded()   && arg_1_Z.isConstant()) {

      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
      DataConstant* tmp_1=dynamic_cast<DataConstant*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_0,dataPointNo_0;
      int numSamples_0 = arg_0_Z.getNumSamples();
      int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
      int offset_1 = tmp_1->getPointOffset(0,0);
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
          double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
          double *ptr_2 = &(res.getDataAtOffset(offset_2));
          tensor_binary_operation(size0, ptr_0, ptr_1[0], ptr_2, operation);
        }
      }


    }
    else if (arg_0_Z.isExpanded()   && arg_1_Z.isTagged()) {

      // After finding a common function space above the two inputs have the same numSamples and num DPPS
      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
      DataTagged*   tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_0,dataPointNo_0;
      int numSamples_0 = arg_0_Z.getNumSamples();
      int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        int offset_1 = tmp_1->getPointOffset(sampleNo_0,0);
        double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
          double *ptr_2 = &(res.getDataAtOffset(offset_2));
          tensor_binary_operation(size0, ptr_0, ptr_1[0], ptr_2, operation);
        }
      }

    }
    else if (arg_0_Z.isExpanded()   && arg_1_Z.isExpanded()) {

      // After finding a common function space above the two inputs have the same numSamples and num DPPS
      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
      DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_0,dataPointNo_0;
      int numSamples_0 = arg_0_Z.getNumSamples();
      int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
          double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
          double *ptr_2 = &(res.getDataAtOffset(offset_2));
          tensor_binary_operation(size0, ptr_0, ptr_1[0], ptr_2, operation);
        }
      }

    }
    else {
      throw DataException("Error - C_TensorBinaryOperation: unknown combination of inputs");
    }

  } else {
    throw DataException("Error - C_TensorBinaryOperation: arguments have incompatible shapes");
  }

  return res;
}

template <typename UnaryFunction>
Data
C_TensorUnaryOperation(Data const &arg_0,
                       UnaryFunction operation)
{
  if (arg_0.isEmpty())	// do this before we attempt to interpolate
  {
     throw DataException("Error - Operations not permitted on instances of DataEmpty.");
  }
  if (arg_0.isLazy())
  {
     throw DataException("Error - Operations not permitted on lazy data.");
  }
  // Interpolate if necessary and find an appropriate function space
  Data arg_0_Z = Data(arg_0);

  // Get rank and shape of inputs
  int rank0 = arg_0_Z.getDataPointRank();
  const DataTypes::ShapeType& shape0 = arg_0_Z.getDataPointShape();
  int size0 = arg_0_Z.getDataPointSize();

  // Declare output Data object
  Data res;

  if (arg_0_Z.isConstant()) {
    res = Data(0.0, shape0, arg_0_Z.getFunctionSpace());      // DataConstant output
//     double *ptr_0 = &((arg_0_Z.getPointDataView().getData())[0]);
//     double *ptr_2 = &((res.getPointDataView().getData())[0]);
    double *ptr_0 = &(arg_0_Z.getDataAtOffset(0));
    double *ptr_2 = &(res.getDataAtOffset(0));
    tensor_unary_operation(size0, ptr_0, ptr_2, operation);
  }
  else if (arg_0_Z.isTagged()) {

    // Borrow DataTagged input from Data object
    DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());

    // Prepare a DataTagged output 2
    res = Data(0.0, shape0, arg_0_Z.getFunctionSpace());   // DataTagged output
    res.tag();
    DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());

//     // Get the views
//     DataArrayView view_0 = tmp_0->getDefaultValue();
//     DataArrayView view_2 = tmp_2->getDefaultValue();
//     // Get the pointers to the actual data
//     double *ptr_0 = &((view_0.getData())[0]);
//     double *ptr_2 = &((view_2.getData())[0]);
    // Get the pointers to the actual data
    double *ptr_0 = &(tmp_0->getDefaultValue(0));
    double *ptr_2 = &(tmp_2->getDefaultValue(0));
    // Compute a result for the default
    tensor_unary_operation(size0, ptr_0, ptr_2, operation);
    // Compute a result for each tag
    const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
    DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
    for (i=lookup_0.begin();i!=lookup_0.end();i++) {
      tmp_2->addTag(i->first);
//       DataArrayView view_0 = tmp_0->getDataPointByTag(i->first);
//       DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
//       double *ptr_0 = &view_0.getData(0);
//       double *ptr_2 = &view_2.getData(0);
      double *ptr_0 = &(tmp_0->getDataByTag(i->first,0));
      double *ptr_2 = &(tmp_2->getDataByTag(i->first,0));
      tensor_unary_operation(size0, ptr_0, ptr_2, operation);
    }

  }
  else if (arg_0_Z.isExpanded()) {

    res = Data(0.0, shape0, arg_0_Z.getFunctionSpace(),true); // DataExpanded output
    DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

    int sampleNo_0,dataPointNo_0;
    int numSamples_0 = arg_0_Z.getNumSamples();
    int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
    #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
    for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
      for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
//         int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
//         int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
//         double *ptr_0 = &((arg_0_Z.getPointDataView().getData())[offset_0]);
//         double *ptr_2 = &((res.getPointDataView().getData())[offset_2]);
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
        double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
        double *ptr_2 = &(res.getDataAtOffset(offset_2));
        tensor_unary_operation(size0, ptr_0, ptr_2, operation);
      }
    }
  }
  else {
    throw DataException("Error - C_TensorUnaryOperation: unknown combination of inputs");
  }

  return res;
}

}
#endif
