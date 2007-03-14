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

extern "C" {
#include "DataC.h"
#include "paso/Paso.h"
}

#ifndef PASO_MPI
#define MPI_Comm long
#endif

#include <string>
#include <algorithm>

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
   Data creates the appropriate Data object for the given construction 
   arguments. 

   Description:
   Data is essentially a factory class which creates the appropriate Data
   object for the given construction arguments. It retains control over
   the object created for the lifetime of the object.
   The type of Data object referred to may change during the lifetime of
   the Data object.
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
     \brief
     Constructor which copies data from a DataArrayView.

     \param value - Input - Data value for a single point.
     \param what - Input - A description of what this data represents.
     \param expanded - Input - Flag, if true fill the entire container with
                       the value. Otherwise a more efficient storage
                       mechanism will be used.
  */
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
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
    ==>*
  */
  ESCRIPT_DLL_API
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
     \brief
     Destructor
  */
  ESCRIPT_DLL_API
  ~Data();

  /**
     \brief
     Perform a deep copy.
  */
  ESCRIPT_DLL_API
  void
  copy(const Data& other);

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

     The data-point number here corresponds to the data-point number in the
     numarray returned by convertToNumArray.
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

  /**
     \brief
     Return the C wrapper for the Data object - const version.
  */
  ESCRIPT_DLL_API
  escriptDataC
  getDataC() const;

  /**
     \brief
     Write the data as a string.
  */
  ESCRIPT_DLL_API
  inline
  std::string
  toString() const
  {
    return m_data->toString();
  }

  /**
     \brief
     Return the DataArrayView of the point data. This essentially contains
     the shape information for each data point although it also may be used
     to manipulate the point data.
  */
  ESCRIPT_DLL_API
  inline
  const DataArrayView&
  getPointDataView() const
  {
     return m_data->getPointDataView();
  }

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
    ==>*
  */
  ESCRIPT_DLL_API
  void
  tag();

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
     \brief
     Return true if this Data is empty.
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
  const AbstractDomain&
  getDomain() const
  {
     return getFunctionSpace().getDomain();
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
    return m_data->getPointDataView().getRank();
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
     dumps the object into a netCDF file 
  */
  ESCRIPT_DLL_API
  inline
  void
  dump(const std::string fileName) const
  {
    return m_data->dump(fileName);
  }

  /**
     \brief
     Return the sample data for the given sample no. This is not the
     preferred interface but is provided for use by C code.
     \param sampleNo - Input - the given sample no.
  */
  ESCRIPT_DLL_API
  inline
  DataAbstract::ValueType::value_type*
  getSampleData(DataAbstract::ValueType::size_type sampleNo)
  {
    return m_data->getSampleData(sampleNo);
  }

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
  ESCRIPT_DLL_API
  const DataArrayView::ShapeType&
  getDataPointShape() const;

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
  DataArrayView::ValueType::size_type
  getLength() const;

  /**
     \brief
     Assign the given value to the tag. Implicitly converts this
     object to type DataTagged. Throws an exception if this object
     cannot be converted to a DataTagged object.
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
     object to type DataTagged. Throws an exception if this object
     cannot be converted to a DataTagged object.
     \param tagKey - Input - Integer key.
     \param value - Input - Value to associate with given key.
    ==>*
  */
  ESCRIPT_DLL_API
  void
  setTaggedValueFromCPP(int tagKey,
                        const DataArrayView& value);

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
  integrate() const;

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
     *
  */
  ESCRIPT_DLL_API
  double
  Lsup() const;

  /**
     \brief
     Return the minimum absolute value of this Data object.
     *
  */
  ESCRIPT_DLL_API
  double
  Linf() const;

  /**
     \brief
     Return the maximum value of this Data object.
     *
  */
  ESCRIPT_DLL_API
  double
  sup() const;

  /**
     \brief
     Return the minimum value of this Data object.
     *
  */
  ESCRIPT_DLL_API
  double
  inf() const;

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
  unaryOp(UnaryFunction operation);

  /**
     \brief
     Return a Data object containing the specified slice of
     this Data object.
     \param region - Input - Region to copy.
     *
  */
  ESCRIPT_DLL_API
  Data
  getSlice(const DataArrayView::RegionType& region) const;

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
           const DataArrayView::RegionType& region);

  /**
     \brief
     Archive the current Data object to the given file.
     \param fileName - Input - file to archive to.
  */
  ESCRIPT_DLL_API
  void
  archiveData(const std::string fileName);

  /**
     \brief
     Extract the Data object archived in the given file, overwriting
     the current Data object.
     Note - the current object must be of type DataEmpty.
     \param fileName - Input - file to extract from.
     \param fspace - Input - a suitable FunctionSpace descibing the data.
  */
  ESCRIPT_DLL_API
  void
  extractData(const std::string fileName,
              const FunctionSpace& fspace);


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
  */
  ESCRIPT_DLL_API
	DataAbstract*
	borrowData(void) const;

 protected:

 private:

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
  template <class IValueType>
  void
  initialise(const IValueType& value,
             const FunctionSpace& what,
             bool expanded);

  //
  // flag to protect the data object against any update
  bool m_protected;

  //
  // pointer to the actual data object
  boost::shared_ptr<DataAbstract> m_data;

  //
  // pointer to the internal profiling data
  struct profDataEntry *profData;

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
    boost::shared_ptr<DataAbstract> temp_data(temp);
    m_data=temp_data;
  } else {
    DataAbstract* temp=new DataConstant(value,what);
    boost::shared_ptr<DataAbstract> temp_data(temp);
    m_data=temp_data;
  }
}

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
  Return true if operands are equivalent, else return false.
  NB: this operator does very little at this point, and isn't to 
  be relied on. Requires further implementation.
*/
//ESCRIPT_DLL_API bool operator==(const Data& left, const Data& right);

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
   if (getPointDataView().getRank()==0 && right.getPointDataView().getRank()!=0) {
     throw DataException("Error - attempt to update rank zero object with object with rank bigger than zero.");
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
   } else if (isConstant()) {
     DataConstant* leftC=dynamic_cast<DataConstant*>(m_data.get());
     DataConstant* rightC=dynamic_cast<DataConstant*>(tempRight.m_data.get());
     EsysAssert((leftC!=0 && rightC!=0), "Programming error - casting to DataConstant.");
     escript::binaryOp(*leftC,*rightC,operation);
   }
   #if defined DOPROF
   profData->binary++;
   #endif
}

/**
  \brief
  Perform the given unary operation on other and return the result.
  Given operation is performed on each element of each data point, thus
  argument object is a rank n Data object, and returned object is a rank n
  Data object.
  Calls Data::unaryOp.
*/
template <class UnaryFunction>
inline
Data
unaryOp(const Data& other,
        UnaryFunction operation)
{
  Data result;
  result.copy(other);
  result.unaryOp(operation);
  return result;
}

/**
  \brief
  Perform the given unary operation on this.
  Given operation is performed on each element of each data point.
  Calls escript::unaryOp.
*/
template <class UnaryFunction>
inline
void
Data::unaryOp(UnaryFunction operation)
{
  if (isExpanded()) {
    DataExpanded* leftC=dynamic_cast<DataExpanded*>(m_data.get());
    EsysAssert((leftC!=0), "Programming error - casting to DataExpanded.");
    escript::unaryOp(*leftC,operation);
  } else if (isTagged()) {
    DataTagged* leftC=dynamic_cast<DataTagged*>(m_data.get());
    EsysAssert((leftC!=0), "Programming error - casting to DataTagged.");
    escript::unaryOp(*leftC,operation);
  } else if (isConstant()) {
    DataConstant* leftC=dynamic_cast<DataConstant*>(m_data.get());
    EsysAssert((leftC!=0), "Programming error - casting to DataConstant.");
    escript::unaryOp(*leftC,operation);
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
  }
  return 0;
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
  if (isExpanded()) {
    Data result(0,DataArrayView::ShapeType(),getFunctionSpace(),isExpanded());
    DataExpanded* dataE=dynamic_cast<DataExpanded*>(m_data.get());
    DataExpanded* resultE=dynamic_cast<DataExpanded*>(result.m_data.get());
    EsysAssert((dataE!=0), "Programming error - casting data to DataExpanded.");
    EsysAssert((resultE!=0), "Programming error - casting result to DataExpanded.");
    escript::dp_algorithm(*dataE,*resultE,operation,initial_value);
    return result;
  } else if (isTagged()) {
    DataTagged* dataT=dynamic_cast<DataTagged*>(m_data.get());
    DataArrayView::ShapeType viewShape;
    DataArrayView::ValueType viewData(1);
    viewData[0]=0;
    DataArrayView defaultValue(viewData,viewShape);
    DataTagged::TagListType keys;
    DataTagged::ValueListType values;
    DataTagged::DataMapType::const_iterator i;
    for (i=dataT->getTagLookup().begin();i!=dataT->getTagLookup().end();i++) {
      keys.push_back(i->first);
      values.push_back(defaultValue);
    }
    Data result(keys,values,defaultValue,getFunctionSpace());
    DataTagged* resultT=dynamic_cast<DataTagged*>(result.m_data.get());
    EsysAssert((dataT!=0), "Programming error - casting data to DataTagged.");
    EsysAssert((resultT!=0), "Programming error - casting result to DataTagged.");
    escript::dp_algorithm(*dataT,*resultT,operation,initial_value);
    return result;
  } else if (isConstant()) {
    Data result(0,DataArrayView::ShapeType(),getFunctionSpace(),isExpanded());
    DataConstant* dataC=dynamic_cast<DataConstant*>(m_data.get());
    DataConstant* resultC=dynamic_cast<DataConstant*>(result.m_data.get());
    EsysAssert((dataC!=0), "Programming error - casting data to DataConstant.");
    EsysAssert((resultC!=0), "Programming error - casting result to DataConstant.");
    escript::dp_algorithm(*dataC,*resultC,operation,initial_value);
    return result;
  }
  Data falseRetVal; // to keep compiler quiet
  return falseRetVal;
}

}
#endif
