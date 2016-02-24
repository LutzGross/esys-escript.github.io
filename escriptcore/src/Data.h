
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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


/** \file Data.h */

#ifndef DATA_H
#define DATA_H
#include "system_dep.h"

#include "BinaryOp.h"
#include "DataAbstract.h"
#include "DataAlgorithm.h"
#include "DataException.h"
#include "DataTypes.h"
#include "FunctionSpace.h"
#include "UnaryOp.h"

#include "esysUtils/Esys_MPI.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <sstream>
#include <string>

#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/math/special_functions/bessel.hpp>

namespace escript {

//
// Forward declaration for various implementations of Data.
class DataConstant;
class DataTagged;
class DataExpanded;
class DataLazy;

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
  Data(const DataTypes::RealVectorType& value,
		 const DataTypes::ShapeType& shape,
                 const FunctionSpace& what=FunctionSpace(),
                 bool expanded=false);

  /**
     \brief
     Constructor which creates a Data with points having the specified shape.

     \param value - Input - Single value applied to all Data.
     \param dataPointShape - Input - The shape of each data point.
     \param what - Input - A description of what this data represents.
     \param expanded - Input - Flag, if true fill the entire container with
                       the given value. Otherwise a more efficient storage
                       mechanism will be used.
  */
  ESCRIPT_DLL_API
  Data(DataTypes::real_t value,
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
     Constructor which copies data from any object that can be treated like a python array/sequence.

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
     Constructor which copies data from a wrapped array.

     \param w - Input - Input data.
     \param what - Input - A description of what this data represents.
     \param expanded - Input - Flag, if true fill the entire container with
                       the value. Otherwise a more efficient storage
                       mechanism will be used.
  */       
  ESCRIPT_DLL_API     
  Data(const WrappedArray& w, const FunctionSpace& what,
           bool expanded=false);       
       

  /**
     \brief
     Constructor which creates a DataConstant.
     Copies data from any object that can be treated like a python array/sequence.
     All other parameters are copied from other.

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
  Data(DataTypes::real_t value,
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
  Data
  copySelf() const;


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
     Returns true, if the data object is protected against update

  */
  ESCRIPT_DLL_API
  bool
  isProtected() const;


/**
   \brief 
   Return the value of a data point as a python tuple.
*/
  ESCRIPT_DLL_API
  const boost::python::object
  getValueOfDataPointAsTuple(int dataPointNo);

  /**
     \brief
     sets the values of a data-point from a python object on this process
  */
  ESCRIPT_DLL_API
  void
  setValueOfDataPointToPyObject(int dataPointNo, const boost::python::object& py_object);

  /**
     \brief
     sets the values of a data-point from a array-like object on this process
  */
  ESCRIPT_DLL_API
  void
  setValueOfDataPointToArray(int dataPointNo, const boost::python::object&);

  /**
     \brief
     sets the values of a data-point on this process
  */
  ESCRIPT_DLL_API
  void
  setValueOfDataPoint(int dataPointNo, const DataTypes::real_t);
  
  ESCRIPT_DLL_API
  void
  setValueOfDataPointC(int dataPointNo, const DataTypes::cplx_t);  
  

  /**
     \brief Return a data point across all processors as a python tuple.
  */
  ESCRIPT_DLL_API
  const boost::python::object
  getValueOfGlobalDataPointAsTuple(int procNo, int dataPointNo);

  
  /**
     \brief Set the value of a global data point
  */
  ESCRIPT_DLL_API
  void
  setTupleForGlobalDataPoint(int id, int proc, boost::python::object);
  
  /**
     \brief
     Return the tag number associated with the given data-point.

  */
  ESCRIPT_DLL_API
  int
  getTagNumber(int dpno);


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
    \brief returns return true if data contains NaN.
  \warning This is dependent on the ability to reliably detect NaNs on your compiler.
   See the nancheck function in LocalOps for details.
  */
  ESCRIPT_DLL_API
  bool
  hasNaN();

  /**
  \brief replaces all NaN values with value 
  */
  ESCRIPT_DLL_API
  void
  replaceNaN(DataTypes::real_t value);

  /**
   \brief Ensures data is ready for write access.
  This means that the data will be resolved if lazy and will be copied if shared with another Data object.
  \warning This method should only be called in single threaded sections of code. (It modifies m_data).
  Do not create any Data objects from this one between calling requireWrite and getSampleDataRW.
  Doing so might introduce additional sharing.
  */
  ESCRIPT_DLL_API
  void
  requireWrite();

  /**
     \brief
     Return true if this Data is expanded.
     \note To determine if a sample will contain separate values for each datapoint. Use actsExpanded instead.
  */
  ESCRIPT_DLL_API
  bool
  isExpanded() const;

  /**
     \brief
     Return true if this Data is expanded or resolves to expanded.
     That is, if it has a separate value for each datapoint in the sample.
  */
  ESCRIPT_DLL_API
  bool
  actsExpanded() const;
  

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
    True if components of this data are stored as complex
  */
  ESCRIPT_DLL_API
  bool
  isComplex() const;

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
     Return the rank of the point data.
  */
  ESCRIPT_DLL_API
  inline
  unsigned int
  getDataPointRank() const
  {
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
     Returns true if the number of data points per sample and the number of
     samples match the respective argument. DataEmpty always returns true.
  */
  ESCRIPT_DLL_API
  inline
  bool numSamplesEqual(int numDataPointsPerSample, int numSamples) const
  {
    return (isEmpty() ||
            (numDataPointsPerSample==getNumDataPointsPerSample() && numSamples==getNumSamples()));
  }

  /**
     \brief
     Returns true if the shape matches the vector (dimensions[0],...,
     dimensions[rank-1]). DataEmpty always returns true.
  */
  inline
  bool isDataPointShapeEqual(int rank, const int* dimensions) const
  {
    if (isEmpty())
      return true;
    const DataTypes::ShapeType givenShape(&dimensions[0],&dimensions[rank]);
    return (getDataPointShape()==givenShape);
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
  \brief returns the values of the object as a list of tuples (one for each datapoint).

  \param scalarastuple If true, scalar data will produce single valued tuples [(1,) (2,) ...]
If false, the result is a list of scalars [1, 2, ...]
 */
  ESCRIPT_DLL_API
  const boost::python::object
  toListOfTuples(bool scalarastuple=true);


 /**
    \brief
    Return the sample data for the given sample no.
    Please do not use this unless you NEED to access samples individually
    \param sampleNo - Input - the given sample no.
    \return pointer to the sample data.
*/
  ESCRIPT_DLL_API
  inline
  const DataTypes::real_t*
  getSampleDataRO(DataTypes::RealVectorType::size_type sampleNo) const;


  /**
     \brief
     Return the sample data for the given sample no.
    Please do not use this unless you NEED to access samples individually
     \param sampleNo - Input - the given sample no.
     \return pointer to the sample data.
  */
  ESCRIPT_DLL_API
  inline
  DataTypes::real_t*
  getSampleDataRW(DataTypes::RealVectorType::size_type sampleNo);


 /**
    \brief
    Return a pointer to the beginning of the underlying data
    \warning please avoid using this method since it by-passes possible lazy improvements. May be removed without notice.
    \return pointer to the data.
*/
  ESCRIPT_DLL_API
  const DataTypes::real_t*
  getDataRO(DataTypes::real_t dummy=0) const;  
  
  ESCRIPT_DLL_API
  const DataTypes::cplx_t*
  getDataRO(DataTypes::cplx_t dummy) const;    
  
  
  
  /**
     \brief
     Return the sample data for the given tag. If an attempt is made to
     access data that isn't tagged an exception will be thrown.
     \param tag - Input - the tag key.
  */
  ESCRIPT_DLL_API
  inline
  DataTypes::real_t*
  getSampleDataByTag(int tag)
  {
    return m_data->getSampleDataByTag(tag);
  }

  /**
     \brief
     Return a reference into the DataVector which points to the specified data point.
     \param sampleNo - Input -
     \param dataPointNo - Input -
  */
  ESCRIPT_DLL_API
  DataTypes::RealVectorType::const_reference
  getDataPointRO(int sampleNo, int dataPointNo);

  /**
     \brief
     Return a reference into the DataVector which points to the specified data point.
     \param sampleNo - Input -
     \param dataPointNo - Input -
  */
  ESCRIPT_DLL_API
  DataTypes::RealVectorType::reference
  getDataPointRW(int sampleNo, int dataPointNo);



  /**
     \brief 
     Return the offset for the given sample and point within the sample
  */
  ESCRIPT_DLL_API
  inline
  DataTypes::RealVectorType::size_type
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
  DataTypes::RealVectorType::size_type
  getLength() const;

  /**
  \brief Return true if this object contains no samples.
  This is not the same as isEmpty() 
  */
  ESCRIPT_DLL_API
  bool
  hasNoSamples() const
  {
	return getLength()==0;
  }

  /**
     \brief
     Assign the given value to the tag assocciated with name. Implicitly converts this
     object to type DataTagged. Throws an exception if this object
     cannot be converted to a DataTagged object or name cannot be mapped onto a tag key.
     \param name - Input - name of tag.
     \param value - Input - Value to associate with given key.
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
                        const DataTypes::RealVectorType& value,
			int dataOffset=0);


  ESCRIPT_DLL_API
  void
  setTaggedValueFromCPP(int tagKey,
			const DataTypes::ShapeType& pointshape,
                        const DataTypes::CplxVectorType& value,
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

  ESCRIPT_DLL_API
  Data
  interpolateFromTable3D(const WrappedArray& table, DataTypes::real_t Amin, DataTypes::real_t Astep,
                       DataTypes::real_t undef, Data& B, DataTypes::real_t Bmin, DataTypes::real_t Bstep, Data& C, 
			DataTypes::real_t Cmin, DataTypes::real_t Cstep, bool check_boundaries);

  ESCRIPT_DLL_API
  Data
  interpolateFromTable2D(const WrappedArray& table, DataTypes::real_t Amin, DataTypes::real_t Astep,
                       DataTypes::real_t undef, Data& B, DataTypes::real_t Bmin, DataTypes::real_t Bstep,bool check_boundaries);

  ESCRIPT_DLL_API
  Data
  interpolateFromTable1D(const WrappedArray& table, DataTypes::real_t Amin, DataTypes::real_t Astep,
                       DataTypes::real_t undef,bool check_boundaries);


  ESCRIPT_DLL_API
  Data
  interpolateFromTable3DP(boost::python::object table, DataTypes::real_t Amin, DataTypes::real_t Astep,
                        Data& B, DataTypes::real_t Bmin, DataTypes::real_t Bstep, Data& C, DataTypes::real_t Cmin, DataTypes::real_t Cstep, DataTypes::real_t undef,bool check_boundaries);


  ESCRIPT_DLL_API
  Data
  interpolateFromTable2DP(boost::python::object table, DataTypes::real_t Amin, DataTypes::real_t Astep,
                        Data& B, DataTypes::real_t Bmin, DataTypes::real_t Bstep, DataTypes::real_t undef,bool check_boundaries);

  ESCRIPT_DLL_API
  Data
  interpolateFromTable1DP(boost::python::object table, DataTypes::real_t Amin, DataTypes::real_t Astep,
                        DataTypes::real_t undef,bool check_boundaries);
  
  ESCRIPT_DLL_API
  Data
  nonuniforminterp(boost::python::object in, boost::python::object out, bool check_boundaries);

  ESCRIPT_DLL_API
  Data
  nonuniformslope(boost::python::object in, boost::python::object out, bool check_boundaries);  
  
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
     Calculate the integral over the function space domain as a python tuple.
  */
  ESCRIPT_DLL_API
  boost::python::object
  integrateToTuple_const() const;


  /**
    \brief
     Calculate the integral over the function space domain as a python tuple.
  */
  ESCRIPT_DLL_API
  boost::python::object
  integrateToTuple();



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
  whereZero(DataTypes::real_t tol=0.0) const;

  /**
     \brief
     Return a Data with a 0 for 0 values and a 1 for +ive or -ive values.
     *
  */
  ESCRIPT_DLL_API
  Data
  whereNonZero(DataTypes::real_t tol=0.0) const;

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
  DataTypes::real_t
  Lsup();

  ESCRIPT_DLL_API
  DataTypes::real_t
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
  DataTypes::real_t
  sup();

  ESCRIPT_DLL_API
  DataTypes::real_t
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
  DataTypes::real_t
  inf();

  ESCRIPT_DLL_API
  DataTypes::real_t
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
     the minimum component value in this Data object.
     \note If you are working in python, please consider using Locator 
instead of manually manipulating process and point IDs.
  */
  ESCRIPT_DLL_API
  const boost::python::tuple
  minGlobalDataPoint() const;

  /**
     \brief
     Return the (sample number, data-point number) of the data point with
     the minimum component value in this Data object.
     \note If you are working in python, please consider using Locator 
instead of manually manipulating process and point IDs.
  */
  ESCRIPT_DLL_API
  const boost::python::tuple
  maxGlobalDataPoint() const;



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
  eigenvalues_and_eigenvectors(const DataTypes::real_t tol=1.e-12) const;

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
     For complex values return the conjugate values.
     For non-complex data return a copy
  */
  ESCRIPT_DLL_API
  Data
  conjugate() const;
  
  ESCRIPT_DLL_API
  Data
  real() const;  
  
  ESCRIPT_DLL_API
  Data
  imag() const;  

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
     Bessel worker function.
     *
  */
  ESCRIPT_DLL_API
  Data
  bessel(int order, DataTypes::real_t (*besselfunc) (int,DataTypes::real_t) );
  

  /**
     \brief
     Return the Bessel function of the first kind for each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  besselFirstKind(int order);

  /**
     \brief
     Return the Bessel function of the second kind for each data point of this Data object.
     *
  */
  ESCRIPT_DLL_API
  Data
  besselSecondKind(int order);


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
    Newer style division operator for python
  */
  ESCRIPT_DLL_API
  Data truedivD(const Data& right);

  /**
    \brief
    Newer style division operator for python
  */
  ESCRIPT_DLL_API
  Data truedivO(const boost::python::object& right);

  /**
    \brief
    Newer style division operator for python
  */
  ESCRIPT_DLL_API
  Data rtruedivO(const boost::python::object& left);

  /**
    \brief
    wrapper for python add operation
  */
  ESCRIPT_DLL_API
  boost::python::object __add__(const boost::python::object& right);
  

  /**
    \brief
    wrapper for python subtract operation
  */
  ESCRIPT_DLL_API
  boost::python::object __sub__(const boost::python::object& right);
  
  /**
    \brief
    wrapper for python reverse subtract operation
  */
  ESCRIPT_DLL_API
  boost::python::object __rsub__(const boost::python::object& right);  

  /**
    \brief
    wrapper for python multiply operation
  */
  ESCRIPT_DLL_API
  boost::python::object __mul__(const boost::python::object& right);
    
  /**
    \brief
    wrapper for python divide operation
  */
  ESCRIPT_DLL_API
  boost::python::object __div__(const boost::python::object& right);
  
  /**
    \brief
    wrapper for python reverse divide operation
  */
  ESCRIPT_DLL_API
  boost::python::object __rdiv__(const boost::python::object& right);    
  
  /**
	\brief return inverse of matricies.
  */
  ESCRIPT_DLL_API
  Data
  matrixInverse() const;

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
        DataTypes::RealVectorType::const_reference
        getDataAtOffsetRO(DataTypes::RealVectorType::size_type i, DataTypes::real_t dummy=0);


  ESCRIPT_DLL_API
        DataTypes::RealVectorType::reference
        getDataAtOffsetRW(DataTypes::RealVectorType::size_type i, DataTypes::real_t dummy=0);
	
  ESCRIPT_DLL_API
        DataTypes::CplxVectorType::const_reference
        getDataAtOffsetRO(DataTypes::CplxVectorType::size_type i, DataTypes::cplx_t dummy);


  ESCRIPT_DLL_API
        DataTypes::CplxVectorType::reference
        getDataAtOffsetRW(DataTypes::CplxVectorType::size_type i, DataTypes::cplx_t dummy);	
	

  /**
    \brief Ensures that the Data is expanded and returns its underlying vector
    Does not check for exclusive write so do that before calling if sharing
    Is a posibility.
    \warning For domain implementors only. Using this function will
    avoid using optimisations like lazy evaluation. It is intended
    to allow quick initialisation of Data by domain; not as a bypass around 
    escript's other mechanisms.
  */
  ESCRIPT_DLL_API
  DataTypes::RealVectorType&
  getExpandedVectorReference(DataTypes::real_t dummy=0);

  ESCRIPT_DLL_API
  DataTypes::CplxVectorType&
  getExpandedVectorReference(DataTypes::cplx_t dummy);
  
  
  /**
   * \brief For tagged Data returns the number of tags with values.
   * For non-tagged data will return 0 (even Data which has been expanded from tagged).
  */ 
  ESCRIPT_DLL_API
  size_t
  getNumberOfTaggedValues() const;

  /*
  * \brief make the data complex
  */
  ESCRIPT_DLL_API
  void complicate();
 
 protected:

 private:

template <class BinaryOp>
  DataTypes::real_t 
#ifdef ESYS_MPI
  lazyAlgWorker(DataTypes::real_t init, MPI_Op mpiop_type);
#else
  lazyAlgWorker(DataTypes::real_t init);
#endif

  DataTypes::real_t
  LsupWorker() const;

  DataTypes::real_t
  supWorker() const;

  DataTypes::real_t
  infWorker() const;

  boost::python::object
  integrateWorker() const;

  void
  calc_minGlobalDataPoint(int& ProcNo,  int& DataPointNo) const;

  void
  calc_maxGlobalDataPoint(int& ProcNo,  int& DataPointNo) const;

  // For internal use in Data.cpp only!
  // other uses should call the main entry points and allow laziness
  Data
  minval_nonlazy() const;

  // For internal use in Data.cpp only!
  Data
  maxval_nonlazy() const;


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
  DataTypes::real_t
  algorithm(BinaryFunction operation,
            DataTypes::real_t initial_value) const;

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
               DataTypes::real_t initial_value) const;

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
  
  void
  binaryDataOp(const Data& right,
		   escript::ESFunction operation); 
  

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
  initialise(const DataTypes::RealVectorType& value,
	     const DataTypes::ShapeType& shape,
             const FunctionSpace& what,
             bool expanded);

  void
  initialise(const WrappedArray& value,
                 const FunctionSpace& what,
                 bool expanded);

  void
  initialise(const DataTypes::real_t value,
	     const DataTypes::ShapeType& shape,
             const FunctionSpace& what,
             bool expanded);

  //
  // flag to protect the data object against any update
  bool m_protected;
  mutable bool m_shared;
  bool m_lazy;

  //
  // pointer to the actual data object
//   boost::shared_ptr<DataAbstract> m_data;
  DataAbstract_ptr m_data;

// If possible please use getReadyPtr instead.
// But see warning below.
  const DataReady*
  getReady() const;

  DataReady*
  getReady();


// Be wary of using this for local operations since it (temporarily) increases reference count.
// If you are just using this to call a method on DataReady instead of DataAbstract consider using 
// getReady() instead
  DataReady_ptr
  getReadyPtr();

  const_DataReady_ptr
  getReadyPtr() const;


  /**
   \brief Update the Data's shared flag
   This indicates that the DataAbstract used by this object is now shared (or no longer shared).
   For internal use only.
  */
  void updateShareStatus(bool nowshared) const
  {
	m_shared=nowshared;		// m_shared is mutable
  }

  // In the isShared() method below:
  // A problem would occur if m_data (the address pointed to) were being modified 
  // while the call m_data->is_shared is being executed.
  // 
  // Q: So why do I think this code can be thread safe/correct?
  // A: We need to make some assumptions.
  // 1. We assume it is acceptable to return true under some conditions when we aren't shared.
  // 2. We assume that no constructions or assignments which will share previously unshared
  //    will occur while this call is executing. This is consistent with the way Data:: and C are written.
  //
  // This means that the only transition we need to consider, is when a previously shared object is
  // not shared anymore. ie. the other objects have been destroyed or a deep copy has been made.
  // In those cases the m_shared flag changes to false after m_data has completed changing.
  // For any threads executing before the flag switches they will assume the object is still shared.
  bool isShared() const
  {
	return m_shared;
/*	if (m_shared) return true;
	if (m_data->isShared())			
	{					
		updateShareStatus(true);
		return true;
	}
	return false;*/
  }

  void forceResolve()
  {
	if (isLazy())
	{
	    #ifdef _OPENMP
	    if (omp_in_parallel())
	    {	// Yes this is throwing an exception out of an omp thread which is forbidden.
		throw DataException("Please do not call forceResolve() in a parallel region.");
	    }
	    #endif
	    resolve();
	}
  }

  /**
  \brief if another object is sharing out member data make a copy to work with instead. 
  This code should only be called from single threaded sections of code.
  */
  void exclusiveWrite()
  {
#ifdef _OPENMP
  if (omp_in_parallel())
  {
	throw DataException("Programming error. Please do not run exclusiveWrite() in multi-threaded sections.");
  }
#endif
	forceResolve();
	if (isShared())
	{
		DataAbstract* t=m_data->deepCopy();
   		set_m_data(DataAbstract_ptr(t));
	}
#ifdef EXWRITECHK		
	m_data->exclusivewritecalled=true;
#endif	
  }

  /**
  \brief checks if caller can have exclusive write to the object
  */
  void checkExclusiveWrite()
  {
	if  (isLazy() || isShared())
	{
		throw DataException("Programming error. ExclusiveWrite required - please call requireWrite()");
	}
  }

  /**
  \brief Modify the data abstract hosted by this Data object
  For internal use only.
  Passing a pointer to null is permitted (do this in the destructor)
  \warning Only to be called in single threaded code or inside a single/critical section. This method needs to be atomic.
  */
  void set_m_data(DataAbstract_ptr p);

  friend class DataAbstract;		// To allow calls to updateShareStatus
  friend class TestDomain;		// so its getX will work quickly
#ifdef IKNOWWHATIMDOING
  friend ESCRIPT_DLL_API Data applyBinaryCFunction(boost::python::object cfunc, boost::python::tuple shape, escript::Data& d, escript::Data& e);
#endif
  friend ESCRIPT_DLL_API Data condEval(escript::Data& mask, escript::Data& trueval, escript::Data& falseval);
  friend ESCRIPT_DLL_API Data randomData(const boost::python::tuple& shape, const FunctionSpace& what, long seed, const boost::python::tuple& filter);

};


#ifdef IKNOWWHATIMDOING
ESCRIPT_DLL_API
Data
applyBinaryCFunction(boost::python::object func, boost::python::tuple shape, escript::Data& d, escript::Data& e);
#endif


ESCRIPT_DLL_API
Data
condEval(escript::Data& mask, escript::Data& trueval, escript::Data& falseval);



/**
 \brief Create a new Expanded Data object filled with pseudo-random data.
*/
ESCRIPT_DLL_API
Data randomData(const boost::python::tuple& shape,
       const FunctionSpace& what,
       long seed, const boost::python::tuple& filter);


}   // end namespace escript


// No, this is not supposed to be at the top of the file
// DataAbstact needs to be declared first, then DataReady needs to be fully declared
// so that I can dynamic cast between them below.
#include "DataReady.h"
#include "DataLazy.h"

namespace escript
{

inline
const DataReady*
Data::getReady() const
{
   const DataReady* dr=dynamic_cast<const DataReady*>(m_data.get());
   ESYS_ASSERT(dr!=0, "error casting to DataReady.");
   return dr;
}

inline
DataReady*
Data::getReady()
{
   DataReady* dr=dynamic_cast<DataReady*>(m_data.get());
   ESYS_ASSERT(dr!=0, "error casting to DataReady.");
   return dr;
}

// Be wary of using this for local operations since it (temporarily) increases reference count.
// If you are just using this to call a method on DataReady instead of DataAbstract consider using 
// getReady() instead
inline
DataReady_ptr
Data::getReadyPtr()
{
   DataReady_ptr dr=boost::dynamic_pointer_cast<DataReady>(m_data);
   ESYS_ASSERT(dr.get()!=0, "error casting to DataReady.");
   return dr;
}


inline
const_DataReady_ptr
Data::getReadyPtr() const
{
   const_DataReady_ptr dr=boost::dynamic_pointer_cast<const DataReady>(m_data);
   ESYS_ASSERT(dr.get()!=0, "error casting to DataReady.");
   return dr;
}

inline
DataTypes::real_t*
Data::getSampleDataRW(DataTypes::RealVectorType::size_type sampleNo)
{
   if (isLazy())
   {
	throw DataException("Error, attempt to acquire RW access to lazy data. Please call requireWrite() first.");
   }
#ifdef EXWRITECHK
   if (!getReady()->exclusivewritecalled)
   {
        throw DataException("Error, call to Data::getSampleDataRW without a preceeding call to requireWrite/exclusiveWrite.");
   }
#endif
   return getReady()->getSampleDataRW(sampleNo);
}

inline
const DataTypes::real_t*
Data::getSampleDataRO(DataTypes::RealVectorType::size_type sampleNo) const
{
   DataLazy* l=dynamic_cast<DataLazy*>(m_data.get());
   if (l!=0)
   {
	size_t offset=0;
	const DataTypes::RealVectorType* res=l->resolveSample(sampleNo,offset);
	return &((*res)[offset]);
   }
   return getReady()->getSampleDataRO(sampleNo);
}

inline
const DataTypes::real_t*
Data::getDataRO(DataTypes::real_t dummy) const
{
    if (isLazy())
    {
        throw DataException("Programmer error - getDataRO must not be called on Lazy Data.");
    }
    if (getNumSamples()==0)
    {
	return 0;
    }
    else
    {
	return &(getReady()->getTypedVectorRO(0)[0]);
    }
}

inline
const DataTypes::cplx_t*
Data::getDataRO(DataTypes::cplx_t dummy) const
{
    if (isLazy())
    {
        throw DataException("Programmer error - getDataRO must not be called on Lazy Data.");
    }
    if (getNumSamples()==0)
    {
	return 0;
    }
    else
    {
	return &(getReady()->getTypedVectorRO(dummy)[0]);
    }
}


/**
   Binary Data object operators.
*/
inline DataTypes::real_t rpow(DataTypes::real_t x,DataTypes::real_t y)
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
  \param arg_0 - Input - Data object
  \param arg_1 - Input - Data object
  \param axis_offset - Input - axis offset
  \param transpose - Input - 0: transpose neither, 1: transpose arg0, 2: transpose arg1
*/
ESCRIPT_DLL_API
Data
C_GeneralTensorProduct(Data& arg_0,
                     Data& arg_1,
                     int axis_offset=0,
                     int transpose=0);

/**
  \brief
  Operator/
  Takes RHS Data object.
*/
inline
Data
Data::truedivD(const Data& right)
{
    return *this / right;
}

/**
  \brief
  Operator/
  Takes RHS python::object.
*/
inline
Data
Data::truedivO(const boost::python::object& right)
{
    Data tmp(right, getFunctionSpace(), false);
    return truedivD(tmp);
}

/**
  \brief
  Operator/
  Takes LHS python::object.
*/
inline
Data
Data::rtruedivO(const boost::python::object& left)
{
    Data tmp(left, getFunctionSpace(), false);
    return tmp.truedivD(*this);
}


inline 
void
Data::binaryDataOp(const Data& right,
		   escript::ESFunction operation)
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
   FunctionSpace fsl=getFunctionSpace();
   FunctionSpace fsr=right.getFunctionSpace();
   if (fsl!=fsr) {
     signed char intres=fsl.getDomain()->preferredInterpolationOnDomain(fsr.getTypeCode(), fsl.getTypeCode());
     if (intres==0)
     {
         std::string msg="Error - attempt to combine incompatible FunctionSpaces.";
	 msg+=fsl.toString();
	 msg+="  ";
	 msg+=fsr.toString();
         throw DataException(msg.c_str());
     } 
     else if (intres==1)
     {
       // an interpolation is required so create a new Data
       tempRight=Data(right,fsl);
     }
     else	// reverse interpolation preferred
     {
        // interpolate onto the RHS function space
       Data tempLeft(*this,fsr);
       set_m_data(tempLeft.m_data);
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
     ESYS_ASSERT(leftC!=0, "Programming error - casting to DataExpanded.");
     escript::binaryOpDataReady(*leftC,*(tempRight.getReady()),operation);
   } else if (isTagged()) {
     //
     // Tagged data is operated on serially, the right hand side can be
     // either DataConstant or DataTagged
     DataTagged* leftC=dynamic_cast<DataTagged*>(m_data.get());
     ESYS_ASSERT(leftC!=0, "Programming error - casting to DataTagged.");
     if (right.isTagged()) {
       DataTagged* rightC=dynamic_cast<DataTagged*>(tempRight.m_data.get());
       ESYS_ASSERT(rightC!=0, "Programming error - casting to DataTagged.");
       escript::binaryOpDataReady(*leftC,*rightC,operation);
     } else {
       DataConstant* rightC=dynamic_cast<DataConstant*>(tempRight.m_data.get());
       ESYS_ASSERT(rightC!=0, "Programming error - casting to DataConstant.");
       escript::binaryOpDataReady(*leftC,*rightC,operation);
     }
   } else if (isConstant()) {
     DataConstant* leftC=dynamic_cast<DataConstant*>(m_data.get());
     DataConstant* rightC=dynamic_cast<DataConstant*>(tempRight.m_data.get());
     ESYS_ASSERT(leftC!=0 && rightC!=0, "Programming error - casting to DataConstant.");
     escript::binaryOpDataReady(*leftC,*rightC,operation);
   }  
}

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
   //Data tempRight(right);				// unfortunately, this adds an owner and breaks d+=d
   
   DataAbstract_ptr rp=right.m_data;
   
   FunctionSpace fsl=getFunctionSpace();
   FunctionSpace fsr=right.getFunctionSpace();
   if (fsl!=fsr) {
     signed char intres=fsl.getDomain()->preferredInterpolationOnDomain(fsr.getTypeCode(), fsl.getTypeCode());
     if (intres==0)
     {
         std::string msg="Error - attempt to combine incompatible FunctionSpaces.";
	 msg+=fsl.toString();
	 msg+="  ";
	 msg+=fsr.toString();
         throw DataException(msg.c_str());
     } 
     else if (intres==1)
     {
       // an interpolation is required so create a new Data
       rp=Data(right,fsl).m_data;
     }
     else	// reverse interpolation preferred
     {
        // interpolate onto the RHS function space
       Data tempLeft(*this,fsr);
       set_m_data(tempLeft.m_data);
     }
   }
   m_data->operandCheck(*rp.get());
   //operandCheck(tempRight);
   //
   // ensure this has the right type for the RHS
   typeMatchRight(right);		// yes we are actually processing rp 
			      // (which may be pointing at a different object) not right.
			      // but rp and right will be referring to the same type of data
   //
   // Need to cast to the concrete types so that the correct binaryOp
   // is called.
   if (isExpanded()) {
     //
     // Expanded data will be done in parallel, the right hand side can be
     // of any data type
     DataExpanded* leftC=dynamic_cast<DataExpanded*>(m_data.get());
     ESYS_ASSERT(leftC!=0, "Programming error - casting to DataExpanded.");
     DataReady* dr=dynamic_cast<DataReady*>(rp.get());
     escript::binaryOp(*leftC,*dr,operation);
   } else if (isTagged()) {
     //
     // Tagged data is operated on serially, the right hand side can be
     // either DataConstant or DataTagged
     DataTagged* leftC=dynamic_cast<DataTagged*>(m_data.get());
     ESYS_ASSERT(leftC!=0, "Programming error - casting to DataTagged.");
     if (right.isTagged()) {
       DataTagged* rightC=dynamic_cast<DataTagged*>(rp.get());
       ESYS_ASSERT(rightC!=0, "Programming error - casting to DataTagged.");
       escript::binaryOp(*leftC,*rightC,operation);
     } else {
       DataConstant* rightC=dynamic_cast<DataConstant*>(rp.get());
       ESYS_ASSERT(rightC!=0, "Programming error - casting to DataConstant.");
       escript::binaryOp(*leftC,*rightC,operation);
     }
   } else if (isConstant()) {
     DataConstant* leftC=dynamic_cast<DataConstant*>(m_data.get());
     DataConstant* rightC=dynamic_cast<DataConstant*>(rp.get());
     ESYS_ASSERT(leftC!=0 && rightC!=0, "Programming error - casting to DataConstant.");
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
DataTypes::real_t
Data::algorithm(BinaryFunction operation, DataTypes::real_t initial_value) const
{
  if (isExpanded()) {
    DataExpanded* leftC=dynamic_cast<DataExpanded*>(m_data.get());
    ESYS_ASSERT(leftC!=0, "Programming error - casting to DataExpanded.");
    return escript::algorithm(*leftC,operation,initial_value);
  } else if (isTagged()) {
    DataTagged* leftC=dynamic_cast<DataTagged*>(m_data.get());
    ESYS_ASSERT(leftC!=0, "Programming error - casting to DataTagged.");
    return escript::algorithm(*leftC,operation,initial_value);
  } else if (isConstant()) {
    DataConstant* leftC=dynamic_cast<DataConstant*>(m_data.get());
    ESYS_ASSERT(leftC!=0, "Programming error - casting to DataConstant.");
    return escript::algorithm(*leftC,operation,initial_value);
  } else if (isEmpty()) {
    throw DataException("Error - Operations (algorithm) not permitted on instances of DataEmpty.");
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
Data::dp_algorithm(BinaryFunction operation, DataTypes::real_t initial_value) const
{
  if (isEmpty()) {
    throw DataException("Error - Operations (dp_algorithm) not permitted on instances of DataEmpty.");
  } 
  else if (isExpanded()) {
    Data result(0,DataTypes::ShapeType(),getFunctionSpace(),isExpanded());
    DataExpanded* dataE=dynamic_cast<DataExpanded*>(m_data.get());
    DataExpanded* resultE=dynamic_cast<DataExpanded*>(result.m_data.get());
    ESYS_ASSERT(dataE!=0, "Programming error - casting data to DataExpanded.");
    ESYS_ASSERT(resultE!=0, "Programming error - casting result to DataExpanded.");
    escript::dp_algorithm(*dataE,*resultE,operation,initial_value);
    return result;
  }
  else if (isTagged()) {
    DataTagged* dataT=dynamic_cast<DataTagged*>(m_data.get());
    ESYS_ASSERT(dataT!=0, "Programming error - casting data to DataTagged.");
    DataTypes::RealVectorType defval(1);
    defval[0]=0;
    DataTagged* resultT=new DataTagged(getFunctionSpace(), DataTypes::scalarShape, defval, dataT);
    escript::dp_algorithm(*dataT,*resultT,operation,initial_value);
    return Data(resultT);   // note: the Data object now owns the resultT pointer
  } 
  else if (isConstant()) {
    Data result(0,DataTypes::ShapeType(),getFunctionSpace(),isExpanded());
    DataConstant* dataC=dynamic_cast<DataConstant*>(m_data.get());
    DataConstant* resultC=dynamic_cast<DataConstant*>(result.m_data.get());
    ESYS_ASSERT(dataC!=0, "Programming error - casting data to DataConstant.");
    ESYS_ASSERT(resultC!=0, "Programming error - casting result to DataConstant.");
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
  \param arg_0 - Input - Data object
  \param arg_1 - Input - Data object
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
     throw DataException("Error - Operations (C_TensorBinaryOperation) not permitted on instances of DataEmpty.");
  }
  if (arg_0.isLazy() || arg_1.isLazy())
  {
     throw DataException("Error - Operations not permitted on lazy data.");
  }
  
  // for disambiguation of some methods later
  typename decltype(operation)::first_argument_type dummy0=0;
  typename decltype(operation)::second_argument_type dummy1=0;
  typename decltype(operation)::result_type dummyr=0;  
  
  // Interpolate if necessary and find an appropriate function space
  Data arg_0_Z, arg_1_Z;
  FunctionSpace fsl=arg_0.getFunctionSpace();
  FunctionSpace fsr=arg_1.getFunctionSpace();
  if (fsl!=fsr) {
     signed char intres=fsl.getDomain()->preferredInterpolationOnDomain(fsr.getTypeCode(), fsl.getTypeCode());
     if (intres==0)
     {
         std::string msg="Error - C_TensorBinaryOperation: arguments have incompatible function spaces.";
         msg+=fsl.toString();
         msg+=" ";
         msg+=fsr.toString();
         throw DataException(msg.c_str());
     } 
     else if (intres==1)
     {
      arg_1_Z=arg_1.interpolate(arg_0.getFunctionSpace());
      arg_0_Z =Data(arg_0);      
     }
     else	// reverse interpolation preferred
     {
      arg_0_Z = arg_0.interpolate(arg_1.getFunctionSpace());
      arg_1_Z = Data(arg_1);
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
      typename decltype(operation)::first_argument_type dummy=0;
      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace());      // DataConstant output
      if (arg_0_Z.isComplex() || arg_1_Z.isComplex())
      {
	  res.complicate();	// It would be much better to create the Data object as complex to start with
      }				// But that would require more work so let's just get this case working first
      const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(0, dummy));
      const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(0, dummy1));
      typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(0, dummyr));

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
      const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));

      // Get the pointers to the actual data
      const typename decltype(operation)::second_argument_type *ptr_1 = &(tmp_1->getDefaultValueRO(0, dummy1));
      typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyr));

      // Compute a result for the default
      tensor_binary_operation(size0, ptr_0, ptr_1, ptr_2, operation);
      // Compute a result for each tag
      const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_1.begin();i!=lookup_1.end();i++) {
        tmp_2->addTag(i->first);
        const typename decltype(operation)::second_argument_type *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0, dummy1));
        typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummyr));

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
      res.requireWrite();
      #pragma omp parallel for private(sampleNo_1,dataPointNo_1) schedule(static)
      for (sampleNo_1 = 0; sampleNo_1 < numSamples_1; sampleNo_1++) {
        for (dataPointNo_1 = 0; dataPointNo_1 < numDataPointsPerSample_1; dataPointNo_1++) {
          int offset_1 = tmp_1->getPointOffset(sampleNo_1,dataPointNo_1);
          int offset_2 = tmp_2->getPointOffset(sampleNo_1,dataPointNo_1);
          const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));
          const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1)); 
          typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr)); 
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

      const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
      // Get the pointers to the actual data
      const typename decltype(operation)::first_argument_type *ptr_0 = &(tmp_0->getDefaultValueRO(0, dummy0));
      typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyr));
      // Compute a result for the default
      tensor_binary_operation(size0, ptr_0, ptr_1, ptr_2, operation);
      // Compute a result for each tag
      const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_0.begin();i!=lookup_0.end();i++) {
        tmp_2->addTag(i->first);
        const typename decltype(operation)::first_argument_type *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0, dummy0));
        typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0,dummyr));
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

      // Get the pointers to the actual data
      const typename decltype(operation)::first_argument_type *ptr_0 = &(tmp_0->getDefaultValueRO(0, dummy0));
      const typename decltype(operation)::second_argument_type *ptr_1 = &(tmp_1->getDefaultValueRO(0, dummy1));
      typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyr));

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

        const typename decltype(operation)::first_argument_type *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0, dummy0));
        const typename decltype(operation)::second_argument_type *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0, dummy1));
        typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummyr));

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
      res.requireWrite();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,0); // They're all the same, so just use #0
        const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
          typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
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
      res.requireWrite();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);

          const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));
          const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
          typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));


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
      res.requireWrite();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        int offset_1 = tmp_1->getPointOffset(sampleNo_0,0);
        const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));
          typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
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
      res.requireWrite();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
	  dataPointNo_0=0;
//        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));
          const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
          typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
          tensor_binary_operation(size0*numDataPointsPerSample_0, ptr_0, ptr_1, ptr_2, operation);
//       }
      }

    }
    else {
      throw DataException("Error - C_TensorBinaryOperation: unknown combination of inputs");
    }

  } else if (0 == rank0) {
    if (arg_0_Z.isConstant()   && arg_1_Z.isConstant()) {
      res = Data(0.0, shape1, arg_1_Z.getFunctionSpace());      // DataConstant output
      const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(0, dummy0));
      const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(0, dummy1));
      typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(0, dummyr));
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
      const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));

      const typename decltype(operation)::second_argument_type *ptr_1 = &(tmp_1->getDefaultValueRO(0, dummy1));
      typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyr));

      // Compute a result for the default
      tensor_binary_operation(size1, ptr_0[0], ptr_1, ptr_2, operation);
      // Compute a result for each tag
      const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_1.begin();i!=lookup_1.end();i++) {
        tmp_2->addTag(i->first);
        const typename decltype(operation)::second_argument_type *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0, dummy1));
        typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummyr));
        tensor_binary_operation(size1, ptr_0[0], ptr_1, ptr_2, operation);
      }

    }
    else if (arg_0_Z.isConstant()   && arg_1_Z.isExpanded()) {

      res = Data(0.0, shape1, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataConstant* tmp_0=dynamic_cast<DataConstant*>(arg_0_Z.borrowData());
      DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_1;
      int numSamples_1 = arg_1_Z.getNumSamples();
      int numDataPointsPerSample_1 = arg_1_Z.getNumDataPointsPerSample();
      int offset_0 = tmp_0->getPointOffset(0,0);
      const typename decltype(operation)::first_argument_type *ptr_src = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));
      typename decltype(operation)::first_argument_type ptr_0 = ptr_src[0];
      int size = size1*numDataPointsPerSample_1;
      res.requireWrite();
      #pragma omp parallel for private(sampleNo_1) schedule(static)
      for (sampleNo_1 = 0; sampleNo_1 < numSamples_1; sampleNo_1++) {
          int offset_1 = tmp_1->getPointOffset(sampleNo_1,0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_1,0);
          const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
          typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
          tensor_binary_operation(size, ptr_0, ptr_1, ptr_2, operation);
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
      const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));

      // Get the pointers to the actual data
      const typename decltype(operation)::first_argument_type *ptr_0 = &(tmp_0->getDefaultValueRO(0, dummy0));
      typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyr));


      // Compute a result for the default
      tensor_binary_operation(size1, ptr_0[0], ptr_1, ptr_2, operation);
      // Compute a result for each tag
      const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_0.begin();i!=lookup_0.end();i++) {
        tmp_2->addTag(i->first);
        const typename decltype(operation)::first_argument_type *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0, dummy0));
        typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummyr));

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

      // Get the pointers to the actual data
      const typename decltype(operation)::first_argument_type *ptr_0 = &(tmp_0->getDefaultValueRO(0, dummy0));
      const typename decltype(operation)::second_argument_type *ptr_1 = &(tmp_1->getDefaultValueRO(0, dummy1));
      typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyr));

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
        const typename decltype(operation)::first_argument_type *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0, dummy0));
        const typename decltype(operation)::second_argument_type *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0, dummy1));
        typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummyr));

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
      res.requireWrite();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,0); // They're all the same, so just use #0
        const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
          typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
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
      res.requireWrite();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));
          const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
          typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
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
      res.requireWrite();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        int offset_1 = tmp_1->getPointOffset(sampleNo_0,0);
        const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));
          typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
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
      res.requireWrite();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));
          const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
          typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
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
      const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(0, dummy0));
      const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(0, dummy1));
      typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(0, dummyr));
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
      const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));

      //Get the pointers to the actual data
      const typename decltype(operation)::second_argument_type *ptr_1 = &(tmp_1->getDefaultValueRO(0, dummy1));
      typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyr));

      // Compute a result for the default
      tensor_binary_operation(size0, ptr_0, ptr_1[0], ptr_2, operation);
      // Compute a result for each tag
      const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_1.begin();i!=lookup_1.end();i++) {
        tmp_2->addTag(i->first);
        const typename decltype(operation)::second_argument_type *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0, dummy1));
        typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummyr));
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
      res.requireWrite();
      #pragma omp parallel for private(sampleNo_1,dataPointNo_1) schedule(static)
      for (sampleNo_1 = 0; sampleNo_1 < numSamples_1; sampleNo_1++) {
        for (dataPointNo_1 = 0; dataPointNo_1 < numDataPointsPerSample_1; dataPointNo_1++) {
          int offset_1 = tmp_1->getPointOffset(sampleNo_1,dataPointNo_1);
          int offset_2 = tmp_2->getPointOffset(sampleNo_1,dataPointNo_1);
          const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));
          const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
          typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
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
      const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
      // Get the pointers to the actual data
      const typename decltype(operation)::first_argument_type *ptr_0 = &(tmp_0->getDefaultValueRO(0, dummy0));
      typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyr));
      // Compute a result for the default
      tensor_binary_operation(size0, ptr_0, ptr_1[0], ptr_2, operation);
      // Compute a result for each tag
      const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_0.begin();i!=lookup_0.end();i++) {
        tmp_2->addTag(i->first);
        const typename decltype(operation)::first_argument_type *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0, dummy0));
        typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummyr));
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

      // Get the pointers to the actual data
      const typename decltype(operation)::first_argument_type *ptr_0 = &(tmp_0->getDefaultValueRO(0, dummy0));
      const typename decltype(operation)::second_argument_type *ptr_1 = &(tmp_1->getDefaultValueRO(0, dummy1));
      typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDefaultValueRW(0, dummyr));

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
        const typename decltype(operation)::first_argument_type *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0, dummy0));
        const typename decltype(operation)::second_argument_type *ptr_1 = &(tmp_1->getDataByTagRO(i->first,0, dummy1));
        typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0, dummyr));
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
      res.requireWrite();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,0); // They're all the same, so just use #0
        const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
          typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
          tensor_binary_operation(size0, ptr_0, ptr_1[0], ptr_2, operation);
        }
      }

    }
    else if (arg_0_Z.isExpanded()   && arg_1_Z.isConstant()) {
      res = Data(0.0, shape0, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
      DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
      DataConstant* tmp_1=dynamic_cast<DataConstant*>(arg_1_Z.borrowData());
      DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());

      int sampleNo_0;
      int numSamples_0 = arg_0_Z.getNumSamples();
      int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
      int offset_1 = tmp_1->getPointOffset(0,0);
      const typename decltype(operation)::second_argument_type *ptr_src = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
      typename decltype(operation)::second_argument_type ptr_1 = ptr_src[0];
      int size = size0 * numDataPointsPerSample_0;
      res.requireWrite();
      #pragma omp parallel for private(sampleNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,0);
          const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));
          typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
          tensor_binary_operation(size, ptr_0, ptr_1, ptr_2, operation);
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
      res.requireWrite();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        int offset_1 = tmp_1->getPointOffset(sampleNo_0,0);
        const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));
          typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
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
      res.requireWrite();
      #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
      for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
        for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
          int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
          int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
          const typename decltype(operation)::first_argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0, dummy0));
          const typename decltype(operation)::second_argument_type *ptr_1 = &(arg_1_Z.getDataAtOffsetRO(offset_1, dummy1));
          typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2, dummyr));
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
     throw DataException("Error - Operations (C_TensorUnaryOperation) not permitted on instances of DataEmpty.");
  }
  if (arg_0.isLazy())
  {
     throw DataException("Error - Operations not permitted on lazy data.");
  }
  // Interpolate if necessary and find an appropriate function space
  Data arg_0_Z = Data(arg_0);

  // Get rank and shape of inputs
  const DataTypes::ShapeType& shape0 = arg_0_Z.getDataPointShape();
  int size0 = arg_0_Z.getDataPointSize();

  // Sanity check on the types of the function
  if (typeid(typename UnaryFunction::argument_type)!=typeid(typename UnaryFunction::result_type))
  {
      throw DataException("Error - Types for C_TensorUnaryOperation must be consistent");
  }
  
  // Declare output Data object
  Data res;

  if (arg_0_Z.isConstant()) {
    res = Data(0.0, shape0, arg_0_Z.getFunctionSpace());      // DataConstant output
    const typename UnaryFunction::argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(0));
    typename UnaryFunction::result_type *ptr_2 = &(res.getDataAtOffsetRW(0));
    tensor_unary_operation(size0, ptr_0, ptr_2, operation);
  }
  else if (arg_0_Z.isTagged()) {

    // Borrow DataTagged input from Data object
    DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());

    // Prepare a DataTagged output 2
    res = Data(0.0, shape0, arg_0_Z.getFunctionSpace());   // DataTagged output
    res.tag();
    DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());

    // Get the pointers to the actual data
    const typename decltype(operation)::argument_type *ptr_0 = &(tmp_0->getDefaultValueRO(0));
    typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDefaultValueRW(0));
    // Compute a result for the default
    tensor_unary_operation(size0, ptr_0, ptr_2, operation);
    // Compute a result for each tag
    const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
    DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
    for (i=lookup_0.begin();i!=lookup_0.end();i++) {
      tmp_2->addTag(i->first);
      const typename decltype(operation)::argument_type *ptr_0 = &(tmp_0->getDataByTagRO(i->first,0));
      typename decltype(operation)::result_type *ptr_2 = &(tmp_2->getDataByTagRW(i->first,0));
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
	dataPointNo_0=0;
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
        const typename decltype(operation)::argument_type *ptr_0 = &(arg_0_Z.getDataAtOffsetRO(offset_0));
        typename decltype(operation)::result_type *ptr_2 = &(res.getDataAtOffsetRW(offset_2));
        tensor_unary_operation(size0*numDataPointsPerSample_0, ptr_0, ptr_2, operation);
    }
  }
  else {
    throw DataException("Error - C_TensorUnaryOperation: unknown combination of inputs");
  }

  return res;
}


Data
C_TensorUnaryOperation(Data const &arg_0,
                       escript::ESFunction operation,
		       DataTypes::real_t tol=0
		      );





}
#endif
