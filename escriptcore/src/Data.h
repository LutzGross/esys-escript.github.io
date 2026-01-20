
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

/** \file Data.h */

#ifndef __ESCRIPT_DATA_H__
#define __ESCRIPT_DATA_H__

#include "system_dep.h"
#include "DataAbstract.h"
#include "DataException.h"
#include "DataTypes.h"
#include "EsysMPI.h"
#include "FunctionSpace.h"
#include "DataVectorOps.h"
#include <algorithm>
#include <string>
#include <sstream>

#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/math/special_functions/bessel.hpp>

#ifndef ESCRIPT_MAX_DATA_RANK
#define ESCRIPT_MAX_DATA_RANK 4
#endif

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
class ESCRIPT_DLL_API Data {

  public:

  /**
     Constructors.
  */

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
     function space of inData the inData are tried to be interpolated to what,
     otherwise a shallow copy of inData is returned.
  */
  Data(const Data& inData,
       const FunctionSpace& what);

  /**
     \brief Copy Data from an existing vector
  */ 
  Data(const DataTypes::RealVectorType& value,
                 const DataTypes::ShapeType& shape,
                 const FunctionSpace& what,
                 bool expanded);

  /**
     \brief
     Constructor which creates a Data with points having the specified shape.

     \param value - Input - Single real value applied to all Data.
     \param dataPointShape - Input - The shape of each data point.
     \param what - Input - A description of what this data represents.
     \param expanded - Input - Flag, if true fill the entire container with
                       the given value. Otherwise a more efficient storage
                       mechanism will be used.
  */
  Data(DataTypes::real_t value,
       const DataTypes::ShapeType& dataPointShape,
       const FunctionSpace& what,
       bool expanded);

  /**
     \brief
     Constructor which creates a Data with points having the specified shape.

     \param value - Input - Single complex value applied to all Data.
     \param dataPointShape - Input - The shape of each data point.
     \param what - Input - A description of what this data represents.
     \param expanded - Input - Flag, if true fill the entire container with
                       the given value. Otherwise a more efficient storage
                       mechanism will be used.
  */
  explicit
  Data(DataTypes::cplx_t value,
       const DataTypes::ShapeType& dataPointShape,
       const FunctionSpace& what,
       bool expanded);

  /**
     \brief
     Constructor which performs a deep copy of a region from another Data object.

     \param inData - Input - Input Data object.
     \param region - Input - Region to copy.
  */
  Data(const Data& inData,
       const DataTypes::RegionType& region);


  /**
     \brief
     Constructor which copies data from a wrapped array.

     \param w - Input - Input data.
     \param what - Input - A description of what this data represents.
     \param expanded - Input - Flag, if true fill the entire container with
                       the value. Otherwise a more efficient storage
                       mechanism will be used.
  */       
  Data(const WrappedArray& w, const FunctionSpace& what,
           bool expanded);       
       

  /**
     \brief
     Constructor which creates a DataConstant.
     Copies data from any object that can be treated like a python array/sequence.
     All other parameters are copied from other.

     \param value - Input - Input data.
     \param other - Input - contains all other parameters.
  */
  Data(const boost::python::object& value,
       const Data& other);
  
  /**
     This constructor subsumes a number of previous python ones.
     
  Data(const boost::python::object& value,
       const FunctionSpace& what=FunctionSpace(),
       bool expanded=false);
       
  Data(DataTypes::real_t value,
       const boost::python::tuple& shape=boost::python::make_tuple(),
       const FunctionSpace& what=FunctionSpace(),
       bool expanded=false);
       
  and a new 
  
  Data(cplx_t value,
       const boost::python::tuple& shape=boost::python::make_tuple(),
       const FunctionSpace& what=FunctionSpace(),
       bool expanded=false);  
  
  */
  Data(boost::python::object value,
       boost::python::object par1=boost::python::object(),
       boost::python::object par2=boost::python::object(),
       boost::python::object par3=boost::python::object());  
  
  
  
  

  /**
        \brief Create a Data using an existing DataAbstract. Warning: The new object assumes ownership of the pointer!
        Once you have passed the pointer, do not delete it.
  */
  explicit Data(DataAbstract* underlyingdata);

  /**
        \brief Create a Data based on the supplied DataAbstract
  */
  explicit Data(DataAbstract_ptr underlyingdata);

  /**
     \brief
     Destructor
  */
  ~Data();

  /**
     \brief Make this object a deep copy of "other".
  */
  void
  copy(const Data& other);

  /**
     \brief Return a pointer to a deep copy of this object.
  */
  Data
  copySelf() const;


  /**
     \brief produce a delayed evaluation version of this Data.
  */
  Data
  delay();

  /**
     \brief convert the current data into lazy data.
  */
  void 
  delaySelf();


  /**
     Member access methods.
  */

  /**
     \brief
     switches on update protection

  */
  void
  setProtection();

  /**
     \brief
     Returns true, if the data object is protected against update

  */
  bool
  isProtected() const;


/**
   \brief 
   Return the value of a data point as a python tuple.
*/
  const boost::python::object
  getValueOfDataPointAsTuple(int dataPointNo);

  /**
     \brief
     sets the values of a data-point from a python object on this process
  */
  void
  setValueOfDataPointToPyObject(int dataPointNo, const boost::python::object& py_object);

  /**
     \brief
     sets the values of a data-point from a array-like object on this process
  */
  void
  setValueOfDataPointToArray(int dataPointNo, const boost::python::object&);

  /**
     \brief
     sets the values of a data-point on this process
  */
  void
  setValueOfDataPoint(int dataPointNo, const DataTypes::real_t);
  
  void
  setValueOfDataPointC(int dataPointNo, const DataTypes::cplx_t);  
  

  /**
     \brief Return a data point across all processors as a python tuple.
  */
  const boost::python::object
  getValueOfGlobalDataPointAsTuple(int procNo, int dataPointNo);

  
  /**
     \brief Set the value of a global data point
  */
  void
  setTupleForGlobalDataPoint(int id, int proc, boost::python::object);
  
  /**
     \brief
     Return the tag number associated with the given data-point.

  */
  int
  getTagNumber(int dpno);


  /**
     \brief
     Write the data as a string. For large amounts of data, a summary is printed.
  */
  std::string
  toString() const;

  /**
     \brief
     Whatever the current Data type make this into a DataExpanded.
  */
  void
  expand();

  /**
     \brief
     If possible convert this Data to DataTagged. This will only allow
     Constant data to be converted to tagged. An attempt to convert
     Expanded data to tagged will throw an exception.
  */
  void
  tag();

  /**
    \brief If this data is lazy, then convert it to ready data.
    What type of ready data depends on the expression. For example, Constant+Tagged==Tagged.
  */
  void
  resolve();

  /**
    \brief returns return true if data contains NaN.
  \warning This is dependent on the ability to reliably detect NaNs on your compiler.
   See the nancheck function in LocalOps for details.
  */
  bool
  hasNaN();

  /**
  \brief replaces all NaN values with value 
  */
  void
  replaceNaN(DataTypes::real_t value);
  
  /**
  \brief replaces all NaN values with value 
  */
  void
  replaceNaN(DataTypes::cplx_t value);  
  
  /**
  \brief replaces all NaN values with value 
  */
  void
  replaceNaNPython(boost::python::object obj);  

  bool
  hasInf();

  void
  replaceInf(DataTypes::real_t value);

  void
  replaceInf(DataTypes::cplx_t value);

  void
  replaceInfPython(boost::python::object obj);
  

  /**
   \brief Ensures data is ready for write access.
  This means that the data will be resolved if lazy and will be copied if shared with another Data object.
  \warning This method should only be called in single threaded sections of code. (It modifies m_data).
  Do not create any Data objects from this one between calling requireWrite and getSampleDataRW.
  Doing so might introduce additional sharing.
  */
  void
  requireWrite();

  /**
     \brief
     Return true if this Data is expanded.
     \note To determine if a sample will contain separate values for each datapoint. Use actsExpanded instead.
  */
  bool
  isExpanded() const;

  /**
     \brief
     Return true if this Data is expanded or resolves to expanded.
     That is, if it has a separate value for each datapoint in the sample.
  */
  bool
  actsExpanded() const;
  

  /**
     \brief
     Return true if this Data is tagged.
  */
  bool
  isTagged() const;

  /**
     \brief
     Return true if this Data is constant.
  */
  bool
  isConstant() const;

  /**
     \brief Return true if this Data is lazy.
  */
  bool
  isLazy() const;

  /**
     \brief Return true if this data is ready.
  */
  bool
  isReady() const;

  /**
     \brief
     Return true if this Data holds an instance of DataEmpty. This is _not_ the same as asking if the object 
contains datapoints.
  */
  bool
  isEmpty() const;

  /**
    \brief
    True if components of this data are stored as complex
  */
  bool
  isComplex() const;

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
     Returns the spatial locations of the data points.
  */
  inline
  escript::Data
  getXFromFunctionSpace() const
  {
      // This is exposed to Python as [Data object].getX()
      return m_data->getFunctionSpace().getX();
  }

  /**
     \brief
     Return the domain.
  */
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
  inline
//   const AbstractDomain&
  Domain_ptr
  getDomainPython() const
  {
     return getFunctionSpace().getDomainPython();
  }

  /**
     \brief
     Returns the MPI communicator for this data's domain
     as a Python mpi4py.MPI.Comm object (or None if MPI/mpi4py not enabled)
  */
  inline
  boost::python::object
  getMPIComm() const
  {
     return getFunctionSpace().getMPIComm();
  }

  /**
     \brief
     Return the rank of the point data.
  */
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
  int
  getNoValues() const
  {
    return m_data->getNoValues();
  }

  /**
     \brief
     dumps the object into an HDF5 file
  */
  void
  dump(const std::string fileName) const;

 /**
  \brief returns the values of the object as a list of tuples (one for each datapoint).

  \param scalarastuple If true, scalar data will produce single valued tuples [(1,) (2,) ...]
If false, the result is a list of scalars [1, 2, ...]
 */
  const boost::python::object
  toListOfTuples(bool scalarastuple=true);


 /**
    \brief
    Return the sample data for the given sample no.
    Please do not use this unless you NEED to access samples individually
    \param sampleNo - Input - the given sample no.
    \return pointer to the sample data.
*/
  const DataTypes::real_t*
  getSampleDataRO(DataTypes::RealVectorType::size_type sampleNo, DataTypes::real_t dummy=0) const;

  const DataTypes::cplx_t*
  getSampleDataRO(DataTypes::CplxVectorType::size_type sampleNo, DataTypes::cplx_t dummy) const;
  

  /**
     \brief
     Return the sample data for the given sample no.
    Please do not use this unless you NEED to access samples individually
     \param sampleNo - Input - the given sample no.
     \return pointer to the sample data.
  */
  DataTypes::real_t*
  getSampleDataRW(DataTypes::RealVectorType::size_type sampleNo, DataTypes::real_t dummy=0);

  DataTypes::cplx_t*
  getSampleDataRW(DataTypes::RealVectorType::size_type sampleNo, DataTypes::cplx_t dummy);  
  
  

 /**
    \brief
    Return a pointer to the beginning of the underlying data
    \warning please avoid using this method since it by-passes possible lazy improvements. May be removed without notice.
    \return pointer to the data.
*/
  const DataTypes::real_t*
  getDataRO(DataTypes::real_t dummy=0) const;  
  
  const DataTypes::cplx_t*
  getDataRO(DataTypes::cplx_t dummy) const;    
  
  
  
  /**
     \brief
     Return the sample data for the given tag. If an attempt is made to
     access data that isn't tagged an exception will be thrown.
     \param tag - Input - the tag key.
  */
  inline
  DataTypes::real_t*
  getSampleDataByTag(int tag, DataTypes::real_t dummy=0)
  {
    return m_data->getSampleDataByTag(tag, dummy);
  }
  
  inline
  DataTypes::cplx_t*
  getSampleDataByTag(int tag, DataTypes::cplx_t dummy)
  {
    return m_data->getSampleDataByTag(tag, dummy);
  }  
  

  /**
     \brief
     Return a reference into the DataVector which points to the specified data point.
     \param sampleNo - Input -
     \param dataPointNo - Input -
  */
  DataTypes::RealVectorType::const_reference
  getDataPointRO(int sampleNo, int dataPointNo);

  /**
     \brief
     Return a reference into the DataVector which points to the specified data point.
     \param sampleNo - Input -
     \param dataPointNo - Input -
  */
  DataTypes::RealVectorType::reference
  getDataPointRW(int sampleNo, int dataPointNo);



  /**
     \brief 
     Return the offset for the given sample and point within the sample
  */
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
  const boost::python::tuple
  getShapeTuple() const;

  /**
     \brief
     Returns the product of the data point shapes
  */
  long
  getShapeProduct() const;


  /**
     \brief
     Return the size of the data point. It is the product of the
     data point shape dimensions.
  */
  int
  getDataPointSize() const;

  /**
     \brief
     Return the number of doubles stored for this Data.
  */
  DataTypes::RealVectorType::size_type
  getLength() const;

  /**
  \brief Return true if this object contains no samples.
  This is not the same as isEmpty() 
  */
  bool
  hasNoSamples() const
  {
        return m_data->getNumSamples()==0;
  }

  /**
     \brief
     Assign the given value to the tag assocciated with name. Implicitly converts this
     object to type DataTagged. Throws an exception if this object
     cannot be converted to a DataTagged object or name cannot be mapped onto a tag key.
     \param name - Input - name of tag.
     \param value - Input - Value to associate with given key.
  */
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
  void
  setTaggedValueFromCPP(int tagKey,
                        const DataTypes::ShapeType& pointshape,
                        const DataTypes::RealVectorType& value,
                        int dataOffset=0);


  void
  setTaggedValueFromCPP(int tagKey,
                        const DataTypes::ShapeType& pointshape,
                        const DataTypes::CplxVectorType& value,
                        int dataOffset=0);  

  /**
    \brief
    Copy other Data object into this Data object where mask is positive.
  */
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
  void
  setToZero();

  /**
     \brief
     Interpolates this onto the given functionspace and returns
     the result as a Data object.
     *
  */
  Data
  interpolate(const FunctionSpace& functionspace) const;

  Data
  interpolateFromTable3D(const WrappedArray& table, DataTypes::real_t Amin, DataTypes::real_t Astep,
                       DataTypes::real_t undef, Data& B, DataTypes::real_t Bmin, DataTypes::real_t Bstep, Data& C, 
                        DataTypes::real_t Cmin, DataTypes::real_t Cstep, bool check_boundaries);

  Data
  interpolateFromTable2D(const WrappedArray& table, DataTypes::real_t Amin, DataTypes::real_t Astep,
                       DataTypes::real_t undef, Data& B, DataTypes::real_t Bmin, DataTypes::real_t Bstep,bool check_boundaries);

  Data
  interpolateFromTable1D(const WrappedArray& table, DataTypes::real_t Amin, DataTypes::real_t Astep,
                       DataTypes::real_t undef,bool check_boundaries);


  Data
  interpolateFromTable3DP(boost::python::object table, DataTypes::real_t Amin, DataTypes::real_t Astep,
                        Data& B, DataTypes::real_t Bmin, DataTypes::real_t Bstep, Data& C, DataTypes::real_t Cmin, DataTypes::real_t Cstep, DataTypes::real_t undef,bool check_boundaries);


  Data
  interpolateFromTable2DP(boost::python::object table, DataTypes::real_t Amin, DataTypes::real_t Astep,
                        Data& B, DataTypes::real_t Bmin, DataTypes::real_t Bstep, DataTypes::real_t undef,bool check_boundaries);

  Data
  interpolateFromTable1DP(boost::python::object table, DataTypes::real_t Amin, DataTypes::real_t Astep,
                        DataTypes::real_t undef,bool check_boundaries);
  
  Data
  nonuniforminterp(boost::python::object in, boost::python::object out, bool check_boundaries);

  Data
  nonuniformslope(boost::python::object in, boost::python::object out, bool check_boundaries);  
  
  /**
     \brief
     Calculates the gradient of the data at the data points of functionspace.
     If functionspace is not present the function space of Function(getDomain()) is used.
     *
  */
  Data
  gradOn(const FunctionSpace& functionspace) const;

  Data
  grad() const;

  /**
    \brief
     Calculate the integral over the function space domain as a python tuple.
  */
  boost::python::object
  integrateToTuple_const() const;


  /**
    \brief
     Calculate the integral over the function space domain as a python tuple.
  */
  boost::python::object
  integrateToTuple();



  /**
     \brief
     Returns 1./ Data object
     *
  */
  Data
  oneOver() const;
  /**
     \brief
     Return a Data with a 1 for +ive values and a 0 for 0 or -ive values.
     *
  */
  Data
  wherePositive() const;

  /**
     \brief
     Return a Data with a 1 for -ive values and a 0 for +ive or 0 values.
     *
  */
  Data
  whereNegative() const;

  /**
     \brief
     Return a Data with a 1 for +ive or 0 values and a 0 for -ive values.
     *
  */
  Data
  whereNonNegative() const;

  /**
     \brief
     Return a Data with a 1 for -ive or 0 values and a 0 for +ive values.
     *
  */
  Data
  whereNonPositive() const;

  /**
     \brief
     Return a Data with a 1 for 0 values and a 0 for +ive or -ive values.
     *
  */
  Data
  whereZero(DataTypes::real_t tol=0.0) const;

  /**
     \brief
     Return a Data with a 0 for 0 values and a 1 for +ive or -ive values.
     *
  */
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
  DataTypes::real_t
  Lsup();

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
  DataTypes::real_t
  sup();

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
  DataTypes::real_t
  inf();

  DataTypes::real_t
  inf_const() const;



  /**
     \brief
     Return the absolute value of each data point of this Data object.
     *
  */
  Data
  abs() const;
  
  /**
     \brief
     Return the phase/arg/angular-part of complex values.
     *
  */
  Data
  phase() const;
  

  /**
     \brief
     Return the maximum value of each data point of this Data object.
     *
  */
  Data
  maxval() const;

  /**
     \brief
     Return the minimum value of each data point of this Data object.
     *
  */
  Data
  minval() const;

  /**
     \brief
     Return the (sample number, data-point number) of the data point with
     the minimum component value in this Data object.
     \note If you are working in python, please consider using Locator 
instead of manually manipulating process and point IDs.
  */
  const boost::python::tuple
  minGlobalDataPoint() const;

  /**
     \brief
     Return the (sample number, data-point number) of the data point with
     the minimum component value in this Data object.
     \note If you are working in python, please consider using Locator 
instead of manually manipulating process and point IDs.
  */
  const boost::python::tuple
  maxGlobalDataPoint() const;



  /**
     \brief
     Return the sign of each data point of this Data object.
     -1 for negative values, zero for zero values, 1 for positive values.
     *
  */
  Data
  sign() const;

  /**
     \brief
     Return the symmetric part of a matrix which is half the matrix plus its transpose.
     *
  */
  Data
  symmetric() const;

  /**
     \brief
     Return the antisymmetric part of a matrix which is half the matrix minus its transpose.
     *
  */
  Data
  antisymmetric() const;


  /**
     \brief
     Return the hermitian part of a matrix which is half the matrix plus its adjoint.
     *
  */
  Data
  hermitian() const;

  /**
     \brief
     Return the anti-hermitian part of a matrix which is half the matrix minus its hermitian.
     *
  */
  Data
  antihermitian() const;

  /**
     \brief
     Return the trace of a matrix
     *
  */
  Data
  trace(int axis_offset) const;

  /**
     \brief
     Transpose each data point of this Data object around the given axis.
     *
  */
  Data
  transpose(int axis_offset) const;

  /**
     \brief
     Return the eigenvalues of the symmetric part at each data point of this Data object in increasing values.
     Currently this function is restricted to rank 2, square shape, and dimension 3.
     *
  */
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
  const boost::python::tuple
  eigenvalues_and_eigenvectors(const DataTypes::real_t tol=1.e-12) const;

  /**
     \brief
     swaps the components axis0 and axis1
     *
  */
  Data
  swapaxes(const int axis0, const int axis1) const;

  /**
     \brief
     Return the error function erf of each data point of this Data object.
     *
  */
  Data
  erf() const;


  /**
     \brief
     For complex values return the conjugate values.
     For non-complex data return a copy
  */
  Data
  conjugate() const;
  
  Data
  real() const;  
  
  Data
  imag() const;  

  /**
     \brief
     Return the sin of each data point of this Data object.
     *
  */
  Data
  sin() const;

  /**
     \brief
     Return the cos of each data point of this Data object.
     *
  */
  Data
  cos() const;

  /**
     \brief
     Bessel worker function.
     *
  */
  Data
  bessel(int order, DataTypes::real_t (*besselfunc) (int,DataTypes::real_t) );
  

  /**
     \brief
     Return the Bessel function of the first kind for each data point of this Data object.
     *
  */
  Data
  besselFirstKind(int order);

  /**
     \brief
     Return the Bessel function of the second kind for each data point of this Data object.
     *
  */
  Data
  besselSecondKind(int order);


  /**
     \brief
     Return the tan of each data point of this Data object.
     *
  */
  Data
  tan() const;

  /**
     \brief
     Return the asin of each data point of this Data object.
     *
  */
  Data
  asin() const;

  /**
     \brief
     Return the acos of each data point of this Data object.
     *
  */
  Data
  acos() const;

  /**
     \brief
     Return the atan of each data point of this Data object.
     *
  */
  Data
  atan() const;

  /**
     \brief
     Return the sinh of each data point of this Data object.
     *
  */
  Data
  sinh() const;

  /**
     \brief
     Return the cosh of each data point of this Data object.
     *
  */
  Data
  cosh() const;

  /**
     \brief
     Return the tanh of each data point of this Data object.
     *
  */
  Data
  tanh() const;

  /**
     \brief
     Return the asinh of each data point of this Data object.
     *
  */
  Data
  asinh() const;

  /**
     \brief
     Return the acosh of each data point of this Data object.
     *
  */
  Data
  acosh() const;

  /**
     \brief
     Return the atanh of each data point of this Data object.
     *
  */
  Data
  atanh() const;

  /**
     \brief
     Return the log to base 10 of each data point of this Data object.
     *
  */
  Data
  log10() const;

  /**
     \brief
     Return the natural log of each data point of this Data object.
     *
  */
  Data
  log() const;

  /**
     \brief
     Return the exponential function of each data point of this Data object.
     *
  */
  Data
  exp() const;

  /**
     \brief
     Return the square root of each data point of this Data object.
     *
  */
  Data
  sqrt() const;

  /**
     \brief
     Return the negation of each data point of this Data object.
     *
  */
  Data
  neg() const;

  /**
     \brief
     Return the identity of each data point of this Data object.
     Simply returns this object unmodified.
     *
  */
  Data
  pos() const;

  /**
     \brief
     Return the given power of each data point of this Data object.

     \param right Input - the power to raise the object to.
     *
  */
  Data
  powD(const Data& right) const;

  /**
     \brief
     Return the given power of each data point of this boost python object.

     \param right Input - the power to raise the object to.
     *
   */
  Data
  powO(const boost::python::object& right) const;

  /**
     \brief
     Return the given power of each data point of this boost python object.

     \param left Input - the bases
     *
   */

  Data
  rpowO(const boost::python::object& left) const;

  /**
     \brief
     Overloaded operator +=
     \param right - Input - The right hand side.
     *
  */
  Data& operator+=(const Data& right);
  Data& operator+=(const boost::python::object& right);

  Data& operator=(const Data& other);

  /**
     \brief
     Overloaded operator -=
     \param right - Input - The right hand side.
     *
  */
  Data& operator-=(const Data& right);
  Data& operator-=(const boost::python::object& right);

 /**
     \brief
     Overloaded operator *=
     \param right - Input - The right hand side.
     *
  */
  Data& operator*=(const Data& right);
  Data& operator*=(const boost::python::object& right);

 /**
     \brief
     Overloaded operator /=
     \param right - Input - The right hand side.
     *
  */
  Data& operator/=(const Data& right);
  Data& operator/=(const boost::python::object& right);

  /**
    \brief
    Newer style division operator for python
  */
  Data truedivD(const Data& right);

  /**
    \brief
    Newer style division operator for python
  */
  Data truedivO(const boost::python::object& right);

  /**
    \brief
    Newer style division operator for python
  */
  Data rtruedivO(const boost::python::object& left);

  /**
    \brief
    wrapper for python add operation
  */
  boost::python::object __add__(const boost::python::object& right);
  

  /**
    \brief
    wrapper for python subtract operation
  */
  boost::python::object __sub__(const boost::python::object& right);
  
  /**
    \brief
    wrapper for python reverse subtract operation
  */
  boost::python::object __rsub__(const boost::python::object& right);  

  /**
    \brief
    wrapper for python multiply operation
  */
  boost::python::object __mul__(const boost::python::object& right);
    
  /**
    \brief
    wrapper for python divide operation
  */
  boost::python::object __div__(const boost::python::object& right);
  
  /**
    \brief
    wrapper for python reverse divide operation
  */
  boost::python::object __rdiv__(const boost::python::object& right);    
  
  /**
        \brief return inverse of matricies.
  */
  Data
  matrixInverse() const;

  /**
     \brief
     Returns true if this can be interpolated to functionspace.
  */
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
  void
  setItemD(const boost::python::object& key,
           const Data& value);

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
  void
  setSlice(const Data& value,
           const DataTypes::RegionType& region);

  /**
     \brief
     print the data values to stdout. Used for debugging
  */
  void
  print(void);

  /**
     \brief
     return the MPI rank number of the local data.
  */
        int
        get_MPIRank(void) const;

  /**
     \brief
     return the number of MPI ranks of the data which is the size of communicator group of the associated domain.
  */
        int
        get_MPISize(void) const;

  /**
     \brief
     return the MPI communicator of the data which is the communicator of the associated domain.
  */
        MPI_Comm
        get_MPIComm(void) const;

  /**
     \brief
     return the object produced by the factory, which is a DataConstant or DataExpanded
        TODO Ownership of this object should be explained in doco.
  */
        DataAbstract*
        borrowData(void) const;

        DataAbstract_ptr
        borrowDataPtr(void) const;

        DataReady_ptr
        borrowReadyPtr(void) const;



  /**
     \brief
     Return a pointer to the beginning of the datapoint at the specified offset.
     TODO Eventually these should be inlined.
     \param i - position(offset) in the underlying datastructure
  */

        DataTypes::RealVectorType::const_reference
        getDataAtOffsetRO(DataTypes::RealVectorType::size_type i, DataTypes::real_t dummy);

        DataTypes::RealVectorType::reference
        getDataAtOffsetRW(DataTypes::RealVectorType::size_type i, DataTypes::real_t dummy);
        
        DataTypes::CplxVectorType::const_reference
        getDataAtOffsetRO(DataTypes::CplxVectorType::size_type i, DataTypes::cplx_t dummy);

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
  DataTypes::RealVectorType&
  getExpandedVectorReference(DataTypes::real_t dummy=0);

  DataTypes::CplxVectorType&
  getExpandedVectorReference(DataTypes::cplx_t dummy);
  
  
  /**
   * \brief For tagged Data returns the number of tags with values.
   * For non-tagged data will return 0 (even Data which has been expanded from tagged).
  */ 
  size_t
  getNumberOfTaggedValues() const;

  /*
  * \brief make the data complex
  */
  void complicate();
 
 protected:

 private:
   void init_from_data_and_fs(const Data& inData,
           const FunctionSpace& functionspace);   
   
   template <typename S>
   void 
   maskWorker(Data& other2, Data& mask2, S sentinel);

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

  template<typename Scalar>
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
  reduction(BinaryFunction operation,
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
  initialise(const DataTypes::CplxVectorType& value,
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

  void
  initialise(const DataTypes::cplx_t value,
             const DataTypes::ShapeType& shape,
             const FunctionSpace& what,
             bool expanded);  
  //
  // flag to protect the data object against any update
  bool m_protected;
  bool m_lazy;

  //
  // pointer to the actual data object
//   boost::shared_ptr<DataAbstract> m_data;
  DataAbstract_ptr m_data;

// If possible please use getReadyPtr instead.
// But see warning below.
  const DataReady*
  getReady() const
{
   const DataReady* dr=dynamic_cast<const DataReady*>(m_data.get());
   ESYS_ASSERT(dr!=0, "error casting to DataReady.");
   return dr;
}  

  DataReady*
  getReady()
{
   DataReady* dr=dynamic_cast<DataReady*>(m_data.get());
   ESYS_ASSERT(dr!=0, "error casting to DataReady.");
   return dr;
}  


// Be wary of using this for local operations since it (temporarily) increases reference count.
// If you are just using this to call a method on DataReady instead of DataAbstract consider using 
// getReady() instead
  DataReady_ptr
  getReadyPtr()
{
   DataReady_ptr dr=REFCOUNTNS::dynamic_pointer_cast<DataReady>(m_data);
   ESYS_ASSERT(dr.get()!=0, "error casting to DataReady.");
   return dr;
}  

  const_DataReady_ptr
  getReadyPtr() const
{
   const_DataReady_ptr dr=REFCOUNTNS::dynamic_pointer_cast<const DataReady>(m_data);
   ESYS_ASSERT(dr.get()!=0, "error casting to DataReady.");
   return dr;
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
#ifdef SLOWSHARECHECK        
	return m_data->isShared();      // single threadsafe check for this
#else
	return !m_data.unique();
#endif	
  }

  void forceResolve()
  {
        if (isLazy())
        {
            #ifdef _OPENMP
            if (omp_in_parallel())
            {   // Yes this is throwing an exception out of an omp thread which is forbidden.
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
                std::ostringstream oss;
                oss << "Programming error. ExclusiveWrite required - please call requireWrite() isLazy=" << isLazy() << " isShared()=" << isShared(); 
                throw DataException(oss.str());
        }
  }

  /**
  \brief Modify the data abstract hosted by this Data object
  For internal use only.
  Passing a pointer to null is permitted (do this in the destructor)
  \warning Only to be called in single threaded code or inside a single/critical section. This method needs to be atomic.
  */
  void set_m_data(DataAbstract_ptr p);

  
  void TensorSelfUpdateBinaryOperation(const Data& right, escript::ES_optype operation);  
  
  friend class DataAbstract;            // To allow calls to updateShareStatus
  friend class TestDomain;              // so its getX will work quickly
#ifdef IKNOWWHATIMDOING
  friend Data applyBinaryCFunction(boost::python::object cfunc, boost::python::tuple shape, escript::Data& d, escript::Data& e);
#endif
  template <typename S>
  friend Data condEvalWorker(escript::Data& mask, escript::Data& trueval, escript::Data& falseval, S sentinel);
  friend ESCRIPT_DLL_API Data randomData(const boost::python::tuple& shape, const FunctionSpace& what, long seed, const boost::python::tuple& filter);

};


#ifdef IKNOWWHATIMDOING
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
#include "DataExpanded.h"
#include "DataConstant.h"
#include "DataTagged.h"

namespace escript
{



inline
DataTypes::real_t*
Data::getSampleDataRW(DataTypes::RealVectorType::size_type sampleNo, DataTypes::real_t dummy)
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
   return getReady()->getSampleDataRW(sampleNo, dummy);
}

inline
DataTypes::cplx_t*
Data::getSampleDataRW(DataTypes::CplxVectorType::size_type sampleNo, DataTypes::cplx_t dummy)
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
   return getReady()->getSampleDataRW(sampleNo, dummy);
}


inline
const DataTypes::real_t*
Data::getSampleDataRO(DataTypes::RealVectorType::size_type sampleNo,DataTypes::real_t dummy) const
{
   DataLazy* l=dynamic_cast<DataLazy*>(m_data.get());
   if (l!=0)
   {
        size_t offset=0;
        const DataTypes::RealVectorType* res=l->resolveSample(sampleNo,offset);
        return &((*res)[offset]);
   }
   return getReady()->getSampleDataRO(sampleNo, dummy);
}

inline
const DataTypes::cplx_t*
Data::getSampleDataRO(DataTypes::RealVectorType::size_type sampleNo, DataTypes::cplx_t dummy) const
{
   DataLazy* l=dynamic_cast<DataLazy*>(m_data.get());
   if (l!=0)
   {
	throw DataException("Programming error: complex lazy objects are not supported.");	
   }
   return getReady()->getSampleDataRO(sampleNo, dummy);
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
ESCRIPT_DLL_API
Data operator+(const Data& left, const Data& right);

/**
  \brief
  Operator-
  Takes two Data objects.
*/
ESCRIPT_DLL_API
Data operator-(const Data& left, const Data& right);

/**
  \brief
  Operator*
  Takes two Data objects.
*/
ESCRIPT_DLL_API
Data operator*(const Data& left, const Data& right);

/**
  \brief
  Operator/
  Takes two Data objects.
*/
ESCRIPT_DLL_API
Data operator/(const Data& left, const Data& right);

/**
  \brief
  Operator+
  Takes LHS Data object and RHS python::object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API
Data operator+(const Data& left, const boost::python::object& right);

/**
  \brief
  Operator-
  Takes LHS Data object and RHS python::object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API
Data operator-(const Data& left, const boost::python::object& right);

/**
  \brief
  Operator*
  Takes LHS Data object and RHS python::object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API
Data operator*(const Data& left, const boost::python::object& right);

/**
  \brief
  Operator/
  Takes LHS Data object and RHS python::object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API
Data operator/(const Data& left, const boost::python::object& right);

/**
  \brief
  Operator+
  Takes LHS python::object and RHS Data object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API
Data operator+(const boost::python::object& left, const Data& right);

/**
  \brief
  Operator-
  Takes LHS python::object and RHS Data object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API
Data operator-(const boost::python::object& left, const Data& right);

/**
  \brief
  Operator*
  Takes LHS python::object and RHS Data object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API
Data operator*(const boost::python::object& left, const Data& right);

/**
  \brief
  Operator/
  Takes LHS python::object and RHS Data object.
  python::object must be convertable to Data type.
*/
ESCRIPT_DLL_API
Data operator/(const boost::python::object& left, const Data& right);



/**
  \brief
  Output operator
*/
ESCRIPT_DLL_API
std::ostream& operator<<(std::ostream& o, const Data& data);

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
Data::reduction(BinaryFunction operation, DataTypes::real_t initial_value) const
{
  if (isExpanded()) {
    DataExpanded* leftC=dynamic_cast<DataExpanded*>(m_data.get());
    ESYS_ASSERT(leftC!=0, "Programming error - casting to DataExpanded.");

    DataExpanded& data=*leftC;
    int i,j;
    int numDPPSample=data.getNumDPPSample();
    int numSamples=data.getNumSamples();
    DataTypes::real_t global_current_value=initial_value;
    DataTypes::real_t local_current_value;
    const auto& vec=data.getTypedVectorRO(typename BinaryFunction::first_argument_type(0));
    const DataTypes::ShapeType& shape=data.getShape();
    // calculate the reduction operation value for each data point
    // reducing the result for each data-point into the current_value variables
    #pragma omp parallel private(local_current_value)
    {
	local_current_value=initial_value;
	#pragma omp for private(i,j) schedule(static)
	for (i=0;i<numSamples;i++) {
	  for (j=0;j<numDPPSample;j++) {
	    local_current_value=operation(local_current_value,escript::reductionOpVector(vec,shape,data.getPointOffset(i,j),operation,initial_value));

	  }
	}
	#pragma omp critical
	global_current_value=operation(global_current_value,local_current_value);
    }
    return global_current_value;
  } else if (isTagged()) {
    DataTagged* leftC=dynamic_cast<DataTagged*>(m_data.get());
    ESYS_ASSERT(leftC!=0, "Programming error - casting to DataTagged.");
    
    DataTagged& data=*leftC;
    DataTypes::real_t current_value=initial_value;

    const auto& vec=data.getTypedVectorRO(typename BinaryFunction::first_argument_type(0));
    const DataTypes::ShapeType& shape=data.getShape();
    const DataTagged::DataMapType& lookup=data.getTagLookup();
    const std::list<int> used=data.getFunctionSpace().getListOfTagsSTL();
    for (std::list<int>::const_iterator i=used.begin();i!=used.end();++i)
    {
      int tag=*i;
      DataTagged::DataMapType::const_iterator it=lookup.find(tag);
      if ((tag==0) || (it==lookup.end()))	// check for the default tag
      {
        current_value=operation(current_value,escript::reductionOpVector(vec,shape,data.getDefaultOffset(),operation,initial_value));
      }
      else
      {
        current_value=operation(current_value,escript::reductionOpVector(vec,shape,it->second,operation,initial_value));
      }
    }
    return current_value;    
  } else if (isConstant()) {
    DataConstant* leftC=dynamic_cast<DataConstant*>(m_data.get());
    ESYS_ASSERT(leftC!=0, "Programming error - casting to DataConstant.");
    return escript::reductionOpVector(leftC->getTypedVectorRO(typename BinaryFunction::first_argument_type(0)),leftC->getShape(),0,operation,initial_value);    
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
    
    
    
    int i,j;
    int numSamples=dataE->getNumSamples();
    int numDPPSample=dataE->getNumDPPSample();
  //  DataArrayView dataView=data.getPointDataView();
  //  DataArrayView resultView=result.getPointDataView();
    const auto& dataVec=dataE->getTypedVectorRO(initial_value);
    const DataTypes::ShapeType& shape=dataE->getShape();
    auto& resultVec=resultE->getTypedVectorRW(initial_value);
    // perform the operation on each data-point and assign
    // this to the corresponding element in result
    #pragma omp parallel for private(i,j) schedule(static)
    for (i=0;i<numSamples;i++) {
      for (j=0;j<numDPPSample;j++) {
	resultVec[resultE->getPointOffset(i,j)] =
	  escript::reductionOpVector(dataVec, shape, dataE->getPointOffset(i,j),operation,initial_value);

      }
    }    
    //escript::dp_algorithm(*dataE,*resultE,operation,initial_value);
    return result;
  }
  else if (isTagged()) {
    DataTagged* dataT=dynamic_cast<DataTagged*>(m_data.get());
    ESYS_ASSERT(dataT!=0, "Programming error - casting data to DataTagged.");
    DataTypes::RealVectorType defval(1);
    defval[0]=0;
    DataTagged* resultT=new DataTagged(getFunctionSpace(), DataTypes::scalarShape, defval, dataT);
    
    
    const DataTypes::ShapeType& shape=dataT->getShape();
    const auto& vec=dataT->getTypedVectorRO(initial_value);
    const DataTagged::DataMapType& lookup=dataT->getTagLookup();
    for (DataTagged::DataMapType::const_iterator i=lookup.begin(); i!=lookup.end(); i++) {
      resultT->getDataByTagRW(i->first,0) =
	  escript::reductionOpVector(vec,shape,dataT->getOffsetForTag(i->first),operation,initial_value);
    }    
    resultT->getTypedVectorRW(initial_value)[resultT->getDefaultOffset()] = escript::reductionOpVector(dataT->getTypedVectorRO(initial_value),dataT->getShape(),dataT->getDefaultOffset(),operation,initial_value);
    
    
    
    
    //escript::dp_algorithm(*dataT,*resultT,operation,initial_value);
    return Data(resultT);   // note: the Data object now owns the resultT pointer
  } 
  else if (isConstant()) {
    Data result(0,DataTypes::ShapeType(),getFunctionSpace(),isExpanded());
    DataConstant* dataC=dynamic_cast<DataConstant*>(m_data.get());
    DataConstant* resultC=dynamic_cast<DataConstant*>(result.m_data.get());
    ESYS_ASSERT(dataC!=0, "Programming error - casting data to DataConstant.");
    ESYS_ASSERT(resultC!=0, "Programming error - casting result to DataConstant.");
    
    DataConstant& data=*dataC;
    resultC->getTypedVectorRW(initial_value)[0] =
	escript::reductionOpVector(data.getTypedVectorRO(initial_value),data.getShape(),0,operation,initial_value);    
    
    //escript::dp_algorithm(*dataC,*resultC,operation,initial_value);
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
Data
C_TensorBinaryOperation(Data const &arg_0,
                        Data const &arg_1,
                        ES_optype operation);


Data
C_TensorUnaryOperation(Data const &arg_0,
                       escript::ES_optype operation,
                       DataTypes::real_t tol=0);

} // namespace escript

#endif // __ESCRIPT_DATA_H__
