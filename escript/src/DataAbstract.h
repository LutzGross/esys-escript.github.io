/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#if !defined escript_DataAbstract_20040315_H
#define escript_DataAbstract_20040315_H
#include "system_dep.h"

#include "DataTypes.h"
#include "FunctionSpace.h"

#include <boost/scoped_ptr.hpp>
#include <boost/python/numeric.hpp>

#include <string>
#include <fstream>

namespace escript {

/**
   \brief
   DataAbstract provides an abstract interface for the class of containers
   which hold ESyS data.

   Description:
   DataAbstract provides an abstract interface for the class of containers
   which hold ESyS data. The container may be thought of as a 2 dimensional
   array of data points where one dimension corresponds to the number of samples
   and the other to the number of data points per sample as defined by the function
   space associated with each Data object. The data points themselves are arrays of
   doubles of rank 0-4.
*/

class DataAbstract {

 public:

  typedef DataTypes::ValueType ValueType;
  typedef DataTypes::ShapeType ShapeType;

  /**
     \brief
     Constructor for DataAbstract.

     Description:
     Constructor for DataAbstract.

     \param what - Input - A description of what this data represents.
  */
  ESCRIPT_DLL_API
  DataAbstract(const FunctionSpace& what, const ShapeType& shape);

  /**
    \brief
    Destructor for DataAbstract.
  */
  ESCRIPT_DLL_API
  virtual
  ~DataAbstract();

  /**
     \brief
     Write the data as a string.
  */
  ESCRIPT_DLL_API
  virtual
  std::string
  toString() const = 0;

  /**
     \brief Return a deep copy of the current object.
  */
  ESCRIPT_DLL_API
  virtual
  DataAbstract*
  deepCopy()=0;

 /**
     \brief
     dumps the object into a netCDF file
  */
  ESCRIPT_DLL_API
  virtual
  void
  dump(const std::string fileName) const;

  /**
     \brief
     Return the number of data points per sample.
  */
  ESCRIPT_DLL_API
  int
  getNumDPPSample() const;

  /**
     \brief
     Return the number of samples.
  */
  ESCRIPT_DLL_API
  int
  getNumSamples() const;

//  /**
//      \brief
//      Return the DataArrayView of the point data. This essentially contains
//      the shape information for each data point although it also may be used
//      to manipulate the point data.
//  */
//   ESCRIPT_DLL_API
//   DataArrayView&
//   getPointDataView();

//   ESCRIPT_DLL_API
//   const DataArrayView&
//   getPointDataView() const;

  /**
     \brief 
     Return the shape information for the point data.

     The omission of a non-constant form is deliberate.
  */
  ESCRIPT_DLL_API
  const DataTypes::ShapeType& 
  getShape() const;

  /**
     \brief 
     Return the rank information for the point data.
  */
  ESCRIPT_DLL_API
  int 
  getRank() const;



  /**
     \brief
     Return the offset for the given sample. This returns the offset for the given
     point into the container holding the point data. Only really necessary to
     avoid creating many DataArrayView objects.

     \param sampleNo - Input - sample number.
     \param dataPointNo - Input - data point number.
   */
  ESCRIPT_DLL_API
  virtual
  ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const = 0;

//  /**
//	return the container storing points for this object
//  */
//   ESCRIPT_DLL_API
//   virtual 
//   ValueType&
//   getVector();
// 
//   ESCRIPT_DLL_API
//   virtual 
//   const ValueType&
//   getVector() const;


  /**
     \brief
     Return the sample data for the given sample number.
  */
  ESCRIPT_DLL_API
  double*
  getSampleData(ValueType::size_type sampleNo);

  /**
     \brief
     Return the number of doubles stored for this Data object.
  */
  ESCRIPT_DLL_API
  virtual
  ValueType::size_type
  getLength() const = 0;

  /**
     \brief
     Return the sample data for the given tag key.
     NB: If the data isn't tagged an exception will be thrown.
  */
  ESCRIPT_DLL_API
  virtual
  double*
  getSampleDataByTag(int tag);


  /**
     \brief
     Check this and the given RHS operands are compatible. Throws
     an exception if they aren't.

     \param right - Input - The right hand side.
  */
  ESCRIPT_DLL_API
  void
  operandCheck(const DataAbstract& right) const;

  /**
     \brief
     Return true if a valid sample point number.
  */
  ESCRIPT_DLL_API
  bool
  validSamplePointNo(int samplePointNo) const;

  /**
     \brief
     Return true if a valid sample number.
  */
  ESCRIPT_DLL_API
  bool
  validSampleNo(int sampleNo) const;

//  /**
//      \brief
//      Return a view into the data for the data point specified.
//      NOTE: Construction of the DataArrayView is a relatively expensive
//      operation.
// 
//      \param sampleNo - Input - the sample number.
//      \param dataPointNo - Input - the data point number.
//  */
//   ESCRIPT_DLL_API
//   virtual
//   DataArrayView
//   getDataPoint(int sampleNo,
//                int dataPointNo) = 0;

  /**
     \brief
     Return the function space associated with this Data object.
  */
  ESCRIPT_DLL_API
  const
  FunctionSpace&
  getFunctionSpace() const;

  /**
     \brief
     Return the given slice from this object.

     NB: The caller is responsible for managing the object created.
  */
  ESCRIPT_DLL_API
  virtual
  DataAbstract*
  getSlice(const DataTypes::RegionType& region) const = 0;

  /**
     \brief
     Copy the specified region from the given object.

     \param value - Input - Data to copy from
     \param region - Input - Region to copy.
  */
  ESCRIPT_DLL_API
  virtual
  void
  setSlice(const DataAbstract* value,
           const DataTypes::RegionType& region) = 0;


//  /**
//      \brief
//      setTaggedValue
// 
//      Description:
//      Assign the given value to the given tag.
// 
//      NB: If the data isn't tagged an exception will be thrown.
// 
//      \param tagKey - Input - Integer key.
//      \param value - Input - Single DataArrayView value to be assigned to the tag.
//  */
//   ESCRIPT_DLL_API
//   virtual
//   void
//   setTaggedValue(int tagKey,
//                  const DataArrayView& value);


  /**
     \brief
     setTaggedValue

     Description:
     Assign the given value to the given tag.

     NB: If the data isn't tagged an exception will be thrown.

     \param tagKey - Input - Integer key.
     \param pointshape - Input - the shape of the value parameter.
     \param value - Input - 
  */
  ESCRIPT_DLL_API
  virtual
  void
  setTaggedValue(int tagKey,
		 const DataTypes::ShapeType& pointshape,
                 const DataTypes::ValueType& value,
		 int dataOffset=0);




  /**
     \brief
     Copy the numarray object to the data points in this object.

     Description:
     Copy the numarray object to the data points in this object.

     \param value Input - new values for the data points
  */
  ESCRIPT_DLL_API
  virtual void
  copyAll(const boost::python::numeric::array& value);

  /**
     \brief
     Copy a double value to the data point dataPointNo of sample sampleNo in this object.

     Description:
     Copy a double value to the data point dataPointNo of sample sampleNo in this object.

     \param sampleNo Input - sample number
     \param dataPointNo Input - data point of the sample
     \param value Input - new values for the data point
  */
  ESCRIPT_DLL_API
  virtual void
  copyToDataPoint(const int sampleNo, const int dataPointNo, const double value);

  /**
     \brief
     Copy the numarray object to the data point dataPointNo of sample sampleNo in this object.

     Description:
     Copy the numarray object to the data point dataPointNo of sample sampleNo in this object.

     \param sampleNo Input - sample number
     \param dataPointNo Input - data point of the sample
     \param value Input - new values for the data point
  */
  ESCRIPT_DLL_API
  virtual void
  copyToDataPoint(const int sampleNo, const int dataPointNo, const boost::python::numeric::array& value);


  /**
     \brief
     Return the tag number associated with the given data-point number.

     If the object cannot be referenced by tag numbers, an exception
     will be thrown.
  */
  ESCRIPT_DLL_API
  virtual
  int
  getTagNumber(int dpno);

  /**
     \brief
     Computes a symmetric matrix (A + AT) / 2

     \param ev - Output - a symmetric matrix

  */
  ESCRIPT_DLL_API
  virtual void
  symmetric(DataAbstract* ev);

  /**
     \brief
     Computes a nonsymmetric matrix (A - AT) / 2

     \param ev - Output - a nonsymmetric matrix

  */
  ESCRIPT_DLL_API
  virtual void
  nonsymmetric(DataAbstract* ev);

  /**
     \brief
     Computes the trace of a matrix

     \param ev - Output - the trace of a matrix

  */
  ESCRIPT_DLL_API
  virtual void
  trace(DataAbstract* ev, int axis_offset);

  /**
     \brief
     Transpose each data point of this Data object around the given axis.

     \param ev - Output - the transpose of a matrix

  */
  ESCRIPT_DLL_API
  virtual void
  transpose(DataAbstract* ev, int axis_offset);

  /**
     \brief
     swaps components axis0 and axis1

     \param ev - Output - swapped components

  */
  ESCRIPT_DLL_API
  virtual void
  swapaxes(DataAbstract* ev, int axis0, int axis1);
  /**
     \brief
     solves the eigenvalue problem this*V=ev*V for the eigenvalues ev

     \param ev - Output - eigenvalues in increasing order at each data point

  */
  ESCRIPT_DLL_API
  virtual void
  eigenvalues(DataAbstract* ev);

  /**
     \brief
     sets values to zero

  */
  ESCRIPT_DLL_API
  virtual void
  setToZero();

  /**
     \brief
     solves the eigenvalue problem this*V=ev*V for the eigenvalues ev and eigenvectors V

     \param ev - Output - eigenvalues in increasing order at each data point
     \param V - Output - corresponding eigenvectors. They are normalized such that their length is one
                         and the first nonzero component is positive.
     \param tol - Input - eigenvalue with relative distance tol are treated as equal.

  */

  ESCRIPT_DLL_API
  virtual void
  eigenvalues_and_eigenvectors(DataAbstract* ev,DataAbstract* V,const double tol=1.e-13);

  /**
     \brief
     reorders data sample ordered by reference_ids to the ordering of the functions space

     \param reference_ids - Input - reference_ids used for current ordering
  */
  ESCRIPT_DLL_API
  virtual void
  reorderByReferenceIDs(int *reference_ids);



  /**
	\brief
	Return the number of values in the shape for this object.
  */
  ESCRIPT_DLL_API
  int
  getNoValues() const;



  /**
      \brief get a reference to the beginning of a data point
  */
  ESCRIPT_DLL_API
  DataTypes::ValueType::const_reference
  getDataAtOffset(DataTypes::ValueType::size_type i) const;


  ESCRIPT_DLL_API
  DataTypes::ValueType::reference
  getDataAtOffset(DataTypes::ValueType::size_type i);


  /**
	\brief Provide access to underlying storage. Internal use only!
  */
  ESCRIPT_DLL_API
  virtual DataTypes::ValueType&
  getVector()=0;

  ESCRIPT_DLL_API
  virtual const DataTypes::ValueType&
  getVector() const=0;

 protected:

//  /**
//     \brief
//     Set the pointDataView DataArrayView associated with this object.
//
//     \param input - Input - The point data view. DataAbstract takes ownership
//     of the DataArrayView provided. It will delete it when it is destructed.
//  */
//   ESCRIPT_DLL_API
//   void
//   setPointDataView(const DataArrayView& input);
// 
//   ESCRIPT_DLL_API
//   void
//   resetPointDataView();






 private:

  //
  // The number of samples in this Data object.
  // This is derived directly from the FunctionSpace.
  int m_noSamples;

  //
  // The number of data points per sample in this Data object.
  // This is derived directly from the FunctionSpace.
  int m_noDataPointsPerSample;

  //
  // The DataArrayView of the data array associated with this object.
  // The data array is defined only in child classes of this class, it
  // is not defined in this abstract parent class.
//  boost::scoped_ptr<DataArrayView> m_pointDataView;

  //
  // A FunctionSpace which provides a description of the data associated
  // with this Data object.
  FunctionSpace m_functionSpace;

  //
  // The shape of the points stored in this view
  DataTypes::ShapeType m_shape;

  //
  // The number of values in each point
  int m_novalues;

  //
  // The rank of the points stored in this view
  int m_rank;

};


inline
DataTypes::ValueType::const_reference
DataAbstract::getDataAtOffset(DataTypes::ValueType::size_type i) const
{
	return getVector()[i];
}

inline
DataTypes::ValueType::reference
DataAbstract::getDataAtOffset(DataTypes::ValueType::size_type i)
{
	return getVector()[i];
}


inline
bool
DataAbstract::validSamplePointNo(int samplePointNo) const
{
  return ((0 <= samplePointNo) && (samplePointNo < m_noDataPointsPerSample));
}

inline
bool
DataAbstract::validSampleNo(int sampleNo) const
{
  return ((0 <= sampleNo) && (sampleNo < m_noSamples));
}

inline
DataAbstract::ValueType::value_type*
DataAbstract::getSampleData(ValueType::size_type sampleNo)
{
//   return &(m_pointDataView->getData(getPointOffset(sampleNo,0)));
  return &(getVector()[getPointOffset(sampleNo,0)]);
}

inline
int
DataAbstract::getNumDPPSample() const
{
  return m_noDataPointsPerSample;
}

inline
int
DataAbstract::getNumSamples() const
{
  return m_noSamples;
}

inline
const
FunctionSpace&
DataAbstract::getFunctionSpace() const
{
  return m_functionSpace;
}

// inline
// const
// DataArrayView&
// DataAbstract::getPointDataView() const
// {
//   return *(m_pointDataView.get());
// }

// inline
// DataArrayView&
// DataAbstract::getPointDataView()
// {
//   return *(m_pointDataView.get());
// }

inline
const DataTypes::ShapeType& 
DataAbstract::getShape() const
{
	return m_shape;
}

inline 
int
DataAbstract::getRank() const
{
	return m_rank;
}

inline
int
DataAbstract::getNoValues() const
{
	return m_novalues;
}

} // end of namespace

#endif
