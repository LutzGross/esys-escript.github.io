// $Id$
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
                                                                           
#if !defined escript_DataAbstract_20040315_H
#define escript_DataAbstract_20040315_H

#include "escript/Data/DataException.h"
#include "escript/Data/DataArrayView.h"
#include "escript/Data/DataArray.h"
#include "escript/Data/FunctionSpace.h"

#include <iostream>
#include <boost/scoped_ptr.hpp>
#include <functional>
#include <string>

namespace escript {

/**
   \brief
   DataAbstract provides an interface for the class of containers
   which hold ESyS data. 

   Description:
   DataAbstract provides an interface for the class of containers
   which hold ESyS data. The container may be thought of as a 2 dimensional
   array of data points. The data points themselves are arrays of rank 0-4.
*/

class DataAbstract {

 public:

  typedef DataArrayView::ValueType ValueType;
  typedef DataArrayView::ShapeType ShapeType;

  /**
     \brief
     Constructor for DataAbstract.

     Description:
     Constructor for DataAbstract.
     \param what - Input - A description of what this data represents.
  */
  DataAbstract(const FunctionSpace& what);

  /**
    \brief
    Destructor for DataAbstract.
  */
  virtual
  ~DataAbstract();

  /**
     \brief
     Write the data as a string.
  */
  virtual
  std::string
  toString() const = 0;

  /**
     \brief
     Return the number of data points per sample.
  */
  int
  getNumDPPSample() const;

  /**
     \brief
     Return the number of samples.
  */
  int
  getNumSamples() const;

  /**
     \brief
     Return the DataArrayView of the point data. This essentially contains 
     the shape information for each data point although it also may be used
     to manipulate the point data.
  */
  const DataArrayView&
  getPointDataView() const;

  DataArrayView&
  getPointDataView();

  /**
     \brief
     Return the offset for the given sample. This is somewhat artificial notion
     but returns the offset for the given point into the container
     holding the point data. Only really necessary to avoid many DataArrayView
     objects.
     \param sampleNo - Input - sample number.
     \param dataPointNo - Input - Input.
   */
  virtual
  ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const = 0;

  /**
     \brief
     Return the sample data for the given sample no.
  */
  double*
  getSampleData(ValueType::size_type sampleNo);

  /**
     \brief
     Return the number of doubles stored for the Data.
  */
  virtual
  ValueType::size_type
  getLength() const = 0;

  /**
     \brief
     Return the sample data for the given tag key.
     NB: If the data isn't tagged an exception will be thrown.
  */
  virtual
  double*
  getSampleDataByTag(int tag);

  /**
     \brief
     Assign the given value to the data-points(s) referenced by the given
     reference number.

     If this Data object cannot be accessed by reference numbers an
     exception will be thrown.

     \param ref - Input - reference number.
     \param value - Input - value to assign to data-points associated with
                            the given reference number.
  */
  virtual
  void
  setRefValue(int ref,
              const DataArray& value);

  /**
     \brief
     Return the values associated with the data-point(s) referenced by the given
     reference number.

     If this Data object cannot be accessed by reference numbers an
     exception will be thrown.

     \param ref - Input - reference number.
     \param value - Output - object to receive data-points associated with
                             the given reference number.
  */
  virtual
  void
  getRefValue(int ref,
              DataArray& value);

  /**
     \brief
     Check this and the right operands are compatible. Throws
     an exception if they aren't.
     \param right - Input - The right hand side.
  */
  void
  operandCheck(const DataAbstract& right) const;

  /**
     \brief
     Return true if a valid sample point number.
  */
  bool
  validSamplePointNo(int samplePointNo) const;

  /**
     \brief
     Return true if a valid number.
  */
  bool
  validSampleNo(int sampleNo) const;
 
  /**
     \brief
     Return a view into the data for the data point specified.
     NOTE: Construction of the DataArrayView is a relatively expensive 
     operation.
     \param samplesNo Input
     \param dataPointNo Input
  */
  virtual
  DataArrayView
  getDataPoint(int samplesNo,
               int dataPointNo) = 0;

  /**
     \brief
     Return the function space.
  */
  const FunctionSpace&
  getFunctionSpace() const;

  /**
     \brief
     Return a newly constructed DataAbstract. The caller is responsible for 
     managing the object created.
  */
  virtual
  DataAbstract*
  getSlice(const DataArrayView::RegionType& region) const = 0;

  /**
     \brief
     Copy the specified region from the given value.
     \param value - Input - Data to copy from
     \param region - Input - Region to copy.
  */
  virtual
  void
  setSlice(const DataAbstract* value,
           const DataArrayView::RegionType& region) = 0;

  /**
     \brief
     Reshape the data point if the data point is currently rank 0.
     Will throw an exception if the data points are not rank 0.
     The original data point value is used for all values of the new
     data point.
  */
  virtual
  void
  reshapeDataPoint(const DataArrayView::ShapeType& shape) = 0;

  /**
     \brief
     setTaggedValue
                                                                                                                                   
     Description:
     Assign the given value to the given tag.
     \param tagKey - Input - Integer key.
     \param value - Input - Single DataArrayView value to be assigned to the tag.
     NB: If the data isn't tagged an exception will be thrown.
  */
  virtual
  void
  setTaggedValue(int tagKey,
                 const DataArrayView& value);

 protected:

  /**
     \brief
     Set the pointDataView
     \param right - Input - The point data view. DataAbstract takes ownership
     of the DataArrayView provided. It will delete it when it is destructed.
  */
  void
  setPointDataView(const DataArrayView& input);

 private:

  int m_noDataPointsPerSample;

  int m_noSamples;

  //
  // Provides a view of the data as point data
  boost::scoped_ptr<DataArrayView> m_pointDataView;

  //
  // function space
  FunctionSpace m_functionSpace;

};

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
  return &(m_pointDataView->getData(getPointOffset(sampleNo,0)));
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

inline
const
DataArrayView&
DataAbstract::getPointDataView() const
{
  return *(m_pointDataView.get());
}

inline
DataArrayView&
DataAbstract::getPointDataView()
{
  return *(m_pointDataView.get());
}

} // end of namespace
#endif
