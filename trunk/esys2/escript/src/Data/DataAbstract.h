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
#include <fstream>
#include <boost/scoped_ptr.hpp>
#include <functional>
#include <string>

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
  DataArrayView&
  getPointDataView();

  const DataArrayView&
  getPointDataView() const;

  /**
     \brief
     Return the offset for the given sample. This returns the offset for the given
     point into the container holding the point data. Only really necessary to
     avoid creating many DataArrayView objects.

     \param sampleNo - Input - sample number.
     \param dataPointNo - Input - data point number.
   */
  virtual
  ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const = 0;

  /**
     \brief
     Return the sample data for the given sample number.
  */
  double*
  getSampleData(ValueType::size_type sampleNo);

  /**
     \brief
     Return the number of doubles stored for this Data object.
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
     Check this and the given RHS operands are compatible. Throws
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
     Return true if a valid sample number.
  */
  bool
  validSampleNo(int sampleNo) const;
 
  /**
     \brief
     Return a view into the data for the data point specified.
     NOTE: Construction of the DataArrayView is a relatively expensive 
     operation.

     \param sampleNo - Input - the sample number.
     \param dataPointNo - Input - the data point number.
  */
  virtual
  DataArrayView
  getDataPoint(int sampleNo,
               int dataPointNo) = 0;

  /**
     \brief
     Return the function space associated with this Data object.
  */
  const
  FunctionSpace&
  getFunctionSpace() const;

  /**
     \brief
     Return the given slice from this object.

     NB: The caller is responsible for managing the object created.
  */
  virtual
  DataAbstract*
  getSlice(const DataArrayView::RegionType& region) const = 0;

  /**
     \brief
     Copy the specified region from the given object.

     \param value - Input - Data to copy from
     \param region - Input - Region to copy.
  */
  virtual
  void
  setSlice(const DataAbstract* value,
           const DataArrayView::RegionType& region) = 0;

  /**
     \brief
     Reshape the data points if they are currently rank 0.
     Will throw an exception if the data points are not rank 0.
     The original data point value is used for all values of the new
     data point.
  */
  virtual
  void
  reshapeDataPoint(const ShapeType& shape) = 0;

  /**
     \brief
     setTaggedValue
                                                                                                                                   
     Description:
     Assign the given value to the given tag.

     NB: If the data isn't tagged an exception will be thrown.

     \param tagKey - Input - Integer key.
     \param value - Input - Single DataArrayView value to be assigned to the tag.
  */
  virtual
  void
  setTaggedValue(int tagKey,
                 const DataArrayView& value);

  /**
    \brief
    Archive the underlying data values to the file referenced
    by ofstream. A count of the number of values expected to be written
    is provided as a cross-check.

    The return value indicates success (0) or otherwise (1).
  */
  virtual
  int
  archiveData(std::ofstream& archiveFile,
              const ValueType::size_type noValues) const;

  /**
    \brief
    Extract the number of values specified by noValues from the file
    referenced by ifstream to the underlying data structure.

    The return value indicates success (0) or otherwise (1).
  */
  virtual
  int
  extractData(std::ifstream& archiveFile,
              const ValueType::size_type noValues);

 protected:

  /**
     \brief
     Set the pointDataView DataArrayView associated with this object.

     \param input - Input - The point data view. DataAbstract takes ownership
     of the DataArrayView provided. It will delete it when it is destructed.
  */
  void
  setPointDataView(const DataArrayView& input);

  void
  resetPointDataView();

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
  boost::scoped_ptr<DataArrayView> m_pointDataView;

  //
  // A FunctionSpace which provides a description of the data associated
  // with this Data object.
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
