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
                                                                           
#if !defined  escript_DataConstant_20040323_H
#define escript_DataConstant_20040323_H

#include "escript/Data/DataAbstract.h"
#include "escript/Data/DataArray.h"
#include "escript/Data/DataArrayView.h"

#include <boost/python/numeric.hpp>

namespace escript {
/**
   \brief
   DataConstant stores single data point which represents the entire
   function space.

   Description:
   DataConstant stores single data point which represents the entire
   function space.
*/
class DataConstant:public DataAbstract  {

 public:

  /**
     \brief
     Default constructor for DataConstant

     Description:
     Default constructor for DataConstant

     \param value Input - Data value for a single point.
     \param noSamples Input - number of samples.
     \param noDataPointsPerSample Input - Input.
     \param what Input - A description of what this data represents.

  */
  DataConstant(const boost::python::numeric::array& value, const FunctionSpace& what);
  /**
     \brief
     Copy constructor. Performs a deep copy.
  */
  DataConstant(const DataConstant& other);
  /**
     \brief
     Alternative constructor for DataConstant

     Description:
     Alternative Constructor for DataConstant
     \param value Input - Data value for a single point.
     \param noSamples Input - number of samples.
     \param noDataPointsPerSample Input - Input.
     \param what Input - A description of what this data represents.
  */
  DataConstant(const DataArrayView& value, const FunctionSpace& what);
  /**
     \brief
     Alternative constructor for DataConstant

     Description:
     Alternative Constructor for DataConstant
     \param other Input - Other DataConstant
     \param region Input - region to copy
  */
  DataConstant(const DataConstant& other, const DataArrayView::RegionType& region);
  /**
     \brief
     Write the data as a string.
  */
  std::string toString() const;
  /**
     \brief
     Return the offset for the given sample. This is somewhat artificial notion
     but returns the offset in bytes for the given point into the container
     holding the point data. Only really necessary to avoid many DataArrayView
     objects.
     \param sampleNo - Input, sample number.
     \param dataPointNo - Input, data point number for the sample.
   */
  virtual DataArrayView::ValueType::size_type getPointOffset(int sampleNo, int dataPointNo) const;
  /**
     \brief
     Return a view into the data for the data point specified.
     NOTE: Construction of the DataArrayView is a relatively expensive 
     operation
     \param sampleNo Input
     \param dataPointNo Input
  */
  virtual DataArrayView getDataPoint(int sampleNo, int dataPointNo);
  /**
     \brief
     Return the number of doubles stored for the Data
  */
  virtual ValueType::size_type getLength() const;
  /**
     \brief
     Factory method that returns a newly created DataConstant.
     The caller is reponsible for managing the object created.
  */
  virtual DataAbstract* getSlice(const DataArrayView::RegionType& region) const;  /**
     \brief
     Copy the specified region from the given value.
     \param value Input - Data to copy from
     \param region Input - Region to copy.
  */
  virtual void setSlice(const DataAbstract* value, const DataArrayView::RegionType& region);
  /**
     \brief
     Reshape the data point if the data point is currently rank 0.
     The original data point value is used for all values of the new
     data point.
  */
  void reshapeDataPoint(const DataArrayView::ShapeType& shape);
 protected:

 private:
  //
  // data
  DataArrayView::ValueType m_data;
};

} // end of namespace
#endif

