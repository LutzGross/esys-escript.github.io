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
                                                                           
#if !defined  escript_DataEmpty_20040726_H
#define escript_DataEmpty_20040726_H

#include "DataAbstract.h"

namespace escript {

/**
   \brief
   Impliments the DataAbstract interface for an empty or null Data.

   Description:
   Impliments the DataAbstract interface for an empty or null Data.
*/
class DataEmpty:public DataAbstract {

 public:

  /**
     \brief
     Default constructor for DataEmpty

     Description:
     Default constructor for DataEmpty

  */
  DataEmpty();
  /**
     \brief
     Destructor
  */
  virtual ~DataEmpty();
  /**
     \brief
     Return a textual representation of the data
  */
  virtual std::string toString() const;
  /**
     \brief
     Return the offset for the given sample. This is somewhat artificial notion
     but returns the offset in bytes for the given point into the container
     holding the point data.
     \param sampleNo - Input, number of samples.
     \param dataPointNo - Input.
   */
  virtual DataArrayView::ValueType::size_type getPointOffset(int sampleNo,int dataPointNo) const;
  /**
     \brief
     Return a view into the data for the data point specified.
     NOTE: Construction of the DataArrayView is a relatively expensive 
     operation
     \param samplesNo - Input
     \param dataPointNo - Input
  */
  virtual DataArrayView getDataPoint(int samplesNo, int dataPointNo);
 /**
     \brief
     Return the sample data for the given tag key. If the data isn't tagged
     an exception will be thrown.
  */
  virtual double* getSampleDataByTag(int tag);
  /**
     \brief
     Throw a standard exception. This function is called if an attempt
     is made to use functions of DataEmpty that are not allowed.
  */
  void throwStandardException(const std::string& functionName) const;
  /**
     \brief
     Return the number of doubles stored for the Data
  */
  virtual ValueType::size_type getLength() const;
  /**
     \brief
     Factory method that returns a newly created DataEmpty.
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
  virtual void reshapeDataPoint(const DataArrayView::ShapeType& shape);
 protected:

 private:
};

} // end of namespace
#endif
