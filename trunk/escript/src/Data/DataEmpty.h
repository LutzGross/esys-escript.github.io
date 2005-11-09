//$Id$
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

#if !defined escript_DataEmpty_20040726_H
#define escript_DataEmpty_20040726_H

#include "DataAbstract.h"

namespace escript {

/**
   \brief
   Implements the DataAbstract interface for an empty Data object.

   Description:
   Implements the DataAbstract interface for an empty Data object.
*/

class DataEmpty : public DataAbstract {

 public:

  /**
     \brief
     Default constructor for DataEmpty.

     Description:
     Default constructor for DataEmpty.

  */
  DataEmpty();

  /**
     \brief
     Destructor for DataEmpty.
  */
  virtual
  ~DataEmpty();

  /**
     \brief
     Return a textual representation of the Data object.
  */
  virtual
  std::string
  toString() const;

  /**
     \brief
     Return the offset for the given sample.
     NB: This will throw an exception as obviously an empty Data object contains no
     samples. An implementation is required by parent DataAbstract class.
     \param sampleNo - Input - Sample number.
     \param dataPointNo - Input - data-point number.
   */
  virtual
  DataArrayView::ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

  /**
     \brief
     Return a view into the data for the data point specified.
     NB: This will throw an exception as obviously an empty Data object contains no
     data points. An implementation is required by parent DataAbstract class.
     \param sampleNo - Input - Sample number.
     \param dataPointNo - Input - data-point number.
  */
  virtual
  DataArrayView
  getDataPoint(int sampleNo,
               int dataPointNo);

  /**
     \brief
     Return the number of doubles stored for the Data object.
     As this is an empty Data object, this method will always return 0.
  */
  virtual
  ValueType::size_type
  getLength() const;

  /**
     \brief
     Factory method that returns a newly created DataEmpty sliced from the
     current Data object according to the specified region.
     NB: This will throw an exception as obviously an empty Data object contains no
     data to slice from. An implementation is required by parent DataAbstract class.
  */
  virtual
  DataAbstract*
  getSlice(const DataArrayView::RegionType& region) const;

  /**
     \brief
     Set the current Data object according to the specified slice from the
     given input value.
     NB: This will throw an exception as obviously an empty Data object contains no
     data to slice to. An implementation is required by parent DataAbstract class.
     \param value Input - Data to copy from
     \param region Input - Region to copy.
  */
  virtual
  void
  setSlice(const DataAbstract* value,
           const DataArrayView::RegionType& region);

  /**
     \brief
     Reshape the data point if the data point is currently rank 0.
     NB: This will throw an exception as obviously an empty Data object contains no
     data to reshape. An implementation is required by parent DataAbstract class.
  */
  virtual
  void
  reshapeDataPoint(const DataArrayView::ShapeType& shape);

 protected:

 private:

  /**
     \brief
     Throw a standard exception. This function is called if an attempt
     is made to use functions of DataEmpty that are not valid.
  */
  void
  throwStandardException(const std::string& functionName) const;

};

} // end of namespace

#endif
