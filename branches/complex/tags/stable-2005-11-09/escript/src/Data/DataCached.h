/* $Id$ */

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

#if !defined escript_DataCached_20050414_H
#define escript_DataCached_20050414_H

#include "DataAbstract.h"

namespace escript {

/**
   \brief
   Implements the DataAbstract interface for a cached Data object.

   Description:
   Implements the DataAbstract interface for a cached Data object.
*/

class DataCached : public DataAbstract {

 public:

  /**
     \brief
     Default constructor for DataCached.

     Description:
     Default constructor for DataCached.

  */
  DataCached();

  /**
     \brief
     Destructor for DataCached.
  */
  virtual
  ~DataCached();

  /**
     \brief
     Return a textual representation of the Data object.
  */
  virtual
  std::string
  toString() const;

  /**
     \brief
     Return the offset into the data array for the data point specified.
     \param sampleNo - Input - Sample number.
     \param dataPointNo - Input - data-point number.
   */
  virtual
  DataArrayView::ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

  /**
     \brief
     Return a view into the data array for the data point specified.
     \param sampleNo - Input - Sample number.
     \param dataPointNo - Input - data-point number.
  */
  virtual
  DataArrayView
  getDataPoint(int sampleNo,
               int dataPointNo);

  /**
     \brief
     Return the number of doubles stored for this object.
  */
  virtual
  ValueType::size_type
  getLength() const;

  /**
     \brief
     Factory method that returns a newly created DataCached sliced from the
     current Data object according to the specified region.
     The caller is reponsible for managing the object created.
  */
  virtual
  DataAbstract*
  getSlice(const DataArrayView::RegionType& region) const;

  /**
     \brief
     Set the current Data object according to the specified slice from the
     given input value.
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
     The original data point value is used for all values of the new
     data point.
  */
  virtual
  void
  reshapeDataPoint(const DataArrayView::ShapeType& shape);

 protected:

 private:

  /**
     \brief
     Throw a standard exception. This function is called if an attempt
     is made to use functions of DataCached that are not valid.
  */
  void
  throwStandardException(const std::string& functionName) const;

  // data members will need to include:

  // 0. flag to turn caching on, or pass everything through to act like normal Data object (?)

  // 1. a stack of the operations applied

  // 2. a coresponding stack of the Data/DataVector objects which are the operands to the operations

  // 3. a DataVector representing the initial values of this object

};

} // end of namespace

#endif
