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

#if !defined escript_DataConstant_20040323_H
#define escript_DataConstant_20040323_H

#include "DataAbstract.h"
#include "DataArrayView.h"

#include <boost/python/numeric.hpp>

namespace escript {

/**
   \brief
   DataConstant stores a single data point which represents the entire
   function space.

   Description:
   DataConstant stores a single data point which represents the entire
   function space.
*/
class DataConstant : public DataAbstract  {

 public:

  /**
     \brief
     Constructor for DataConstant objects.

     Description:
     Constructor for DataConstant objects.

     \param value - Input - Data value for a single point.
     \param what - Input - A description of what this data object represents.
  */
  DataConstant(const boost::python::numeric::array& value,
               const FunctionSpace& what);

  /**
     \brief
     Copy constructor. Performs a deep copy.
  */
  DataConstant(const DataConstant& other);

  /**
     \brief
     Alternative constructor for DataConstant objects.

     Description:
     Alternative Constructor for DataConstant objects.
     \param value - Input - Data value for a single point.
     \param what - Input - A description of what this data object represents.
  */
  DataConstant(const DataArrayView& value,
               const FunctionSpace& what);

  /**
     \brief
     Alternative constructor for DataConstant objects.

     Description:
     Alternative Constructor for DataConstant objects.
     \param other - Input - Data object to copy from.
     \param region - Input - region to copy.
  */
  DataConstant(const DataConstant& other,
               const DataArrayView::RegionType& region);

  /**
     \brief
     Alternative constructor for DataConstant objects.

     Description:
     Alternative Constructor for DataConstant objects.
     \param what - Input - A description of what this data object represents.
     \param shape - Input - the shape of each data-point.
     \param data - the data values for each data-point.
  */
  DataConstant(const FunctionSpace& what,
               const DataArrayView::ShapeType &shape,
               const DataArrayView::ValueType &data);

  /**
     \brief
     Write the data as a string.
  */
  std::string
  toString() const;

  /**
     \brief
     Return the offset for the given sample. This is a somewhat artificial notion
     but returns the offset in bytes for the given point into the container
     holding the point data. Only really necessary to avoid many DataArrayView
     objects.
     \param sampleNo - Input - sample number.
     \param dataPointNo - Input - data point number for the sample.
   */
  virtual
  DataArrayView::ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

  /**
     \brief
     Return a view into the data for the data point specified.
     \param sampleNo - Input - sample number.
     \param dataPointNo - Input - data point number for the sample.
  */
  virtual
  DataArrayView
  getDataPoint(int sampleNo,
               int dataPointNo);

  /**
     \brief
     Return the number of doubles stored for the Data object.
  */
  virtual
  DataArrayView::ValueType::size_type
  getLength() const;

  /**
     \brief
     Factory method that returns a newly created DataConstant object
     sliced from the specified region of this object.
     The caller is reponsible for managing the object created.
     \param region - Input - region to slice from this object.
  */
  virtual
  DataAbstract*
  getSlice(const DataArrayView::RegionType& region) const;

  /**
     \brief
     Copy the specified region from the given value.
     \param value - Input - Data object to copy from.
     \param region - Input - Region to copy.
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
  void
  reshapeDataPoint(const DataArrayView::ShapeType& shape);

  /**
    \brief
    Archive the underlying data values to the file referenced
    by ofstream. A count of the number of values expected to be written
    is provided as a cross-check.

    The return value indicates success (0) or otherwise (1).
  */
  int
  archiveData(std::ofstream& archiveFile,
              const DataArrayView::ValueType::size_type noValues) const;

  /**
    \brief
    Extract the number of values specified by noValues from the file
    referenced by ifstream to the underlying data structure.

    The return value indicates success (0) or otherwise (1).
  */
  int
  extractData(std::ifstream& archiveFile,
              const DataArrayView::ValueType::size_type noValues);

 protected:

 private:
  //
  // the actual data
  DataArrayView::ValueType m_data;

};

} // end of namespace
#endif
