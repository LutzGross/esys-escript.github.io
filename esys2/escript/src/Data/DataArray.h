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

#if !defined escript_DataArray_20040421_H
#define escript_DataArray_20040421_H

#include "DataArrayView.h"

#include <boost/scoped_ptr.hpp>
#include <boost/python/object.hpp>

namespace escript {

/**
   \brief
   DataArray contains a DataArrayView plus the vector of data values
   associated with the View.

   Description:
   DataArray implements the management of the underlying data values contained in
   an escript Data object. It consists of a vector (m_data) which holds all the individual
   data values, plus a DataArrayView (m_dataView) which defines the shape of the data
   points contained in the Data object.

*/

class DataArray {

 public:

  /**
     \brief
     Default constructor for DataArray.

     Description:
     Default constructor for DataArray.
     Creates a data vector containing a single value, plus a DataArrayView
     which presents this data value as a scalar Data object.
  */
  DataArray(double value=0.0);

  /**
     \brief
     Constructor for DataArray.

     Description:
     Constructor for DataArray of given shape.
     Assigns each element of the shape the given value.
  */
  DataArray(const DataArrayView::ShapeType& shape,
            double value=0.0);

  /**
     \brief
     Copy constructor for DataArray.
     
     Description:
     Copy constructor for DataArray.
     Takes a DataArray and performs a deep copy.
  */
  DataArray(const DataArray& value);

  /**
     \brief
     Constructor for DataArray.
     
     Description:
     Constructor for DataArray.
     Takes a DataArrayView and performs a deep copy.
  */
  DataArray(const DataArrayView& value);

  /**
     \brief
     Constructor for DataArray.

     Description:
     Constructor for DataArray.
     Takes a boost::python::object.

     Throws:
     A DataException if a DataArray cannot be created from the python object.
  */
  DataArray(const boost::python::object& value);

  /**
     \brief
     Constructor for DataArray.

     Description:
     Constructor for DataArray.
     Takes a boost::python::numeric::array.
  */
  DataArray(const boost::python::numeric::array& value);

  /**
     \brief
     Return a reference to the DataArrayView.
  */
  const DataArrayView&
  getView() const;

  DataArrayView&
  getView();

  /**
     \brief
     Return a reference to the the data vector.
  */
  const DataArrayView::ValueType&
  getData() const;

  DataArrayView::ValueType&
  getData();

 protected:

 private:

  /**
     \brief
     Performs initialisation common to DataArray.
  */
  void
  initialise(const boost::python::numeric::array& value);

  //
  // data vector
  // this is a simple STL vector of floats
  DataArrayView::ValueType m_data;

  //
  // pointer to view of the data vector
  // this is a DataArrayView
  boost::scoped_ptr<DataArrayView> m_dataView;

};

} // end of namespace
#endif
