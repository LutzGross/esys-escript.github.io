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
   DataArray contains a DataArrayView plus the Vector data associated with the view.

   Description:

*/

class DataArray {

 public:
  /**
     \brief
     Default constructor for DataArray.

     Description:
     Default constructor for DataArray. Creates a scalar.
  */
  DataArray(double value=0.0);

 /**
     \brief
     Constructor for DataArray.

     Description:
     Constructor for DataArray of shape "shape".
     Assigns each element the given value.
  */
  DataArray(const DataArrayView::ShapeType& shape,
            double value=0.0);

  /**
     \brief
     Copy constructor for DataArray.
     
     Description:
     Copy constructor for DataArray. Takes a DataArray.
  */
  DataArray(const DataArray& value);

  /**
     \brief
     Constructor for DataArray.
     
     Description:
     Constructor for DataArray.
     Takes a DataArrayView.
  */
  DataArray(const DataArrayView& value);

  /**
     \brief
     Constructor for DataArray.

     Description:
     Constructor for DataArray.
     Takes a boost::python::object.

     Throws:
     A DataException if a DataArray cannot be created from the python object
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
     Return the DataArrayView of the data.
  */
  const DataArrayView&
  getView() const;

  /**
     \brief
     Return the DataArrayView of the data.
     Non-const version.
  */
  DataArrayView&
  getView();

  /**
     \brief
     Return the data.
  */
  const DataArrayView::ValueType&
  getData() const;

  /**
     \brief
     Return the data, non-const version.
  */
  DataArrayView::ValueType&
  getData();

 protected:

 private:

  /**
     \brief
     Performs initialisation common of DataArray.
  */
  void
  initialise(const boost::python::numeric::array& value);

  //
  // data 
  DataArrayView::ValueType m_data;

  //
  // view of the data
  boost::scoped_ptr<DataArrayView> m_dataView;

};

} // end of namespace
#endif
