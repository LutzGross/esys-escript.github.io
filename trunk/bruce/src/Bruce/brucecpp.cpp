// $Id$
/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT ACcESS 2005 -  All Rights Reserved                         *
 *                                                                            *
 * This software is the property of ACcESS.  No part of this code             *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that                          *
 * person has a software license agreement with ACcESS.                       *
 *                                                                            *
 ******************************************************************************
*/

#include "bruce/Bruce/Bruce.h"
#include "bruce/Bruce/BruceFactory.h"

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/detail/defaults_gen.hpp>

using namespace boost::python;

/**
   \page bruce Bruce
   Bruce is the python module name that contains the interfaces
   to the C++ wrapper to bruce.

   \version 1.0.0 

   \section class_desc Class Description:
   Bruce

   \section class_limits Class Limitations:
   None

   \section class_conds Class Conditions of Use:
   None

   \section throws Throws:
   None

*/

BOOST_PYTHON_MODULE(brucecpp)
{

  //
  // NOTE: The return_value_policy is necessary for functions that
  // return pointers.

  def ("Brick",bruce::brick,
      (arg("n0")=2,arg("n1")=2,arg("n2")=2,
      arg("l0")=1.0,arg("l1")=1.0,arg("l2")=1.0),
      return_value_policy<manage_new_object>());

  def ("Rectangle",bruce::rectangle,
      (arg("n0")=2,arg("n1")=2,
      arg("l0")=1.0,arg("l1")=1.0),
      return_value_policy<manage_new_object>());

  class_<bruce::Bruce, bases<escript::AbstractContinuousDomain> >("Bruce",init<>())
      .def(init<const bruce::Bruce&>())
      .def("getDescription",&bruce::Bruce::getDescription)
      .def("isValidFunctionSpaceType",&bruce::Bruce::isValidFunctionSpaceType)
      .def("functionSpaceTypeAsString",&bruce::Bruce::functionSpaceTypeAsString)
      .def("getContinuousFunctionCode",&bruce::Bruce::getContinuousFunctionCode)
      .def("getFunctionCode",&bruce::Bruce::getFunctionCode)
      .def("getDim",&bruce::Bruce::getDim)
      .def("getDataShape",&bruce::Bruce::getDataShape)
      .def("getX",&bruce::Bruce::getX)
      .def("setToX",&bruce::Bruce::setToX)
      .def("getSize",&bruce::Bruce::getSize)
      .def("setToSize",&bruce::Bruce::setToSize)
      .def(self == other<object>())
      .def(self == self)
      .def(self != other<object>())
      .def(self != self);

}
