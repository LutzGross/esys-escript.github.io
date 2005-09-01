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

  class_<bruce::Bruce, bases<escript::AbstractContinuousDomain> >("Bruce",init<>())
      .def(init<const bruce::Bruce&>())
      .def("getDim",&bruce::Bruce::getDim)
      .def("getX",&bruce::Bruce::getX)
      .def("getSize",&bruce::Bruce::getSize);

}
