
/*****************************************************************************
*
* Copyright (c) 2003-2015 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"


#include <paso/Paso.h>

#include "SystemMatrixAdapter.h"
#include "TransportProblemAdapter.h"

#include "PasoException.h"   
#include "esysUtils/esysExceptionTranslator.h"

#include "escript/AbstractContinuousDomain.h"

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/detail/defaults_gen.hpp>
#include <boost/version.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(pasowrapcpp)
{
// This feature was added in boost v1.34
#if ((BOOST_VERSION/100)%1000 > 34) || (BOOST_VERSION/100000 >1)
  // params are: bool show_user_defined, bool show_py_signatures, bool show_cpp_signatures
  docstring_options docopt(true, true, false);
#endif


  register_exception_translator<paso::PasoException>(&(esysUtils::RuntimeErrorTranslator));

  class_<paso::SystemMatrixAdapter, bases<escript::AbstractSystemMatrix> >
      ("OperatorAdapter","A concrete class representing an operator. For more details, please see the c++ documentation.", no_init)
      .def("print_matrix_info",&paso::SystemMatrixAdapter::Print_Matrix_Info,(arg("full")=false),"prints information about a system matrix")
      .def("nullifyRowsAndCols",&paso::SystemMatrixAdapter::nullifyRowsAndCols)
      .def("resetValues",&paso::SystemMatrixAdapter::resetValues, "resets the matrix entries")
      .def("saveMM",&paso::SystemMatrixAdapter::saveMM,args("fileName"), 
"writes the matrix to a file using the Matrix Market file format")
      .def("saveHB",&paso::SystemMatrixAdapter::saveHB, args("filename"),
"writes the matrix to a file using the Harwell-Boeing file format");

  class_<paso::TransportProblemAdapter, bases<escript::AbstractTransportProblem> >
      ("TransportProblemAdapter","",no_init)
      .def("getSafeTimeStepSize",&paso::TransportProblemAdapter::getSafeTimeStepSize)
      .def("getUnlimitedTimeStepSize",&paso::TransportProblemAdapter::getUnlimitedTimeStepSize)
      .def("resetTransport",&paso::TransportProblemAdapter::resetTransport,
"resets the transport operator typically as they have been updated");
}
