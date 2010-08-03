
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include <weipa/FileSavers.h>

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/object.hpp>
#include <boost/version.hpp>

using namespace boost::python;

/*! \mainpage Esys Documentation
 *
 * \version 3.0.0
 *
 * - \ref escript
 *
 * - \ref esys_exception "Esys Exception"
 *
 * - \ref finley
 *
 * - <a href="../../epydoc/index.html">Python module documentation (epydoc generated)</a>
 *
 */

/*! \page weipa Weipa
 * Weipa is the python module that contains the interfaces
 * to the C++ side of the escript data exporter.
 *
 * 
 *
 * \section class_desc Class Description:
 * None
 *
 * \section class_limits Class Limitations:
 * None
 *
 * \section class_conds Class Conditions of Use:
 * None
 *
 * \section class_throws Throws:
 * None
 *
 */

BOOST_PYTHON_MODULE(weipacpp)
{
#if BOOST_VERSION >= 103500
// params are: bool show_user_defined, bool show_py_signatures, bool show_cpp_signatures
  docstring_options docopt(true,true,false);
#endif

  def("_saveSilo", weipa::saveSilo, (args("filename", "cycle", "time", "domain", "datavars")));

}

