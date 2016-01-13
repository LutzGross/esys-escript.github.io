/*****************************************************************************
*
* Copyright (c) 2015 by University of Queensland
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
#include "first.h"
#include "pyerr.h"

#include <sstream>
#include <boost/python.hpp>

// Function factored out of SubWorld code

void getStringFromPyException(boost::python::error_already_set e, std::string& errormsg)
{
	using namespace boost::python;
  	PyObject* ptype=0;
 	PyObject* pvalue=0;
 	PyObject* ptraceback=0;
 	PyErr_Fetch(&ptype, &pvalue, &ptraceback);
	PyErr_NormalizeException(&ptype, &pvalue, &ptraceback);
 
	PyObject* errobj=PyObject_Str(pvalue);

#ifdef ESPYTHON3	
	PyObject* rr=PyUnicode_AsASCIIString(errobj);
	errormsg=PyBytes_AsString(rr);
	Py_XDECREF(rr);
#else
	errormsg=PyString_AsString(errobj);
#endif
	Py_XDECREF(errobj);

	Py_XDECREF(ptype);
	Py_XDECREF(pvalue);
	Py_XDECREF(ptraceback);
}

void getTraceStringFromPyException(boost::python::error_already_set e, std::string& errormsg)
{
	using namespace boost::python;
  	PyObject* ptype=0;
 	PyObject* pvalue=0;
 	PyObject* ptraceback=0;
 	PyErr_Fetch(&ptype, &pvalue, &ptraceback);
	PyErr_NormalizeException(&ptype, &pvalue, &ptraceback);
 
	PyObject* errobj=PyObject_Str(pvalue);
	
	// Need to invoke the equivalent of traceback.format_exception(ptype, pvalue, ptraceback)
	try
	{
	    std::ostringstream oss;
	    boost::python::object tb(boost::python::handle<>(boost::python::borrowed(ptraceback)));
	    
	    while (tb!=object())
	    {
		int lineno=boost::python::extract<int>(tb.attr("tb_lineno"))();
		boost::python::object code=tb.attr("tb_frame").attr("f_code");
		std::string fname=boost::python::extract<std::string>(code.attr("co_filename"))();
		std::string funct=boost::python::extract<std::string>(code.attr("co_name"))();
		oss << "File \"" << fname << "\", line " << lineno << ", in " << funct << std::endl;
		tb=tb.attr("tb_next");
	    }
	    boost::python::object eo(boost::python::handle<>(boost::python::borrowed(errobj)));
	    boost::python::object typ(boost::python::handle<>(boost::python::borrowed(ptype)));
	    boost::python::object typs=typ.attr("__name__").attr("__str__")();
	    oss << extract<std::string>(typs)() << ": ";
	    oss << extract<std::string>(eo)();
	    oss << std::endl;
	    errormsg=oss.str();
	}
	catch (boost::python::error_already_set zz)
	{
	    getStringFromPyException(zz, errormsg);
	    errormsg="Encountered a problem extracting an exception trace:"+errormsg;
	}
	Py_XDECREF(errobj);

	Py_XDECREF(ptype);
	Py_XDECREF(pvalue);
	Py_XDECREF(ptraceback);
}