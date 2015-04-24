/*****************************************************************************
*
* Copyright (c) 2015 by The University of Queensland
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
