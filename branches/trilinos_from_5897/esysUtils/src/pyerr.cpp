/*****************************************************************************
*
* Copyright (c)2015-2016 by The University of Queensland
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
#include <boost/python/object.hpp>
#include <boost/python/import.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>

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
	object tb = import("traceback"); 
	object trace(handle<>(borrowed(ptraceback)));
	object li=tb.attr("extract_tb")(trace);
	object li2=tb.attr("format_list")(li);
	list l=extract<list>(li2)();
	


#ifdef ESPYTHON3	
	std::string ss;
	for (int i=0;i<len(l);++i) {
	    object o=l[i];
	    PyObject* rr=PyUnicode_AsASCIIString(o.ptr());
	    ss+=PyBytes_AsString(rr);
	    Py_XDECREF(rr);
	}
	
	PyObject* errobj=PyObject_Str(pvalue);	
	
	PyObject* rr=PyUnicode_AsASCIIString(errobj);
	errormsg=PyBytes_AsString(rr);
	errormsg+="\n";
	Py_XDECREF(rr);
	errormsg+=ss;
#else
	
	std::string ss;
	for (int i=0;i<len(l);++i) {
	    ss+=extract<std::string>(l[i])();
	}
	
	PyObject* errobj=PyObject_Str(pvalue);	
	
	errormsg=PyString_AsString(errobj);
	errormsg+="\n";
	errormsg+=ss;
#endif
	Py_XDECREF(errobj);

	Py_XDECREF(ptype);
	Py_XDECREF(pvalue);
	Py_XDECREF(ptraceback);
}
