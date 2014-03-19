
/*****************************************************************************
*
* Copyright (c) 2014 by University of Queensland
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

#include "boost/python/import.hpp"
#include "SubWorld.h"


#include <iostream>

using namespace escript;
namespace bp=boost::python;
using namespace esysUtils;

SubWorld::SubWorld(JMPI& comm)
    :mpiinfo(comm), domain((AbstractDomain*)0)
{


}

SubWorld::~SubWorld()
{
}

JMPI& SubWorld::getMPI()
{
    return mpiinfo;
}

void SubWorld::setDomain(Domain_ptr d)
{
    domain=d;
}

Domain_ptr SubWorld::getDomain()
{
    return domain;
}

void SubWorld::addJob(boost::python::object j)
{
    jobvec.push_back(j);
}

void SubWorld::clearJobs()
{
    jobvec.clear();
}


// if 3, a Job threw an exception 
// if 2, a Job did not return a bool
// if 1, at least one Job returned False
// if 0, all jobs in this world returned True
char SubWorld::runJobs(std::string& errormsg)
{
    errormsg.clear();
    bp::object gettrace=bp::import("traceback").attr("format_exc");
    int ret=0;
    try
    {
	for (size_t i=0;i<jobvec.size();++i)
	{
	    boost::python::object result=jobvec[i].attr("work")();
	    boost::python::extract<bool> ex(result);
	    if (!ex.check() || (result.is_none()))
	    {
		return 2;	
	    }
	    if (!ex())
	    {
		ret=1;
	    }
	}
    } 
    catch (boost::python::error_already_set e)
    {
	using namespace boost::python;
  	PyObject* ptype=0;
 	PyObject* pvalue=0;
 	PyObject* ptraceback=0;
 	PyErr_Fetch(&ptype, &pvalue, &ptraceback);
	PyErr_NormalizeException(&ptype, &pvalue, &ptraceback);
 
	PyObject* errobj=PyObject_Str(pvalue);

	errormsg=PyString_AsString(errobj);
	Py_XDECREF(errobj);

	Py_XDECREF(ptype);
	Py_XDECREF(pvalue);
	Py_XDECREF(ptraceback);
 	
	return 3;
    }  
    return ret;
}
