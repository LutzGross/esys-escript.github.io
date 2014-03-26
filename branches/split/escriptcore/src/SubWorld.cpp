
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
#include "SplitWorldException.h"

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

// Deal with the i-th job's exports
bool SubWorld::processExports(size_t i, std::string& errmsg)
{
    bp::dict expmap=bp::extract<bp::dict>(jobvec[i].attr("exportedvalues"))();	
    bp::list items=expmap.items();
    size_t l=bp::len(items);
    for (int j=0;j<l;++j)
    {
	bp::object o1=items[j][0];
	bp::object o2=items[j][1];
	bp::extract<std::string> ex1(o1);
	if (!ex1.check())
	{
	    errmsg="Job attempted export using a name which was not a string.";
	    return false;
	}
	std::string name=ex1();
	std::map<std::string, Reducer_ptr>::iterator it=reducemap.find(name);
	if (it==reducemap.end())
	{
	    errmsg="Attempt to export variable \""+name+"\". SplitWorld was not told about this variable.";
	    return false;
	}
	// so now we know it is a known name, we check that it is not None and that it is compatible
	if (o2.is_none())
	{
	    errmsg="Attempt to export variable \""+name+"\" with value of None, this is not permitted.";
	    return false;
	}
	if (!(it->second)->valueCompatible(o2))
	{
	    errmsg="Attempt to export variable \""+name+"\" with an incompatible value.";
	    return false;
	}
	if (!(it->second)->reduceLocalValue(o2, errmsg))
	{
	    return false;	// the error string will be set by the reduceLocalValue
	}
    }
    return true;
}

// clears out all the old import and export values from _Jobs_
// does not clear values out of reducers
void SubWorld::clearImportExports()
{
    for (size_t i=0;i<jobvec.size();++i)
    {
	jobvec[i].attr("clearImports")();
	jobvec[i].attr("clearExports")();
    }
}

// if 4, a Job performed an invalid export
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
	    // now we check the exports from that Job,
	    // Do names and type match?
	    if (!processExports(i, errormsg))
	    {    
		return 4;
	    }
	    // check to see if we need to keep running
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

void SubWorld::addVariable(std::string& name, Reducer_ptr& rp)
{
    if (reducemap.find(name)!=reducemap.end())
    {
	std::ostringstream oss;
	oss << "There is already a variable called " << name;
	throw SplitWorldException(oss.str());    
    }
    reducemap[name]=rp;
}

void SubWorld::removeVariable(std::string& s)
{
    reducemap.erase(s);
}

