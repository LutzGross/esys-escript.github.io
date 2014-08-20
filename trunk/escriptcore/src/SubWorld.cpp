
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

#include "SubWorld.h"
#include "SplitWorldException.h"

#include <boost/python/import.hpp>
#include <boost/python/dict.hpp>

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

// takes an vector of bools of size 2*number of variables
// Each entry in the first group is true if this subworld has at least one job which exports that variable.
// Each entry in the second group is true if there is at least one Job in this world which wishes
// to import that variable.
// The order of the variables is determined by the order of keys in the reducemap
// Q: Why does it take chars then?
// A: Because I want raw storage that I can pass via MPI and bool is special cased.
bool SubWorld::localTransport(std::vector<char>& vb, std::string& errmsg)
{
    for (size_t i=0;i<jobvec.size();++i)
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
    }

	// If we get here, all of the (local) exports worked
	// Now, lets make a record for distrubtion    
    size_t l=reducemap.size();
    vb.resize(l*2);
    size_t i=0;
    for (str2reduce::iterator it=reducemap.begin();it!=reducemap.end();++it, ++i)
    {
	if (it->second->hasValue())
	{
	    vb[i]=1;
	}
	else
	{
	    vb[i]=0;
	}
	if (importmap[it->first])
	{
	    vb[i+l]=1;
	}
	else
	{
	    vb[i+l]=0;	  
	}
    }     
    return true;      
}
/*
// Deal with the i-th job's exports
bool SubWorld::processExportsLocal(size_t i, std::string& errmsg)
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
}*/
/*

void SubWorld::populateVarMoveInfo(vector<char>& vb)
{
    size_t l=reducemap.size();
    vb.resize(l*2);
    for (str2reduce::iterator it=reducemap.begin(), int i=0;it!=reducemap.end();++it, ++i)
    {
	if (it->second->hasValue())
	{
	    vb[i]=1;
	}
	else
	{
	    vb[i]=0;
	}
	if (importmap[it->first])
	{
	    vb[i+l]=1;
	}
	else
	{
	    vb[i+l]=0;	  
	}
    }
}*/

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
 	
	return 3;
    }  
    return ret;
}

// if manual import is false, add this new variable to all the Jobs in this world
void SubWorld::addVariable(std::string& name, Reducer_ptr& rp, bool manualimport)
{
    if (reducemap.find(name)!=reducemap.end())
    {
	std::ostringstream oss;
	oss << "There is already a variable called " << name;
	throw SplitWorldException(oss.str());    
    }
    reducemap[name]=rp;
    if (!manualimport)
    {
	for (size_t i=0;i<jobvec.size();++i)
	{
	    jobvec[i].attr("requestImport")(name);
	}
    }
}

void SubWorld::removeVariable(std::string& s)
{
    reducemap.erase(s);
}

