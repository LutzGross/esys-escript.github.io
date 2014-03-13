
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


SubWorld::SubWorld(MPI_Comm comm)
    :communicator(comm), domain((AbstractDomain*)0)
{


}

SubWorld::~SubWorld()
{
}

MPI_Comm SubWorld::getComm()
{
    return communicator;
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
std::cerr << "Here" << __LINE__ << std::endl;	
// 	std::cerr << bp::extract<std::string>(e.attr("__str__")())() << std::cerr;

// 	PyObject* xx=PyErr_Occurred();
// 	handle<> xy(borrowed(xx));
// 	object zzz(xy); 

	
// 	Still trying to work out how to get the exception text from the values below
// 	It seems the best approach at the moment is to work out how do do it using the C API
// 	Should probably have a #def to switch this stuff off in case we need to give people instructions
// 	about how to disable it if it goes nasty with a python or boost update
	

//	object mess=zzz.attr("message");
	
//	std::cerr << (extract<std::string>(mess)()) << std::endl;

std::cerr << "Here" << __LINE__ << std::endl;	

  	PyObject* ptype=0;
 	PyObject* pvalue=0;
 	PyObject* ptraceback=0;
 	PyErr_Fetch(&ptype, &pvalue, &ptraceback);
	PyErr_NormalizeException(&ptype, &pvalue, &ptraceback);
 
std::cerr << "Here" << __LINE__ << std::endl;	
 
 
std::cerr << PyString_Check(pvalue) << std::endl;

	PyObject* errobj=PyObject_Str(pvalue);

	errormsg=PyString_AsString(errobj);
	Py_XDECREF(errobj);

	Py_XDECREF(ptype);
	Py_XDECREF(pvalue);
	Py_XDECREF(ptraceback);
 	
std::cerr << errormsg << std::endl;


	
	
//  	handle<> htype(ptype);
//  	handle<> hvalue(ptype);
//  	handle<> htraceback(ptype);
//  	
// 	object otype(htype);
// 	object ovalue(hvalue);
// 	
// 
// 	
//  	object otraceback(htraceback);
 	
//  	std::cerr << (extract<std::string>(otype.attr("__str__")())()) << std::endl;
//  	std::cerr << (extract<std::string>(ovalue.attr("__str__")())()) << std::endl;
// 	std::cerr << (extract<std::string>(otraceback.attr("__str__")())()) << std::endl;
 	
 	
//  	object btraceback(handle<>(ptraceback));
//  	
//  	std::cerr << (extract<std::string>(btype.attr("__str__")())()) << std::endl;
 	
// 	bp::object eo(p);
//  	bp::object zz=eo.attr("__str__")();
 	//std::cerr << (bp::extract<std::string>(zz)()) << std::endl;
// 	std::cerr<< "[[[" << (bp::extract<std::string>(gettrace())()) << "]]]\n";
	return 3;
    }  
    return ret;
}