
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
namespace rs=escript::reducerstatus;

SubWorld::SubWorld(JMPI& comm, JMPI& corr)
    :swmpi(comm), corrmpi(corr), domain((AbstractDomain*)0)
{


}

SubWorld::~SubWorld()
{
}

JMPI& SubWorld::getMPI()
{
    return swmpi;
}

JMPI& SubWorld::getCorrMPI()
{
    return corrmpi;
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

  // query the jobs to find the names of requested imports
bool SubWorld::findImports(bool manualimports, std::string& errmsg)
{
      // set all imports to false (if manual import is true) else set all to true
      // query each job:
      //		1) get its wantedvalues field
      //		2) check to ensure each value is actually a known variable
      //		3) set its value in import map to true
    if (!manualimports)
    {
	  // import everything
	for (str2bool::iterator it=importmap.begin();it!=importmap.end();++it)
	{
	    it->second=true;
	}
    }
    else		// manual imports
    {
	  // import nothing
 	for (str2bool::iterator it=importmap.begin();it!=importmap.end();++it)
	{
	    it->second=false;
	}     
	
	for (size_t i=0;i<jobvec.size();++i)
	{
	    bp::list wanted=bp::extract<bp::list>(jobvec[i].attr("wantedvalues"))();		  
	    for (size_t j=0;j<len(wanted);++j)
	    {
		bp::extract<std::string> exs(wanted[j]);
		if (!exs.check()) 
		{
		    errmsg="names in wantedvalues must be strings";
		    return false;
		}
		std::string n=exs();
		  // now we need to check to see if this value is known
		str2reduce::iterator it=reducemap.find(n);
		if (it==reducemap.end())
		{
		    errmsg="Attempt to import variable \""+name+"\". SplitWorld was not told about this variable.";
		    return false;
		}
		importmap[n]=true;
	    }
	}
    }
    return true;
}

// this will give the imported values to interested jobs
bool SubWorld::deliverImports(std::string& errmsg)
{
    if (!manualimports)	// import all available vars into all jobs
    {
	for (size_t s=0;s<jobvec.size();++s)
	{
	    bp::object f=jobvec[s].attr("setImportValue");
	    
	    for (str2reduce::iterator it=reducemap.begin();it!=reducemap.end();++it)
	    {
		f(it->first, it->second->getValue());
	    }
	}
    }
    else		// manual imports
    {
      	for (size_t i=0;i<jobvec.size();++i)
	{
	    bp::object f=jobvec[s].attr("setImportValue");	    
	    bp::list wanted=bp::extract<bp::list>(jobvec[i].attr("wantedvalues"))();		  
	    for (size_t j=0;j<len(wanted);++j)
	    {
		bp::extract<std::string> exs(wanted[j]);
		if (!exs.check()) 
		{
		    errmsg="names in wantedvalues must be strings";
		    return false;
		}
		std::string n=exs();
		  // now we need to check to see if this value is known
		str2reduce::iterator it=reducemap.find(n);
		if (it==reducemap.end())
		{
		    errmsg="Attempt to import variable \""+name+"\". SplitWorld was not told about this variable.";
		    return false;
		}
		f(wanted[j], it->second->getValue());
	    }
	}
    }
    
    
    // TODO: This function needs error checking
    
    
    return true;
}


// takes a vector of bools of size 2*number of variables
// Each entry in the first group is true if this subworld has at least one job which exports that variable.
// Each entry in the second group is true if there is at least one Job in this world which wishes
// to import that variable.
// The order of the variables is determined by the order of keys in the reducemap
// Q: Why does it take chars then?
// A: Because I want raw storage that I can pass via MPI and vector<bool> is special cased.
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

bool SubWorld::checkRemoteCompatibility(std::string& errmsg)
{
    for (str2reduce::iterator it=reducemap.begin();it!=reducemap.end();++it)
    {
	if (! it->second->checkRemoteCompatibility(errmsg))
	{
	    return false;
	}
    }
    return true;
}


bool SubWorld::reduceRemoteValues()
{
    for (str2reduce::iterator it=reducemap.begin();it!=reducemap.end();++it)
    {
	if (! it->second->reduceRemoteValues(errmsg))
	{
	    return false;
	}
    }
    return true;    
}


  // Work out which variables we need to make available to python
bool SubWorld::deliverImports(std::vector<char>& vb, std::string& errmsg)
{
    errmsg="Haven't implemented this yet";
    return false;
}

bool SubWorld::amLeader()
{
    return swmpi->rank==0;
}


// Look at the vector of interest and send around variables to where they need to be
// To try to make error reporting consistant, every variable transfer will be attempted,
// even if an early one fails.
bool SubWorld::deliverGlobalImports(std::vector<char>& vb, string& errmsg)
{
#ifdef ESYS_MPI  
    size_t maxv=std::numeric_limits<size_t>::max();   // largest possible value in size_t
    size_t nv=getNumVars();
    str2reduce::iterator rmi=reducemap.begin();
    char error=0;
    for (size_t i=0;i<nv;++i, ++rmi)
    {
	size_t earliest=maxv;	
	for (int s=0;s<swcount;++s)
	{
	    if (vb[s*nv+i]==rs::HAVE)	// This world has a copy of the value
	    {
		earliest=s;	      
	    }
	}
	if (earliest==maxv)
	{
	      // noone has a value for the requested variable
	      // We could make that an error condition but for now, I'll let
	      // any job scripts deal with that
	    continue;
	}
	  // do I need that variable (but don't have a value)?
	if (vb[localid*nv+i]==rs::INTERESTED)
	{
	    if (!rmi->second->recvFrom(localid, earliest, corrmpi->comm)) 
	    {
		error=1;
		errmsg="Error receving value for variable.";
	    }
	} 
	else if (localid==earliest)	// we will be the one to ship the value to others
	{
	    for (uint target=0;target<swcount;++target)	// look through all worlds and see who we need to send to
	    {
		if (vb[target*nv+i]==rs::INTERESTED)	// nobody who has the value will be marked as INTERESTED
		{
		    if (!rmi->second->sendTo(localid, target, corrmpi))
		    {
			error=1;
			errmsg="Error sending value for variable.";
		    }
		}
	    }
	}
    }
      // It is possible that one or more subworlds were not participating in a failed transfer to we need to check
    char nerr;  
    if (MPI_Allreduce(&error, &nerr, 1, MPI_CHAR, MPI_MAX, globalcom->comm)!=MPI_SUCCESS)
    {
	if (errmsg.empty())
	{
	    errmsg="Error checking for other errors.";
	}
	return false;
    }
    if (nerr)
    {
	errmsg="Another process experienced an error delivering variables.";
	return false;
    }
    return !error;
#else
    return true;
#endif    
}



// Resize the vector so it holds the current number of variables
// walk along all variables in the maps and determine, what our interest is in each one
void SubWorld::getVariableStatus(std::vector<char>& vb)
{
    vb.resize(reducemap.size());
    str2reduce::const_iterator rm=reducemap.begin();
    str2bool::const_iterator im=importmap.begin();
    for (uint i=0;rm!=reducemap.end();++rm,++it,++i)
    {
	if (rm->second->hasValue())
	{
	    vb[i]=rs::HAVE;	// we have a value to contribute
	}
	else if (im->second)
	{
	    vb[i]=rs::INTERESTED;	// we are interested in this variable
	}
	else
	{
	    vb[i]=rs::NONE;	// we have no interest in this variable
	}
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

size_t SubWorld::getNumVars()
{
    return reducemap.size();
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

