/*****************************************************************************
*
* Copyright (c) 2014-2015 by University of Queensland
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

#include "esysUtils/Esys_MPI.h"
#include "SplitWorld.h"
#include "AbstractDomain.h"
#include "SplitWorldException.h"
#include "SplitWorldException.h"
#include "esysUtils/pyerr.h"

#include <iostream>
#include <sstream>

using namespace boost::python;
using namespace escript;
namespace rs=escript::reducerstatus;


double SplitWorld::getScalarVariable(const std::string& name)
{
    // do we have a variable of that name?
    return localworld->getScalarVariable(name);
}

SplitWorld::SplitWorld(unsigned int numgroups, MPI_Comm global)
    :localworld((SubWorld*)0), swcount(numgroups>0?numgroups:1), jobcounter(1), manualimport(false)
{
    globalcom=esysUtils::makeInfo(global);
    
    int grank=0;
    int wsize=1;		// each world has this many processes
    esysUtils::JMPI subcom;	// communicator linking other processes in this subworld
    esysUtils::JMPI corrcom;	// communicator linking corresponding processes in different subworlds
    #ifdef ESYS_MPI
	int gsize=globalcom->size;
	grank=globalcom->rank;
	if (gsize%swcount!=0)
	{
	    throw SplitWorldException("SplitWorld error: requested number of groups is not a factor of global communicator size.");
	}
	wsize=gsize/swcount;	// each world has this many processes
	MPI_Comm sub;
	int res=MPI_Comm_split(global, grank/wsize, grank%wsize, &sub);
	if (res!=MPI_SUCCESS)
	{
	    throw SplitWorldException("SplitWorld error: Unable to form communicator.");
	}
	subcom=esysUtils::makeInfo(sub,true);
	
	
	MPI_Comm corrsub;
	res=MPI_Comm_split(global, grank%wsize, grank/wsize, &corrsub);
	if (res!=MPI_SUCCESS)
	{
	    throw SplitWorldException("SplitWorld error: Unable to form communicator.");
	}
	corrcom=esysUtils::makeInfo(corrsub,true);
	
    #else
	subcom=esysUtils::makeInfo(0);
	corrcom=esysUtils::makeInfo(0);
    #endif
    localworld=SubWorld_ptr(new SubWorld(globalcom, subcom,corrcom, swcount, grank%wsize,manualimport));
    localid=grank/wsize;
}

SplitWorld::~SplitWorld()
{
    // communicator cleanup handled by the MPI_Info
}


// The boost wrapper will ensure that there is at least one entry in the tuple
object SplitWorld::buildDomains(tuple t, dict kwargs)
{
    int tsize=len(t);
    // get the callable that we will invoke in a sec
    object tocall=t[0];
    // make a new tuple without the first element
    tuple ntup=tuple(t.slice(1,tsize));
    // now add the subworld to the kwargs
    kwargs["escriptworld"]=localworld;

    // pass the whole package to the python call
    object dobj=tocall(*ntup, **kwargs);
    extract<Domain_ptr> ex1(dobj);
    Domain_ptr dptr=ex1();
    
    // now do a sanity check to see if the domain has respected the communicator info we passed it.
    if (dptr->getMPIComm()!=localworld->getMPI()->comm)
    {
	throw SplitWorldException("The newly constructed domain is not using the correct communicator.");
    }
    localworld->setDomain(dptr);
    return object();	// return None
}


// Executes all pending jobs on all subworlds
void SplitWorld::runJobs()
{
    esysUtils::NoCOMM_WORLD ncw;	// it's destructor will unset the flag
    localworld->resetInterest();    
    try 
    {
	distributeJobs();
	int mres=0;
	std::string err;
	std::vector<char> variableinterest(localworld->getNumVars(), reducerstatus::NONE);
	do
	{  
		// find out which variables each world wants
	    if (!localworld->synchVariableInfo(err))
	    {
		mres=4;
		break;
	    }
		// distribute values to worlds as needed
	    if (!localworld->synchVariableValues(err))
	    {
		mres=4;
		break;
	    }
		// give values to jobs
	    if (!localworld->deliverImports(err))
	    {
		mres=4;
		break;
	    }
	    // now we actually need to run the jobs
	    // everybody will be executing their localworld's jobs
	    int res=localworld->runJobs(err);	
	    // now we find out about the other worlds
	    if (!esysUtils::checkResult(res, mres, globalcom))
	    {
		throw SplitWorldException("MPI appears to have failed.");
	    }
	    if (mres>1)	// 1 and 0 are normal returns, >1 is some sort of error
	    {
	      break; 
	    }
	    if (!localworld->localTransport(err))
	    {
		mres=4;
		break;
	    }
	} while (mres==1);
	
	localworld->clearJobs();
	  // at this point, the remote world has all the reductions done
	  // now we need to do the global merges
	if (!localworld->checkRemoteCompatibility(err))
	{
	    mres=4;
	    err="Error in checkRemoteCompatibility.";
	}
	if (mres==0)	
	{
	    return;
	}
	else if (mres==2)
	{
	    throw SplitWorldException("At least one Job's work() function did not return True/False.");
	}
	else if (mres==3)
	{
	    char* resultstr=0;
	    // now we ship around the error message - This should be safe since
	    // eveyone must have finished their Jobs to get here
	    if (!esysUtils::shipString(err.c_str(), &resultstr, globalcom->comm))
	    {
		throw SplitWorldException("MPI appears to have failed.");
	    }
	    //throw SplitWorldException("At least one Job's work() function raised an exception.");
	    std::string s("At least one Job's work() function raised the following exception:\n");
	    s+=resultstr;
	    throw SplitWorldException(s);
	}
	else if (mres==4)
	{
	    throw SplitWorldException("While processing exports: "+err);
	}
	else
	{ 
	    throw SplitWorldException("Unexpected return value from runJobs.");
	}

    }
    catch (SplitWorldException e )
    {
	clearAllJobs();
	throw e;
    }
}

/**
  stores the constructor/factory to make Jobs and the parameters.
*/
void SplitWorld::addJob(boost::python::object creator, boost::python::tuple tup, boost::python::dict kw)
{
    create.push_back(creator);
    tupargs.push_back(tup);
    kwargs.push_back(kw);  
}

// At some point, we may need there to be more isolation here
// and trap any python exceptions etc, but for now I'll just call the constructor
void SplitWorld::addVariable(std::string name, boost::python::object creator, boost::python::tuple ntup, boost::python::dict kwargs)
{
    object red=creator(*ntup, **kwargs);
    extract<Reducer_ptr> ex(red);
    if (!ex.check())
    {
	throw SplitWorldException("Creator function did not produce a reducer.");
    }
    Reducer_ptr rp=ex();
    localworld->addVariable(name, rp);
}


void SplitWorld::removeVariable(std::string name)
{
    localworld->removeVariable(name);
}

void SplitWorld::clearAllJobs()
{
    clearPendingJobs();
    localworld->clearJobs();
}

void SplitWorld::clearPendingJobs()
{
    create.clear();
    tupargs.clear();
    kwargs.clear();
}

// All the job params are known on all the ranks.
void SplitWorld::distributeJobs()
{
    std::string errmsg;
    unsigned int numjobs=create.size()/swcount;
    unsigned int start=create.size()/swcount*localid;
    if (localid<create.size()%swcount)
    {
	numjobs++;
	start+=localid;
    }
    else
    {
	start+=create.size()%swcount;
    }
    int errstat=0;
    try
    {
	// No other subworld will be looking at this portion of the array
	// so jobs will only be created on one subworld
	for (unsigned int i=start;i<start+numjobs;++i)
	{
	    // we need to add some things to the kw map
	    kwargs[i]["domain"]=localworld->getDomain();
	    kwargs[i]["jobid"]=object(jobcounter+i);
	    object job=create[i](*(tupargs[i]), **(kwargs[i]));
	    localworld->addJob(job);
	}
    }
    catch (boost::python::error_already_set e)
    {
	errstat=1;
	getStringFromPyException(e, errmsg);
    }
    jobcounter+=create.size();
    clearPendingJobs();
    
    // MPI check to ensure that it worked for everybody
    int mstat=0;
    if (!esysUtils::checkResult(errstat, mstat, globalcom))
    {
	throw SplitWorldException("MPI appears to have failed.");
    }
    
      // Now we need to find out if anyone else had an error      
    if (!esysUtils::checkResult(errstat, mstat, globalcom))
    {
	throw SplitWorldException("MPI appears to have failed.");
    }
    errstat=mstat;  
      
    if (errstat==1)
    {
	char* resultstr=0;
	// now we ship around the error message - This should be safe since
	// eveyone must have finished their Jobs to get here
	if (!esysUtils::shipString(errmsg.c_str(), &resultstr, globalcom->comm))
	{
	    throw SplitWorldException("MPI appears to have failed.");
	}      
	throw SplitWorldException(std::string("(During Job creation/distribution) ")+resultstr);
    }
}


namespace escript
{

boost::python::object raw_buildDomains(boost::python::tuple t, boost::python::dict kwargs)
{
    int l=len(t);
    if (l<2)
    {
	throw SplitWorldException("Insufficient parameters to buildDomains.");
    }
    extract<SplitWorld&> exw(t[0]);
    if (!exw.check())
    {
	throw SplitWorldException("First parameter to buildDomains must be a SplitWorld.");
    }
    SplitWorld& ws=exw();
    tuple ntup=tuple(t.slice(1,l));	// strip off the object param
    return ws.buildDomains(ntup, kwargs);
}

boost::python::object raw_addJob(boost::python::tuple t, boost::python::dict kwargs)
{
    int l=len(t);
    if (l<2)
    {
	throw SplitWorldException("Insufficient parameters to addJob.");
    }
    extract<SplitWorld&> exw(t[0]);
    if (!exw.check())
    {
	throw SplitWorldException("First parameter to addJob must be a SplitWorld.");
    }
    SplitWorld& ws=exw();
    object creator=t[1]; 
    tuple ntup=tuple(t.slice(2,l));	// strip off the object param
    ws.addJob(creator, ntup, kwargs);
    return object();
}

// expects, splitworld, name of var, constructor function for the reducer, any constructor params
boost::python::object raw_addVariable(boost::python::tuple t, boost::python::dict kwargs)
{
    int l=len(t);
    if (l<3)
    {
	throw SplitWorldException("Insufficient parameters to addReducer.");
    }
    extract<SplitWorld&> exw(t[0]);
    if (!exw.check())
    {
	throw SplitWorldException("First parameter to addVariable must be a SplitWorld.");
    }
    SplitWorld& ws=exw();
    object pname=t[1];
    extract<std::string> ex2(pname);
    if (!ex2.check())
    {
	throw SplitWorldException("Second parameter to addVariable must be a string");
    }
    std::string name=ex2();
    object creator=t[2]; 
    tuple ntup=tuple(t.slice(3,l));	// strip off the object param
    ws.addVariable(name, creator, ntup, kwargs);
    return object();
}



}
