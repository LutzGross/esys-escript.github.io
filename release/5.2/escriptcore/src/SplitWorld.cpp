/*****************************************************************************
*
* Copyright (c) 2014-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include "SplitWorld.h"

#include "AbstractDomain.h"
#include "SplitWorldException.h"
#include "pyerr.h"

#include <sstream>

using namespace boost::python;
using namespace escript;
namespace rs=escript::reducerstatus;


DataTypes::real_t SplitWorld::getScalarVariable(const std::string& name)
{
    // do we have a variable of that name?
    return localworld->getScalarVariable(name);
}

boost::python::object SplitWorld::getLocalObjectVariable(const std::string& name)
{
    // do we have a variable of that name?
    return localworld->getLocalObjectVariable(name);  
}


SplitWorld::SplitWorld(unsigned int numgroups, MPI_Comm global)
    : localworld((SubWorld*)0), swcount(numgroups > 0 ? numgroups : 1),
      jobcounter(1), manualimport(false)
{
    globalcom = makeInfo(global);
    
    int grank=0;
    int wsize=1;		// each world has this many processes
    JMPI subcom;	// communicator linking other processes in this subworld
    JMPI corrcom;	// communicator linking corresponding processes in different subworlds
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
	subcom=makeInfo(sub,true);
	
	
	MPI_Comm corrsub;
	res=MPI_Comm_split(global, grank%wsize, grank/wsize, &corrsub);
	if (res!=MPI_SUCCESS)
	{
	    throw SplitWorldException("SplitWorld error: Unable to form communicator.");
	}
	corrcom=makeInfo(corrsub,true);
	
    #else
	if (numgroups!=1)
	{
	    throw SplitWorldException("SplitWorld error: non-MPI builds can only create 1 subworld.");
	  
	}
	subcom=makeInfo(0);
	corrcom=makeInfo(0);
    #endif
    localworld=SubWorld_ptr(new SubWorld(globalcom, subcom,corrcom, swcount, grank/wsize,manualimport));
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
    NoCOMM_WORLD ncw;	// it's destructor will unset the flag
    localworld->resetInterest();  
    localworld->newRunJobs();
    try 
    {
	distributeJobs();
	int mres=0;
	std::string err;
	do	// only here so I can "break" to the end
	{
		// find out which variables each world has and wants  (global op)
	    if (!localworld->synchVariableInfo(err))
	    {
		mres=4;
		break;
	    }
		// distribute values to worlds as needed  (global op)
	    if (!localworld->synchVariableValues(err))
	    {
		mres=4;
		break;
	    }
		// give values to jobs (local op)	
	    if (!localworld->deliverImports(err))
	    {
		mres=4;	// can't jump to the end because this is a local op
	    }
	    else	// import delivery was successful
	    {
	        // now we actually need to run the jobs
	        // everybody will be executing their localworld's jobs
	        mres=localworld->runJobs(err);	

                if (mres<2)
	        {
                    if (!localworld->localTransport(err))
		    {
		        mres=4;		// both running jobs and local reduction are local ops
		    }
	        }
	    }
	} while (false);
        int res=mres;
        // now we find out about the other worlds
        if (!checkResult(res, mres, globalcom))
        {
	    throw SplitWorldException("MPI appears to have failed.");
        }
	localworld->clearJobs();
	  // at this point, the remote world has all the reductions done
	  // now we need to do the global merges
	if (!localworld->checkRemoteCompatibility(err))
	{
	    mres=4;
	    err=std::string("Error in checkRemoteCompatibility. ")+err;
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
	    if (!shipString(err.c_str(), &resultstr, globalcom->comm))
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

/**
  creates exactly one instance of the job on each world.
  This bypasses the normal job allocation method (since that does not guarantee one job per world)
*/
void SplitWorld::addJobPerWorld(boost::python::object creator, boost::python::tuple tup, boost::python::dict kw)
{
    std::string errmsg;
    int errstat=0;
    try
    {
	// we need to add some things to the kw map
	kw["domain"]=localworld->getDomain();
	kw["jobid"]=object(jobcounter+localid);
	kw["swcount"]=object(swcount);
	kw["swid"]=object(localid);
	object job=creator(*tup, **kw);
	localworld->addJob(job);
    }
    catch (boost::python::error_already_set e)
    {
	errstat=1;
	getStringFromPyException(e, errmsg);
    }
    jobcounter+=swcount;
    clearPendingJobs();
    
    // MPI check to ensure that it worked for everybody
    int mstat=0;
    if (!checkResult(errstat, mstat, globalcom))
    {
	throw SplitWorldException("MPI appears to have failed.");
    }
    
      // Now we need to find out if anyone else had an error      
    if (!checkResult(errstat, mstat, globalcom))
    {
	throw SplitWorldException("MPI appears to have failed.");
    }
    errstat=mstat;  
      
    if (errstat==1)
    {
	char* resultstr=0;
	// now we ship around the error message - This should be safe since
	// eveyone must have finished their Jobs to get here
	if (!shipString(errmsg.c_str(), &resultstr, globalcom->comm))
	{
	    throw SplitWorldException("MPI appears to have failed.");
	}      
	throw SplitWorldException(std::string("(During Job creation/distribution) ")+resultstr);
    }  
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


void SplitWorld::clearVariable(std::string name)
{
    localworld->clearVariable(name);
}

std::list<std::pair<std::string, bool> > SplitWorld::getVarList()
{
    return localworld->getVarList();
}


boost::python::object SplitWorld::getVarPyList()
{
    std::list<std::pair<std::string, bool> > r=localworld->getVarList();
    boost::python::list l;
    for (std::list<std::pair<std::string, bool> >::iterator it=r.begin();it!=r.end();++it)
    {
        boost::python::list t;
	t.append(it->first);
	t.append(it->second);
        l.append(t);
    }
    return l;
}

// We make no committments that the format of the output will remain the same
boost::python::object SplitWorld::getVarPyInfo()
{
    std::list<std::pair<std::string, std::string> > r=localworld->getVarInfo();
    boost::python::list l;
    for (std::list<std::pair<std::string, std::string> >::iterator it=r.begin();it!=r.end();++it)
    {
        boost::python::list t;
	t.append(it->first);
	t.append(it->second);
        l.append(t);
    }
    return l;
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
// note that perWorld jobs will already be loaded directly to the worlds
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
	    kwargs[i]["swcount"]=object(swcount);
	    kwargs[i]["swid"]=object(localid);    
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
    if (!checkResult(errstat, mstat, globalcom))
    {
	throw SplitWorldException("MPI appears to have failed.");
    }
    
      // Now we need to find out if anyone else had an error      
    if (!checkResult(errstat, mstat, globalcom))
    {
	throw SplitWorldException("MPI appears to have failed.");
    }
    errstat=mstat;  
      
    if (errstat==1)
    {
	char* resultstr=0;
	// now we ship around the error message - This should be safe since
	// eveyone must have finished their Jobs to get here
	if (!shipString(errmsg.c_str(), &resultstr, globalcom->comm))
	{
	    throw SplitWorldException("MPI appears to have failed.");
	}      
	throw SplitWorldException(std::string("(During Job creation/distribution) ")+resultstr);
    }
}

int SplitWorld::getSubWorldCount()
{
    return swcount;
}

int SplitWorld::getSubWorldID()
{
    return localid;
}

void SplitWorld::copyVariable(const std::string& src, const std::string& dest)
{
    if (manualimport)
    {
	throw SplitWorldException("copyVariable is not yet supported for manualimport.");
	  // we would need to make sure that the value is on this world
	  // so we need to do a value transport here
    }
    else
    {
        localworld->copyVariable(src, dest);
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

boost::python::object raw_addJobPerWorld(boost::python::tuple t, boost::python::dict kwargs)
{
    int l=len(t);
    if (l<2)
    {
	throw SplitWorldException("Insufficient parameters to addJobPerWorld.");
    }
    extract<SplitWorld&> exw(t[0]);
    if (!exw.check())
    {
	throw SplitWorldException("First parameter to addJobPerWorld must be a SplitWorld.");
    }
    SplitWorld& ws=exw();
    object creator=t[1]; 
    tuple ntup=tuple(t.slice(2,l));	// strip off the object param
    ws.addJobPerWorld(creator, ntup, kwargs);
    return object();
}


// expects, splitworld, name of var, constructor function for the reducer, any constructor params
boost::python::object raw_addVariable(boost::python::tuple t, boost::python::dict kwargs)
{
    int l=len(t);
    if (l<3)
    {
	throw SplitWorldException("Insufficient parameters to addVariable.");
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

