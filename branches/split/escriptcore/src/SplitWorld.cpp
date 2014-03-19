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

#include "esysUtils/Esys_MPI.h"
#include "SplitWorld.h"
#include "AbstractDomain.h"
#include "SplitWorldException.h"
#include "SplitWorldException.h"

#include <iostream>

using namespace boost::python;
using namespace escript;

SplitWorld::SplitWorld(unsigned int numgroups, MPI_Comm global)
    :localworld((SubWorld*)0), swcount(numgroups>0?numgroups:1), jobcounter(1)
{
    globalcom=esysUtils::makeInfo(global);
    
    int grank=0;
    int wsize=1;		// each world has this many processes
    #ifdef ESYS_MPI
	int gsize=globalcom->size;
	grank=globalcom->rank;
	if (gsize%swcount!=0)
	{
	    throw SplitWorldException("SplitWorld error: requested number of groups is not a factor of global communicator size.");
	}
	wsize=gsize/swcount;	// each world has this many processes
	MPI_Comm sub;
	int res=MPI_Comm_split(MPI_COMM_WORLD, grank/wsize, grank%wsize, &sub);
	if (res!=MPI_SUCCESS)
	{
	    throw SplitWorldException("SplitWorld error: Unable to form communicator.");
	}
	subcom=esysUtils::makeInfo(sub,true);
    #else
	subcom=esysUtils::makeInfo(0);
    #endif
    localworld=SubWorld_ptr(new SubWorld(subcom));
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

// std::cerr << "About to call function with:\n";
// //extract<std::string> ex(ntup.attr("__str__")());
//  std::cerr << extract<std::string>(ntup.attr("__str__")())() << std::endl;
// // for (int i=0;i<tsize-1;++i)
// // {
// //     std::cout << extract<const char*>(ntup[i])() << " ";
// // }
// std::cerr << std::endl;
    
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


namespace
{

// Throw all values in and get the maximum
// This is not an AllReduce because we need to use Tagged Communication so as not to interfere with 
// other messages / collective ops which the SubWorld's Jobs might be using
// This is implemented by sending to rank 0
bool checkResultInt(int res, int& mres, esysUtils::JMPI& info)
{
    const int leader=0;
    const int BIGTAG=esysUtils::getSubWorldTag();
    if (info->size==1)
    {
	mres=res;
	return true;
    }
    else
    {
	if (info->rank!=leader)
	{
	    if (MPI_Send(&res, 1, MPI_INT, leader, BIGTAG, info->comm)!=MPI_SUCCESS)
	    {
		return false;
	    }
	    if (MPI_Recv(&mres, 1, MPI_INT, leader, BIGTAG, info->comm,0)!=MPI_SUCCESS)
	    {
		return false;
	    }
	}
	else
	{
	    MPI_Request* reqs=new MPI_Request[info->size-1];
	    int* eres=new int[info->size-1];
	    for (int i=0;i<info->size-1;++i)
	    {
		MPI_Irecv(eres+i, 1, MPI_INT, i+1, BIGTAG, info->comm, reqs+i);	  
	    }
	    if (MPI_Waitall(info->size-1, reqs, 0)!=MPI_SUCCESS)
	    {
		delete[] reqs;
		return false;
	    }
	    // now we have them all, find the max
	    mres=res;
	    for (int i=0;i<info->size-1;++i)
	    {
		if (mres<eres[i])
		{
		    mres=eres[i];
		}
	    }
	    // now we know what the result should be
	    // send it to the others
	    for (int i=0;i<info->size-1;++i)
	    {
		MPI_Isend(&mres, 1, MPI_INT, i+1, BIGTAG, info->comm, reqs+i);	  
	    }
	    if (MPI_Waitall(info->size-1, reqs,0)!=MPI_SUCCESS)
	    {
		return false;
	    }
	}
      
    }
    return true;
}


}

// Executes all pending jobs on all subworlds
void SplitWorld::runJobs()
{
    distributeJobs();
    int mres=0;
    std::string err;
    do
    {
	// now we actually need to run the jobs
	// everybody will be executing their localworld's jobs
	int res=localworld->runJobs(err);	
std::cerr << "Done local jobs" << std::endl;	
	// now we find out about the other worlds
	if (!checkResultInt(res, mres, globalcom))
	{
	    throw SplitWorldException("MPI appears to have failed.");
	}
    } while (mres==1);
    if (mres==2)
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
}

void SplitWorld::registerCrate(escript::crate_ptr c)
{
    protocrates.push_back(c);
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

void SplitWorld::clearPendingJobs()
{
    create.clear();
    tupargs.clear();
    kwargs.clear();
}

void SplitWorld::clearActiveJobs()
{
    localworld->clearJobs();
}

// All the job params are known on all the ranks.
void SplitWorld::distributeJobs()
{
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
std::cerr << "Numjobs=" << numjobs << " start=" << start << std::endl;    
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
    }
    jobcounter+=create.size();
    clearPendingJobs();
    
    // MPI check to ensure that it worked for everybody
    int mstat=0;
    if (!esysUtils::checkResult(errstat, mstat, globalcom->comm))
    {
	throw SplitWorldException("MPI appears to have failed.");
    }
    
    if (errstat==1)
    {
	throw SplitWorldException("distributeJobs: Job creation failed.");
	clearActiveJobs();
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


}
