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
    :globalcom(global), subcom(MPI_COMM_NULL), localworld((SubWorld*)0), swcount(numgroups>0?numgroups:1), jobcounter(1)
{
    int grank=0;
    int wsize=1;		// each world has this many processes
    #ifdef ESYS_MPI
	int gsize=1;
	if ((MPI_Comm_size(global, &gsize)!=MPI_SUCCESS) || (MPI_Comm_rank(global, &grank)!=MPI_SUCCESS))
	{
	    throw SplitWorldException("MPI appears to be inoperative.");
	}
	if (gsize%swcount!=0)
	{
	    throw SplitWorldException("SplitWorld error: requested number of groups is not a factor of global communicator size.");
	}
	wsize=gsize/swcount;	// each world has this many processes
	int res=MPI_Comm_split(MPI_COMM_WORLD, grank/wsize, grank%wsize, &subcom);
	if (res!=MPI_SUCCESS)
	{
	    throw SplitWorldException("SplitWorld error: Unable to form communicator.");
	}
    #endif
    localworld=SubWorld_ptr(new SubWorld(subcom));
    localid=grank/wsize;
}

// We may need to look into this more closely.
// What if the domain lives longer than the world splitter?
SplitWorld::~SplitWorld()
{
#ifdef ESYS_MPI  
    if (subcom!=MPI_COMM_NULL)
    {
	MPI_Comm_free(&subcom);
    }
#endif    
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
    if (dptr->getMPIComm()!=localworld->getComm())
    {
	throw SplitWorldException("The newly constructed domain is not using the correct communicator.");
    }
    localworld->setDomain(dptr);
    return object();	// return None
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
	// now we find out about the other worlds
	if (!esysUtils::checkResult(res, mres, globalcom))
	{
	    throw SplitWorldException("MPI appears to have failed.");
	}
std::cerr << "I got a res of " << mres << std::endl;	
	
    } while (mres==1);
    if (mres==2)
    {
	throw SplitWorldException("At least one Job's work() function did not return True/False.");
    }
    else if (mres==3)
    {
std::cerr << "My err string is [" << err <<"]\n";	
      
	char* resultstr=0;
	// now we ship around the error message
	if (!esysUtils::shipString(err.c_str(), &resultstr, globalcom))
	{
	    throw SplitWorldException("MPI appears to have failed.");
	}
	//throw SplitWorldException("At least one Job's work() function raised an exception.");
	std::string s("At least one Job's work() function raised the following exception:\n");
	s+=resultstr;
// std::cerr << "My combined [" << s.c_str() << std::endl;
// 	char* testing=new char[s.size()+1];
// 	strcpy(testing, s.c_str());
std::cerr << "Pre-throw [[[" << s << "]]]\n";	
	throw SplitWorldException(s);
    }
}


/** a list of tuples/sequences:  (Job class, number of instances)*/
// void SplitWorld::runJobs(boost::python::list l)
// {
//     // first count up how many jobs we have in total
//     unsigned int numjobs=0;
//     std::vector<object> classvec;
//     std::vector<unsigned int> countvec;
//     std::vector<unsigned int> partialcounts;
//     for (int i=0;i<len(l);++i)
//     {
// 	extract<tuple> ex(l[i]);
// 	if (!ex.check())
// 	{
// 	    throw SplitWorldException("runJobs takes a list of tuples (jobclass, number).");
// 	}
// 	tuple t=ex();
// 	if (len(t)!=2)
// 	{
// 	    throw SplitWorldException("runJobs takes a list of tuples (jobclass, number).");
// 	}
// 	extract<unsigned int> ex2(t[1]);
// 	unsigned int c=0;
// 	if (!ex2.check() || ((c=ex2())==0))
// 	{
// 	    throw SplitWorldException("Number of jobs must be a strictly positive integer.");
// 	  
// 	}
// 	classvec.push_back(t[0]);
// 	countvec.push_back(c);
// 	numjobs+=c;
// 	partialcounts.push_back(numjobs);
//     }
//     unsigned int classnum=0;
//     unsigned int lowend=1;
//     unsigned int highend=lowend+numjobs/swcount+(numjobs%swcount);
// // std::cout << localid << std::endl;
//     for (int i=1;i<=localid;++i)
//     {
// 	lowend=highend;
// 	highend=lowend+numjobs/swcount;
// 	if (i<numjobs%swcount)
// 	{
// 	    highend++;
// 	}
//     }
// // std::cout << "There are " << numjobs << " jobs with range [" << lowend << ", " << highend << ")\n";    
//     // We could do something more clever about trying to fit Jobs to subworlds
//     // to ensure that instances sharing the same Job class would share the same 
//     // world as much as possible but for now we'll do this:
//     for (unsigned int j=1;j<=numjobs;++j)		// job #0 is a sentinel
//     {
// 	if (j>partialcounts[classnum])
// 	{
// 	    classnum++;	// we dont' need to loop this because each count >0
// 	}
// 	// now if this is one of the job numbers in our local range,
// 	// create an instance of the appropriate class
// 	if (j>=lowend and j<highend)
// 	{  
// 	    object o=classvec[classnum](localworld->getDomain(), object(j));
// 	    localworld->addJob(o);
// 	}
//     }
//     int mres=0;
//     do
//     {
// 	// now we actually need to run the jobs
// 	// everybody will be executing their localworld's jobs
// 	int res=localworld->runJobs();	
// 	// now we find out about the other worlds
// 	mres=0;
// 	if (MPI_Allreduce(&res, &mres, 1, MPI_INT, MPI_MAX, globalcom)!=MPI_SUCCESS)
// 	{
// 	    throw SplitWorldException("MPI appears to have failed.");
// 	}
//     } while (mres==1);
//     if (mres==2)
//     {
// 	throw SplitWorldException("At least one Job's work() function did not return True/False.");
//     }
//     else if (mres==3)
//     {
// 	throw SplitWorldException("At least one Job's work() function raised an exception.");
//     }
// }


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
    if (!esysUtils::checkResult(errstat, mstat, globalcom))
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
