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
    :globalcom(global), subcom(MPI_COMM_NULL), localworld((SubWorld*)0), groupcount(numgroups)
{
    int gsize;
    int grank;
    if ((MPI_Comm_size(global, &gsize)!=MPI_SUCCESS) || (MPI_Comm_rank(global, &grank)!=MPI_SUCCESS))
    {
	throw SplitWorldException("MPI appears to be inoperative.");
    }
    if (gsize%numgroups!=0)
    {
	throw SplitWorldException("SplitWorld error: requested number of groups is not a factor of global communicator size.");
    }
    int wsize=gsize/numgroups;	// each world has this many processes
    int res=MPI_Comm_split(MPI_COMM_WORLD, grank/wsize, grank%wsize, &subcom);
    if (res!=MPI_SUCCESS)
    {
	throw SplitWorldException("SplitWorld error: Unable to form communicator.");
    } 
    localworld=SubWorld_ptr(new SubWorld(subcom));
    localid=grank/wsize;
}

// We may need to look into this more closely.
// What if the domain lives longer than the world splitter?
SplitWorld::~SplitWorld()
{
    if (subcom!=MPI_COMM_NULL)
    {
	MPI_Comm_free(&subcom);
    }
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

/** a list of tuples/sequences:  (Job class, number of instances)*/
void SplitWorld::runJobs(boost::python::list l)
{
    // first count up how many jobs we have in total
    unsigned int numjobs=0;
    std::vector<object> classvec;
    std::vector<unsigned int> countvec;
    std::vector<unsigned int> partialcounts;
    for (int i=0;i<len(l);++i)
    {
	extract<tuple> ex(l[i]);
	if (!ex.check())
	{
	    throw SplitWorldException("runJobs takes a list of tuples (jobclass, number).");
	}
	tuple t=ex();
	if (len(t)!=2)
	{
	    throw SplitWorldException("runJobs takes a list of tuples (jobclass, number).");
	}
	extract<unsigned int> ex2(t[1]);
	unsigned int c=0;
	if (!ex2.check() || ((c=ex2())==0))
	{
	    throw SplitWorldException("Number of jobs must be a strictly positive integer.");
	  
	}
	classvec.push_back(t[0]);
	countvec.push_back(c);
	numjobs+=c;
	partialcounts.push_back(numjobs);
    }
    unsigned int classnum=0;
    unsigned int lowend=1;
    unsigned int highend=lowend+numjobs/groupcount+(numjobs%groupcount);
// std::cout << localid << std::endl;
    for (int i=1;i<=localid;++i)
    {
	lowend=highend;
	highend=lowend+numjobs/groupcount;
	if (i<numjobs%groupcount)
	{
	    highend++;
	}
    }
// std::cout << "There are " << numjobs << " jobs with range [" << lowend << ", " << highend << ")\n";    
    // We could do something more clever about trying to fit Jobs to subworlds
    // to ensure that instances sharing the same Job class would share the same 
    // world as much as possible but for now we'll do this:
    for (unsigned int j=1;j<=numjobs;++j)		// job #0 is a sentinel
    {
	if (j>partialcounts[classnum])
	{
	    classnum++;	// we dont' need to loop this because each count >0
	}
	// now if this is one of the job numbers in our local range,
	// create an instance of the appropriate class
	if (j>=lowend and j<highend)
	{  
	    object o=classvec[classnum](localworld->getDomain(), object(j));
	    localworld->addJob(o);
	}
    }
    int mres=0;
    do
    {
	// now we actually need to run the jobs
	// everybody will be executing their localworld's jobs
	int res=localworld->runJobs();	
	// now we find out about the other worlds
	mres=0;
	if (MPI_Allreduce(&res, &mres, 1, MPI_INT, MPI_MAX, globalcom)!=MPI_SUCCESS)
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
	throw SplitWorldException("At least one Job's work() function raised an exception.");
    }
}


void SplitWorld::registerCrate(escript::crate_ptr c)
{
    protocrates.push_back(c);
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

}
