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
#include <sstream>

using namespace boost::python;
using namespace escript;
namespace rs=escript::reducerstatus;

SplitWorld::SplitWorld(unsigned int numgroups, MPI_Comm global)
    :localworld((SubWorld*)0), swcount(numgroups>0?numgroups:1), jobcounter(1), manualimport(false)
{
    globalcom=esysUtils::makeInfo(global);
    
    int grank=0;
    int wsize=1;		// each world has this many processes
    esysUtils::JMPI subcom;	// all procs in the local subworld
    esysUtils::JMPI corcom;	// all procs in corresponding ranks in each subworld
				// eg: all procs having rank 0 in their world
    #ifdef ESYS_MPI
	int gsize=globalcom->size;
	grank=globalcom->rank;
	if (gsize%swcount!=0)
	{
	    throw SplitWorldException("SplitWorld error: requested number of groups is not a factor of global communicator size.");
	}
	wsize=gsize/swcount;	// each world has this many processes
	int lrank=grank%gsize;
	MPI_Comm sub;
	int res=MPI_Comm_split(MPI_COMM_WORLD, grank/wsize, grank%wsize, &sub);
	if (res!=MPI_SUCCESS)
	{
	    throw SplitWorldException("SplitWorld error: Unable to form communicator.");
	}
	subcom=esysUtils::makeInfo(sub,true);
	
	  // Now we need to make a communicator which links all the world leaders
	  // first we need to make a group with the correct procs in it
	MPI_Group globalg;	// group for global communicator
	MPI_Group corrg;
	MPI_Comm lcom;
	int* mem=new int[numgroups];
	for (int i=0;i<numgroups;++i)
	{
	    mem[i]=i*wsize+lrank;	  
	}
	scoped_array<int> memg(mem);	// so we get it cleaned up if an exception is thrown
	if (MPI_Comm_group(global, &globalg)!=MPI_SUCCESS) 
	{
	    throw SplitWorldException("SplitWorld error: Unable to form group.");
	}
	if (MPI_Group_incl(globalg, numgroups, memg.get(), &corrg)!=MPI_SUCCESS) 
	{
	    throw SplitWorldException("SplitWorld error: Unable to form group.");
	}
	  // now we have the group we can build the corresponding communicator
	if (MPI_Comm_create(global, corrg, &lcom)!=MPI_SUCCESS)
	{
	    MPI_Group_free(corrg);
	    throw SplitWorldException("SplitWorld error: Unable to form leader communicator.");
	}
	MPI_Group_free(corrg);	
	  // we have the communicator, 
	corcom=esysUtils::makeInfo(lcom, true);	
    #else
	subcom=esysUtils::makeInfo(0);
	corcom=esysUtils::makeInfo(0);
    #endif
    localworld=SubWorld_ptr(new SubWorld(subcom, corcom));
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


namespace
{

// Throw all values in and get the maximum
// This is not an AllReduce because we need to use Tagged Communication so as not to interfere with 
// other messages / collective ops which the SubWorld's Jobs might be using
// This is implemented by sending to rank 0
bool checkResultInt(int res, int& mres, esysUtils::JMPI& info)
{
#ifdef ESYS_MPI
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
#else
    mres=res;
    return true;
#endif
}


}


// Update the vector to indicate the interest of all subworlds in each variable
bool SplitWorld::getVariableInterest(std::vector<char>& vb)
{

	std::vector<char> lvec;
	localworld->getVariableStatus(lvec);	// this tells us what our local world wants  
	
	vb.resize(localworld->getNumVars()*swcount);

      // now we need to ask all the world leaders about their involvement
#ifdef ESYS_MPI
    if (localworld->amLeader())
    {
	// The leaders of each world, send their variable information to the proc "0" in
	// the global world (which will be the leader of subworld "0").
	//    There is an issue here if this operation fails
	if (MPI_Gather(&lvec[0], localworld->getNumVars(), MPI_CHAR, &vb[0], vb.size(), 
		   MPI_CHAR, 0, localworld->getCorrMPI()->comm)!=MPI_SUCCESS) 
	{
	    for (size_t i=0;i<vb.size;++i)
	    {
		vb[i]=rs::ERROR;
	    }
	}
    }
    // now share the combined info with all processes
    if ((MPI_Bcast(&vb[0], vb.size, MPI_CHAR, 0, global->getMPI()->comm)!=MPI_SUCCESS)
	  || (vb[0]==rs::ERROR))
    {
	return false;
    }
#else
      // since there is only one world, we just copy lvec to global
    vb=lvec;
#endif  
    return true;
}

// Executes all pending jobs on all subworlds
void SplitWorld::runJobs()
{
    esysUtils::NoCOMM_WORLD ncw;	// it's destructor will unset the flag
    try 
    {
	distributeJobs();
	int mres=0;
	std::string errmsg;
	std::vector<char> impexpdetail;
	std::vector<char> variableinterest;
	do
	{
	      // make sure that any jobs which register as needing imports get them
	      // first check local jobs to find out what they need
	    if (!localworld->findImports(manualimport, errmsg))
	    {
		mres=4;
		break;
	    }
	      // Now we find out what the other worlds want 
	    if (!getVariableInterest(variableinterest))
	    {
		mres=4;
		errmsg="Error while gathering variable use information.";
		break;
	    }
	      // If there are any variables which are needed but not present on local
	    if (!localworld->deliverGlobalImports(variableinterest, errmsg))
	    {
		mres=4;
		break;
	    }
	      
	    if (!localworld->deliverImports(errmsg))
	    {
		mres=4;
		break;
	    }	  
	    // now we actually need to run the jobs
	    // everybody will be executing their localworld's jobs
	    int res=localworld->runJobs(errmsg);	

// 	    // take this opportunity to clean up (the python side)
// 	    localworld->clearImportExports();
	    // now we find out about the other worlds
	    if (!checkResultInt(res, mres, globalcom))
	    {
		throw SplitWorldException("MPI appears to have failed.");
	    }
	    if (mres>1)	// 1 and 0 are normal returns, >1 is some sort of error
	    {
	      break; 
	    }
	    if (!localworld->localTransport(impexpdetail, errmsg))
	    {
		mres=4;
		break;
	    }
	    
	      // at this point, the remote world has all the reductions done
	      // now we need to do the global merges
	      
	    if (!localworld->checkRemoteCompatibility(errmsg))
	    {
		mres=4;
		break;
	    }
	    if (!localworld->reduceRemoteValue(errmsg))
	    {
		mres=4;
		break;	      
	    }	    
	} while (mres==1);
	if (mres==0)
	{
	    clearAllJobs();
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
	    if (!esysUtils::shipString(errmsg.c_str(), &resultstr, globalcom->comm))
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
	    throw SplitWorldException("While processing exports: "+errmsg);
	
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
    localworld->addVariable(name, rp, manualimport);
}


void SplitWorld::removeVariable(std::string name)
{
    localworld->removeVariable(name);
}

void SplitWorld::clearAllJobs()
{
    clearPendingJobs();
    clearActiveJobs();
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
