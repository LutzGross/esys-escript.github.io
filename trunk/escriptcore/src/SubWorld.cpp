
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

#include "SubWorld.h"
#include "SplitWorldException.h"

#include <boost/python/import.hpp>
#include <boost/python/dict.hpp>

#include <iostream>

using namespace escript;
namespace bp=boost::python;
using namespace esysUtils;
namespace rs=escript::reducerstatus;

SubWorld::SubWorld(JMPI& global, JMPI& comm, JMPI& corr, unsigned int subworldcount, unsigned int local_id, bool manualimport)
    :everyone(global), swmpi(comm), corrmpi(corr), domain((AbstractDomain*)0), swcount(subworldcount), localid(local_id), manualimports(manualimport)
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
bool SubWorld::findImports(std::string& errmsg)
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
		    errmsg="Attempt to import variable \""+n+"\". SplitWorld was not told about this variable.";
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
    for (size_t i=0;i<jobvec.size();++i)
    {
	bp::list wanted=bp::extract<bp::list>(jobvec[i].attr("wantedvalues"))();		  
	for (size_t j=0;j<len(wanted);++j)
	{
	    bp::extract<std::string> exs(wanted[j]);	// must have been checked by now
	    std::string n=exs();
	      // now we need to check to see if this value is known
	    str2reduce::iterator it=reducemap.find(n);
	    if (it==reducemap.end())
	    {
		errmsg="Attempt to import variable \""+n+"\". SplitWorld was not told about this variable.";
		return false;
	    }
	    try
	    {
	        jobvec[i].attr("setImportValue")(it->first, reducemap[it->first]->getPyObj());
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
		errmsg=PyBytes_AsString(rr);
		Py_XDECREF(rr);
	#else
		errmsg=PyString_AsString(errobj);
	#endif
		Py_XDECREF(errobj);

		Py_XDECREF(ptype);
		Py_XDECREF(pvalue);
		Py_XDECREF(ptraceback);
		
		return false;
	  } 

	}
    }
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
		errmsg="Attempt to export variable \""+name+"\" with an incompatible value. Using ";
		errmsg+=(it->second)->description();
		return false;
	    }
	    if (!(it->second)->reduceLocalValue(o2, errmsg))
	    {
		return false;	// the error string will be set by the reduceLocalValue
	    }
	    varstate[name]=rs::NEW;
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

double SubWorld::getScalarVariable(const std::string& name)
{
      str2reduce::iterator it=reducemap.find(name);
      if (it==reducemap.end())
      {
	  throw SplitWorldException("No variable of that name.");
      }
      if (dynamic_cast<MPIScalarReducer*>(it->second.get())==0)
      {
	  throw SplitWorldException("Variable is not scalar.");
      }
      return dynamic_cast<MPIScalarReducer*>(it->second.get())->getDouble();
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
	if (! it->second->checkRemoteCompatibility(corrmpi, errmsg))
	{
	    return false;
	}
    }
    return true;
}

#ifdef ESYS_MPI  
  
bool SubWorld::makeComm(MPI_Comm& sourcecom, MPI_Comm& subcom,std::vector<int>& members)
{
      MPI_Group sourceg, g;
      if (MPI_Comm_group(sourcecom, &sourceg)!=MPI_SUCCESS) {return false;}
      if (MPI_Group_incl(sourceg, members.size(), &members[0], &g)!=MPI_SUCCESS) {return false;}
      // then create a communicator with that group
      if (MPI_Comm_create(sourcecom, g, &subcom)) {return false;}
      return true;
}

// a group with NEW nodes at the front and INT and OLDINT at the back
// NONE worlds get an empty communicator
bool SubWorld::makeGroupComm1(MPI_Comm& srccom, int vnum, char mystate, MPI_Comm& com)
{
      if ((mystate==rs::NEW)
	    || (mystate==rs::INTERESTED)
	    || (mystate==rs::OLDINTERESTED))
      {
      // first create a group with [updates, interested and oldinterested in it]
	  std::vector<int> members;
	  for (int i=0+vnum;i<globalvarinfo.size();i+=getNumVars())
	  {
	      // make a vector of the involved procs with New at the front
	      switch (globalvarinfo[i])
	      {
		case rs::NEW:   members.insert(members.begin(), i%getNumVars()); break;
		case rs::INTERESTED:     
		case rs::OLDINTERESTED:
			  members.push_back(i%getNumVars());
			  break;
	      }
	  }		
	  
	  return makeComm(srccom, com, members)==MPI_SUCCESS;
      }
      else	// for people not in involved in the value shipping
      {		// This would be a nice time to use MPI_Comm_create_group
		// but it does not exist in MPI2.1
	  return MPI_Comm_create(srccom, MPI_GROUP_EMPTY, &com)==MPI_SUCCESS;
      }
}

// A group with a single OLD or OLDINT at the front and all the INT worlds 
// following it
bool SubWorld::makeGroupComm2(MPI_Comm& srccom, int vnum, char mystate, MPI_Comm& com)
{
      if ((mystate==rs::OLD)
	    || (mystate==rs::INTERESTED)
	    || (mystate==rs::OLDINTERESTED))
      {
	  // first create a group with [old, interested and oldinterested in it]
	  std::vector<int> members;
	  bool havesrc=false;
	  for (int i=0+vnum;i<globalvarinfo.size();i+=getNumVars())
	  {
	      // make a vector of the involved procs with New at the front
	      switch (globalvarinfo[i])
	      {
		case rs::NEW:   return false;  break;
		case rs::INTERESTED: members.push_back(i%getNumVars());  break;     
		case rs::OLD: 
		case rs::OLDINTERESTED:
			  if (!havesrc)
			  {
			      members.insert(members.begin(), i%getNumVars());
			      havesrc=true;
			  }
			  break;
	      }
	  }		
	  
	  return makeComm(srccom, com, members)==MPI_SUCCESS;
      }
      else	// for people not in involved in the value shipping
      {		// This would be a nice time to use MPI_Comm_create_group
		// but it does not exist in MPI2.1
	  return MPI_Comm_create(srccom, MPI_GROUP_EMPTY, &com)==MPI_SUCCESS;
      }
}

#endif  


bool SubWorld::reduceRemoteValues(std::string& err)
{
#ifdef ESYS_MPI  
  
  
    // There are three possibilities here but since all worlds have the same knowledge
    // we can be sure that they will all make the same choice
    // 1) No updates are required
    // 2) There is a single world with a new value so it can broadcast it
    // 3) There are multiple worlds with updates
  

    
    // need to keep track of which vars have updates
    std::vector<std::string> varswithupdates;
    
    int vnum=0;    
    for (str2reduce::iterator it=reducemap.begin();it!=reducemap.end();++it, ++vnum)
    {
         // check to see if anyone needs it
	int needcount=0; // who wants a new value
	int newcount=0;	// who has a new version
	int oldcount=0;	// who has an old version
	int oldintcount=0;
	
	    // loop over all the instructions about this variable
	for (int i=0+vnum;i<globalvarinfo.size();i+=getNumVars())
	{
	    switch (globalvarinfo[i])
	    {
	      case rs::NEW: newcount++; break;
	      case rs::INTERESTED: needcount++; break;
	      case rs::OLD: oldcount++;	      
	      case rs::OLDINTERESTED: 
			needcount++;
			oldintcount++;
			break;
	    }
	}
	if (newcount>0)
	{
	    varswithupdates.push_back(it->first);
	}
	if (needcount+newcount+oldcount==0)
	{
	    continue;		// noone cares about this variable
	}
	if (needcount>0 && (oldcount+newcount)==0)
	{
	    err="Job has asked to import a variable with no value.";
	    return false;
	}
	    // worlds have the variable but noone is interested in it
	if (needcount==0)
	{
	    continue;
	}
	if (swcount==1)
	{		// nobody else to communicate with
	    continue;
	}
	    // to reach this point, there must be >=1 source and >=1 sink and multiple worlds
	    // first deal updates as source(s)
	if (newcount==1)	// only one update so send from that
	{
	    MPI_Comm com;
	    if (!makeGroupComm1(corrmpi->comm, vnum, varstate[it->first],com))
	    {
		err="Error creating group for sharing values,";
		return false;
	    }
	    it->second->groupSend(com);
	    // disolve the group
	    MPI_Comm_free(&com);
	    continue;
	}
	if (newcount>1)
	{
	    // form a group to send to [updates and interested and oldinterested]
	    MPI_Comm com;
	    if (!makeGroupComm1(corrmpi->comm, vnum, varstate[it->first],com))
	    {
		err="Error creating group for sharing values,";
		return false;
	    }
	    it->second->groupReduce(com,varstate[it->first]);
	    MPI_Comm_free(&com);	    
	    continue;
	}
	    // at this point, we need to ship info around but there are no updates
	    // that is, we are shipping an old copy
	    // picking a source arbitarily (the last one)
	    
	    // but first, eliminate the special case where the only interested ones
	    // already have a copy
	if (oldintcount==needcount)
	{
	    continue;
	}
	MPI_Comm com;
	if (!makeGroupComm2(corrmpi->comm, vnum, varstate[it->first],com))
	{
	    err="Error creating group for sharing values";
	    return false;
	}
	// form group to send to [latestsource and interested]
	it->second->groupSend(com);
	// dissolve the group	
	MPI_Comm_free(&com);
    }
	// now we need to age any out of date copies of vars
    for (size_t i=0;i<varswithupdates.size();++i)
    {
        std::string vname=varswithupdates[i];
        	if (varstate[vname]==rs::NEW)
	{
	    varstate[vname]=rs::OLD;
	}
	else if (varstate[vname]==rs::OLD)
	{
	    varstate[vname]=rs::NONE;
	    reducemap[vname]->clear();
	} 
    }
#endif    
    return true;    
}

bool SubWorld::amLeader()
{
    return swmpi->rank==0;
}


// Look at the vector of interest and send around variables to where they need to be
// To try to make error reporting consistant, every variable transfer will be attempted,
// even if an early one fails.
bool SubWorld::deliverGlobalImports(std::vector<char>& vb, std::string& errmsg)
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
	    if (vb[s*nv+i]==rs::NEW)	// This world has a copy of the value
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
	    if (!rmi->second->recvFrom(localid, earliest, corrmpi)) 
	    {
		error=1;
		errmsg="Error receving value for variable.";
	    }
	} 
	else if (localid==earliest)	// we will be the one to ship the value to others
	{
	    for (unsigned int target=0;target<swcount;++target)	// look through all worlds and see who we need to send to
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
      // It is possible that one or more subworlds were not participating in a failed transfer so we need to check
    char nerr;
    int mpierr=MPI_Allreduce(&error, &nerr, 1, MPI_UNSIGNED_CHAR, MPI_SUM, everyone->comm);
    if (mpierr!=MPI_SUCCESS)
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


// Find out which variables the local queued jobs are interested in
// share that info around
bool SubWorld::synchVariableInfo(std::string& err)
{

    if (!manualimports)
    {
	for (str2char::iterator it=varstate.begin();it!=varstate.end();++it)
	{
	    if (it->second==rs::NONE)
	    {
	        it->second=rs::INTERESTED;
	    }
	    else if (it->second==rs::OLD)
	    {
	        it->second=rs::OLDINTERESTED; 
	    }
	}      
    }
    else		// manual control over imports
    {
	for (size_t i=0;i<jobvec.size();++i)
	{
	    bp::list wanted=bp::extract<bp::list>(jobvec[i].attr("wantedvalues"))();		  
	    for (size_t j=0;j<len(wanted);++j)
	    {
		bp::extract<std::string> exs(wanted[j]);
		if (!exs.check()) 
		{
		    err="names in wantedvalues must be strings";
		    return false;
		}
		std::string n=exs();
		  // now we need to check to see if this value is known
		str2char::iterator it=varstate.find(n);
		if (it==varstate.end())
		{
		    err="Attempt to import variable \""+n+"\". SplitWorld was not told about this variable.";
		    return false;
		}
		// So at least one job wants this variable
		switch (it->second)
		{
		  case rs::NONE: it->second=rs::INTERESTED; break;
		  case rs::INTERESTED: break;
		  case rs::OLD: it->second=rs::OLDINTERESTED; break;
		  case rs::NEW: break;
		  default:
		    err="Unknown variable state";
		    return false;
		}
	    }
	}
    }
      
        // Make a vector to hold the info from the map (so we can send it around)
    std::vector<char> lb(getNumVars(), rs::NONE);
    size_t i=0;
    for (str2char::iterator it=varstate.begin();it!=varstate.end();++it,++i)
    {
        lb[i]=it->second;
    }


#ifdef ESYS_MPI
        // Vector to hold the result
    globalvarinfo.resize(getNumVars()*swcount, rs::NONE);
    if (amLeader())	// we only need on representative from each world to send
    {
	// The leaders of each world, send their variable information to the proc "0" in
	// the global world (which will be the leader of subworld "0").
	//    There is an issue here if this operation fails
	if (MPI_Gather(&lb[0], getNumVars(), MPI_CHAR, &globalvarinfo[0], getNumVars(), 
		   MPI_CHAR, 0, getCorrMPI()->comm)!=MPI_SUCCESS) 
	{
	    for (size_t i=0;i<globalvarinfo.size();++i)
	    {
		globalvarinfo[i]=rs::ERROR;
	    }
	}      
    }
    // now share the combined info with all processes
    if ((MPI_Bcast(&globalvarinfo[0], globalvarinfo.size(), MPI_CHAR, 0, everyone->comm)!=MPI_SUCCESS)
	  || (globalvarinfo[0]==rs::ERROR))
    {
	err="Error while gathering variable use information.";
	return false;	
    }
#endif    
    return true;
}

// merge / ship values as required
bool SubWorld::synchVariableValues(std::string& err)
{
#ifdef ESYS_MPI  
    // so now, everybody knows what everyone else wants
    size_t maxv=std::numeric_limits<size_t>::max();   // largest possible value in size_t
    size_t nv=getNumVars();
    str2reduce::iterator rmi=reducemap.begin();
    char error=0;
    for (size_t i=0;i<nv;++i, ++rmi)
    {
	size_t earliest=maxv;
	char lookingfor=rs::OLD;
	    // check to see if anyone has an updated copy of this variable
	for (int s=0;s<swcount;++s)
	{
	    if (globalvarinfo[s*nv+i]==rs::NEW)	// This world has a copy of the value
	    {
		lookingfor=rs::NEW;    
	    }
	}
	
	for (int s=0;s<swcount;++s)
	{
	    if (globalvarinfo[s*nv+i]==lookingfor)	// This world has a copy of the value
	    {
		earliest=s;
		break;
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
	if (globalvarinfo[localid*nv+i]==rs::INTERESTED)
	{
	    if (!rmi->second->recvFrom(localid, earliest, corrmpi)) 
	    {
		err="Error receving value for variable.";
		return false;
	    }
	} 
	else if (localid==earliest)	// we will be the one to ship the value to others
	{
	    for (unsigned int target=0;target<swcount;++target)	// look through all worlds and see who we need to send to
	    {
		if (globalvarinfo[target*nv+i]==rs::INTERESTED)	// nobody who has the value will be marked as INTERESTED
		{
		    if (!rmi->second->sendTo(localid, target, corrmpi))
		    {
			err="Error sending value for variable.";
			return false;
		    }
		}
	    }
	}
    }
      // It is possible that one or more subworlds were not participating in a failed transfer so we need to check
    char nerr;
    int mpierr=MPI_Allreduce(&error, &nerr, 1, MPI_UNSIGNED_CHAR, MPI_SUM, everyone->comm);
    if (mpierr!=MPI_SUCCESS)
    {
	if (err.empty())
	{
	    err="Error checking for other errors.";
	}
	return false;
    }
    if (nerr)
    {
	err="Another process experienced an error delivering variables.";
	return false;
    }
#endif

        // This is a simpler version of the above, since we already know it will work
    for (size_t i=0;i<jobvec.size();++i)
    {
	bp::list wanted=bp::extract<bp::list>(jobvec[i].attr("wantedvalues"))();		  
	for (size_t j=0;j<len(wanted);++j)
	{
	    bp::extract<std::string> exs(wanted[j]);
	    std::string n=exs();
	      // now we need to check to see if this value is known
	    str2reduce::iterator it=reducemap.find(n);
	    if (it!=reducemap.end())	// for sanity's sake
	    {
                jobvec[i].attr("setImportValue")(it->first, reducemap[it->first]->getPyObj());
	    }
	}
    }
    return true;
  
}

// Resize the vector so it holds the current number of variables
// walk along all variables in the maps and determine, what our interest is in each one
void SubWorld::getVariableStatus(std::vector<char>& vb)
{
    vb.resize(reducemap.size());
    str2reduce::const_iterator rm=reducemap.begin();
    str2bool::const_iterator im=importmap.begin();
    for (unsigned int i=0;rm!=reducemap.end();++rm,++im,++i)
    {
	if (rm->second->hasValue())
	{
	    vb[i]=rs::NEW;	// we have a value to contribute
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
void SubWorld::addVariable(std::string& name, Reducer_ptr& rp)
{
    if (reducemap.find(name)!=reducemap.end())
    {
	std::ostringstream oss;
	oss << "There is already a variable called " << name;
	throw SplitWorldException(oss.str());    
    }
    reducemap[name]=rp;
    importmap[name]=false;
    varstate[name]=reducerstatus::NONE;
    if (!manualimports)
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
    importmap.erase(s);
    varstate.erase(s);
#ifdef ESYS_MPI    
    globalvarinfo.resize(0);	// a size of 0 indicates the info is invalid
#endif    
}

void SubWorld::resetInterest()
{
    for (str2char::iterator it=varstate.begin();it!=varstate.end();++it)
    {
	if (it->second==rs::INTERESTED)
	{
	    it->second=rs::NONE;
	}
	else if (it->second==rs::OLDINTERESTED)
	{
	    it->second=rs::OLD;
	}
    }
}

