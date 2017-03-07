
/*****************************************************************************
*
* Copyright (c) 2014-2017 by The University of Queensland
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

#include "SubWorld.h"
#include "MPIDataReducer.h"
#include "MPIScalarReducer.h"
#include "NonReducedVariable.h"
#include "SplitWorldException.h"
#include "pyerr.h"

#include <boost/python/import.hpp>
#include <boost/python/dict.hpp>

#include <iostream>

using namespace escript;
namespace bp = boost::python;
namespace rs = escript::reducerstatus;

using namespace std;

SubWorld::SubWorld(JMPI& global, JMPI& comm, JMPI& corr,
                   unsigned int subworldcount, unsigned int local_id,
                   bool manualimport)
    : everyone(global), swmpi(comm), corrmpi(corr), domain((AbstractDomain*)0),
    swcount(subworldcount), localid(local_id), manualimports(manualimport)
#ifdef ESYS_MPI
    ,globalinfoinvalid(true)
#endif
{
        swcount=subworldcount;  // redundant to keep clang happy
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

void SubWorld::addJob(bp::object j)
{
    jobvec.push_back(j);
}

void SubWorld::clearJobs()
{
    jobvec.clear();
}

void SubWorld::setMyVarState(const std::string& vname, char state)
{
    setVarState(vname, state, localid);
}

void SubWorld::setAllVarsState(const std::string& vname, char state)
{
#ifdef ESYS_MPI
      // we need to know where the variable is in the sequence
    str2char::iterator it=varstate.find(vname);
    size_t c=0;
    for (;it!=varstate.end();++it,++c)
    {
        if (it->first==vname)
        {
            break;
        }
    }
    if (it==varstate.end())
    {
        return;
    }
    it->second=state;
    c--;                // we now have the sequence position of the variable
    for (char z=rs::NONE; z<=rs::NEW;++z)
    {
        globalvarcounts[vname][z]=0;
    }
    globalvarcounts[vname][state]=swcount;
    if (!globalinfoinvalid)     // it will be updated in the next synch
    {
        for (size_t p=c;p<globalvarinfo.size();p+=getNumVars())
        {
            globalvarinfo[p]=state;
        }
    }
#else
    varstate[vname]=state;
#endif
}


void SubWorld::setVarState(const std::string& vname, char state, int swid)
{
#ifdef ESYS_MPI
      // we need to know where the variable is in thbe sequence
    str2char::iterator it;
    size_t c=0;
    for (it=varstate.begin();it!=varstate.end();++it,++c)
    {
        if (it->first==vname)
        {
            break;
        }
    }
    if (it==varstate.end())
    {
        return;
    }
        // we now have the sequence position of the variable
    if (!globalinfoinvalid)     // it will be updated in the next synch
    {
        unsigned char ostate=globalvarinfo[c+getNumVars()*swid];
        globalvarinfo[c+getNumVars()*swid]=state;
        globalvarcounts[vname][ostate]--;
        globalvarcounts[vname][state]++;
    }
    if (swid==localid)  // we are updating our own state so we need to change "varstate"
    {
        it->second=state;
    }
#else
    varstate[vname]=state;
#endif
}


// this will give the imported values to interested jobs
bool SubWorld::deliverImports(std::string& errmsg)
{
    for (size_t i=0;i<jobvec.size();++i)
    {
        if (manualimports)
        {
            bp::list wanted=bp::extract<bp::list>(jobvec[i].attr("wantedvalues"))();
            for (size_t j=0;j<len(wanted);++j)
            {
                bp::extract<std::string> exs(wanted[j]);        // must have been checked by now
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
                catch (bp::error_already_set e)
                {
                    getStringFromPyException(e, errmsg);
                    return false;
                }
            }
        }
        else
        {
              // For automatic imports, we want to import "Everything" into every job.
              // However, we don't want to import things with no value yet
            for (str2reduce::iterator it=reducemap.begin();it!=reducemap.end();++it)
            {
                if (it->second->hasValue())
                {
                    try
                    {
                        jobvec[i].attr("setImportValue")(it->first, it->second->getPyObj());
                    }
                    catch (bp::error_already_set e)
                    {
                        getStringFromPyException(e, errmsg);
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

// Gather exported values from jobs and merge them in the reducer
bool SubWorld::localTransport(std::string& errmsg)
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
                return false;   // the error string will be set by the reduceLocalValue
            }
            setMyVarState(name, rs::NEW);
        }
    }
    return true;
}

void SubWorld::debug()
{
    using namespace std;
    using namespace escript::reducerstatus;
    std::cout << "Variables:";
#ifdef ESYS_MPI
        if (!globalinfoinvalid)
        {
            cout << "{ NONE INTR OLD OINT NEW }";
        }
        else
        {
            cout << "(no valid global info)";
        }
#endif
    std::cout << std::endl;
    int i=0;
    for (str2char::iterator it=varstate.begin();it!=varstate.end();++it,++i)
    {
        std::cout << it->first << ": ";
        std::cout << reducemap[it->first]->description() << " ";
        switch (it->second)
        {
          case NONE: cout << "NONE "; break;
          case INTERESTED: cout << "INTR "; break;
          case OLDINTERESTED: cout << "OINT "; break;
          case OLD: cout << "OLD  "; break;
          case NEW: cout << "NEW  "; break;
        }
#ifdef ESYS_MPI
        if (!globalinfoinvalid)
        {
            cout << "{ ";
            for (unsigned char z=rs::NONE;z<=rs::NEW;++z)
            {
                cout << globalvarcounts[it->first][z] << ' ';
            }
            cout << " } ";
        }
        else
        {
            cout << "(no valid global info)";
        }
#endif  
        cout << endl;
    }

#ifdef ESYS_MPI
    if (!globalinfoinvalid)
    {
        cout << "[";
        for (size_t i=0;i<globalvarinfo.size();++i)
        {
            if (i%getNumVars()==0)
            {
                cout << " ";
            }
            cout << (short)globalvarinfo[i];
        }
        cout << " ] ";
        
    }


#endif
    std::cout << "Debug end\n";
    std::cout.flush();
}


// not to be called while running jobs
// The tricky bit, is that this could be be called between job runs
// this means that the values of variables may not have been synched yet
DataTypes::real_t SubWorld::getScalarVariable(const std::string& name)
{
    str2reduce::iterator it=reducemap.find(name);
    if (it==reducemap.end())
    {
        throw SplitWorldException("No variable of that name.");
    }
        // need to indicate we are interested in the variable
    if (varstate[name]==rs::NONE)
    {
        setMyVarState(name, rs::INTERESTED);
    }
    else if (varstate[name]==rs::OLD)
    {
        setMyVarState(name, rs::OLDINTERESTED);
    }
        // anything else, indicates interest anyway
#ifdef ESYS_MPI
    std::string errmsg;
    if (!synchVariableInfo(errmsg))
    {
        throw SplitWorldException(std::string("(Getting scalar --- Variable information) ")+errmsg);
    }
    if (!synchVariableValues(errmsg))
    {
        throw SplitWorldException(std::string("(Getting scalar --- Variable value) ")+errmsg);
    }
#endif
    if (dynamic_cast<MPIScalarReducer*>(it->second.get()))
    {
        return dynamic_cast<MPIScalarReducer*>(it->second.get())->getDouble();
    }
    if (dynamic_cast<NonReducedVariable*>(it->second.get()))
    {
        bp::extract<DataTypes::real_t> ex(it->second->getPyObj());
        if (!ex.check())
        {
            throw SplitWorldException("Variable is not scalar.");
        }
        return ex();
    }
    throw SplitWorldException("Variable is not scalar.");
}


// not to be called while running jobs
// The tricky bit, is that this could be be called between job runs
// this means that the values of variables may not have been synched yet
bp::object SubWorld::getLocalObjectVariable(const std::string& name)
{
    str2reduce::iterator it=reducemap.find(name);
    if (it==reducemap.end())
    {
        throw SplitWorldException("No variable of that name.");
    }
        // need to indicate we are interested in the variable
    if (varstate[name]==rs::NONE)
    {
        setMyVarState(name, rs::INTERESTED);
    }
    else if (varstate[name]==rs::OLD)
    {
        setMyVarState(name, rs::OLDINTERESTED);
    }
        // anything else, indicates interest anyway
#ifdef ESYS_MPI
    std::string errmsg;
    if (!synchVariableInfo(errmsg))
    {
        throw SplitWorldException(std::string("(Getting local object --- Variable information) ")+errmsg);
    }
    if (!synchVariableValues(errmsg))
    {
        throw SplitWorldException(std::string("(Getting local object --- Variable value) ")+errmsg);
    }
#endif
    if (dynamic_cast<NonReducedVariable*>(it->second.get()))
    {
        return dynamic_cast<NonReducedVariable*>(it->second.get())->getPyObj();
    }
    throw SplitWorldException("Variable is not a local object.");
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
bool SubWorld::makeComm(MPI_Comm& sourcecom, JMPI& ncom,std::vector<int>& members)
{
      MPI_Comm subcom;
      MPI_Group sourceg, g;
      if (MPI_Comm_group(sourcecom, &sourceg)!=MPI_SUCCESS) {return false;}
      if (MPI_Group_incl(sourceg, members.size(), &members[0], &g)!=MPI_SUCCESS) {return false;}
      // then create a communicator with that group
      if (MPI_Comm_create(sourcecom, g, &subcom)!=MPI_SUCCESS)
      {
          return false; 
      }
      ncom=makeInfo(subcom, true);
      return true;
}


// The mystate, could be computed from vnum, this is just to shortcut
// creates two groups, the first contains procs which need to reduce
// the second group contains a single process with the new value and
// all other interested parties
bool SubWorld::makeGroupReduceGroups(MPI_Comm& srccom, int vnum, char mystate, JMPI& red, JMPI& cop, bool& incopy)
{
    incopy=false;
    if ((mystate==rs::NEW)
            || (mystate==rs::INTERESTED)
            || (mystate==rs::OLDINTERESTED))
    {
        // first create a group with all the updates in it
        std::vector<int> redmembers;
        std::vector<int> copmembers;
        for (int i=0+vnum;i<globalvarinfo.size();i+=getNumVars())
        {
            bool havesrc=false;
            int world=i/getNumVars();
            // make a vector of the involved procs with New at the front
            switch (globalvarinfo[i])
            {
                case rs::NEW:
                    if (!havesrc)
                    {
                        copmembers.insert(copmembers.begin(), world);
                        havesrc=true;
                        if (world==localid)
                        {
                            incopy=true;                        
                        }
                    }
                    redmembers.push_back(world);
                    break;
                case rs::INTERESTED:
                case rs::OLDINTERESTED:
                          copmembers.push_back(world);
                          if (world==localid)
                          {
                              incopy=true;
                          }
                          break;
            }
        }
        if (!makeComm(srccom, red, redmembers))
        {
            return false;
        }
        if (!makeComm(srccom, cop, copmembers))
        {
            return false;
        }
        return true;

    }
    else  // for people not in involved in the value shipping
    {     // This would be a nice time to use MPI_Comm_create_group
          // but it does not exist in MPI2.1
        MPI_Comm temp;
        if (MPI_Comm_create(srccom, MPI_GROUP_EMPTY, &temp)!=MPI_SUCCESS)
        {
            return false;
        }
        red=makeInfo(temp, true);
        if (MPI_Comm_create(srccom, MPI_GROUP_EMPTY, &temp)!=MPI_SUCCESS)
        {
            return false;
        }
        cop=makeInfo(temp, true);
        return true;
    }

}


// a group with NEW nodes at the front and INT and OLDINT at the back
// NONE worlds get an empty communicator
bool SubWorld::makeGroupComm1(MPI_Comm& srccom, int vnum, char mystate, JMPI& com)
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
                case rs::NEW:   members.insert(members.begin(), i/getNumVars()); break;
                case rs::INTERESTED:
                case rs::OLDINTERESTED:
                          members.push_back(i/getNumVars());
                          break;
              }
          }
          return makeComm(srccom, com, members);
      }
      else      // for people not in involved in the value shipping
      {         // This would be a nice time to use MPI_Comm_create_group
                // but it does not exist in MPI2.1
          MPI_Comm temp;
          MPI_Comm_create(srccom, MPI_GROUP_EMPTY, &temp);
          com=makeInfo(temp, true);
          return true;
      }
}

// A group with a single OLD or OLDINT at the front and all the INT worlds
// following it
bool SubWorld::makeGroupComm2(MPI_Comm& srccom, int vnum, char mystate, JMPI& com, bool& ingroup)
{
      ingroup=false;
      if ((mystate==rs::OLD)
            || (mystate==rs::INTERESTED)
            || (mystate==rs::OLDINTERESTED))
      {
          // first create a group with [old, interested and oldinterested in it]
          std::vector<int> members;
          bool havesrc=false;
          for (int i=0+vnum;i<globalvarinfo.size();i+=getNumVars())
          {
              int world=i/getNumVars();
              // make a vector of the involved procs with OLD/OLDINTERESTED at the front
              switch (globalvarinfo[i])
              {
                case rs::NEW:   return false;  break;
                case rs::INTERESTED: members.push_back(world);
                          if (world==localid)
                          {
                              ingroup=true;
                          }
                          break;
                case rs::OLD:
                case rs::OLDINTERESTED:
                          if (!havesrc)
                          {
                              members.insert(members.begin(), world);
                              havesrc=true;
                              if (world==localid)
                              {
                                ingroup=true;
                              }
                          }
                          break;
              }
          }             
          return makeComm(srccom, com, members);
      }
      else      // for people not in involved in the value shipping
      {         // This would be a nice time to use MPI_Comm_create_group
                // but it does not exist in MPI2.1      
          MPI_Comm temp;
          MPI_Comm_create(srccom, MPI_GROUP_EMPTY, &temp);
          com=makeInfo(temp, true);
          return true;
      }
}
#endif


bool SubWorld::synchVariableValues(std::string& err)
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
        int newcount=0; // who has a new version
        int oldcount=0; // who has an old version
        int oldintcount=0;
        newcount=globalvarcounts[it->first][rs::NEW];
        oldcount=globalvarcounts[it->first][rs::OLD];
        oldintcount=globalvarcounts[it->first][rs::OLDINTERESTED];
        needcount=globalvarcounts[it->first][rs::INTERESTED]+oldintcount;
        if (newcount>0)
        {
            varswithupdates.push_back(it->first);
        }
        if (needcount+newcount+oldcount==0)
        {
            continue;           // noone cares about this variable
        }
        if (needcount>0 && (oldcount+oldintcount+newcount)==0)
        {
            err="Import attempted for a variable \""+(it->first)+"\" with no value.";
            return false;
        }
            // worlds have the variable but noone is interested in it
            // note that if there are multiple new values, we still want to merge them
        if ((needcount==0) && (newcount<=1))
        {
            continue;   
        }
        if (swcount==1)
        {               // nobody else to communicate with
            continue;
        }
            // to reach this point, there must be >=1 source and >=1 sink and multiple worlds
            // first deal updates as source(s)
        if (newcount==1)        // only one update so send from that
        {
            JMPI com;
            if (!makeGroupComm1(corrmpi->comm, vnum, varstate[it->first],com))
            {
                err="Error creating group for sharing values,";
                return false;
            }
            if (varstate[it->first]!=rs::NONE && varstate[it->first]!=rs::OLD)
            {
                it->second->groupSend(com->comm, (varstate[it->first]==rs::NEW));
                  // Now record the fact that we have the variable now
                if (varstate[it->first]==rs::INTERESTED)
                {
                    setMyVarState(it->first, rs::OLDINTERESTED);
                }
            }
            continue;
        }
        if (newcount==swcount)          // everybody is in on this
        {
            if (!it->second->reduceRemoteValues(corrmpi->comm))
            {
                it->second->reset();
                setAllVarsState(it->first, rs::NONE);
                //setMyVarState(it->first, rs::NONE);
                err=it->first+"Either MPI failed, or there were multiple simultaneous updates to a variable with the SET operation.";
                return false;
            }
                // Now record the fact that we have the variable now
            if (varstate[it->first]==rs::INTERESTED)
            {
                setMyVarState(it->first, rs::OLDINTERESTED);
            }
            continue;
        }
        if (newcount>1)
        {
            // make groups to reduce and then copy
            JMPI red;
            JMPI cop;
            bool incopy;
            if (!makeGroupReduceGroups(corrmpi->comm, vnum, varstate[it->first], red, cop, incopy))
            {
                err="Error creating groups for sharing values,";
                return false;
            }
            char reduceresult=0;
                // only new values get reduced
            if (varstate[it->first]==rs::NEW)
            {
                if (!it->second->reduceRemoteValues(red->comm))
                {
                    char s=1;
                    MPI_Allreduce(&s, &reduceresult, 1, MPI_CHAR, MPI_MAX, corrmpi->comm);
                    reduceresult=1;

                }
                else
                {
                    if (it->second->canClash())
                    {
                        char s=0;
                        MPI_Allreduce(&s, &reduceresult, 1, MPI_CHAR, MPI_MAX, corrmpi->comm);
                    }
                }
            }
            else
            {
                if (it->second->canClash())
                {
                    char s=0;
                    MPI_Allreduce(&s, &reduceresult, 1, MPI_CHAR, MPI_MAX, corrmpi->comm);
                }
            }
                // if there was a clash somewhere
            if (reduceresult!=0)
            {
                it->second->reset();
                setAllVarsState(it->first, rs::NONE);
                err="Either MPI failed, or there were multiple simultaneous updates to a variable with the SET operation.";
                return false;
            }

                // if we are involved in copying the new value around
            if (incopy)
            {
                it->second->groupSend(cop->comm, (varstate[it->first]==rs::NEW));
                if (varstate[it->first]==rs::INTERESTED)
                {
                    setMyVarState(it->first, rs::OLDINTERESTED);
                }
            }
            if (varstate[it->first]==rs::NEW)
            {
                setMyVarState(it->first, rs::OLDINTERESTED);
            }
            continue;
        }
            // at this point, we need to ship info around but there are no updates
            // that is, we are shipping an old copy
            // picking a source arbitarily (the first one in the array)

            // but first, eliminate the special case where the only interested ones
            // already have a copy
        if (oldintcount==needcount)
        {
            continue;
        }
        JMPI com;
        bool ingroup=false;
        if (!makeGroupComm2(corrmpi->comm, vnum, varstate[it->first],com, ingroup))
        {
            err="Error creating group for sharing values";
            return false;
        }
        // form group to send to [latestsource and interested]
        
        if (ingroup)            // since only one holder needs to send
        {
            bool imsending=(varstate[it->first]==rs::NEW);
            it->second->groupSend(com->comm, imsending);
        }
    }
        // now we need to age any out of date copies of vars
    for (size_t i=0;i<varswithupdates.size();++i)
    {
        std::string vname=varswithupdates[i];
                if (varstate[vname]==rs::NEW)
        {
            setMyVarState(vname, rs::OLD);
        }
        else if (varstate[vname]==rs::OLD)
        {
            setMyVarState(vname, rs::NONE);
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

// Find out which variables the local queued jobs are interested in
// share that info around
bool SubWorld::synchVariableInfo(std::string& err)
{
    if (getNumVars()==0)
    {
        return true;
    }
    if (manualimports)          // manual control over imports
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
    if (amLeader())     // we only need on representative from each world to send
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
      // now we convert that info into a form which is easier to read
    int p=0;
    for (str2reduce::iterator it=reducemap.begin();it!=reducemap.end();++it,++p)
    {
        globalvarcounts[it->first][rs::NONE]=0;
        globalvarcounts[it->first][rs::INTERESTED]=0;
        globalvarcounts[it->first][rs::OLD]=0;
        globalvarcounts[it->first][rs::OLDINTERESTED]=0;
        globalvarcounts[it->first][rs::NEW]=0;
        for (int j=p;j<globalvarinfo.size();j+=getNumVars())
        {
            if (globalvarinfo[j]<=rs::NEW)
            {
                globalvarcounts[it->first][globalvarinfo[j]]++;
            }
        }
    }

#endif
    if (!manualimports) 
    {
            // import all known variables _BUT_ don't import something if noone has a value
            // for it
        int vnum=0;
        for (str2char::iterator it=varstate.begin();it!=varstate.end();++it, ++vnum)
        {
#ifdef ESYS_MPI
              // if at least one world has a value for a variable
            if (globalvarcounts[it->first][rs::OLDINTERESTED]
                 + globalvarcounts[it->first][rs::OLD]
                 + globalvarcounts[it->first][rs::NEW] > 0 )
            {
#endif

                if (it->second==rs::NONE)
                {
                    it->second=rs::INTERESTED;
                }
                else if (it->second==rs::OLD)
                {
                    it->second=rs::OLDINTERESTED;
                }
#ifdef ESYS_MPI
                  // now we need to update the globalvarinfo to record all the extra interest
                for (int j=vnum;j<globalvarinfo.size();j+=getNumVars())
                {
                    if (globalvarinfo[j]==rs::NONE)
                    {
                        globalvarinfo[j]=rs::INTERESTED;
                        globalvarcounts[it->first][rs::NONE]--;
                        globalvarcounts[it->first][rs::INTERESTED]++;                   
                    }
                    else if (globalvarinfo[j]==rs::OLD)
                    {
                        globalvarinfo[j]=rs::OLDINTERESTED;
                        globalvarcounts[it->first][rs::OLD]--;
                        globalvarcounts[it->first][rs::OLDINTERESTED]++;                        
                    }
                }
            }
#endif  
        }
    }
#ifdef ESYS_MPI
    globalinfoinvalid=false;
#endif

    return true;
}

// if 4, a Job performed an invalid export
// if 3, a Job threw an exception
// if 2, a Job did not return a bool
// if 1, at least one Job returned False
// if 0, all jobs in this world returned True
char SubWorld::runJobs(std::string& errormsg)
{
    errormsg.clear();
    int ret=0;
    try
    {
        for (size_t i=0;i<jobvec.size();++i)
        {
            bp::object result=jobvec[i].attr("work")();
            bp::extract<bool> ex(result);
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
    catch (bp::error_already_set e)
    {
        getStringFromPyException(e, errormsg);
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
        throw SplitWorldException(oss.str());
    }
    if (domain.get()==0)
    {
        throw SplitWorldException("No domain has been set yet.");
    }
    rp->setDomain(domain);
    reducemap[name]=rp;
    varstate[name]=reducerstatus::NONE;
    if (!manualimports)
    {
        for (size_t i=0;i<jobvec.size();++i)
        {
            jobvec[i].attr("declareImport")(name);
        }
    }
#ifdef ESYS_MPI
    globalinfoinvalid=true;     // since we need to regenerate globalvarinfo
#endif
}

void SubWorld::removeVariable(std::string& s)
{
    reducemap.erase(s);
    varstate.erase(s);
#ifdef ESYS_MPI
    globalinfoinvalid=true;
    globalvarinfo.resize(0);
    globalvarcounts.erase(s);
#endif
}

void SubWorld::clearVariable(std::string& name)
{
    str2reduce::iterator it=reducemap.find(name);
    if (it==reducemap.end())
    {
        return;
    }
    it->second->reset();
      // if we got here, we must have a valid name so we can change state directly
    setAllVarsState(name, rs::NONE);
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

void SubWorld::newRunJobs()
{
    for (str2reduce::iterator it=reducemap.begin();it!=reducemap.end();++it)
    {
        it->second->newRunJobs();
    }
}

std::list<std::pair<std::string, bool> > SubWorld::getVarList()
{
    std::list<std::pair<std::string,bool> > res;
    for (std::map<std::string, Reducer_ptr>::iterator it=reducemap.begin();it!=reducemap.end();++it)
    {
        res.push_back(std::pair<std::string, bool>(it->first, it->second->hasValue()));
    }
    return res;
}

std::list<std::pair<std::string, std::string> > SubWorld::getVarInfo()
{
    std::list<std::pair<std::string,std::string> > res;
    for (std::map<std::string, Reducer_ptr>::iterator it=reducemap.begin();it!=reducemap.end();++it)
    {
        std::string desc=it->second->description();
        if (!it->second->hasValue())
        {
            desc+=" (No value)";
        }
        res.push_back(std::pair<std::string, std::string>(it->first, desc));
    }
    return res;
}


void SubWorld::copyVariable(const std::string& src, const std::string& dest)
{
        if (reducemap.find(src)==reducemap.end())
        {
            throw SplitWorldException("Source variable name is not known");
        }
        if (reducemap.find(dest)==reducemap.end())
        {
            throw SplitWorldException("Destination variable name is not known");
        }
        Reducer_ptr sptr=reducemap[src];
        Reducer_ptr dptr=reducemap[dest];
        dptr->copyValueFrom(sptr);
}

