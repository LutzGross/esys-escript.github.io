
/*****************************************************************************
*
* Copyright (c) 2014-2015 by The University of Queensland
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

#ifndef escript_SubWorld_H
#define escript_SubWorld_H

#include "esysUtils/Esys_MPI.h"
#include "AbstractDomain.h"
#include "AbstractReducer.h"

namespace escript
{
  
/** This class represents a collection of MPI processes which will execute a number of Jobs (in series).
 * There could be a number of SubWorlds (executing jobs in parallel) in the overall system.
 * All jobs running in a SubWorld will use a common domain object.
 * After each job runs, any values it exports will be merged into local reducer objects.
 * Global (ie with the participation of other SubWorlds) reductions and interworld transfers
 * are handled after the current batch of jobs have completed.
 * That is, variable values should not be considered up to date until the whole batch is complete.
 * Further, after a batch has completed, multiple subworlds may have copies of the variable,
 * if the variable is modified in a later batch, this may result in unwanted double counting.
 * eg: v (reduce:+)
 * Batch 1:
 *   world 1:   v+=1,2,3  --- local v=6
 *   world 2:   v+=1,2    --- local v=3
 * What is the value of v in this split world?   v=9
 * 
 * Batch 2:
 *   world 1:  v+=1   --- local v=1+9
 *   world 2:  v+=1   --- local v=1+9
 * What is the value of v? 20, not 11
*/
class SubWorld : public boost::enable_shared_from_this<SubWorld>
{
public:
    SubWorld(esysUtils::JMPI& globalcom, esysUtils::JMPI& comm, esysUtils::JMPI& corr, unsigned int subworldcount, unsigned int local_id, bool manualimport);
    ~SubWorld();
    void setDomain(Domain_ptr d);
    Domain_ptr getDomain();
    esysUtils::JMPI& getMPI();
    esysUtils::JMPI& getCorrMPI();
    void addJob(boost::python::object j);		// add a Job to the current batch
    char runJobs(std::string& errmsg);		// run all jobs in the current batch
    void clearJobs();				// remove all jobs in the current batch

    void addVariable(std::string&, Reducer_ptr& red);
    void removeVariable(std::string& name);  
    void clearVariable(std::string& name);
    std::list<std::pair<std::string, bool> > getVarList();
    size_t getNumVars();
    
    bool localTransport(std::string& errmsg);	// gather exported values from jobs
    bool checkRemoteCompatibility(std::string& errmsg);	// check to ensure values
						// in all worlds are compatible
    
    bool deliverImports(std::string& errmsg);	// load imports into Job objects
    bool amLeader();	// true if this proc is the leader for its world
    
    double getScalarVariable(const std::string& name);
    
    void debug();	// print out current state information
    
    
    
    bool synchVariableInfo(std::string& err);
    bool synchVariableValues(std::string& err);    
    void resetInterest();    

    void copyVariable(const std::string& src, const std::string& dest);
    
    void newRunJobs();
    
private:
    esysUtils::JMPI everyone;	// communicator linking all procs in all subworlds
    esysUtils::JMPI swmpi;	// communicator linking all procs in this subworld
    esysUtils::JMPI corrmpi;	// communicator linking corresponding procs in all subworlds
				// eg: If this proc is the first in its domain, then corrmpi
				//     links to the other "first in its domain" processes.
				//      (So one in each SubWorld).
    escript::Domain_ptr domain;
    std::vector<boost::python::object> jobvec;	// jobs in the current batch
    
    
    unsigned int swcount;		// number of subwords
    unsigned int localid;    	// position of this subworld in that sequence
    
typedef std::map<std::string, Reducer_ptr> str2reduce;  
typedef std::map<std::string, unsigned char> str2char;
    str2reduce reducemap;		// map: name ->reducer for that variable
    str2char varstate;		// using the state values from AbstractReducer.h

    bool manualimports;
    
#ifdef ESYS_MPI    
    std::vector<unsigned char> globalvarinfo;	// info about which worlds want which vars
				  // [vars on process0][vars on process 1][vars on ...]
typedef std::map<unsigned char, int> countmap;
typedef std::map<std::string, countmap> str2countmap;
    str2countmap globalvarcounts;
    bool globalinfoinvalid;
    
    
    bool makeComm(MPI_Comm& sourcecom, MPI_Comm& subcom,std::vector<int>& members);


    // a group with NEW nodes at the front and INT and OLDINT at the back
    // NONE worlds get an empty communicator
    bool makeGroupComm1(MPI_Comm& srccom, int vnum, char mystate, MPI_Comm& com);

    // A group with a single OLD or OLDINT at the front and all the INT worlds 
    // following it
    bool makeGroupComm2(MPI_Comm& srccom, int vnum, char mystate, MPI_Comm& com, bool& ingroup);    
    
#endif
    
      // change the various views of a variable's state
    void setMyVarState(const std::string& vname, char state);
    void setVarState(const std::string& vname, char state, int swid);
    void setAllVarsState(std::string& name, char state);
};

typedef boost::shared_ptr<SubWorld> SubWorld_ptr;



}
#endif
