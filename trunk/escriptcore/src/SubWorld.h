
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

#ifndef escript_SubWorld_H
#define escript_SubWorld_H

#include "esysUtils/Esys_MPI.h"
#include "AbstractDomain.h"
#include "Reducer.h"

namespace escript
{
  
/** class to hold a collection of MPI processes and a communicator linking them
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
    void addJob(boost::python::object j);
    char runJobs(std::string& errmsg);
    void clearJobs();
    void clearImportExports();
    void addVariable(std::string&, Reducer_ptr& red);
    void removeVariable(std::string& name);  
    size_t getNumVars();
    
    bool localTransport(/*std::vector<char>& vb, */ std::string& errmsg);
    bool checkRemoteCompatibility(std::string& errmsg);
    bool reduceRemoteValues(std::string& errmsg);
    bool deliverImports(std::vector<char>& vb, std::string& errmsg);
    
    
//    bool findImports(std::string& errmsg);
    bool deliverImports(std::string& errmsg);
    bool deliverGlobalImports(std::vector<char>& vb, std::string& errmsg);
    /*void getVariableStatus(std::vector<char>& vb); */
    bool reduceRemoteValues();    
    bool amLeader();	// true if this proc is the leader for its world
    
    double getScalarVariable(const std::string& name);
    
    void debug();
    
    
    
    bool synchVariableInfo(std::string& err);
    bool synchVariableValues(std::string& err);    
    void ageVariables();
    void resetInterest();    
    
private:
    esysUtils::JMPI everyone;	// communicator linking all procs in all subworlds
    esysUtils::JMPI swmpi;	// communicator linking all procs in this subworld
    esysUtils::JMPI corrmpi;	// communicator linking corresponding procs in all subworlds
    escript::Domain_ptr domain;
    std::vector<boost::python::object> jobvec;
    
    
    unsigned int swcount;		// number of subwords
    unsigned int localid;    	// my position within the sequence
    
typedef std::map<std::string, Reducer_ptr> str2reduce;  
typedef std::map<std::string, bool> str2bool;
typedef std::map<std::string, unsigned char> str2char;
    str2reduce reducemap;		// map: name ->reducer for that variable
    str2char varstate;		// using the state values from Reducer.h

    bool manualimports;
    
#ifdef ESYS_MPI    
    std::vector<unsigned char> globalvarinfo;	// info about which worlds want which vars
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
    bool makeGroupComm2(MPI_Comm& srccom, int vnum, char mystate, MPI_Comm& com);    
    
#endif
};

typedef boost::shared_ptr<SubWorld> SubWorld_ptr;



}
#endif
