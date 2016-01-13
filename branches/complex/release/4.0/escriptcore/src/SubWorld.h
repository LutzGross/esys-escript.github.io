
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
    SubWorld(esysUtils::JMPI& globalcom, esysUtils::JMPI& comm, esysUtils::JMPI& corr, unsigned int subworldcount, unsigned int local_id);
    ~SubWorld();
    void setDomain(Domain_ptr d);
    Domain_ptr getDomain();
    esysUtils::JMPI& getMPI();
    esysUtils::JMPI& getCorrMPI();
    void addJob(boost::python::object j);
    char runJobs(std::string& errmsg);
    void clearJobs();
    void clearImportExports();
    void addVariable(std::string&, Reducer_ptr& red, bool manualimport);
    void removeVariable(std::string& name);  
    size_t getNumVars();
    
    bool localTransport(std::vector<char>& vb, std::string& errmsg);
    bool checkRemoteCompatibility(std::string& errmsg);
    bool reduceRemoteValues(std::string& errmsg);
    bool deliverImports(std::vector<char>& vb, std::string& errmsg);
    
    
    bool findImports(bool manualimports, std::string& errmsg);
    bool deliverImports(std::string& errmsg);
    bool deliverGlobalImports(std::vector<char>& vb, std::string& errmsg);
    void getVariableStatus(std::vector<char>& vb);
    bool reduceRemoteValues();    
    bool amLeader();	// true if this proc is the leader for its world
    
    
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
    str2reduce reducemap;		// map: name ->reducer for that variable
    str2bool importmap;

};

typedef boost::shared_ptr<SubWorld> SubWorld_ptr;



}
#endif
