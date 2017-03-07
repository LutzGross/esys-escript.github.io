
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

#ifndef escript_SplitWorld_H
#define escript_SplitWorld_H

#include "AbstractReducer.h"
#include "SubWorld.h"

#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>

namespace escript
{

/** 
 * Provides an interface to a collection of subworlds.
 * Variables are declared and jobs are submitted using this interface.
 * Internally, the work is done by a local subworld instance (and associated communicators which 
 * this process belongs to). The local subworld will communicate with subworlds 
 * in other processes as needed.
 * The main reason for this class, is to insulate users from the MPI type thinking needed for
 * subworlds and instead provide an interface which allows them to think about subworlds as a group.
*/
class SplitWorld
{
public:
    SplitWorld(unsigned int numgroups, MPI_Comm global=MPI_COMM_WORLD);
    ~SplitWorld();
    boost::python::object buildDomains(boost::python::tuple t, boost::python::dict kwargs);	
    
    void runJobs();
    
    void addJob(boost::python::object creator, boost::python::tuple tup, boost::python::dict kw);
    void addJobPerWorld(boost::python::object creator, boost::python::tuple tup, boost::python::dict kw);
    
    void addVariable(std::string name, boost::python::object creator, boost::python::tuple ntup, boost::python::dict kwargs);
    void removeVariable(std::string name); 
    void clearVariable(std::string name); 
    std::list<std::pair<std::string, bool> > getVarList();
    boost::python::object getVarPyList();
    boost::python::object getVarPyInfo();    

    void clearAllJobs();

    DataTypes::real_t getScalarVariable(const std::string& name);
    boost::python::object getLocalObjectVariable(const std::string& name);
    
    int getSubWorldCount();
    int getSubWorldID();

    void copyVariable(const std::string& src, const std::string& dest);     
    
    
private:    
    JMPI globalcom;	// communicator linking all procs used in this splitworld
    JMPI leadercom;	// communicator linking the first proc in each subworld
    escript::SubWorld_ptr localworld;	// subworld which this process belongs to
    unsigned int swcount;		// number of subwords
    unsigned int localid;		// position of localworld in overall world sequence
    
    // details of jobs to be created
    std::vector<boost::python::object> create;
    std::vector<boost::python::tuple> tupargs;
    std::vector<boost::python::dict> kwargs;
    
    unsigned int jobcounter;		// note that the id of the first job is 1 not 0.
    bool manualimport;		// if false, all reduced vars will be shipped to all subworlds    
    void clearPendingJobs();
    void distributeJobs();

};


/**
  used to invoke the SplitWorld version from python (in lieu of a method based equivalent to raw_function) 
*/
boost::python::object raw_buildDomains(boost::python::tuple t, boost::python::dict kwargs);

/**
 used to invoke the SplitWorld version from python (in lieu of a method based equivalent to raw_function)
*/
boost::python::object raw_addJob(boost::python::tuple t, boost::python::dict kwargs);

/**
 used to invoke the SplitWorld version from python (in lieu of a method based equivalent to raw_function)
*/
boost::python::object raw_addJobPerWorld(boost::python::tuple t, boost::python::dict kwargs);

/**
 used to add a reducer for shared values.
*/
boost::python::object raw_addVariable(boost::python::tuple t, boost::python::dict kwargs);
}
#endif

