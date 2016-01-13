
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

#ifndef escript_SplitWorld_H
#define escript_SplitWorld_H
#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>
#include "esysUtils/Esys_MPI.h"
#include "SubWorld.h"
#include "Reducer.h"
namespace escript
{

/** class to hold a collection of MPI ranks and a communicator linking them
*/
class SplitWorld
{
public:
    SplitWorld(unsigned int numgroups, MPI_Comm global=MPI_COMM_WORLD);
    ~SplitWorld();
    boost::python::object buildDomains(boost::python::tuple t, boost::python::dict kwargs);	
    
    void runJobs();
    
    void addJob(boost::python::object creator, boost::python::tuple tup, boost::python::dict kw);
    
    void addVariable(std::string name, boost::python::object creator, boost::python::tuple ntup, boost::python::dict kwargs);
    void removeVariable(std::string name);    
    void clearActiveJobs();
    void clearAllJobs();

    
    
    
private:    
    esysUtils::JMPI globalcom;	// communicator linking all procs used in this splitworld
    esysUtils::JMPI leadercom;	// communicator linking the first proc in each subworld
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
    bool getVariableInterest(std::vector<char>& vb);    
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
 used to add a reducer for shared values.
*/
boost::python::object raw_addVariable(boost::python::tuple t, boost::python::dict kwargs);
}
#endif
