
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

#ifndef escript_WorldSplitter_H
#define escript_WorldSplitter_H
#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>
#include "esysUtils/Esys_MPI.h"
#include "SubWorld.h"
namespace escript
{

/** class to hold a collection of MPI ranks and a communicator linking them
*/
class WorldSplitter
{
public:
    WorldSplitter(unsigned int numgroups, MPI_Comm global=MPI_COMM_WORLD);
    ~WorldSplitter();
    boost::python::object buildDomains(boost::python::tuple t, boost::python::dict kwargs);	
	// first param will be the factory method / constructor
	// return value will always be None (it would be void but apparently boost doesn't like that
	// this method will be def("set....", raw_function(WorldSplitter::setDomainParams, 1)
    
    void runJobs(boost::python::list l);
private:    
    MPI_Comm globalcom;	// don't free this because we don't own it
    MPI_Comm subcom;
    escript::SubWorld_ptr localworld;	// subworld which this process belongs to
    unsigned int groupcount;
    unsigned int localid;		// position of localworld in overall world sequence
};


/**
  used to invoke the WorldSplitter version from python (in lieu of a method based equivalent to raw_function) 
*/
boost::python::object raw_buildDomains(boost::python::tuple t, boost::python::dict kwargs);


}
#endif
