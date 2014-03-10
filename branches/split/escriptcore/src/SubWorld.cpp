
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

#include "SubWorld.h"

using namespace escript;


SubWorld::SubWorld(MPI_Comm comm)
    :communicator(comm), domain((AbstractDomain*)0)
{


}

SubWorld::~SubWorld()
{
}

MPI_Comm SubWorld::getComm()
{
    return communicator;
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

// if 3, a Job threw an exception 
// if 2, a Job did not return a bool
// if 1, at least one Job returned False
// if 0, all jobs in this world returned True
char SubWorld::runJobs()
{
    int ret=0;
    try
    {
	for (size_t i=0;i<jobvec.size();++i)
	{
	    boost::python::object result=jobvec[i].attr("work")();
	    boost::python::extract<bool> ex(result);
	    if (!ex.check())
	    {
		return 2;	
	    }
	    if (!ex())
	    {
		ret=1;
	    }
	}
    } catch (boost::python::error_already_set e)
    {
	return 3;
    }  
    return ret;
}