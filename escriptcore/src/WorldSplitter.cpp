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
#include "WorldSplitter.h"
#include "AbstractDomain.h"
#include "DomainException.h"

using namespace boost::python;
using namespace escript;

WorldSplitter::WorldSplitter(unsigned int numgroups, MPI_Comm global)
    :globalcom(global), subcom(MPI_COMM_NULL), localworld((SubWorld*)0), groupcount(numgroups)
{
    int gsize;
    int grank;
    if ((MPI_Comm_size(global, &gsize)!=MPI_SUCCESS) || (MPI_Comm_size(global, &grank)!=MPI_SUCCESS))
    {
	throw DomainException("MPI appears to be inoperative.");
    }
    if (gsize%numgroups!=0)
    {
	throw DomainException("WorldSplitter error: requested number of groups is not a factor of global communicator size.");
    }
    int res=MPI_Comm_split(MPI_COMM_WORLD, grank/gsize, grank%gsize, &subcom);
    if (res!=MPI_SUCCESS)
    {
	throw DomainException("WorldSplitter error: Unable to form communicator.");
    }
std::cerr << subcom << std::endl;    
    localworld=SubWorld_ptr(new SubWorld(subcom));
}

// We may need to look into this more closely.
// What if the domain lives longer than the world splitter?
WorldSplitter::~WorldSplitter()
{
    if (subcom!=MPI_COMM_NULL)
    {
	MPI_Comm_free(&subcom);
    }
}


// The boost wrapper will ensure that there is at least one entry in the tuple
object WorldSplitter::buildDomains(tuple t, dict kwargs)
{
    int tsize=len(t);
    // get the callable that we will invoke in a sec
    object tocall=t[0];
    // make a new tuple without the first element
    tuple ntup=tuple(t.slice(1,tsize));
    // now add the subworld to the kwargs
    kwargs["escriptworld"]=localworld;

// std::cerr << "About to call function with:\n";
// //extract<std::string> ex(ntup.attr("__str__")());
//  std::cerr << extract<std::string>(ntup.attr("__str__")())() << std::endl;
// // for (int i=0;i<tsize-1;++i)
// // {
// //     std::cout << extract<const char*>(ntup[i])() << " ";
// // }
// std::cerr << std::endl;
    
    // pass the whole package to the python call
    object dobj=tocall(*ntup, **kwargs);
    extract<Domain_ptr> ex1(dobj);
    Domain_ptr dptr=ex1();
    
    // now do a sanity check to see if the domain has respected the communicator info we passed it.
    if (dptr->getMPIComm()!=localworld->getComm())
    {
	throw DomainException("The newly constructed domain is not using the correct communicator.");
    }
    localworld->setDomain(dptr);
    return object();	// return None
}

namespace escript
{

boost::python::object raw_buildDomains(boost::python::tuple t, boost::python::dict kwargs)
{
    int l=len(t);
    if (l<2)
    {
	throw DomainException("Insufficient parameters to buildDomains.");
    }
    extract<WorldSplitter&> exw(t[0]);
    if (!exw.check())
    {
	throw DomainException("First parameter to buildDomains must be a WorldSplitter.");
    }
    WorldSplitter& ws=exw();
    tuple ntup=tuple(t.slice(1,l));	// strip off the object param
    return ws.buildDomains(ntup, kwargs);
}

}
