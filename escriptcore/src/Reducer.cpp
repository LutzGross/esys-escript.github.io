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

#include <sstream>
#include <boost/python/extract.hpp>

#include "Reducer.h"
#include "SplitWorldException.h"

using namespace boost::python;
using namespace escript;


namespace escript
{
Reducer_ptr makeDataReducer(std::string type)
{
    MPI_Op op;
    if (type=="SUM")
    {
	op=MPI_SUM;
    }
    else
    {
	throw SplitWorldException("Unsupported operation for makeDataReducer.");
    }
    MPIDataReducer* m=new MPIDataReducer(op);
    return Reducer_ptr(m);
}

}

namespace
{

void combineData(Data& d1, const Data& d2, MPI_Op op)
{
    if (op==MPI_SUM)
    {
	d1+=d2;
    }
}


}



MPIDataReducer::MPIDataReducer(MPI_Op op)
  : reduceop(op)
{
    valueadded=false;
    if (op==MPI_SUM)
    {
	// deliberately left blank
    }
    else
    {
	throw SplitWorldException("Unsupported MPI_Op");
    }
}

void MPIDataReducer::setDomain(escript::Domain_ptr d)
{
    dom=d;
}


bool MPIDataReducer::valueCompatible(boost::python::object v)
{
    extract<const Data&> ex(v);
    if (!ex.check())
    {
	return false;
    }
    const Data& d=ex();
    if (d.getDomain()!=dom)
    {
	return false;	// the domains don't match
    }
    return true;
}


bool MPIDataReducer::reduceLocalValue(boost::python::object v, std::string& errstring)
{
    extract<Data&> ex(v);
    if (!ex.check())
    {
	errstring="reduceLocalValue: expected Data object. Got something else.";
	return false;
    }
    Data& d=ex();
    if (d.getDomain()!=dom)
    {
	errstring="reduceLocalValue: Got a Data object, but it was not using the SubWorld's domain.";
	return false;
    }
    d.expand();		// because I don't want to mess about with types of Data
    if (!valueadded)	// first value so answer becomes this one
    {
	value=d;
    }
    else
    {
	if (d.getFunctionSpace()!=value.getFunctionSpace())
	{
	    errstring="reduceLocalValue: FunctionSpaces for Data objects being combined must match.";
	    return false;
	}
	combineData(value, d, reduceop);
    }
    return true;
}

void MPIDataReducer::reset()
{
    valueadded=false;
    value=Data();
}

bool MPIDataReducer::checkRemoteCompatibility(esysUtils::JMPI& mpi_info, std::string& errstring)
{
    // since they can't add it unless it is using the proper domain, we need to 
    // check the following:
    //   FunctionSpace,
    //   Shape (for which we can just use dpps*dpsize)
    //   Domain version stamp [ if it exists ]
    
    // So We'll use an INT built up in the following way
    int check=0;
    check+=(value.getFunctionSpace().getTypeCode()%30)+1;				// guarantees we don't get all zeros
    check+=(value.getDataPointSize()*value.getNumDataPointsPerSample()%255)<<5;

    // still need to incorporate domain version into this
    // or are domains not mutable in any way that matters?
    int* rbuff=new int[mpi_info->size];
    for (int i=0;i<mpi_info->size;++i)
    {
	rbuff[i]=0;	// since this won't match any valid value we can use it as a failure check
    }
    if (MPI_Allgather(&check, 1, MPI_INT, rbuff, 1, MPI_INT, mpi_info->comm)!=MPI_SUCCESS)
    {
	errstring="MPI failure in checkRemoteCompatibility.";
	return false;
    }
    for (int i=0;i<mpi_info->size-1;++i)
    {
	if (rbuff[i]!=rbuff[i+1])
	{
	    std::ostringstream oss;
	    oss << "Incompatible value found for SubWorld " << i+1 << '.';
	    errstring=oss.str();
	    return false;
	}
    }
    return true;
}

// By the time this function is called, we know that all the values 
// are compatible
bool MPIDataReducer::reduceRemoteValues(esysUtils::JMPI& mpi_info)
{
    DataTypes::ValueType& vr=value.getExpandedVectorReference();
    Data result(0, value.getDataPointShape(), value.getFunctionSpace(), true);
    DataTypes::ValueType& rr=value.getExpandedVectorReference();
    if (MPI_Allreduce(&(vr[0]), &(rr[0]), vr.size(), MPI_DOUBLE, reduceop, mpi_info->comm)!=MPI_SUCCESS)
    {
	return false;
    }
    return true;
}
