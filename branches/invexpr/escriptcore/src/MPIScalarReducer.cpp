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

#define ESNEEDPYTHON
#include "esysUtils/first.h"


#include <sstream>
#include <limits>
#include <boost/python/extract.hpp>
#include <boost/scoped_array.hpp>

#include "MPIScalarReducer.h"
#include "SplitWorldException.h"

using namespace boost::python;
using namespace escript;


namespace escript
{

Reducer_ptr makeScalarReducer(std::string type)
{
    MPI_Op op;
    if (type=="SUM")
    {
	op=MPI_SUM;
    }
    else if (type=="MAX")
    {
	op=MPI_MAX;
    }
    else if (type=="MIN")
    {
	op=MPI_MIN;
    }
    else
    {
	throw SplitWorldException("Unsupported operation for makeScalarReducer.");
    }
    MPIScalarReducer* m=new MPIScalarReducer(op);
    return Reducer_ptr(m);    
}


}

namespace
{

void combineDouble(double& d1, const double d2, MPI_Op op)
{
    if (op==MPI_SUM)
    {
	d1+=d2;
    }  
    else if (op==MPI_MAX)
    {
	d1=(d2>d1)?d2:d1;
    }
    else if (op==MPI_MIN)
    {
	d1=(d2<d1)?d2:d1;      
    }    
}
}


MPIScalarReducer::MPIScalarReducer(MPI_Op op)
  : reduceop(op)
{
    valueadded=false;
    if (op==MPI_SUM)	// why not switch? because we don't know MPI_Op is scalar
    {
	identity=0;
    }
    else if (op==MPI_MAX)
    {
	identity=std::numeric_limits<double>::min();
    }
    else if (op==MPI_MIN)
    {
	identity=std::numeric_limits<double>::max();
    }
    else      
    {
	throw SplitWorldException("Unsupported MPI_Op");
    }
}

void MPIScalarReducer::setDomain(escript::Domain_ptr d)
{
    // deliberately left blank
}

std::string MPIScalarReducer::description()
{
    std::string op;
    if (reduceop==MPI_SUM)
    {
	op="SUM";
    }
    else if (reduceop==MPI_MAX)
    {
	op="MAX";
    } 
    else if (reduceop==MPI_MIN)
    {
	op="MIN";
    }
    else
    {
	throw SplitWorldException("Unsupported MPI reduction operation");
    }
    return "Reducer("+op+") for double scalars"; 
}

bool MPIScalarReducer::valueCompatible(boost::python::object v)
{
    extract<double> ex(v);
    if (!ex.check())
    {
	return false;
    }
    return true;
}


bool MPIScalarReducer::reduceLocalValue(boost::python::object v, std::string& errstring)
{
    extract<double> ex(v);
    if (!ex.check())
    {
	errstring="reduceLocalValue: expected double value. Got something else.";
	return false;
    }
    if (!valueadded)	// first value so answer becomes this one
    {
	value=ex();
	valueadded=true;
    }
    else
    {
	combineDouble(value, ex(), reduceop);
    }
    return true;
}

void MPIScalarReducer::reset()
{
    valueadded=false;
    value=0;
}

bool MPIScalarReducer::checkRemoteCompatibility(esysUtils::JMPI& mpi_info, std::string& errstring)
{
    return true;
}

// By the time this function is called, we know that all the values 
// are compatible
bool MPIScalarReducer::reduceRemoteValues(esysUtils::JMPI& mpi_info, bool active)
{
#ifdef ESYS_MPI
    if (!active)
    {
        value=identity;
    }
std::cout << "Value in " << value << std::endl;    
    if (MPI_Allreduce(&value, &value, 1, MPI_DOUBLE, reduceop, mpi_info->comm)!=MPI_SUCCESS)
    {
	return false;
    }
std::cout << "Value out " << value << std::endl;    
    return true;
#else
    return true;
#endif
}

// populate a vector of ints with enough information to ensure two values are compatible
// or to construct a container for incomming data
// Format for this:
//  [0]    Type of Data:  {0 : error,  1: DataEmpty, 10: constant, 11:tagged, 12:expanded}
//  [1]    Functionspace type code
//  [2]    Only used for tagged --- gives the number of tags (which exist in the data object)
//  [3..6] Components of the shape  
void MPIScalarReducer::getCompatibilityInfo(std::vector<unsigned>& params)
{
    params.resize(1);	// in case someone tries to do something with it
}


	// Get a value for this variable from another process
	// This is not a reduction and will replace any existing value
bool MPIScalarReducer::recvFrom(Esys_MPI_rank localid, Esys_MPI_rank source, esysUtils::JMPI& mpiinfo)
{
#ifdef ESYS_MPI  
    MPI_Status stat;
    if (MPI_Recv(&value, 1, MPI_DOUBLE, source, PARAMTAG, mpiinfo->comm, &stat)!=MPI_SUCCESS)
    {
	return false;
    }
#endif    
    return true;
}

	// Send a value to this variable to another process
	// This is not a reduction and will replace any existing value    
bool MPIScalarReducer::sendTo(Esys_MPI_rank localid, Esys_MPI_rank target, esysUtils::JMPI& mpiinfo)
{
#ifdef ESYS_MPI  
      if (MPI_Send(&value, 1, MPI_DOUBLE, target, PARAMTAG, mpiinfo->comm)!=MPI_SUCCESS)
      {
	  return false;
      }
#endif      
      return true;
}

double MPIScalarReducer::getDouble()
{
    return value;
}


boost::python::object MPIScalarReducer::getPyObj()
{
    boost::python::object o(value);
    return o;
}

#ifdef ESYS_MPI

	// send from proc 0 in the communicator to all others
bool MPIScalarReducer::groupSend(MPI_Comm& com)
{
    if (MPI_Bcast(&value, 1, MPI_DOUBLE, 0, com)==MPI_SUCCESS)
    {
	valueadded=true;
	return true;
    }
    return false;
}

bool MPIScalarReducer::groupReduce(MPI_Comm& com, char mystate)
{
    double answer=0;
    if (MPI_Allreduce((mystate==reducerstatus::NEW)?&value:&identity, &answer, 1, MPI_DOUBLE, reduceop, com)==MPI_SUCCESS)
    {
	value=answer;
	valueadded=true;
	return true;
    }
    return false;
}

#endif
