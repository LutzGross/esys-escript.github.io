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
#include <boost/scoped_array.hpp>

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

bool AbstractReducer::hasValue()
{
    return valueadded;
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
    extract<Data&> ex(v);
    if (!ex.check())
    {
	return false;
    }
    if (dom.get()!=0)
    {
	const Data& d=ex();
	if (d.getDomain().get()!=dom.get())
	{
	    return false;	// the domains don't match
	}
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
    if ((d.getDomain()!=dom) && (dom.get()!=0))
    {
	errstring="reduceLocalValue: Got a Data object, but it was not using the SubWorld's domain.";
	return false;
    }
    d.expand();		// because I don't want to mess about with types of Data
    if (!valueadded)	// first value so answer becomes this one
    {
	value=d;
	dom=d.getDomain();
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
#ifdef ESYS_MPI    
    // since they can't add it unless it is using the proper domain, we need to check 
    
    std::vector<unsigned> compat(6);
    getCompatibilityInfo(compat);

    // still need to incorporate domain version into this
    // or are domains not mutable in any way that matters?
    int* rbuff=new int[mpi_info->size*compat.size()];
    boost::scoped_array<int> dummy(rbuff);	// to ensure cleanup
    for (int i=0;i<mpi_info->size;++i)
    {
	rbuff[i]=0;	// since this won't match any valid value we can use it as a failure check
    }
    if (MPI_Allgather(&check, compat.size(), MPI_UNSIGNED, rbuff, 
	    compat.size(), MPI_UNSIGNED, mpi_info->comm)!=MPI_SUCCESS)
    {
	errstring="MPI failure in checkRemoteCompatibility.";
	return false;
    }
    for (int i=0;i<mpi_info->size-1;++i)
    {
	for (int j=0;j<compat.size();++i)
	{
	    if (rbuff[i*compat.size()+j]!=rbuff[(i+1)*compat.size()+j])
	    {
		std::ostringstream oss;
		oss << "Incompatible value found for SubWorld " << i+1 << '.';
		errstring=oss.str();
		return false;	      
	    } 
	}
    }
    return true;
#else
    return true;
#endif
}

// By the time this function is called, we know that all the values 
// are compatible
bool MPIDataReducer::reduceRemoteValues(esysUtils::JMPI& mpi_info)
{
#ifdef ESYS_MPI
    DataTypes::ValueType& vr=value.getExpandedVectorReference();
    Data result(0, value.getDataPointShape(), value.getFunctionSpace(), true);
    DataTypes::ValueType& rr=value.getExpandedVectorReference();
    if (MPI_Allreduce(&(vr[0]), &(rr[0]), vr.size(), MPI_DOUBLE, reduceop, mpi_info->comm)!=MPI_SUCCESS)
    {
	return false;
    }
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
void MPIDataReducer::getCompatibilityInfo(std::vector<unsigned>& params)
{
    params.resize(7);
    for (int i=0;i<7;++i)
    {
	params[0]=0;
    }
    if (value.isConstant())
    {
	params[0]=10;
    }
    else if (value.isTagged())
    {
	params[0]=11;
    }
    else if (value.isExpanded())
    {
	params[0]=12;
    }
    else	// This could be DataEmpty or some other weirdness but we won't allow that
    {
	params[0]=0;	// invalid type to send
    }    
    params[1]=value.getFunctionSpace().getTypeCode();
    params[2]=static_cast<unsigned>(value.getNumberOfTaggedValues());    
    const DataTypes::ShapeType& s=value.getDataPointShape();
    for (int i=0;i<s.size();++i)
    {
	params[3+i]=s[i];
    }    
}


	// Get a value for this variable from another process
	// This is not a reduction and will replace any existing value
bool MPIDataReducer::recvFrom(Esys_MPI_rank localid, Esys_MPI_rank source, esysUtils::JMPI& mpiinfo)
{
#ifdef ESYS_MPI	
      // first we need to find out what we are expecting
    unsigned params[6];
    if (MPI_Recv(params, 6, MPI_UNSIGNED, source, PARAMTAG, mpiinfo)!=MPI_SUCCESS)
    {
	return false;
    }
    if (params[0]<10)	// the sender somehow tried to send something invalid
    {
	return false;
    }

      // To avoid making this more complicated than it needs to be for now,
      // we will only allow expanded Data to be sent
    if (params[0]!=12)
    {
	return false;
    }
      // now we put the shape object together
    escript::DataTypes::ShapeType s;
    for (int i=0;i<4;++i)
    {
	if (params[3+i]>0)
	{
	    s.push_back(params[3+i]);
	}
	else
	{
	    break;
	}
    }
      // Now we need the FunctionSpace
    FunctionSpace fs=FunctionSpace(dom, static_cast<int>(params[1]));
    value=Data(0, s, fs, params[0]==12);
    /*
    if (params[0]==11)	// The Data is tagged so we need to work out what tags we need
    {
	// TODO:  Need to ship the tags and names over but for now just make sure there
	// are the same number of tags
	value.tag();
	
	DataVector dv(DataTypes::noValues(s), 0, 1);
	for (unsigned i=0;i<params[2];++i)
	{
	    value.setTaggedValueFromCPP(static_cast<int>(i)+1, s, dv, 0);
	}
    }
    */
    	// if this proc won't have any samples, the other end won't send any
    if (value.getNumSamples()>0)
    {
    	    // since we just created this data object, we no there is no
	    // sharing and it is safe to call this
    	DataTypes::ValueType& vec=value.getExpandedVectorReference();

	    // Now we get the values themselves
    	if (MPI_Recv(&vec[0], vec.size(), MPI_DOUBLE, source, DATATAG, mpiinfo)!=MPI_SUCCESS)
    	{
	    return false;
    	}
    }
#endif    
    return true;

}

	// Send a value to this variable to another process
	// This is not a reduction and will replace any existing value    
bool MPIDataReducer::sendTo(Esys_MPI_rank localid, Esys_MPI_rank target, esysUtils::JMPI& mpiinfo)
{
      // first step is to let the other world know what sort of thing it needs to make

      if (value.isLazy())
      {
	  value.resolve();
      }

      	// as noted in :recvFrom, we'll keep this simple for now and
	// only send expanded
      value.expand();
#ifdef ESYS_MPI      
      std::vector<unsigned> params;
      getCompatibilityInfo(params);
      if (MPI_Send(&params[0], 6, MPI_UNSIGNED, target, PARAMTAG, mpiinfo)!=MPI_SUCCESS)
      {
	  return false;
      }
	// now we have informed the other end of what happened
	// are we done or is there actually data to send
      if (params[0]<10)
      {
	  return false;
      }
	// at this point, we know there is data to send
      const DataAbstract::ValueType::value_type* vect=value.getDataRO();
	// now the receiver knows how much data it should be receive
	// need to make sure that we aren't trying to send data with no local samples
      if (vect>0)
      {
	  if (MPI_Send(vect, value.getLength(), MPI_DOUBLE, target, PARAMTAG, mpiinfo)!=MPI_SUCCESS)
	  {
	      return false;
	  }
      }
#endif      
      return true;
}
