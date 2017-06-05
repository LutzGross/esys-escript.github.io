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

#include "MPIDataReducer.h"
#include "SplitWorldException.h"

#include <limits>
#include <sstream>
#include <boost/python/extract.hpp>
#include <boost/scoped_array.hpp>

using namespace boost::python;
using namespace escript;

namespace {

void combineData(Data& d1, const Data& d2, MPI_Op op)
{
    if (op==MPI_SUM)
    {
        d1+=d2;
    } 
    else if (op==MPI_OP_NULL) 
    {
        throw SplitWorldException("Multiple 'simultaneous' attempts to export a 'SET' variable.");
    }
}

} // anonymous namespace


namespace escript {

Reducer_ptr makeDataReducer(std::string type)
{
    MPI_Op op;
    if (type=="SUM")
    {
        op=MPI_SUM;
    }
    else if (type=="SET")
    {
        op=MPI_OP_NULL;
    }
    else
    {
        throw SplitWorldException("Unsupported operation for makeDataReducer.");
    }
    MPIDataReducer* m=new MPIDataReducer(op);
    return Reducer_ptr(m);    
}

MPIDataReducer::MPIDataReducer(MPI_Op op)
  : reduceop(op), had_an_export_this_round(false)
{
    valueadded=false;
    if ((op==MPI_SUM) || (op==MPI_OP_NULL))
    {
        // deliberately left blank
    }
    else
    {
        throw SplitWorldException("Unsupported MPI_Op");
    }
}

void MPIDataReducer::newRunJobs()
{
    had_an_export_this_round=false;
}

void MPIDataReducer::setDomain(escript::Domain_ptr d)
{
    dom=d;
}

std::string MPIDataReducer::description()
{
    std::string op="SUM";
    if (reduceop==MPI_OP_NULL)
    {
        op="SET";
    }
    return "Reducer("+op+") for Data objects"; 
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
            return false;       // the domains don't match
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
    if (d.isEmpty())
    {
        errstring="reduceLocalValue: Got an empty Data object. Not allowed to reduce those.";
        return false;
    }
    if ((d.getDomain()!=dom) && (dom.get()!=0))
    {
        errstring="reduceLocalValue: Got a Data object, but it was not using the SubWorld's domain.";
        return false;
    }
    d.expand();         // because I don't want to mess about with types of Data
    if (!valueadded || !had_an_export_this_round)       // first value so answer becomes this one
    {
        value=d;
        dom=d.getDomain();
        had_an_export_this_round=true;
        valueadded=true;
    }
    else
    {
        if (reduceop==MPI_OP_NULL)
        {
            if (had_an_export_this_round) 
            {
                reset();
                errstring="reduceLocalValue: Multiple 'simultaneous' attempts to export a 'SET' variable.";
                return false;
            }
            value=d;
            dom=d.getDomain();
            had_an_export_this_round=true;
        }
        else
        { 
            had_an_export_this_round=true;
            if (d.getFunctionSpace()!=value.getFunctionSpace())
            {
                errstring="reduceLocalValue: FunctionSpaces for Data objects being combined must match.";
                return false;
            }
            combineData(value, d, reduceop);
        }
    }
    return true;
}

void MPIDataReducer::reset()
{
    valueadded=false;
    value=Data();
}

bool MPIDataReducer::checkRemoteCompatibility(JMPI& mpi_info, std::string& errstring)
{
#ifdef ESYS_MPI    
    // since they can't add it unless it is using the proper domain, we need to check 
    
    std::vector<unsigned> compat(6);
    getCompatibilityInfo(compat);

    // still need to incorporate domain version into this
    // or are domains not mutable in any way that matters?
    int* rbuff=new int[mpi_info->size*compat.size()];
    boost::scoped_array<int> dummy(rbuff);      // to ensure cleanup
    for (int i=0;i<mpi_info->size;++i)
    {
        rbuff[i]=0;     // since this won't match any valid value we can use it as a failure check
    }
    if (MPI_Allgather(&compat[0], compat.size(), MPI_UNSIGNED, rbuff, 
            compat.size(), MPI_UNSIGNED, mpi_info->comm)!=MPI_SUCCESS)
    {
        errstring="MPI failure in checkRemoteCompatibility.";
        return false;
    }
    for (int i=0;i<(mpi_info->size-1);++i)
    {
        if ((rbuff[i*compat.size()]==1) || (rbuff[(i+1)*compat.size()]==1))     // one of them doesn't have a value
        {
            continue;
        }
        for (int j=0;j<compat.size();++j)
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
bool MPIDataReducer::reduceRemoteValues(MPI_Comm& comm)
{
#ifdef ESYS_MPI
    DataTypes::RealVectorType& vr=value.getExpandedVectorReference();
    Data result(0, value.getDataPointShape(), value.getFunctionSpace(), true);
    DataTypes::RealVectorType& rr=result.getExpandedVectorReference();
    if (reduceop==MPI_OP_NULL)
    {
        reset();        // we can't be sure what the value should be
        return false;           // this will stop bad things happening but won't give an informative error message
    }
    if (MPI_Allreduce(&(vr[0]), &(rr[0]), vr.size(), MPI_DOUBLE, reduceop, comm)!=MPI_SUCCESS)
    {
        return false;
    }
    value=result;
    return true;
#else
    return true;
#endif
}

// populate a vector of ints with enough information to ensure two values are compatible
// or to construct a container for incomming data
// Format for this:
//  [0]    Type of Data:  {0 : error,  1:no value, 10: constant, 11:tagged, 12:expanded}
//  [1]    Functionspace type code
//  [2]    Only used for tagged --- gives the number of tags (which exist in the data object)
//  [3..6] Components of the shape  
//  [7]    Complexity: {0: real, 1:complex}
void MPIDataReducer::getCompatibilityInfo(std::vector<unsigned>& params)
{
    params.resize(8);
    for (int i=0;i<8;++i)
    {
        params[i]=0;
    }
    if (!valueadded)
    {
        params[0]=1;
        return;
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
    else        // This could be DataEmpty or some other weirdness but we won't allow that
    {
        params[0]=0;    // invalid type to send
        return;
    }    
    params[1]=value.getFunctionSpace().getTypeCode();
    params[2]=static_cast<unsigned>(value.getNumberOfTaggedValues());    
    const DataTypes::ShapeType& s=value.getDataPointShape();
    for (int i=0;i<s.size();++i)
    {
        params[3+i]=s[i];
    }
    params[7]=value.isComplex();
}


// Get a value for this variable from another process
// This is not a reduction and will replace any existing value
bool MPIDataReducer::recvFrom(int localid, int source, JMPI& mpiinfo)
{
#ifdef ESYS_MPI 
      // first we need to find out what we are expecting
    unsigned params[7];
    MPI_Status stat;
    if (MPI_Recv(params, 7, MPI_UNSIGNED, source, PARAMTAG, mpiinfo->comm, &stat)!=MPI_SUCCESS)
    {
        return false;
    }
    if (params[0]<10)   // the sender somehow tried to send something invalid
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
    if (params[0]==11)  // The Data is tagged so we need to work out what tags we need
    {
        // TODO:  Need to ship the tags and names over but for now just make sure there
        // are the same number of tags
        value.tag();
        
        DataTypes::RealVectorType dv(DataTypes::noValues(s), 0, 1);
        for (unsigned i=0;i<params[2];++i)
        {
            value.setTaggedValueFromCPP(static_cast<int>(i)+1, s, dv, 0);
        }
        return false;   // because I don't trust this yet
    }
#endif    
    return true;
}

// Send a value to this variable to another process
// This is not a reduction and will replace any existing value    
bool MPIDataReducer::sendTo(int localid, int target, JMPI& mpiinfo)
{
      if (!valueadded)
      {
          return false;         // May be misinterpreted as an MPI failure
      }
#ifdef ESYS_MPI  
      // first step is to let the other world know what sort of thing it needs to make
      if (value.isLazy())
      {
          value.resolve();
      }
      std::vector<unsigned> params;
      getCompatibilityInfo(params);
      if (MPI_Send(&params[0], 6, MPI_UNSIGNED, target, PARAMTAG, mpiinfo->comm)!=MPI_SUCCESS)
      {
          return false;
      }
      // now we have informed the other end of what happened
      // are we done or is there actually data to send
      if (params[0]<10)
      {
          return false;
      }
      
      if (value.isComplex())
      {
          DataTypes::cplx_t dummy=0;
            // at this point, we know there is data to send
          const DataTypes::cplx_t* vect=value.getDataRO(dummy);
            // now the receiver knows how much data it should be receive
            // need to make sure that we aren't trying to send data with no local samples
          if (vect!=0)
          {
              // MPI v3 has this first param as a const void* (as it should be)
              // Version on my machine expects void*
              // we don't require MPIv3 yet ... so we can't use MPI_CXX_DOUBLE_COMPLEX
              // We'll try just sending twice as many doubles
              //if (MPI_Send(const_cast<DataTypes::cplx_t*>(vect), value.getLength(), MPI_CXX_DOUBLE_COMPLEX, target, PARAMTAG, mpiinfo->comm)!=MPI_SUCCESS)
              if (MPI_Send(const_cast<DataTypes::cplx_t*>(vect), 2*value.getLength(), MPI_DOUBLE, target, PARAMTAG, mpiinfo->comm)!=MPI_SUCCESS)
              {
                  return false;
              }
          }
      }
      else
      {
          DataTypes::real_t dummy=0;
            // at this point, we know there is data to send
          const DataTypes::real_t* vect=value.getDataRO(dummy);
            // now the receiver knows how much data it should be receive
            // need to make sure that we aren't trying to send data with no local samples
          if (vect!=0)
          {
              // MPI v3 has this first param as a const void* (as it should be)
              // Version on my machine expects void*
              if (MPI_Send(const_cast<DataTypes::real_t*>(vect), value.getLength(), MPI_DOUBLE, target, PARAMTAG, mpiinfo->comm)!=MPI_SUCCESS)
              {
                  return false;
              }
          }
      }
#endif      
      return true;
}

boost::python::object MPIDataReducer::getPyObj()
{
    boost::python::object o(value);
    return o;
}


// send from proc 0 in the communicator to all others
// second argument is true if this rank is sending
bool MPIDataReducer::groupSend(MPI_Comm& comm, bool imsending)
{
      if (dom.get()==0)
      {
          return 0;     // trying to avoid throwing here
                        // this will still cause a lockup if it happens
      }
#ifdef ESYS_MPI
      if (imsending)
      {
          // first step is to let the other world know what sort of thing it needs to make
          if (value.isLazy())
          {
              value.resolve();
          }
          std::vector<unsigned> params;
          getCompatibilityInfo(params);
          if (MPI_Bcast(&params[0], params.size(), MPI_UNSIGNED, 0,comm)!=MPI_SUCCESS)
          {
              return false;
          }
            // now we have informed the other end of what happened
            // are we done or is there actually data to send
          if (params[0]<10)
          {
              return false;
          }
          
          if (value.isComplex())
          {
              DataTypes::cplx_t dummy=0;
                // at this point, we know there is data to send
              const DataTypes::cplx_t* vect=value.getDataRO(dummy);
                // now the receiver knows how much data it should be receive
                // need to make sure that we aren't trying to send data with no local samples
              if (vect!=0)
              {
                  // we don't require MPIv3 yet ... so we can't use MPI_CXX_DOUBLE_COMPLEX
                  // We'll try just sending twice as many doubles               
                  //if (MPI_Bcast(const_cast<DataTypes::cplx_t*>(vect), value.getLength(), MPI_CXX_DOUBLE_COMPLEX, 0, comm)!=MPI_SUCCESS)
                  if (MPI_Bcast(const_cast<DataTypes::cplx_t*>(vect), value.getLength()*2, MPI_DOUBLE, 0, comm)!=MPI_SUCCESS)
                  {
                      return false;
                  }
              }
          }
          else
          {
              DataTypes::real_t dummy=0;
                // at this point, we know there is data to send
              const DataTypes::real_t* vect=value.getDataRO(dummy);
                // now the receiver knows how much data it should be receive
                // need to make sure that we aren't trying to send data with no local samples
              if (vect!=0)
              {
                  if (MPI_Bcast(const_cast<DataTypes::real_t*>(vect), value.getLength(), MPI_DOUBLE, 0, comm)!=MPI_SUCCESS)
                  {
                      return false;
                  }
              }
          }
      }
      else      // we are receiving
      {
          bool createcplx=false;
            // first we need to find out what we are expecting
          unsigned params[8];
          if (MPI_Bcast(params, 8, MPI_UNSIGNED, 0, comm)!=MPI_SUCCESS)
          {
              return false;
          }
          if (params[0]<10)     // the sender somehow tried to send something invalid
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
          
          if (createcplx)       // we need to make a complex data
          {
              value=Data(0, s, fs, params[0]==12);
              value.complicate();
              if (params[0]==11)        // The Data is tagged so we need to work out what tags we need
              {
                  // TODO:  Need to ship the tags and names over but for now just make sure there
                  // are the same number of tags
                  value.tag();
                  
                  DataTypes::CplxVectorType dv(DataTypes::noValues(s), 0, 1);
                  for (unsigned i=0;i<params[2];++i)
                  {
                      value.setTaggedValueFromCPP(static_cast<int>(i)+1, s, dv, 0);
                  }
                  return false; // because I don't trust this yet
              }
              DataTypes::cplx_t* vect=&(value.getExpandedVectorReference(DataTypes::cplx_t(0))[0]);
              //if (MPI_Bcast(const_cast<DataTypes::cplx_t*>(vect), value.getLength(), MPI_CXX_DOUBLE_COMPLEX, 0, comm)!=MPI_SUCCESS)
              if (MPI_Bcast(const_cast<DataTypes::cplx_t*>(vect), value.getLength()*2, MPI_DOUBLE, 0, comm)!=MPI_SUCCESS)
              {
                  return false;
              }     
          }
          else
          {
              
              value=Data(0, s, fs, params[0]==12);
              if (params[0]==11)        // The Data is tagged so we need to work out what tags we need
              {
                  // TODO:  Need to ship the tags and names over but for now just make sure there
                  // are the same number of tags
                  value.tag();
                  
                  DataTypes::RealVectorType dv(DataTypes::noValues(s), 0, 1);
                  for (unsigned i=0;i<params[2];++i)
                  {
                      value.setTaggedValueFromCPP(static_cast<int>(i)+1, s, dv, 0);
                  }
                  return false; // because I don't trust this yet
              }
              DataTypes::real_t* vect=&(value.getExpandedVectorReference(0)[0]);
              if (MPI_Bcast(const_cast<DataTypes::real_t*>(vect), value.getLength(), MPI_DOUBLE, 0, comm)!=MPI_SUCCESS)
              {
                  return false;
              }
          }
          valueadded=true;
      }
#endif        
    return true;
}

// We assume compatible values at this point
bool MPIDataReducer::groupReduce(MPI_Comm& com, char mystate)
{
    throw SplitWorldException("groupReduce Not implemented yet.");
}

void MPIDataReducer::copyValueFrom(boost::shared_ptr<AbstractReducer>& src)
{
    MPIDataReducer* sr=dynamic_cast<MPIDataReducer*>(src.get());
    if (sr==0)
    {
        throw SplitWorldException("Source and destination need to be the same reducer types.");
    }
    if (sr->value.isEmpty())
    {
        throw SplitWorldException("Attempt to copy DataEmpty.");
    }
    if (sr==this)
    {
        throw SplitWorldException("Source and destination can not be the same variable.");
    }
    value.copy(sr->value);    
    valueadded=true;
}

bool MPIDataReducer::canClash()
{
    return (reduceop==MPI_OP_NULL);
}

} // namespace escript

