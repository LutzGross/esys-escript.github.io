/*****************************************************************************
*
* Copyright (c) 2014-2018 by The University of Queensland
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

#ifndef __ESCRIPT_ABSTRACTREDUCER_H__
#define __ESCRIPT_ABSTRACTREDUCER_H__

#include <escript/Data.h>
#include <escript/EsysMPI.h>

#include <boost/shared_ptr.hpp>

namespace escript
{

  
namespace reducerstatus
{
  
// Because these may be used in loops, the values must form a contiguous block (except ERROR)  
const unsigned char NONE=0;     // I have no value for this var and no interest in it
const unsigned char INTERESTED=1;   // I am interested in this variable but I have no value for it
const unsigned char OLD=2;  // I have a copy from elsewhere but no new values to contribute
const unsigned char OLDINTERESTED=3;    // interested but only have a cached copy (no new values)
const unsigned char NEW=4;  // I have a new value for this variable
const unsigned char ERROR='!';  // Something bad happened  
}
  
// There is currently no way to get a completely generic result out of this
class AbstractReducer
{
public:
    virtual ~AbstractReducer() {}
    // Is the value compatible with this reduction function?
    // does not guarantee the value is compatible with
    // other values added so far
    virtual bool valueCompatible(boost::python::object v)=0;
    // merge the parameter with the answer we already have
    virtual bool reduceLocalValue(boost::python::object v, std::string& errstring)=0;
    // clear previous result ready for a new set of reductions
    virtual void reset()=0;
    
    virtual std::string description()=0;
    
    // converse with other subworlds to ensure subtype information matches
    // The main problem case here would be Data on different function spaces
    // same communicator requirements for reduceRemoteValues
    // Must give the same answer when called on any process in the subworlds
    // Must only be called on 
    virtual bool checkRemoteCompatibility(JMPI& mpi_info, std::string& errstring)=0;
    // Some reducers need to know what domain they are operating in
    virtual void setDomain(Domain_ptr dom) {} 
    

#ifdef ESYS_MPI  
    // send from proc 0 in the communicator to all others
    // second param is true if we have rank o
    virtual bool groupSend(MPI_Comm& com, bool imsending)=0;
    
    // reduction with some procs submitting identity values
    virtual bool groupReduce(MPI_Comm& com, char mystate)=0;  
#endif  
    
    // call to merge with values on other subworlds
    // It does not take a value argument because local values should have 
    // already been added with reduceLocal
    // Must only be called on participating SubWorlds
    // the mpi_info holds a communicator linking corresponding processes
    // in every participating subworld
    virtual bool reduceRemoteValues(MPI_Comm& comm)=0;
    
    // true if at least one localValue has been added
    // used to check if this subworld should participate in remote merges
    bool hasValue();
    
    // true if reductions could fail for some reason other than MPI failure
    // for example SET type variables 
    virtual bool canClash();
    
    // Get a value for this variable from another process
    // This is not a reduction and will replace any existing value
    virtual bool recvFrom(int localid, int source, JMPI& mpiinfo)=0;

    // Send a value to this variable to another process
    // This is not a reduction and will replace any existing value    
    virtual bool sendTo(int localid, int target, JMPI& mpiinfo)=0;
    
    virtual double getDouble();
   
    virtual boost::python::object getPyObj()=0; 
    
    // notify the reducer that a new runJobs() call is being executed
    virtual void newRunJobs();

    virtual void clear();

    virtual void copyValueFrom(boost::shared_ptr<AbstractReducer>& src)=0;

protected:
    bool valueadded;
    bool had_an_export_this_round;
    static const int PARAMTAG;    
};


typedef boost::shared_ptr<AbstractReducer> Reducer_ptr;

}

#endif // __ESCRIPT_ABSTRACTREDUCER_H__

