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

#ifndef __ESCRIPT_NONREDUCEDVARIABLE_H__
#define __ESCRIPT_NONREDUCEDVARIABLE_H__

#include "AbstractReducer.h"
#include "Data.h"

namespace escript
{

// plugs into the import/export mechanism but stays on the 
// subworld it was created by (no actual reduction takes place)
class NonReducedVariable : public AbstractReducer
{
public:
    NonReducedVariable();
    ~NonReducedVariable();

    // This is not a constructor parameter because 
    // if these are created outside the subworld, they won't have
    // access to a domain yet.
    // I also want SplitWorld to be able to set this
    void setDomain(escript::Domain_ptr d);
    bool valueCompatible(boost::python::object v);
    bool reduceLocalValue(boost::python::object v, std::string& errstring);
    void reset();
    bool checkRemoteCompatibility(JMPI& mpi_info, std::string& errstring);
    
    void getCompatibilityInfo(std::vector<unsigned>& params);
    
      // talk to corresponding processes in other subworlds
    bool reduceRemoteValues(MPI_Comm& mpi_info);
    
      // human readable description
    std::string description();
    
    // Get a value for this variable from another process
    // This is not a reduction and will replace any existing value
    bool recvFrom(int localid, int source, JMPI& mpiinfo);

    // Send a value to this variable to another process
    // This is not a reduction and will replace any existing value    
    bool sendTo(int localid, int target, JMPI& mpiinfo);    
    double getDouble();
    virtual boost::python::object getPyObj(); 
    
        // send from proc 0 in the communicator to all others
    bool groupSend(MPI_Comm& com, bool imsending);
    
    // reduction with some procs submitting identity values
    bool groupReduce(MPI_Comm& com, char mystate);
    
    void copyValueFrom(boost::shared_ptr<AbstractReducer>& src);
    
private:    
    boost::python::object value;
    boost::python::object identity;
};

Reducer_ptr makeNonReducedVariable();

} // namespace escript

#endif // __ESCRIPT_NONREDUCEDVARIABLE_H__

