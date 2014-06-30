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

#ifndef __ESCRIPT_REDUCER_H__
#define __ESCRIPT_REDUCER_H__

#include "esysUtils/Esys_MPI.h"
#include "escript/Data.h"
#include <boost/shared_ptr.hpp>

namespace escript
{

// There is currently no way to get a completely generic result out of this
class AbstractReducer
{
public:
    virtual ~AbstractReducer(){};
	// Is the value compatible with this reduction function?
	// does not guarantee the value is compatible with
	// other values added so far
    virtual bool valueCompatible(boost::python::object v)=0;
	// merge the parameter with the answer we already have
    virtual bool reduceLocalValue(boost::python::object v, std::string& errstring)=0;
	// clear previous result ready for a new set of reductions
    virtual void reset()=0;
    
	// converse with other subworlds to ensure subtype information matches
	// The main problem case here would be Data on different function spaces
	// same communicator requirements for reduceRemoteValues
	// Must give the same answer when called on any process in the subworlds
	// Must only be called on 
    virtual bool checkRemoteCompatibility(esysUtils::JMPI& mpi_info, std::string& errstring)=0; 
    
    
	// call to merge with values on other subworlds
	// It does not take a value argument because local values should have 
	// already been added with reduceLocal
	// Must only be called on participating SubWorlds
	// the mpi_info holds a communicator linking corresponding processes
	// in every participating subworld
    virtual bool reduceRemoteValues(esysUtils::JMPI& mpi_info)=0;
    
	// true if at least one localValue has been added
	// used to check if this subworld should participate in remote merges
    bool hasValue();
protected:
    bool valueadded;
    
};


typedef boost::shared_ptr<AbstractReducer> Reducer_ptr;

// Reduces using pointwise MPI operations
class MPIDataReducer : public AbstractReducer
{
public:
    MPIDataReducer(MPI_Op op);
    ~MPIDataReducer(){};
    
        // This is not a constructor parameter because 
        // if these are created outside the subworld, they won't have
        // access to a domain yet.
        // I also want SplitWorld to be able to set this
    void setDomain(escript::Domain_ptr d);
    bool valueCompatible(boost::python::object v);
    bool reduceLocalValue(boost::python::object v, std::string& errstring);
    void reset();
    bool checkRemoteCompatibility(esysUtils::JMPI& mpi_info, std::string& errstring);
    bool reduceRemoteValues(esysUtils::JMPI& mpi_info);
private:    
    escript::Data value;
    escript::const_Domain_ptr dom;
    MPI_Op reduceop;
};

Reducer_ptr makeDataReducer(std::string type);

}

#endif // __ESCRIPT_REDUCER_H__

