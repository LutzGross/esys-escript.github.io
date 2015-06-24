/*****************************************************************************
*
* Copyright (c) 2015 by The University of Queensland
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


#include "NonReducedVariable.h"
#include "SplitWorldException.h"

using namespace escript;

NonReducedVariable::NonReducedVariable()
{
    valueadded=false;
}

NonReducedVariable::~NonReducedVariable()
{
}

void NonReducedVariable::setDomain(escript::Domain_ptr d)
{
	
}


// SInce there is no remote transfer, we don't need to check this
bool NonReducedVariable::valueCompatible(boost::python::object v)
{
    return true;
}

// Any new export, replaces the old value
bool NonReducedVariable::reduceLocalValue(boost::python::object v, std::string& errstring)
{
    value=v;
    valueadded=true;
    return true;
}

void NonReducedVariable::reset()
{
    value=boost::python::object();
    valueadded=false;
}

// Since we aren't actually don't a check here, this call won't function
// as a barrier like other implementations of this method
bool NonReducedVariable::checkRemoteCompatibility(esysUtils::JMPI& mpi_info, std::string& errstring)
{
    return true;
}

void NonReducedVariable::getCompatibilityInfo(std::vector<unsigned>& params)
{
    // empty
}

bool NonReducedVariable::reduceRemoteValues(esysUtils::JMPI& mpi_info, bool active)
{
    return true;
}

std::string NonReducedVariable::description()
{
    return "Non-Reduced Variable.";
}

bool NonReducedVariable::recvFrom(Esys_MPI_rank localid, Esys_MPI_rank source, esysUtils::JMPI& mpiinfo)
{
    return true;
}

bool NonReducedVariable::sendTo(Esys_MPI_rank localid, Esys_MPI_rank source, esysUtils::JMPI& mpiinfo)
{
    return true;
}

double NonReducedVariable::getDouble()
{
    throw SplitWorldException("No double value from this type.");
}

boost::python::object NonReducedVariable::getPyObj()
{
    return value;
}

bool NonReducedVariable::groupSend(MPI_Comm& com, bool imsending)
{
    return true;
}

bool NonReducedVariable::groupReduce(MPI_Comm& com, char mystate)
{
    return true;
}

void NonReducedVariable::copyValueFrom(boost::shared_ptr<AbstractReducer>& src)
{
    NonReducedVariable* sr=dynamic_cast<NonReducedVariable*>(src.get());
    if (sr==0)
    {
	throw SplitWorldException("Source and destination need to be the same reducer types.");
    }
    value=sr->value;
    valueadded=true;
}


namespace escript
{
Reducer_ptr makeNonReducedVariable()
{
    NonReducedVariable* m=new NonReducedVariable();
    return Reducer_ptr(m);

}

}
