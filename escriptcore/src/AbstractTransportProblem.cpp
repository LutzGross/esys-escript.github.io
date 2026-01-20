
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "AbstractTransportProblem.h"
#include "Data.h"
#include "DataTypes.h"
#include "TransportProblemException.h"

namespace bp = boost::python;

namespace escript {

AbstractTransportProblem::AbstractTransportProblem()
{
    m_empty=1;
}

AbstractTransportProblem::AbstractTransportProblem(int blocksize,
                                                   const FunctionSpace& functionspace)
    : m_empty(0),
      m_blocksize(blocksize),
      m_functionspace(functionspace)
{
    ESYS_ASSERT(blocksize>0, "non-positive block size given");
}

AbstractTransportProblem::~AbstractTransportProblem() {
}

int AbstractTransportProblem::isEmpty() const
{
   return m_empty;
}


Data AbstractTransportProblem::solve(Data& u0, Data& source, double dt,
                                     bp::object& options)
{
    if (isEmpty())
        throw TransportProblemException("Error - transport problem is empty.");
    if (dt<=0.)
        throw ValueError("dt needs to be positive.");
    if (source.getFunctionSpace()!=getFunctionSpace())
        throw ValueError("Function space of transport problem and function space of source do not match.");
    if (u0.getFunctionSpace()!=getFunctionSpace())
        throw ValueError("Function space of transport problem and function space of initial value do not match.");
    if (source.getDataPointSize()!=getBlockSize())
        throw ValueError("Block size of transport problem and source do not match.");
    if (u0.getDataPointSize()!=getBlockSize())
        throw ValueError("Block size of transport problem and initial value do not match.");

    DataTypes::ShapeType shape;
    if (getBlockSize()>1) shape.push_back(getBlockSize());
    Data out=Data(0.,shape,getFunctionSpace(),true);
    setToSolution(out, u0, source, dt, options);
    return out;
}

void AbstractTransportProblem::insertConstraint(Data& source, Data& q, Data& r)
{
    source.expand();
    if (isEmpty())
        throw TransportProblemException("insertConstraint(): Transport problem is empty.");
    if (q.isEmpty()) {
        return;
    }
    if (((getBlockSize()==1) && (q.getDataPointRank()>0)) || (q.getDataPointRank()>1))
        throw ValueError("insertConstraint(): illegal rank of constraint location.");
    if (q.getDataPointSize()!=getBlockSize())
        throw ValueError("insertConstraint(): Block size of transport problem and constraint location don't match.");
    Data q2=Data(q,getFunctionSpace());

    if (r.isEmpty()) {
        Data r2=Data(0.,q.getDataPointShape(),getFunctionSpace(), false);
        copyConstraint(source,q2,r2);
    } else {
        if (((getBlockSize()==1) && (r.getDataPointRank()>0)) || (r.getDataPointRank()>1))
            throw ValueError("Illegal rank of constraint value.");
        if (r.getDataPointSize()!=getBlockSize())
            throw ValueError("Block size of transport problem and constraint value don't match.");
        Data r2=Data(r,getFunctionSpace());
        copyConstraint(source,q2,r2);
    }
}

void AbstractTransportProblem::copyConstraint(Data& source, Data& q, Data& r)
{
    throw NotImplementedError("copyConstraint is not available");
}

void AbstractTransportProblem::setToSolution(Data& out, Data &u0, Data& source,
                                             double dt, bp::object& options)
{
    throw NotImplementedError("setToSolution is not available");
}
void AbstractTransportProblem::resetTransport(bool preserveSolverData) const
{
    throw NotImplementedError("resetProblem is not implemented.");
}
double AbstractTransportProblem::getSafeTimeStepSize() const
{
    throw NotImplementedError("getSafeTimeStepSize is not implemented.");
}
double AbstractTransportProblem::getUnlimitedTimeStepSize() const
{
    throw NotImplementedError("getUnlimitedTimeStepSize is not implemented.");
}

}  // end of namespace

