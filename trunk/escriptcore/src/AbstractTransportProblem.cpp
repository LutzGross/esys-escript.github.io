
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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
#include "AbstractTransportProblem.h" 
#include "TransportProblemException.h"
#include "DataTypes.h"
#include "Data.h"
#include <iostream>

namespace bp = boost::python;

namespace escript {

AbstractTransportProblem::AbstractTransportProblem() {
    m_empty=1;
}

AbstractTransportProblem::AbstractTransportProblem(const int blocksize,
                                                   const FunctionSpace& functionspace)
:m_functionspace(functionspace)
{
  if (blocksize<=0) 
     throw TransportProblemException("Error - negative block size of transport problem.");

   m_empty=0;
   m_blocksize=blocksize;
//    m_functionspace=functionspace;
}

AbstractTransportProblem::~AbstractTransportProblem() {
}

int AbstractTransportProblem::isEmpty() const {
   return m_empty;
}


Data AbstractTransportProblem::solve(Data& u0, Data& source, double dt,
                                     bp::object& options)
{
     if (isEmpty())
          throw TransportProblemException("Error - transport problem is empty.");
     if (dt<=0.)
          throw TransportProblemException("Error - dt needs to be positive.");
     if (source.getFunctionSpace()!=getFunctionSpace())
          throw TransportProblemException("Error - function space of transport problem and function space of source do not match.");
     if (u0.getFunctionSpace()!=getFunctionSpace())
          throw TransportProblemException("Error - function space of transport problem and function space of initial value do not match.");
     if (source.getDataPointSize()!=getBlockSize())
          throw TransportProblemException("Error - block size of transport problem and source do not match.");
     if (u0.getDataPointSize()!=getBlockSize())
          throw TransportProblemException("Error - block size of transport problem and initial value do not match.");

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
          throw TransportProblemException("Error - transport problem is empty.");
     if (q.isEmpty()) {
          return;
     }
     if (((getBlockSize()==1) && (q.getDataPointRank()>0)) || (q.getDataPointRank()>1))
          throw TransportProblemException("Error - illegal rank of constraint location.");
     if (q.getDataPointSize()!=getBlockSize())
          throw TransportProblemException("Error - block size of transport problem and constraint location don't match.");
     Data q2=Data(q,getFunctionSpace());

     if (r.isEmpty()) {
          Data r2=Data(0.,q.getDataPointShape(),getFunctionSpace());
          copyConstraint(source,q2,r2);
     } else {
        if (((getBlockSize()==1) && (r.getDataPointRank()>0)) || (r.getDataPointRank()>1))
             throw TransportProblemException("Error - illegal rank of constraint value.");
        if (r.getDataPointSize()!=getBlockSize())
             throw TransportProblemException("Error - block size of transport problem and constraint value don't match.");
        Data r2=Data(r,getFunctionSpace());
        copyConstraint(source,q2,r2);
     }
}

void AbstractTransportProblem::copyConstraint(Data& source, Data& q, Data& r)
{
    throw TransportProblemException("Error - copyConstraint is not available");
}

void AbstractTransportProblem::setToSolution(Data& out, Data &u0, Data& source,
                                             double dt, bp::object& options)
{
    throw TransportProblemException("Error - setToSolution is not available");
}
void AbstractTransportProblem::resetTransport() const
{
    throw TransportProblemException("Error - resetProblem is not implemented.");
}
double AbstractTransportProblem::getSafeTimeStepSize() const
{
    throw TransportProblemException("Error - getSafeTimeStepSize is not implemented.");
}
double AbstractTransportProblem::getUnlimitedTimeStepSize() const
{
    throw TransportProblemException("Error - getUnlimitedTimeStepSize is not implemented.");
}

}  // end of namespace
