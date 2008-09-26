
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "AbstractTransportProblem.h" 
#include "TransportProblemException.h"
#include "DataTypes.h"
#include "Data.h"
#include <iostream>


namespace escript {

AbstractTransportProblem::AbstractTransportProblem() {
    m_empty=1;
}

AbstractTransportProblem::AbstractTransportProblem(const double theta, const int blocksize,
                                                   const FunctionSpace& functionspace)
{
  if (blocksize<=0) 
     throw TransportProblemException("Error - negative block size of transport problem.");
  if ((theta<0.) || (theta>1.))
     throw TransportProblemException("Error - theta needs to be between 0. and 1..");

   m_empty=0;
   m_blocksize=blocksize;
   m_functionspace=functionspace;
   m_theta=theta;
}

AbstractTransportProblem::~AbstractTransportProblem() {
}

int AbstractTransportProblem::isEmpty() const {
   return m_empty;
}


Data AbstractTransportProblem::solve(Data& source, const double dt, const boost::python::dict& options) const
{
     if (isEmpty())
          throw TransportProblemException("Error - transport problem is empty.");
     if (dt<=0.)
          throw TransportProblemException("Error - dt needs to be positive.");
     if (source.getFunctionSpace()!=getFunctionSpace())
          throw TransportProblemException("Error - function space of transport problem and function space of source do not match.");
     if (source.getDataPointSize()!=getBlockSize())
          throw TransportProblemException("Error - block size of transport problem and source do not match.");
     DataTypes::ShapeType shape;
     if (getBlockSize()>1) shape.push_back(getBlockSize());
     Data out=Data(0.,shape,getFunctionSpace(),true);
     setToSolution(out,source,dt,options);
     return out;
}

void AbstractTransportProblem::setInitialValue(Data& u) const
{
     if (isEmpty())
          throw TransportProblemException("Error - transport problem is empty.");
     if (u.isEmpty())
          throw TransportProblemException("Error - empty initial value.");

     if ((getBlockSize()==1) && (u.getDataPointRank()>0) || (u.getDataPointRank()>1))
          throw TransportProblemException("Error - illegal rank of initial value.");

     if (u.getDataPointSize()!=getBlockSize())
          throw TransportProblemException("Error - block size of transport problem and initial value do not match.");

     Data u2=Data(u,getFunctionSpace());
     copyInitialValue(u2);
}
void AbstractTransportProblem::insertConstraint(Data& source, Data& q, Data& r) const
{
     source.expand();
     if (isEmpty())
          throw TransportProblemException("Error - transport problem is empty.");
     if (q.isEmpty()) {
          return;
     }
     if ((getBlockSize()==1) && (q.getDataPointRank()>0) || (q.getDataPointRank()>1))
          throw TransportProblemException("Error - illegal rank of constraint location.");
     if (q.getDataPointSize()!=getBlockSize())
          throw TransportProblemException("Error - block size of transport problem and constraint location don't match.");
     Data q2=Data(q,getFunctionSpace());

     if (r.isEmpty()) {
          Data r2=Data(0.,q.getDataPointShape(),getFunctionSpace());
          copyConstraint(source,q2,r2);
     } else {
        if ((getBlockSize()==1) && (r.getDataPointRank()>0) || (r.getDataPointRank()>1))
             throw TransportProblemException("Error - illegal rank of constraint value.");
        if (r.getDataPointSize()!=getBlockSize())
             throw TransportProblemException("Error - block size of transport problem and constraint value don't match.");
        Data r2=Data(r,getFunctionSpace());
        copyConstraint(source,q2,r2);
     }
}

void AbstractTransportProblem::copyConstraint(Data& source, Data& q, Data& r) const
{
    throw TransportProblemException("Error - copyConstraint is not available");
}
void AbstractTransportProblem::copyInitialValue(Data& u) const
{
    throw TransportProblemException("Error - copyInitialValue is not available");
}
void AbstractTransportProblem::setToSolution(Data& out,Data& source,const double dt, const boost::python::dict& options) const
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

}  // end of namespace
