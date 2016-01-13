
/* $Id$ */

/*******************************************************
 *
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#include "AbstractTransportProblem.h" 
#include "TransportProblemException.h"
#include "DataArrayView.h"
#include "Data.h"

namespace escript {

AbstractTransportProblem::AbstractTransportProblem() {
    m_empty=1;
}

AbstractTransportProblem::AbstractTransportProblem(const double theta,const double dt_max,
                                                   const int blocksize,
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
   m_dt_max=dt_max;
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
     DataArrayView::ShapeType shape;
     if (getBlockSize()>1) shape.push_back(getBlockSize());
     Data out=Data(0.,shape,getFunctionSpace(),true);
     setToSolution(out,source,dt,options);
     return out;
}

void AbstractTransportProblem::setInitialValue(Data& u) const
{
     if (isEmpty())
          throw TransportProblemException("Error - transport problem is empty.");
     if (u.getFunctionSpace()!=getFunctionSpace())
          throw TransportProblemException("Error - function space of transport problem and function space of initial value do not match.");
     if (u.getDataPointSize()!=getBlockSize())
          throw TransportProblemException("Error - block size of transport problem and initial value source do not match.");
     copyInitialValue(u);
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

}  // end of namespace
