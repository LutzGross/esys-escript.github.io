
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
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

#include "FunctionSpaceFactory.h"
#include "AbstractContinuousDomain.h"
#include "FunctionSpaceException.h"

namespace escript {

#define CTS_CHECK const AbstractContinuousDomain* temp=dynamic_cast<const AbstractContinuousDomain*>(&domain);\
if (temp==0)\
{\
   throw FunctionSpaceException("This method will only make FunctionSpaces for ContinuousDomains.");\
}


FunctionSpace continuousFunction(const AbstractDomain& domain) 
{
  CTS_CHECK
  return FunctionSpace(domain.getPtr(),temp->getContinuousFunctionCode());
}

FunctionSpace reducedContinuousFunction(const AbstractDomain& domain) 
{
  CTS_CHECK
  return FunctionSpace(domain.getPtr(),temp->getReducedContinuousFunctionCode());
}
 
FunctionSpace function(const AbstractDomain& domain) 
{
  CTS_CHECK
  return FunctionSpace(domain.getPtr(),temp->getFunctionCode());
}

FunctionSpace reducedFunction(const AbstractDomain& domain) 
{
  CTS_CHECK
  return FunctionSpace(domain.getPtr(),temp->getReducedFunctionCode());
}

FunctionSpace functionOnBoundary(const AbstractDomain& domain) 
{
  CTS_CHECK
  return FunctionSpace(domain.getPtr(),temp->getFunctionOnBoundaryCode());
}

FunctionSpace reducedFunctionOnBoundary(const AbstractDomain& domain) 
{
  CTS_CHECK
  return FunctionSpace(domain.getPtr(),temp->getReducedFunctionOnBoundaryCode());
}

FunctionSpace functionOnContactZero(const AbstractDomain& domain) 
{
  CTS_CHECK
  return FunctionSpace(domain.getPtr(),temp->getFunctionOnContactZeroCode());
}

FunctionSpace reducedFunctionOnContactZero(const AbstractDomain& domain) 
{
  CTS_CHECK
  return FunctionSpace(domain.getPtr(),temp->getReducedFunctionOnContactZeroCode());
}
 
FunctionSpace functionOnContactOne(const AbstractDomain& domain) 
{
  CTS_CHECK
  return FunctionSpace(domain.getPtr(),temp->getFunctionOnContactOneCode());
}

FunctionSpace reducedFunctionOnContactOne(const AbstractDomain& domain) 
{
  CTS_CHECK
  return FunctionSpace(domain.getPtr(),temp->getReducedFunctionOnContactOneCode());
}

FunctionSpace solution(const AbstractDomain& domain) 
{
  CTS_CHECK
  return FunctionSpace(domain.getPtr(),temp->getSolutionCode());
}

FunctionSpace reducedSolution(const AbstractDomain& domain) 
{
  CTS_CHECK
  return FunctionSpace(domain.getPtr(),temp->getReducedSolutionCode());
}

FunctionSpace diracDeltaFunctions(const AbstractDomain& domain)
{
  CTS_CHECK
  return FunctionSpace(domain.getPtr(),temp->getDiracDeltaFunctionsCode());
}

}  // end of namespace
