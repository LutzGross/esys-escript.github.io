
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
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
