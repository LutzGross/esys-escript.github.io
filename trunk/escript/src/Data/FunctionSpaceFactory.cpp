/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/
#include "escript/Data/FunctionSpace.h"
#include "escript/Data/FunctionSpaceFactory.h"
#include "escript/Data/AbstractContinuousDomain.h"

namespace escript {

FunctionSpace continuousFunction(const AbstractDomain& domain) 
{
  const AbstractContinuousDomain& temp=AbstractContinuousDomain::asAbstractContinuousDomain(domain);
  return FunctionSpace(domain,temp.getContinuousFunctionCode());
}
 
FunctionSpace function(const AbstractDomain& domain) 
{
  const AbstractContinuousDomain& temp=AbstractContinuousDomain::asAbstractContinuousDomain(domain);
  return FunctionSpace(domain,temp.getFunctionCode());
}

FunctionSpace functionOnBoundary(const AbstractDomain& domain) 
{
  const AbstractContinuousDomain& temp=AbstractContinuousDomain::asAbstractContinuousDomain(domain);
  return FunctionSpace(domain,temp.getFunctionOnBoundaryCode());
}

FunctionSpace functionOnContactZero(const AbstractDomain& domain) 
{
  const AbstractContinuousDomain& temp=AbstractContinuousDomain::asAbstractContinuousDomain(domain);
  return FunctionSpace(domain,temp.getFunctionOnContactZeroCode());
}
 
FunctionSpace functionOnContactOne(const AbstractDomain& domain) 
{
  const AbstractContinuousDomain& temp=AbstractContinuousDomain::asAbstractContinuousDomain(domain);
  return FunctionSpace(domain,temp.getFunctionOnContactOneCode());
}

FunctionSpace solution(const AbstractDomain& domain) 
{
  const AbstractContinuousDomain& temp=AbstractContinuousDomain::asAbstractContinuousDomain(domain);
  return FunctionSpace(domain,temp.getSolutionCode());
}

FunctionSpace reducedSolution(const AbstractDomain& domain) 
{
  const AbstractContinuousDomain& temp=AbstractContinuousDomain::asAbstractContinuousDomain(domain);
  return FunctionSpace(domain,temp.getReducedSolutionCode());
}

FunctionSpace diracDeltaFunction(const AbstractDomain& domain)
{
  const AbstractContinuousDomain& temp=AbstractContinuousDomain::asAbstractContinuousDomain(domain);
  return FunctionSpace(domain,temp.getDiracDeltaFunctionCode());
}

}  // end of namespace
