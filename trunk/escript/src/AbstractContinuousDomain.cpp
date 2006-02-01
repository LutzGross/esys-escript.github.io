/* $Id$ */
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

#include "AbstractContinuousDomain.h"
#include "Data.h"

using namespace boost::python;

namespace escript {

AbstractContinuousDomain::AbstractContinuousDomain()
{
}

AbstractContinuousDomain::~AbstractContinuousDomain()
{
}

bool AbstractContinuousDomain::isValidFunctionSpaceType(int functionSpaceType) const
{
  throwStandardException("AbstractContinuousDomain::isValidFunctionSpaceType");
  return false;
}

std::string AbstractContinuousDomain::getDescription() const
{
  throwStandardException("AbstractContinuousDomain::getDescription");
  return "";
}

int AbstractContinuousDomain::getContinuousFunctionCode() const
{
  throwStandardException("AbstractContinuousDomain::getContinuousFunctionCode");
  return 0;
}

int AbstractContinuousDomain::getFunctionCode() const
{
  throwStandardException("AbstractContinuousDomain::getFunctionCode");
  return 0;
}

int AbstractContinuousDomain::getFunctionOnBoundaryCode() const
{
  throwStandardException("AbstractContinuousDomain::getFunctionOnBoundaryCode");
  return 0;
}

int AbstractContinuousDomain::getFunctionOnContactZeroCode() const
{
  throwStandardException("AbstractContinuousDomain::getFunctionOnContactZeroCode");
  return 0;
}

int AbstractContinuousDomain::getFunctionOnContactOneCode() const
{
  throwStandardException("AbstractContinuousDomain::getFunctionOnContactOneCode");
  return 0;
}

int AbstractContinuousDomain::getSolutionCode() const
{
  throwStandardException("AbstractContinuousDomain::getSolutionCode");
  return 0;
}

int AbstractContinuousDomain::getReducedSolutionCode() const
{
  throwStandardException("AbstractContinuousDomain::getReducedSolutionCode");
  return 0;
}

int AbstractContinuousDomain::getDiracDeltaFunctionCode() const
{
  throwStandardException("AbstractContinuousDomain::getDiracDeltaFunctionCode");
  return 0;
}

void AbstractContinuousDomain::setToIntegrals(std::vector<double>& integrals,const escript::Data& arg) const
{
  throwStandardException("AbstractContinuousDomain::setToIntegrals");
  return;
}

int AbstractContinuousDomain::getSystemMatrixTypeId(const int solver, const int package, const bool symmetry) const 
{
   return 0;
}

const AbstractContinuousDomain& AbstractContinuousDomain::asAbstractContinuousDomain(const AbstractDomain& domain)
{
  return dynamic_cast<const AbstractContinuousDomain&>(domain);
}

}  // end of namespace
