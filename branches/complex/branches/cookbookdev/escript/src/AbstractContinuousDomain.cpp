
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


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

int AbstractContinuousDomain::getReducedContinuousFunctionCode() const
{
  throwStandardException("AbstractContinuousDomain::getReducedContinuousFunctionCode");
  return 0;
}

int AbstractContinuousDomain::getFunctionCode() const
{
  throwStandardException("AbstractContinuousDomain::getFunctionCode");
  return 0;
}

int AbstractContinuousDomain::getReducedFunctionCode() const
{
  throwStandardException("AbstractContinuousDomain::getReducedFunctionCode");
  return 0;
}

int AbstractContinuousDomain::getFunctionOnBoundaryCode() const
{
  throwStandardException("AbstractContinuousDomain::getFunctionOnBoundaryCode");
  return 0;
}

int AbstractContinuousDomain::getReducedFunctionOnBoundaryCode() const
{
  throwStandardException("AbstractContinuousDomain::getReducedFunctionOnBoundaryCode");
  return 0;
}

int AbstractContinuousDomain::getFunctionOnContactZeroCode() const
{
  throwStandardException("AbstractContinuousDomain::getFunctionOnContactZeroCode");
  return 0;
}

int AbstractContinuousDomain::getReducedFunctionOnContactZeroCode() const
{
  throwStandardException("AbstractContinuousDomain::getReducedFunctionOnContactZeroCode");
  return 0;
}

int AbstractContinuousDomain::getFunctionOnContactOneCode() const
{
  throwStandardException("AbstractContinuousDomain::getFunctionOnContactOneCode");
  return 0;
}

int AbstractContinuousDomain::getReducedFunctionOnContactOneCode() const
{
  throwStandardException("AbstractContinuousDomain::getReducedFunctionOnContactOneCode");
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

int AbstractContinuousDomain::getSystemMatrixTypeId(const int solver, const int precondioner, const int package, const bool symmetry) const 
{
   return 0;
}

int AbstractContinuousDomain::getTransportTypeId(const int solver, const int precondioner, const int package, const bool symmetry) const 
{
   return 0;
}

// const AbstractContinuousDomain& AbstractContinuousDomain::asAbstractContinuousDomain(const AbstractDomain& domain)
// {
//   return dynamic_cast<const AbstractContinuousDomain&>(domain);
// }

}  // end of namespace
