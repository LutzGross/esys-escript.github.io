
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

int AbstractContinuousDomain::getDiracDeltaFunctionsCode() const
{
  throwStandardException("AbstractContinuousDomain::getDiracDeltaFunctionsCode");
  return 0;
}

void AbstractContinuousDomain::setToIntegrals(std::vector<double>& integrals, const escript::Data& arg) const
{
  throwStandardException("AbstractContinuousDomain::setToIntegrals");
  return;
}

int AbstractContinuousDomain::getSystemMatrixTypeId(const boost::python::object& options) const 
{
   return 0;
}

int AbstractContinuousDomain::getTransportTypeId(int solver, int precondioner, int package, bool symmetry) const 
{
   return 0;
}

void AbstractContinuousDomain::addPDEToSystem(
                     AbstractSystemMatrix& mat, escript::Data& rhs,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C, 
                     const escript::Data& D, const escript::Data& X, const escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact, const escript::Data& y_contact, const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
  throwStandardException("AbstractContinuousDomain::addPDEToSystem");
  return;
}

void AbstractContinuousDomain::addPDEToRHS(escript::Data& rhs,
                     const escript::Data& X, const escript::Data& Y,
                     const escript::Data& y, const escript::Data& y_contact, const escript::Data& y_dirac) const
{
  throwStandardException("AbstractContinuousDomain::addPDEToRHS");
  return;
}

void AbstractContinuousDomain::addPDEToTransportProblem(
                     AbstractTransportProblem& tp, escript::Data& source, 
                     const escript::Data& M,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C,const  escript::Data& D,
                     const  escript::Data& X,const  escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact,const escript::Data& y_contact, const escript::Data& d_dirac,const escript::Data& y_dirac) const
{
  throwStandardException("AbstractContinuousDomain::addPDEToTransportProblem");
  return;
}

ASM_ptr AbstractContinuousDomain::newSystemMatrix(
                      const int row_blocksize,
                      const escript::FunctionSpace& row_functionspace,
                      const int column_blocksize,
                      const escript::FunctionSpace& column_functionspace,
                      const int type) const
{
  throwStandardException("AbstractContinuousDomain::newSystemMatrix");
  return ASM_ptr();
}

ATP_ptr AbstractContinuousDomain::newTransportProblem(
                      const int blocksize,
                      const escript::FunctionSpace& functionspace,
                      const int type) const
{
  throwStandardException("AbstractContinuousDomain::newTransportProblem");
  return ATP_ptr();
}

DataTypes::dim_t AbstractContinuousDomain::getNumDataPointsGlobal() const
{
  throwStandardException("AbstractContinuousDomain::getNumDataPointsGlobal");
  return 1;
}

std::pair<int,DataTypes::dim_t> AbstractContinuousDomain::getDataShape(int functionSpaceCode) const
{
  throwStandardException("AbstractContinuousDomain::getDataShape");
  return std::pair<int,DataTypes::dim_t>(0,0);
}

void AbstractContinuousDomain::setNewX(const escript::Data& arg)
{
  throwStandardException("AbstractContinuousDomain::setNewX");
  return;
}

void AbstractContinuousDomain::Print_Mesh_Info(const bool full) const
{
  throwStandardException("AbstractContinuousDomain::Print_Mesh_Info");
  return;
}


}  // end of namespace

