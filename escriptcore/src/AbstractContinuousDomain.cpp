
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

#include "AbstractContinuousDomain.h"
#include "Data.h"

using namespace boost::python;

namespace escript {

AbstractContinuousDomain::AbstractContinuousDomain()
{
}

AbstractContinuousDomain::AbstractContinuousDomain(JMPI mpiInfo)
    : AbstractDomain(mpiInfo)
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

void AbstractContinuousDomain::setToIntegrals(std::vector<DataTypes::real_t>& integrals,
                                              const escript::Data& arg) const
{
  throwStandardException("AbstractContinuousDomain::setToIntegrals<real_t>");
  return;
}

void AbstractContinuousDomain::setToIntegrals(std::vector<DataTypes::cplx_t>& integrals,
                                              const escript::Data& arg) const
{
  throwStandardException("AbstractContinuousDomain::setToIntegrals<cplx_t>");
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

#ifdef ESYS_HAVE_TRILINOS

void AbstractContinuousDomain::finaliseA(AbstractSystemMatrix& mat, 
                        Teuchos::RCP<Tpetra::CrsMatrix<DataTypes::cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> iz) const
{

}

void AbstractContinuousDomain::finaliseRhs(escript::Data& rhs,
                        Teuchos::RCP<Tpetra::CrsMatrix<DataTypes::cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> z) const
{

}

#endif

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

