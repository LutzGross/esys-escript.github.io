#include "BuckleyDomain.h"
#include "BuckleyException.h"

using namespace buckley;

namespace
{
  // Need to put some thought into how many spaces we actually need
    const int initialspace=1;   
    const int invalidspace=-1;
    
    const int notIMPLEMENTED=-1;
}


BuckleyDomain::BuckleyDomain(double x, double y, double z)
:ot(x,y,z)
{
}

BuckleyDomain::~BuckleyDomain()
{
}

bool BuckleyDomain::isValidFunctionSpaceType(int functionSpaceType) const
{
    return functionSpaceType==initialspace;
}

std::string BuckleyDomain::getDescription() const
{
    return "Dummy text for domain";
}

int BuckleyDomain::getContinuousFunctionCode() const
{
    return initialspace;    
}

int BuckleyDomain::getReducedContinuousFunctionCode() const
{
    return initialspace;
}

int BuckleyDomain::getFunctionCode() const
{
    return invalidspace;  
}

int BuckleyDomain::getReducedFunctionCode() const
{
    return invalidspace;  
}

int BuckleyDomain::getFunctionOnBoundaryCode() const
{
    return invalidspace;  
}

int BuckleyDomain::getReducedFunctionOnBoundaryCode() const
{
    return invalidspace;  
}

int BuckleyDomain::getFunctionOnContactZeroCode() const
{
    return invalidspace;  
}

int BuckleyDomain::getReducedFunctionOnContactZeroCode() const
{
    return invalidspace;
}

int BuckleyDomain::getFunctionOnContactOneCode() const
{
    return invalidspace;
}

int BuckleyDomain::getReducedFunctionOnContactOneCode() const
{
    return invalidspace;  
}

int BuckleyDomain::getSolutionCode() const
{
    return initialspace;  
}

int BuckleyDomain::getReducedSolutionCode() const
{
    return initialspace;  
}

// hahaha - no 
int BuckleyDomain::getDiracDeltaFunctionsCode() const
{
    return invalidspace  ;
}

int BuckleyDomain::getSystemMatrixTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const
{
    return notIMPLEMENTED;  
}

int BuckleyDomain::getTransportTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const
{
    return notIMPLEMENTED;  
}

void BuckleyDomain::setToIntegrals(std::vector<double>& integrals,const escript::Data& arg) const
{
    throw BuckleyException("Not Implemented");
}

void BuckleyDomain::addPDEToSystem(
                     escript::AbstractSystemMatrix& mat, escript::Data& rhs,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C, 
                     const escript::Data& D, const escript::Data& X, const escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact, const escript::Data& y_contact, 
                     const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
    throw BuckleyException("Not Implemented");
}

void BuckleyDomain::addPDEToRHS(escript::Data& rhs,
                     const escript::Data& X, const escript::Data& Y,
                     const escript::Data& y, const escript::Data& y_contact, const escript::Data& y_dirac) const
{
    throw BuckleyException("Not Implemented");
}

void BuckleyDomain::addPDEToTransportProblem(
                     escript::AbstractTransportProblem& tp, escript::Data& source, 
                     const escript::Data& M,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C,const  escript::Data& D,
                     const  escript::Data& X,const  escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact,const escript::Data& y_contact,
                     const escript::Data& d_dirac,const escript::Data& y_dirac) const
{
    throw BuckleyException("Not Implemented");
}

ASM_ptr BuckleyDomain::newSystemMatrix(
                      const int row_blocksize,
                      const escript::FunctionSpace& row_functionspace,
                      const int column_blocksize,
                      const escript::FunctionSpace& column_functionspace,
                      const int type) const
{
    throw BuckleyException("Not Implemented");
}

ATP_ptr BuckleyDomain::newTransportProblem(
                      const bool useBackwardEuler,
                      const int blocksize,
                      const escript::FunctionSpace& functionspace,
                      const int type) const
{
    throw BuckleyException("Not Implemented");
}

int BuckleyDomain::getNumDataPointsGlobal() const
{
    throw BuckleyException("Not Implemented");
}


  BUCKLEY_DLL_API
std::pair<int,int> BuckleyDomain::getDataShape(int functionSpaceCode) const
{
   throw BuckleyException("Not Implemented");  
  
}

BUCKLEY_DLL_API
void BuckleyDomain::setNewX(const escript::Data& arg)
{
    throw BuckleyException("This domain does not support changing coordinates");  
}


BUCKLEY_DLL_API
void BuckleyDomain::Print_Mesh_Info(const bool full) const
{
    throw BuckleyException("Not Implemented");
}