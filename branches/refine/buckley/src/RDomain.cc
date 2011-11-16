#include "RDomain.h"
#include "RDomainException.h"

using namespace buckley;

namespace
{
  // Need to put some thought into how many spaces we actually need
    const int initialspace=1;   
    const int invalidspace=-1;
    
    const int notIMPLEMENTED=-1;
}


RDomain::RDomain(double x, double y, double z)
:ot(x,y,z)
{
}

RDomain::~RDomain()
{
}

bool RDomain::isValidFunctionSpaceType(int functionSpaceType) const
{
    return functionSpaceType==initialspace;
}

std::string RDomain::getDescription() const
{
    return "Dummy text for domain";
}

int RDomain::getContinuousFunctionCode() const
{
    return initialspace;    
}

int RDomain::getReducedContinuousFunctionCode() const
{
    return initialspace;
}

int RDomain::getFunctionCode() const
{
    return invalidspace;  
}

int RDomain::getReducedFunctionCode() const
{
    return invalidspace;  
}

int RDomain::getFunctionOnBoundaryCode() const
{
    return invalidspace;  
}

int RDomain::getReducedFunctionOnBoundaryCode() const
{
    return invalidspace;  
}

int RDomain::getFunctionOnContactZeroCode() const
{
    return invalidspace;  
}

int RDomain::getReducedFunctionOnContactZeroCode() const
{
    return invalidspace;
}

int RDomain::getFunctionOnContactOneCode() const
{
    return invalidspace;
}

int RDomain::getReducedFunctionOnContactOneCode() const
{
    return invalidspace;  
}

int RDomain::getSolutionCode() const
{
    return initialspace;  
}

int RDomain::getReducedSolutionCode() const
{
    return initialspace;  
}

// hahaha - no 
int RDomain::getDiracDeltaFunctionsCode() const
{
    return invalidspace  ;
}

int RDomain::getSystemMatrixTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const
{
    return notIMPLEMENTED;  
}

int RDomain::getTransportTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const
{
    return notIMPLEMENTED;  
}

void RDomain::setToIntegrals(std::vector<double>& integrals,const escript::Data& arg) const
{
    throw RDomainException("Not Implemented");
}

void RDomain::addPDEToSystem(
                     escript::AbstractSystemMatrix& mat, escript::Data& rhs,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C, 
                     const escript::Data& D, const escript::Data& X, const escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact, const escript::Data& y_contact, 
                     const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
    throw RDomainException("Not Implemented");
}

void RDomain::addPDEToRHS(escript::Data& rhs,
                     const escript::Data& X, const escript::Data& Y,
                     const escript::Data& y, const escript::Data& y_contact, const escript::Data& y_dirac) const
{
    throw RDomainException("Not Implemented");
}

void RDomain::addPDEToTransportProblem(
                     escript::AbstractTransportProblem& tp, escript::Data& source, 
                     const escript::Data& M,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C,const  escript::Data& D,
                     const  escript::Data& X,const  escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact,const escript::Data& y_contact,
                     const escript::Data& d_dirac,const escript::Data& y_dirac) const
{
    throw RDomainException("Not Implemented");
}

ASM_ptr RDomain::newSystemMatrix(
                      const int row_blocksize,
                      const escript::FunctionSpace& row_functionspace,
                      const int column_blocksize,
                      const escript::FunctionSpace& column_functionspace,
                      const int type) const
{
    throw RDomainException("Not Implemented");
}

ATP_ptr RDomain::newTransportProblem(
                      const bool useBackwardEuler,
                      const int blocksize,
                      const escript::FunctionSpace& functionspace,
                      const int type) const
{
    throw RDomainException("Not Implemented");
}

int RDomain::getNumDataPointsGlobal() const
{
    throw RDomainException("Not Implemented");
}


  BUCKLEY_DLL_API
std::pair<int,int> RDomain::getDataShape(int functionSpaceCode) const
{
   throw RDomainException("Not Implemented");  
  
}

BUCKLEY_DLL_API
void RDomain::setNewX(const escript::Data& arg)
{
    throw RDomainException("This domain does not support changing coordinates");  
}


BUCKLEY_DLL_API
void RDomain::Print_Mesh_Info(const bool full) const
{
    throw RDomainException("Not Implemented");
}