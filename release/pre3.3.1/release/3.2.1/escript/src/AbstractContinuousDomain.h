
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined escript_AbstractContinuousDomain_20040528_H
#define escript_AbstractContinuousDomain_20040528_H

#include "system_dep.h"
#include "AbstractDomain.h"

#include <string>
#include <vector>

#include "AbstractSystemMatrix.h"
#include "AbstractTransportProblem.h"

namespace escript {

//
// Forward declaration
class Data;

/**
   \brief
   AbstractContinuousDomain, base class for continuous domains.

   Description:
   AbstractContinuousDomain, base class for continuous domains.

   NOTE: Most of the virtual functions would be pure virtual except
   boost.python requires a non abstract base class.
*/
class AbstractContinuousDomain : public AbstractDomain {

 public:

  /**
     \brief
     Default constructor for AbstractContinuousDomain

     Description:
     Default constructor for AbstractContinuousDomain

     Preconditions:
     Describe any preconditions

     Throws:
     Describe any exceptions thrown
  */
  ESCRIPT_DLL_API 
  AbstractContinuousDomain();

  /**
     \brief
     Destructor for AbstractContinuousDomain

     Description:
     Destructor for AbstractContinuousDomain

     Preconditions:
     Describe any preconditions

     Throws:
     Describe any exceptions thrown
  */
  ESCRIPT_DLL_API 
  virtual ~AbstractContinuousDomain();

  /**
     \brief
     Returns true if the given integer is a valid function space type
     for this domain.
  */
  ESCRIPT_DLL_API 
  virtual bool isValidFunctionSpaceType(int functionSpaceType) const;

  /**
     \brief
     Return a description for this domain
  */
  ESCRIPT_DLL_API 
  virtual std::string getDescription() const;

  /**
     \brief
     Return a continuous FunctionSpace code
  */
  ESCRIPT_DLL_API 
  virtual int getContinuousFunctionCode() const;

  /**
     \brief
     Return a continuous on reduced order FunctionSpace code
  */
  ESCRIPT_DLL_API 
  virtual int getReducedContinuousFunctionCode() const;

  /**
     \brief
     Return a function FunctionSpace code
  */
  ESCRIPT_DLL_API 
  virtual int getFunctionCode() const;

  /**
     \brief
     Return a function FunctionSpace code with reduced integration order
  */
  ESCRIPT_DLL_API 
  virtual int getReducedFunctionCode() const;

  /**
     \brief
     Return a function on boundary FunctionSpace code
  */
  ESCRIPT_DLL_API 
  virtual int getFunctionOnBoundaryCode() const;

  /**
     \brief
     Return a code for a function on boundary FunctionSpace with reduced integration order
  */
  ESCRIPT_DLL_API 
  virtual int getReducedFunctionOnBoundaryCode() const;


  /**
     \brief
     Return a FunctionOnContactZero code
  */
  ESCRIPT_DLL_API 
  virtual int getFunctionOnContactZeroCode() const;

  /**
     \brief
     Return a FunctionOnContactZero for reduced integration order code 
  */
  ESCRIPT_DLL_API 
  virtual int getReducedFunctionOnContactZeroCode() const;

  /**
     \brief
     Return a FunctionOnContactOne code
  */
  ESCRIPT_DLL_API 
  virtual int getFunctionOnContactOneCode() const;

  /**
     \brief
     Return a FunctionOnContactOne for reduced integration order code
  */
  ESCRIPT_DLL_API 
  virtual int getReducedFunctionOnContactOneCode() const;

  /**
     \brief
     Return a Solution code
  */
  ESCRIPT_DLL_API 
  virtual int getSolutionCode() const;

  /**
     \brief
     Return a ReducedSolution code
  */
  ESCRIPT_DLL_API 
  virtual int getReducedSolutionCode() const;

  /**
     \brief
     Return a DiracDeltaFunctions code
  */
  ESCRIPT_DLL_API 
  virtual int getDiracDeltaFunctionsCode() const;

  /**
     \brief
     return the identifier of the matrix type to be used for the global
     stiffness matrix when a particular solver package
     and symmetric matrix is used.
  */
  ESCRIPT_DLL_API 
  virtual int getSystemMatrixTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const;

  /**
     \brief
     return the identifier of the matrix type to be used for the global
     stiffness matrix when a particular solver package
     and symmetric matrix is used.
  */
  ESCRIPT_DLL_API 
  virtual int getTransportTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const;

  /**
     \brief
     copies the integrals of the function defined by arg into integrals.
     arg has to be defined on this.
     has to be implemented by the Domain Adapter.
  */
  ESCRIPT_DLL_API 
  virtual void setToIntegrals(std::vector<double>& integrals,const escript::Data& arg) const;

//  /**
//     \brief
//     Return the domain as const AbstractContinuousDomain&
//  */
//   ESCRIPT_DLL_API 
//   static const  AbstractContinuousDomain& asAbstractContinuousDomain(const AbstractDomain& domain);




  /**
     \brief
     adds a PDE onto the stiffness matrix mat and a rhs 
  */
  ESCRIPT_DLL_API
  virtual void addPDEToSystem(
                     AbstractSystemMatrix& mat, escript::Data& rhs,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C, 
                     const escript::Data& D, const escript::Data& X, const escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact, const escript::Data& y_contact, 
                     const escript::Data& d_dirac, const escript::Data& y_dirac) const;

// We do not require this method at this level since the python side checks to ensure it exists
// before calling it.

//  /**
//     \brief
//     adds a PDE onto the lumped stiffness matrix matrix
//  */
//  ESCRIPT_DLL_API
//  virtual void addPDEToLumpedSystem(
//                     escript::Data& mat,
//                     const escript::Data& D, 
//                     const escript::Data& d) const;

  /**
     \brief
     adds a PDE onto the stiffness matrix mat and a rhs 
  */
  ESCRIPT_DLL_API
  virtual void addPDEToRHS(escript::Data& rhs,
                     const escript::Data& X, const escript::Data& Y,
                     const escript::Data& y, const escript::Data& y_contact, const escript::Data& y_dirac) const;
  /**
     \brief
     adds a PDE onto a transport problem
  */

  ESCRIPT_DLL_API
  virtual void addPDEToTransportProblem(
                     AbstractTransportProblem& tp, escript::Data& source, 
                     const escript::Data& M,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C,const  escript::Data& D,
                     const  escript::Data& X,const  escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact,const escript::Data& y_contact,
                     const escript::Data& d_dirac,const escript::Data& y_dirac) const;

  /**
     \brief
    creates a SystemMatrixAdapter stiffness matrix and initializes it with zeros:
  */
  ESCRIPT_DLL_API
  virtual ASM_ptr newSystemMatrix(
                      const int row_blocksize,
                      const escript::FunctionSpace& row_functionspace,
                      const int column_blocksize,
                      const escript::FunctionSpace& column_functionspace,
                      const int type) const;
  /**
   \brief 
    creates a TransportProblemAdapter 

  */

  ESCRIPT_DLL_API
  virtual ATP_ptr newTransportProblem(
                      const bool useBackwardEuler,
                      const int blocksize,
                      const escript::FunctionSpace& functionspace,
                      const int type) const;

  /**
     \brief
     Return the number of data points summed across all MPI processes
  */
  ESCRIPT_DLL_API
  virtual int getNumDataPointsGlobal() const;

  /**
     \brief
     Return the number of data points per sample, and the number of samples as a pair.
     \param functionSpaceCode Input -
  */
  ESCRIPT_DLL_API
  virtual std::pair<int,int> getDataShape(int functionSpaceCode) const;

  /**
     \brief
     assigns new location to the domain
  */
  ESCRIPT_DLL_API
  virtual void setNewX(const escript::Data& arg);

  /**
     \brief
     \param full
  */
  ESCRIPT_DLL_API
  virtual void Print_Mesh_Info(const bool full=false) const;


 protected:

 private:

};

} // end of namespace

#endif
