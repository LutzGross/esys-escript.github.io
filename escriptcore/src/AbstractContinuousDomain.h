
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


#ifndef __ESCRIPT_ABSTRACTCONTINUOUSDOMAIN_H__
#define __ESCRIPT_ABSTRACTCONTINUOUSDOMAIN_H__

#include "system_dep.h"
#include "AbstractDomain.h"
#include "AbstractSystemMatrix.h"
#include "AbstractTransportProblem.h"
#include "DataTypes.h"

#ifdef ESYS_HAVE_TRILINOS
#include <trilinoswrap/CrsMatrixWrapper.h>
#include <trilinoswrap/TrilinosMatrixAdapter.h>
#include <trilinoswrap/types.h>
#include <Teuchos_Comm.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#endif

#include <string>
#include <vector>

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
class ESCRIPT_DLL_API AbstractContinuousDomain : public AbstractDomain
{

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
  AbstractContinuousDomain();

 protected:
  /**
     \brief
     Constructor with MPI info for AbstractContinuousDomain
  */
  explicit AbstractContinuousDomain(JMPI mpiInfo);

 public:
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
  virtual ~AbstractContinuousDomain();

  /**
     \brief
     Returns true if the given integer is a valid function space type
     for this domain.
  */
  virtual bool isValidFunctionSpaceType(int functionSpaceType) const;

  /**
     \brief
     Return a description for this domain
  */
  virtual std::string getDescription() const;

  /**
     \brief
     Return a continuous FunctionSpace code
  */
  virtual int getContinuousFunctionCode() const;

  /**
     \brief
     Return a continuous on reduced order FunctionSpace code
  */
  virtual int getReducedContinuousFunctionCode() const;

  /**
     \brief
     Return a function FunctionSpace code
  */
  virtual int getFunctionCode() const;

  /**
     \brief
     Return a function FunctionSpace code with reduced integration order
  */
  virtual int getReducedFunctionCode() const;

  /**
     \brief
     Return a function on boundary FunctionSpace code
  */
  virtual int getFunctionOnBoundaryCode() const;

  /**
     \brief
     Return a code for a function on boundary FunctionSpace with reduced integration order
  */
  virtual int getReducedFunctionOnBoundaryCode() const;


  /**
     \brief
     Return a FunctionOnContactZero code
  */
  virtual int getFunctionOnContactZeroCode() const;

  /**
     \brief
     Return a FunctionOnContactZero for reduced integration order code 
  */
  virtual int getReducedFunctionOnContactZeroCode() const;

  /**
     \brief
     Return a FunctionOnContactOne code
  */
  virtual int getFunctionOnContactOneCode() const;

  /**
     \brief
     Return a FunctionOnContactOne for reduced integration order code
  */
  virtual int getReducedFunctionOnContactOneCode() const;

  /**
     \brief
     Return a Solution code
  */
  virtual int getSolutionCode() const;

  /**
     \brief
     Return a ReducedSolution code
  */
  virtual int getReducedSolutionCode() const;

  /**
     \brief
     Return a DiracDeltaFunctions code
  */
  virtual int getDiracDeltaFunctionsCode() const;

  /**
     \brief
     return the identifier of the matrix type to be used for the global
     stiffness matrix when particular solver options are used.
  */
  virtual int getSystemMatrixTypeId(const boost::python::object& options) const;

  /**
     \brief
     return the identifier of the matrix type to be used for the global
     stiffness matrix when a particular solver package
     and symmetric matrix is used.
  */
  virtual int getTransportTypeId(int solver, int preconditioner, int package, bool symmetry) const;

  /**
     \brief
     copies the integrals of the function defined by arg into integrals.
     arg has to be defined on this.
     has to be implemented by the Domain Adapter.
  */
  virtual void setToIntegrals(std::vector<DataTypes::real_t>& integrals,
                              const escript::Data& arg) const;
  virtual void setToIntegrals(std::vector<DataTypes::cplx_t>& integrals,
                              const escript::Data& arg) const;

//  /**
//     \brief
//     Return the domain as const AbstractContinuousDomain&
//  */
//   static const  AbstractContinuousDomain& asAbstractContinuousDomain(const AbstractDomain& domain);




  /**
     \brief
     adds a PDE onto the stiffness matrix mat and a rhs 
  */
  virtual void addPDEToSystem(
                     AbstractSystemMatrix& mat, escript::Data& rhs,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C, 
                     const escript::Data& D, const escript::Data& X, const escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact, const escript::Data& y_contact, 
                     const escript::Data& d_dirac, const escript::Data& y_dirac) const;

  // Used by Oxley
  #ifdef ESYS_HAVE_TRILINOS
  // template<typename SC>
  // virtual Teuchos::RCP<Tpetra::CrsMatrix<DataTypes::cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> getZ() const;
  // template<typename SC>
  // virtual Teuchos::RCP<Tpetra::CrsMatrix<DataTypes::cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> getIZ() const;
  // template<typename SC>
  // virtual void finalise(AbstractSystemMatrix& mat, 
  //                       escript::Data& rhs,
  //                       Teuchos::RCP<Tpetra::CrsMatrix<DataTypes::real_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> z,
  //                       Teuchos::RCP<Tpetra::CrsMatrix<DataTypes::real_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> iz) const;
  virtual void finaliseA(AbstractSystemMatrix& mat, 
                        Teuchos::RCP<Tpetra::CrsMatrix<DataTypes::cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> iz) const;
  virtual void finaliseRhs(escript::Data& rhs,
                        Teuchos::RCP<Tpetra::CrsMatrix<DataTypes::cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> z) const;
  #endif

// We do not require this method at this level since the python side checks to ensure it exists
// before calling it.

//  /**
//     \brief
//     adds a PDE onto the lumped stiffness matrix matrix
//  */
//  virtual void addPDEToLumpedSystem(
//                     escript::Data& mat,
//                     const escript::Data& D, 
//                     const escript::Data& d) const;

  /**
     \brief
     adds a PDE onto the stiffness matrix mat and a rhs 
  */
  virtual void addPDEToRHS(escript::Data& rhs,
                     const escript::Data& X, const escript::Data& Y,
                     const escript::Data& y, const escript::Data& y_contact, const escript::Data& y_dirac) const;
  /**
     \brief
     adds a PDE onto a transport problem
  */

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

  virtual ATP_ptr newTransportProblem(
                      const int blocksize,
                      const escript::FunctionSpace& functionspace,
                      const int type) const;

  /**
     \brief
     Return the number of data points summed across all MPI processes
  */
  virtual DataTypes::dim_t getNumDataPointsGlobal() const;

  /**
     \brief
     Return the number of data points per sample, and the number of samples as a pair.
     \param functionSpaceCode Input -
  */
  virtual std::pair<int,DataTypes::dim_t> getDataShape(int functionSpaceCode) const;

  /**
     \brief
     assigns new location to the domain
  */
  virtual void setNewX(const escript::Data& arg);

  /**
     \brief
     \param full
  */
  virtual void Print_Mesh_Info(const bool full=false) const;
};

} // end of namespace

#endif // __ESCRIPT_ABSTRACTCONTINUOUSDOMAIN_H__

