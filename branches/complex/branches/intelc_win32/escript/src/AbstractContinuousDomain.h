//$Id$
/* 
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************

*/

#if !defined escript_AbstractContinuousDomain_20040528_H
#define escript_AbstractContinuousDomain_20040528_H

#include "system_dep.h"
#include "AbstractDomain.h"

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
     Return a functon FunctionSpace code
  */
  ESCRIPT_DLL_API 
  virtual int getFunctionCode() const;

  /**
     \brief
     Return a function on boundary FunctionSpace code
  */
  ESCRIPT_DLL_API 
  virtual int getFunctionOnBoundaryCode() const;

  /**
     \brief
     Return a FunctionOnContactZero code
  */
  ESCRIPT_DLL_API 
  virtual int getFunctionOnContactZeroCode() const;

  /**
     \brief
     Return a FunctionOnContactOne code
  */
  ESCRIPT_DLL_API 
  virtual int getFunctionOnContactOneCode() const;

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
     Return a DiracDeltaFunction code
  */
  ESCRIPT_DLL_API 
  virtual int getDiracDeltaFunctionCode() const;

  /**
     \brief
     return the identifier of the matrix type to be used for the global
     stiffness matrix when a particular solver package
     and symmetric matrix is used.
  */
  ESCRIPT_DLL_API 
  virtual int getSystemMatrixTypeId(const int solver, const int package, const bool symmetry) const;

  /**
     \brief
     copies the integrals of the function defined by arg into integrals.
     arg has to be defined on this.
     has to be implemented by the Domain Adapter.
  */
  ESCRIPT_DLL_API 
  virtual void setToIntegrals(std::vector<double>& integrals,const escript::Data& arg) const;

  /**
     \brief
     Return the domain as const AbstractContinuousDomain&
  */
  ESCRIPT_DLL_API 
  static const  AbstractContinuousDomain& asAbstractContinuousDomain(const AbstractDomain& domain);

 protected:

 private:

};

} // end of namespace

#endif
