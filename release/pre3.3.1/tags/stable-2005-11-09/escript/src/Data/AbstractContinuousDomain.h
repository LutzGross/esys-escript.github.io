//$Id$
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

#if !defined escript_AbstractContinuousDomain_20040528_H
#define escript_AbstractContinuousDomain_20040528_H

#include "escript/Data/AbstractDomain.h"
#include "escript/Data/AbstractSystemMatrix.h"
#include "escript/Data/Data.h"

#include <boost/python/tuple.hpp>
#include <boost/python/object.hpp>

#include <string>
#include <vector>

namespace escript {

//
// Forward declaration
class FunctionSpace;

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
     Return a functon FunctionSpace code
  */
  virtual int getFunctionCode() const;

  /**
     \brief
     Return a function on boundary FunctionSpace code
  */
  virtual int getFunctionOnBoundaryCode() const;

  /**
     \brief
     Return a FunctionOnContactZero code
  */
  virtual int getFunctionOnContactZeroCode() const;

  /**
     \brief
     Return a FunctionOnContactOne code
  */
  virtual int getFunctionOnContactOneCode() const;

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
     Return a DiracDeltaFunction code
  */
  virtual int getDiracDeltaFunctionCode() const;

  /**
     \brief
     return the identifier of the matrix type to be used for the global
     stiffness matrix when a particular solver package
     and symmetric matrix is used.
  */
  virtual int getSystemMatrixTypeId(const int solver, const int package, const bool symmetry) const;

  /**
     \brief
     copies the integrals of the function defined by arg into integrals.
     arg has to be defined on this.
     has to be implemented by the Domain Adapter.
  */
  virtual void setToIntegrals(std::vector<double>& integrals,const escript::Data& arg) const;

  /**
     \brief
     Return the domain as const AbstractContinuousDomain&
  */
  static const  AbstractContinuousDomain& asAbstractContinuousDomain(const AbstractDomain& domain);

 protected:

 private:

};

} // end of namespace

#endif
