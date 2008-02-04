
/* $Id$ */

/*******************************************************
 *
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#if !defined  escript_AbstractTransportProblem_H
#define escript_AbstractTransportProblem_H
#include "system_dep.h"

#include "FunctionSpace.h"
#include "TransportProblemException.h"

#include <boost/python/dict.hpp>

//
// Forward declaration
class Data;

namespace escript {

/**
   \brief
   Give a short description of what AbstractTransportProblem does.

   Description:
   Give a detailed description of AbstractTransportProblem

   Template Parameters:
   For templates describe any conditions that the parameters used in the
   template must satisfy
*/
class AbstractTransportProblem {

 public:

  /**
     \brief
     Default constructor for AbstractTransportProblem

     Description:
     Default constructor for AbstractTransportProblem

     Preconditions:
     Describe any preconditions

     Throws:
     Describe any exceptions thrown
  */
  ESCRIPT_DLL_API
  AbstractTransportProblem();

  ESCRIPT_DLL_API
  AbstractTransportProblem(const double theta,
                           const int blocksize,
                           const FunctionSpace& functionspace);

  /**
    \brief
    Destructor.
  */
  ESCRIPT_DLL_API
  virtual ~AbstractTransportProblem();

  ESCRIPT_DLL_API
  int isEmpty() const;

  /**
    \brief
    returns the column function space
  */
  ESCRIPT_DLL_API
  inline FunctionSpace getFunctionSpace() const
  {
       if (isEmpty())
            throw TransportProblemException("Error - Transport Problem is empty.");
       return m_functionspace;
  }

  /**
    \brief
    returns the block size
  */
  ESCRIPT_DLL_API
  inline int getBlockSize() const
  {
       if (isEmpty())
            throw TransportProblemException("Error - Transport Problem is empty.");
       return m_blocksize;
  }

  /**
     \brief
     returns the solution u for a time step dt>0
  */
  ESCRIPT_DLL_API
  Data solve(Data& source, const double dt, const boost::python::dict& options) const;

  /**
     \brief
     sets the value for u at time t=0.
  */
  ESCRIPT_DLL_API
  void setInitialValue(Data& u) const;
  /**
     \brief resets the transport operator typically as they have been updated.
  */
  ESCRIPT_DLL_API
  virtual void resetTransport() const;

  /**
   *      \brief returns a save time step size.
   *        */
  ESCRIPT_DLL_API
  virtual double getSafeTimeStepSize() const;


 protected:

 private:

  /**
     \brief
     sets solution out by time step dt.
  */
  ESCRIPT_DLL_API
  virtual void setToSolution(Data& out,Data& source,const double dt, const boost::python::dict& options) const;

  /**
     \brief
     copies the initial value into the problem
  */
  ESCRIPT_DLL_API
  virtual void copyInitialValue(Data& u) const;

  int m_empty;
  int m_blocksize;
  double m_theta;
  FunctionSpace m_functionspace;

};


} // end of namespace
#endif
