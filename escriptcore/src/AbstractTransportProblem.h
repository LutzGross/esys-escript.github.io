
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESCRIPT_ABSTRACTTRANSPORTPROBLEM_H__
#define __ESCRIPT_ABSTRACTTRANSPORTPROBLEM_H__

#include "system_dep.h"
#include "FunctionSpace.h"
#include "TransportProblemException.h"

#include <boost/python/object.hpp>

namespace escript {

//
// Forward declaration
class Data;

/**
   \brief
   Give a short description of what AbstractTransportProblem does.

   Description:
   Give a detailed description of AbstractTransportProblem

   Template Parameters:
   For templates describe any conditions that the parameters used in the
   template must satisfy
*/
class AbstractTransportProblem
{
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
  AbstractTransportProblem(const int blocksize,
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
     returns the solution u for a time step dt>0 with initial value u0 at time t=0
  */
  ESCRIPT_DLL_API
  Data solve(Data& u0, Data& source, const double dt, boost::python::object& options);


  /**
     \brief resets the transport operator typically as they have been updated.
  */
  ESCRIPT_DLL_API
  virtual void resetTransport(bool preserveSolverData) const;

  /**
     \brief
     inserts constraint u_{,t}=r where q>0  into the problem
  */
  ESCRIPT_DLL_API
  void insertConstraint(Data& source, Data& q, Data& r);
  /*
   *      \brief returns a safe time step size.
   */
  ESCRIPT_DLL_API
  virtual double getSafeTimeStepSize() const;
  /*
   *      \brief returns the value for unlimited time step size.
   */
  ESCRIPT_DLL_API
  virtual double getUnlimitedTimeStepSize() const;


 protected:

 private:

  /**
     \brief
     sets solution out by time step dt.
  */
  ESCRIPT_DLL_API
  virtual void setToSolution(Data& out, Data& u0, Data& source, double dt,
                             boost::python::object& options);

  /**
     \brief
     copy constraint u_{,t}=r where q>0  into the problem
     it can be assumed that q and r are not empty and have
     appropriate shape and function space.
  */
  ESCRIPT_DLL_API
  virtual void copyConstraint(Data& source, Data& q, Data& r);

  int m_empty;
  int m_blocksize;
  FunctionSpace m_functionspace;

};

typedef boost::shared_ptr<AbstractTransportProblem> ATP_ptr;

} // end of namespace

#endif // __ESCRIPT_ABSTRACTTRANSPORTPROBLEM_H__

