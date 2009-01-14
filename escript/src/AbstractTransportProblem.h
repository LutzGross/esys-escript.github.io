
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
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
     \brief
     inserts constraint u_{,t}=r where q>0  into the problem.
  */
  ESCRIPT_DLL_API
  void insertConstraint(Data& source, Data& q, Data& r) const;
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
  virtual void setToSolution(Data& out,Data& source,const double dt, const boost::python::dict& options) const;

  /**
     \brief
     copies the initial value into the problem
  */
  ESCRIPT_DLL_API
  virtual void copyInitialValue(Data& u) const;
  /**
     \brief
     copy constraint u_{,t}=r where q>0  into the problem 
     it can be assumed that q and r are not empty and have  
     appropriate shape and function space.
  */
  ESCRIPT_DLL_API
  virtual void copyConstraint(Data& source, Data& q, Data& r) const;

  int m_empty;
  int m_blocksize;
  double m_theta;
  FunctionSpace m_functionspace;

};


} // end of namespace
#endif
