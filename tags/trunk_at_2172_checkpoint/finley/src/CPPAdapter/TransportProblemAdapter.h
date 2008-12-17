
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


#if !defined  finley_TransportProblemAdapter_H
#define finley_TransportProblemAdapter_H
#include "system_dep.h"

extern "C" {
#include "paso/SolverFCT.h"
#include "paso/Options.h"
}

#include "FinleyAdapterException.h"
#include "FinleyError.h"

#include "escript/AbstractTransportProblem.h"
#include "escript/Data.h"
#include "escript/UtilC.h"

#include <boost/python/dict.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/extract.hpp>

namespace finley {

class TransportProblemAdapter:public escript::AbstractTransportProblem {

/**
   \brief
   Wrapper for Paso_FCTransportProblem. 

   Description:
   Wrapper for Paso_FCTransportProblem.
*/

 public:

  /**
     /brief
     Default Constructor for TransportProblemAdapter.
     NB: Only throws an exception.
  */
  FINLEY_DLL_API
  TransportProblemAdapter();

  /**
     /brief
     Constructor for TransportProblemAdapter.
  */
  FINLEY_DLL_API
  TransportProblemAdapter(Paso_FCTransportProblem* transport_problem,
                          const double theta,
                          const int block_size,
                          const escript::FunctionSpace& functionspace);

  /**
     \brief
     Destructor for TransportProblemAdapter. As specified in the constructor
     this deallocates the pointer given to the constructor.
  */
  FINLEY_DLL_API
  ~TransportProblemAdapter();

  /**
     \brief
     Returns the pointer to the transport problem.
  */
  FINLEY_DLL_API
  Paso_FCTransportProblem* getPaso_FCTransportProblem() const;

  /**
     \brief
     Returns the transport problem as a const AbstractTransportProblem&.
  */
  inline const escript::AbstractTransportProblem& asAbstractTransportProblem() const
  {
     return dynamic_cast<const escript::AbstractTransportProblem&>(*this);
  }

  /**
     \brief
     Returns a transport problem as a const TransportProblemAdapter&.
  */
  inline static const TransportProblemAdapter& asTransportProblemAdapter(const AbstractTransportProblem& transportproblem)
  {
     return dynamic_cast<const TransportProblemAdapter&>(transportproblem);
  }

  /**
  *      \brief resets the transport operator typically as they have been updated.
  *        */
  FINLEY_DLL_API
  virtual void resetTransport() const;

  /**
  *      \brief returns a save time step size.
  */
  FINLEY_DLL_API
  virtual double getSafeTimeStepSize() const;

  /**
  *      \brief \brief returns the value for unlimited time step size.
  */
  FINLEY_DLL_API
  virtual double getUnlimitedTimeStepSize() const;

 protected:

 private:

    /**
    *      \brief
    *           sets solution out by time step dt.
    *             */
    FINLEY_DLL_API
    virtual void setToSolution(escript::Data& out,escript::Data& source,const double dt, const boost::python::dict& options) const;
   
   /**
   *      \brief
   *           copies the initial value into the problem
   *             */
    FINLEY_DLL_API
    virtual void copyInitialValue(escript::Data& u) const;

   /**
    *      \brief
    *           copy constraint u_{,t}=r where q>0  into the problem 
    *            it is assumed that q and r are not empty and has appropriate shape and function space.
    *                       */
    FINLEY_DLL_API
    virtual void copyConstraint(escript::Data& source, escript::Data& q, escript::Data& r) const;


   //
   // pointer to the externally created finley mesh - transport_problem.
   //
   boost::shared_ptr<Paso_FCTransportProblem> m_transport_problem;

};

} // end of namespace
#endif
