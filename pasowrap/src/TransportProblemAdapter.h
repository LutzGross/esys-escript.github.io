
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

/* This file was extracted from finley's CPPAdapter then modified */

#if !defined  TransportProblemAdapter_H
#define TransportProblemAdapter_H
#include "system_dep.h"

extern "C" {
#include "paso/Transport.h"
#include "paso/Options.h"
}

#include "PasoException.h"

#include "escript/AbstractTransportProblem.h"
#include "escript/Data.h"
#include "escript/UtilC.h"

#include <boost/python/object.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/extract.hpp>

namespace paso {

class TransportProblemAdapter:public escript::AbstractTransportProblem {

/**
   \brief
   Wrapper for Paso_TransportProblem. 

   Description:
   Wrapper for Paso_TransportProblem.
*/

 public:

  /**
     /brief
     Default Constructor for TransportProblemAdapter.
     NB: Only throws an exception.
  */
  PASOWRAP_DLL_API
  TransportProblemAdapter();

  /**
     /brief
     Constructor for TransportProblemAdapter.
  */
  PASOWRAP_DLL_API
  TransportProblemAdapter(Paso_TransportProblem* transport_problem,
                          const bool useBackwardEuler,
                          const int block_size,
                          const escript::FunctionSpace& functionspace);

  /**
     \brief
     Destructor for TransportProblemAdapter. As specified in the constructor
     this deallocates the pointer given to the constructor.
  */
  PASOWRAP_DLL_API
  ~TransportProblemAdapter();

  /**
     \brief
     Returns the pointer to the transport problem.
  */
  PASOWRAP_DLL_API
  Paso_TransportProblem* getPaso_TransportProblem() const;

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
  PASOWRAP_DLL_API
  virtual void resetTransport() const;

  /**
  *      \brief returns a save time step size.
  */
  PASOWRAP_DLL_API
  virtual double getSafeTimeStepSize() const;

  /**
  *      \brief \brief returns the value for unlimited time step size.
  */
  PASOWRAP_DLL_API
  virtual double getUnlimitedTimeStepSize() const;

  /**
     \brief
     returns the identifier of the transport problem type to be used
     when a particular solver, preconditioner and package is used
  */
  PASOWRAP_DLL_API
  static int getTransportTypeId(const int solver, const int preconditioner,
          const int package, const bool symmetry, Esys_MPIInfo* mpiInfo);

 protected:

 private:

    /**
    *      \brief
    *           sets solution out by time step dt.
    *             */
    PASOWRAP_DLL_API
    virtual void setToSolution(escript::Data& out,escript::Data& u0, escript::Data& source,const double dt, boost::python::object& options) const;
   

   /**
    *      \brief
    *           copy constraint u_{,t}=r where q>0  into the problem 
    *            it is assumed that q and r are not empty and has appropriate shape and function space.
    *                       */
    PASOWRAP_DLL_API
    virtual void copyConstraint(escript::Data& source, escript::Data& q, escript::Data& r, const double factor) const;


   //
   // pointer to the externally created transport_problem.
   //
   boost::shared_ptr<Paso_TransportProblem> m_transport_problem;

};

} // end of namespace
#endif
