
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

/* This file was extracted from finley's CPPAdapter then modified */

#ifndef __PASOWRAP_TRANSPORTPROBLEMADAPTER_H__
#define __PASOWRAP_TRANSPORTPROBLEMADAPTER_H__

#include "system_dep.h"

#include "paso/Transport.h"
#include "paso/Options.h"

#include "PasoException.h"

#include "escript/AbstractTransportProblem.h"
#include "escript/Data.h"
#include "escript/UtilC.h"

#include <boost/python/object.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/extract.hpp>

namespace paso {

PASOWRAP_DLL_API
class TransportProblemAdapter : public escript::AbstractTransportProblem
{

/**
   \brief
   Wrapper for paso::TransportProblem. 

   Description:
   Wrapper for paso::TransportProblem.
*/

public:

  /**
     /brief
     Default Constructor for TransportProblemAdapter.
     NB: Only throws an exception.
  */
  TransportProblemAdapter();

  /**
     /brief
     Constructor for TransportProblemAdapter.
  */
  TransportProblemAdapter(TransportProblem_ptr transport_problem,
                          int block_size,
                          const escript::FunctionSpace& functionspace);

  /**
     \brief
     Empty destructor for TransportProblemAdapter.
  */
  ~TransportProblemAdapter() {}

  /**
     \brief
     Returns the pointer to the transport problem.
  */
  TransportProblem_ptr getPaso_TransportProblem() const;

  /**
  *  \brief resets the transport operator typically as they have been updated.
  */
  virtual void resetTransport() const;

  /**
  *      \brief returns a save time step size.
  */
  virtual double getSafeTimeStepSize() const;

  /**
  *      \brief \brief returns the value for unlimited time step size.
  */
  virtual double getUnlimitedTimeStepSize() const;

  /**
     \brief
     returns the identifier of the transport problem type to be used
     when a particular solver, preconditioner and package is used
  */
  static int getTransportTypeId(const int solver, const int preconditioner,
          const int package, const bool symmetry, const esysUtils::JMPI& mpiInfo);

 protected:

 private:

    /**
    * \brief
    * sets solution out by time step dt.
    */
    virtual void setToSolution(escript::Data& out, escript::Data& u0,
                               escript::Data& source, double dt,
                               boost::python::object& options) const;
   

   /**
    * \brief
    * copy constraint u_{,t}=r where q>0  into the problem 
    * it is assumed that q and r are not empty and has appropriate shape
    * and function space.
    */
    virtual void copyConstraint(escript::Data& source, escript::Data& q,
                                escript::Data& r) const;


   //
   // shared pointer to the externally created transport_problem.
   //
   TransportProblem_ptr m_transport_problem;

};

} // end of namespace

#endif // __PASOWRAP_TRANSPORTPROBLEMADAPTER_H__

