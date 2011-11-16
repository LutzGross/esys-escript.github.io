
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

#ifndef __RIPLEY_TRANSPORTPROBLEMADAPTER_H__
#define __RIPLEY_TRANSPORTPROBLEMADAPTER_H__

#include <ripley/system_dep.h>

extern "C" {
#include <paso/Transport.h>
#include <paso/Options.h>
}

#include <escript/AbstractTransportProblem.h>
#include <escript/Data.h>
#include <escript/UtilC.h>

#include <boost/python/object.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/extract.hpp>

namespace ripley {

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
    */
    RIPLEY_DLL_API
    TransportProblemAdapter();

    /**
       /brief
       Constructor for TransportProblemAdapter.
    */
    RIPLEY_DLL_API
    TransportProblemAdapter(Paso_TransportProblem* transport_problem,
                            const bool useBackwardEuler,
                            const int block_size,
                            const escript::FunctionSpace& functionspace);

    /**
       \brief
       Destructor for TransportProblemAdapter. As specified in the constructor
       this deallocates the pointer given to the constructor.
    */
    RIPLEY_DLL_API
    ~TransportProblemAdapter();

    /**
       \brief
       returns the pointer to the transport problem.
    */
    RIPLEY_DLL_API
    Paso_TransportProblem* getPaso_TransportProblem() const;

    /**
       \brief
       returns the transport problem as a const AbstractTransportProblem&.
    */
    inline const escript::AbstractTransportProblem& asAbstractTransportProblem() const
    {
        return dynamic_cast<const escript::AbstractTransportProblem&>(*this);
    }

    /**
       \brief
       returns a transport problem as a const TransportProblemAdapter&.
    */
    inline static const TransportProblemAdapter& asTransportProblemAdapter(const AbstractTransportProblem& transportproblem)
    {
        return dynamic_cast<const TransportProblemAdapter&>(transportproblem);
    }

    /**
    *  \brief
    *  resets the transport operator typically as they have been updated.
    */
    RIPLEY_DLL_API
    virtual void resetTransport() const;

    /**
    *  \brief returns a save time step size.
    */
    RIPLEY_DLL_API
    virtual double getSafeTimeStepSize() const;

    /**
    *  \brief returns the value for unlimited time step size.
    */
    RIPLEY_DLL_API
    virtual double getUnlimitedTimeStepSize() const;

 private:
    /**
    *  \brief
    *  sets solution out by time step dt.
    */
    RIPLEY_DLL_API
    virtual void setToSolution(escript::Data& out,escript::Data& u0, escript::Data& source,const double dt, boost::python::object& options) const;
   

    /**
    *  \brief
    *  copies constraint u_{,t}=r where q>0 into the problem.
    *  It is assumed that q and r are not empty and have appropriate shape
    *  and function space.
    */
    RIPLEY_DLL_API
    virtual void copyConstraint(escript::Data& source, escript::Data& q, escript::Data& r, const double factor) const;

    //
    // pointer to the externally created ripley mesh - transport_problem.
    boost::shared_ptr<Paso_TransportProblem> m_transport_problem;
};

} // end of namespace ripley

#endif // __RIPLEY_TRANSPORTPROBLEMADAPTER_H__

