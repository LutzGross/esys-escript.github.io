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

#ifdef PASO_MPI
#include <mpi.h>
#endif
#include "TransportProblemAdapter.h" 
#include "SystemMatrixAdapter.h" 

using namespace std;

namespace finley {

struct null_deleter
{
  void operator()(void const *ptr) const
  {
  }
};


TransportProblemAdapter::TransportProblemAdapter()
{
   throw FinleyAdapterException("Error - Illegal to generate default TransportProblemAdapter.");
}


TransportProblemAdapter::TransportProblemAdapter(Paso_FCTransportProblem* transport_problem,
                                                 const double theta,
                                                 const double dt_max,
                                                 const int block_size,
                                                 const escript::FunctionSpace& functionspace):
AbstractTransportProblem(theta, dt_max, block_size, functionspace)
{
    m_transport_problem.reset(transport_problem,null_deleter());
}

TransportProblemAdapter::~TransportProblemAdapter()
{ 
    if (m_transport_problem.unique()) {
        Paso_FCTransportProblem* transp=m_transport_problem.get();
        Paso_FCTransportProblem_free(transp);
    }
}

Paso_FCTransportProblem* TransportProblemAdapter::getPaso_FCTransportProblem() const 
{
   return m_transport_problem.get();
}


void TransportProblemAdapter::setToSolution(escript::Data& out, escript::Data& source, const double dt, const boost::python::dict& options) const
{
    Paso_FCTransportProblem* transp=getPaso_FCTransportProblem();
    Paso_Options paso_options;
    SystemMatrixAdapter::dictToPasoOptions(&paso_options,options);
    if ( out.getDataPointSize()  != getBlockSize()) {
     throw FinleyAdapterException("solve : block size of solution does not match block size of transport problems.");
    } else if ( source.getDataPointSize() != getBlockSize()) {
     throw FinleyAdapterException("solve : block size of source term does not match block size of transport problems.");
    } else if ( out.getFunctionSpace()  != getFunctionSpace()) {
     throw FinleyAdapterException("solve : function spaces of solution and of transport problem don't match.");
    } else if (source.getFunctionSpace() != getFunctionSpace()) {
     throw FinleyAdapterException("solve : function spaces of source term and of transport problem don't match.");
    } else if (dt<=0.) {
     throw FinleyAdapterException("solve : time increment dt needs to be positive.");
    }
    out.expand();
    source.expand();
    double* out_dp=out.getSampleData(0);
    double* source_dp=source.getSampleData(0);
    Paso_SolverFCT_solve(transp,out_dp,dt,source_dp,&paso_options);
    checkPasoError();
}

void TransportProblemAdapter::resetTransport() const
{
   Paso_FCTransportProblem* transp = getPaso_FCTransportProblem();
   throw FinleyAdapterException("resetTransport() not implemented yet.");
   // 
   //
   // Paso_FCTransportProblem_setValues(transp,0.);
   // Paso_solve_free(transp);
   checkPasoError();
}

void TransportProblemAdapter::copyInitialValue(escript::Data& u) const
{
    Paso_FCTransportProblem* transp=getPaso_FCTransportProblem();
    if ( u.getDataPointSize()  != getBlockSize()) {
     throw FinleyAdapterException("copyInitialValue : block size of solution does not match block size of transport problems.");
    } else if ( u.getFunctionSpace()  != getFunctionSpace()) {
     throw FinleyAdapterException("copyInitialValue : function spaces of solution and of transport problem don't match.");
    }
    u.expand();
    double* u_dp=u.getSampleData(0);
    Paso_FCTransportProblem_checkinSolution( transp,u_dp);
    checkPasoError();
}


}  // end of namespace
