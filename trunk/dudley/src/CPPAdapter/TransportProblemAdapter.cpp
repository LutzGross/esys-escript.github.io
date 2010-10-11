
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


#ifdef ESYS_MPI
#include <mpi.h>
#endif
#include "TransportProblemAdapter.h" 
#include "SystemMatrixAdapter.h" 

using namespace std;

namespace dudley {

struct null_deleter
{
  void operator()(void const *ptr) const
  {
  }
};


TransportProblemAdapter::TransportProblemAdapter()
{
   throw DudleyAdapterException("Error - Illegal to generate default TransportProblemAdapter.");
}


TransportProblemAdapter::TransportProblemAdapter(Paso_TransportProblem* transport_problem,
                                                 const bool useBackwardEuler,
                                                 const int block_size,
                                                 const escript::FunctionSpace& functionspace):
AbstractTransportProblem(useBackwardEuler, block_size, functionspace)
{
    m_transport_problem.reset(transport_problem,null_deleter());
}

TransportProblemAdapter::~TransportProblemAdapter()
{ 
    if (m_transport_problem.unique()) {
        Paso_TransportProblem* transp=m_transport_problem.get();
        Paso_TransportProblem_free(transp);
    }
}

Paso_TransportProblem* TransportProblemAdapter::getPaso_TransportProblem() const 
{
   return m_transport_problem.get();
}


void TransportProblemAdapter::setToSolution(escript::Data& out, escript::Data& u0, escript::Data& source, const double dt, boost::python::object& options) const
{
    Paso_TransportProblem* transp=getPaso_TransportProblem();
    Paso_Options paso_options;
    SystemMatrixAdapter::escriptToPasoOptions(&paso_options,options);
    options.attr("resetDiagnostics")();
    if ( out.getDataPointSize()  != getBlockSize()) {
     throw DudleyAdapterException("solve : block size of solution does not match block size of transport problems.");
    } else if ( source.getDataPointSize() != getBlockSize()) {
     throw DudleyAdapterException("solve : block size of source term does not match block size of transport problems.");
    } else if ( out.getFunctionSpace()  != getFunctionSpace()) {
     throw DudleyAdapterException("solve : function spaces of solution and of transport problem don't match.");
    } else if (source.getFunctionSpace() != getFunctionSpace()) {
     throw DudleyAdapterException("solve : function spaces of source term and of transport problem don't match.");
    } else if (dt<=0.) {
     throw DudleyAdapterException("solve : time increment dt needs to be positive.");
    }
    out.expand();
    source.expand();
    out.requireWrite();
    source.requireWrite();
    double* out_dp=out.getSampleDataRW(0);
    double* u0_dp=u0.getSampleDataRW(0);
    double* source_dp=source.getSampleDataRW(0);
    Paso_TransportProblem_solve(transp,out_dp,dt,u0_dp,source_dp,&paso_options);
    SystemMatrixAdapter::pasoToEscriptOptions(&paso_options,options);
    checkPasoError();
}

void TransportProblemAdapter::resetTransport() const
{
   Paso_TransportProblem* transp = getPaso_TransportProblem();
   Paso_TransportProblem_reset(transp);
   checkPasoError();
}
void TransportProblemAdapter::copyConstraint(escript::Data& source, escript::Data& q, escript::Data& r, const double factor) const
{
    if ( q.getDataPointSize()  != getBlockSize()) {
     throw DudleyAdapterException("copyConstraint : block size does not match the number of components of constraint mask.");
    } else if ( q.getFunctionSpace()  != getFunctionSpace()) {
     throw DudleyAdapterException("copyConstraint : function spaces of transport problem and constraint mask don't match.");
    } else if ( r.getDataPointSize()  != getBlockSize()) {
     throw DudleyAdapterException("copyConstraint : block size does not match the number of components of constraint values.");
    } else if ( r.getFunctionSpace()  != getFunctionSpace()) {
     throw DudleyAdapterException("copyConstraint : function spaces of transport problem and constraint values don't match.");
    } else if ( source.getDataPointSize()  != getBlockSize()) {
     throw DudleyAdapterException("copyConstraint : block size does not match the number of components of source.");
    } else if ( source.getFunctionSpace()  != getFunctionSpace()) {
     throw DudleyAdapterException("copyConstraint : function spaces of transport problem and source don't match.");
    }
    Paso_TransportProblem* transp=getPaso_TransportProblem();

    /* r2=r where q>0, 0 elsewhere */
    escript::Data r2(0.,q.getDataPointShape(),q.getFunctionSpace());
    r2.copyWithMask(r,q);

    /* source-=transp->mass_matrix*r2 */
    r2.expand();
    source.expand();
    q.expand();
    r2.requireWrite();
    source.requireWrite();
    q.requireWrite();
    double* r2_dp=r2.getSampleDataRW(0);
    double* source_dp=source.getSampleDataRW(0);
    double* q_dp=q.getSampleDataRW(0);

    if (false) {
       cout << "v1\n";
       Paso_SystemMatrix_MatrixVector(-1., transp->mass_matrix, r2_dp, 1., source_dp);
       checkPasoError();

       /* insert 0 rows into transport matrix */
       Paso_SystemMatrix_nullifyRows(transp->transport_matrix,q_dp, 0.);
       checkPasoError();

       /* insert 0 rows amd 1 in main diagonal into mass matrix */
       Paso_SystemMatrix_nullifyRowsAndCols(transp->mass_matrix,q_dp,q_dp,1.);
       checkPasoError();

       source.copyWithMask(escript::Data(0.,q.getDataPointShape(),q.getFunctionSpace()),q);
   } else {
       Paso_TransportProblem_setUpConstraint(transp, q_dp, factor);
       checkPasoError();
       Paso_TransportProblem_insertConstraint(transp,r2_dp, source_dp);
       checkPasoError();
   }
}

double TransportProblemAdapter::getSafeTimeStepSize() const
{
    Paso_TransportProblem* transp=getPaso_TransportProblem();
    double dt=Paso_TransportProblem_getSafeTimeStepSize(transp);
    checkPasoError();
    return dt;
}

double TransportProblemAdapter::getUnlimitedTimeStepSize() const
{
    return LARGE_POSITIVE_FLOAT;
}



}  // end of namespace
