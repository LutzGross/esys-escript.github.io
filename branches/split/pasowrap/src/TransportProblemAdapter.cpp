
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

#include "TransportProblemAdapter.h" 
#include "SystemMatrixAdapter.h" 

using namespace std;

namespace paso {

struct null_deleter
{
  void operator()(void const *ptr) const
  {
  }
};

PASOWRAP_DLL_API
TransportProblemAdapter::TransportProblemAdapter()
{
   throw PasoException("Error - Illegal to generate default TransportProblemAdapter.");
}


PASOWRAP_DLL_API
TransportProblemAdapter::TransportProblemAdapter(Paso_TransportProblem* transport_problem,
                                                 const int block_size,
                                                 const escript::FunctionSpace& functionspace):
AbstractTransportProblem(block_size, functionspace)
{
    m_transport_problem.reset(transport_problem,null_deleter());
}

PASOWRAP_DLL_API
TransportProblemAdapter::~TransportProblemAdapter()
{ 
    if (m_transport_problem.unique()) {
        Paso_TransportProblem* transp=m_transport_problem.get();
        Paso_TransportProblem_free(transp);
    }
}

PASOWRAP_DLL_API
Paso_TransportProblem* TransportProblemAdapter::getPaso_TransportProblem() const 
{
   return m_transport_problem.get();
}

PASOWRAP_DLL_API
void TransportProblemAdapter::setToSolution(escript::Data& out, escript::Data& u0, escript::Data& source, const double dt, boost::python::object& options) const
{
    Paso_TransportProblem* transp=getPaso_TransportProblem();
    Paso_Options paso_options;
    SystemMatrixAdapter::escriptToPasoOptions(&paso_options,options);
    options.attr("resetDiagnostics")();
    if ( out.getDataPointSize()  != getBlockSize()) {
     throw PasoException("solve : block size of solution does not match block size of transport problems.");
    } else if ( source.getDataPointSize() != getBlockSize()) {
     throw PasoException("solve : block size of source term does not match block size of transport problems.");
    } else if ( out.getFunctionSpace()  != getFunctionSpace()) {
     throw PasoException("solve : function spaces of solution and of transport problem don't match.");
    } else if (source.getFunctionSpace() != getFunctionSpace()) {
     throw PasoException("solve : function spaces of source term and of transport problem don't match.");
    } else if (dt<=0.) {
     throw PasoException("solve : time increment dt needs to be positive.");
    }
    out.expand();
    source.expand();
    u0.expand();
    out.requireWrite();
    source.requireWrite();
    double* out_dp=out.getSampleDataRW(0);
    double* u0_dp=u0.getSampleDataRW(0);
    double* source_dp=source.getSampleDataRW(0);
    SystemMatrixAdapter::pasoToEscriptOptions(&paso_options,options);
    Paso_TransportProblem_solve(transp,out_dp,dt,u0_dp,source_dp,&paso_options);

    checkPasoError();
}

PASOWRAP_DLL_API
void TransportProblemAdapter::resetTransport() const
{
   Paso_TransportProblem* transp = getPaso_TransportProblem();
   Paso_TransportProblem_reset(transp);
   checkPasoError();
}

PASOWRAP_DLL_API
void TransportProblemAdapter::copyConstraint(escript::Data& source, escript::Data& q, escript::Data& r) const
{
    if ( q.getDataPointSize()  != getBlockSize()) {
     throw PasoException("copyConstraint : block size does not match the number of components of constraint mask.");
    } else if ( q.getFunctionSpace()  != getFunctionSpace()) {
     throw PasoException("copyConstraint : function spaces of transport problem and constraint mask don't match.");
    } else if ( r.getDataPointSize()  != getBlockSize()) {
     throw PasoException("copyConstraint : block size does not match the number of components of constraint values.");
    } else if ( r.getFunctionSpace()  != getFunctionSpace()) {
     throw PasoException("copyConstraint : function spaces of transport problem and constraint values don't match.");
    } else if ( source.getDataPointSize()  != getBlockSize()) {
     throw PasoException("copyConstraint : block size does not match the number of components of source.");
    } else if ( source.getFunctionSpace()  != getFunctionSpace()) {
     throw PasoException("copyConstraint : function spaces of transport problem and source don't match.");
    }
    Paso_TransportProblem* transp=getPaso_TransportProblem();


    if (false) {
      
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
       r.expand();
       source.expand();
       q.expand();
       r.requireWrite();
       source.requireWrite();
       q.requireWrite();
       double* r_dp=r.getSampleDataRW(0);
       double* source_dp=source.getSampleDataRW(0);
       double* q_dp=q.getSampleDataRW(0);
       Paso_TransportProblem_setUpConstraint(transp, q_dp);
       checkPasoError();
       Paso_TransportProblem_insertConstraint(transp,r_dp, source_dp);
       checkPasoError();
   }
}

PASOWRAP_DLL_API
double TransportProblemAdapter::getSafeTimeStepSize() const
{
    Paso_TransportProblem* transp=getPaso_TransportProblem();
    double dt=Paso_TransportProblem_getSafeTimeStepSize(transp);
    checkPasoError();
    return dt;
}

PASOWRAP_DLL_API
double TransportProblemAdapter::getUnlimitedTimeStepSize() const
{
    return LARGE_POSITIVE_FLOAT;
}

PASOWRAP_DLL_API
int TransportProblemAdapter::getTransportTypeId(const int solver,
        const int preconditioner, const int package, const bool symmetry,
        const esysUtils::JMPI& mpiInfo)
{
    int out=Paso_TransportProblem_getTypeId(
            SystemMatrixAdapter::mapOptionToPaso(solver),
            SystemMatrixAdapter::mapOptionToPaso(preconditioner),
            SystemMatrixAdapter::mapOptionToPaso(package),
            symmetry?1:0, mpiInfo);
    checkPasoError();
    return out;
}


}  // end of namespace

