
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

TransportProblemAdapter::TransportProblemAdapter()
{
   throw PasoException("Error - Illegal to generate default TransportProblemAdapter.");
}

TransportProblemAdapter::TransportProblemAdapter(TransportProblem_ptr tp,
                int block_size, const escript::FunctionSpace& functionspace) :
    AbstractTransportProblem(block_size, functionspace),
    m_transport_problem(tp)
{
}

TransportProblem_ptr TransportProblemAdapter::getPaso_TransportProblem() const 
{
   return m_transport_problem;
}

void TransportProblemAdapter::setToSolution(escript::Data& out,
        escript::Data& u0, escript::Data& source, double dt,
        boost::python::object& options) const
{
    Options paso_options;
    SystemMatrixAdapter::escriptToPasoOptions(&paso_options, options);
    options.attr("resetDiagnostics")();
    if ( out.getDataPointSize() != getBlockSize()) {
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
    double* out_dp = out.getSampleDataRW(0);
    double* u0_dp = u0.getSampleDataRW(0);
    double* source_dp = source.getSampleDataRW(0);
    SystemMatrixAdapter::pasoToEscriptOptions(&paso_options, options);
    m_transport_problem->solve(out_dp, dt, u0_dp, source_dp, &paso_options);

    checkPasoError();
}

void TransportProblemAdapter::resetTransport() const
{
    m_transport_problem->reset();
    checkPasoError();
}

void TransportProblemAdapter::copyConstraint(escript::Data& source,
                                             escript::Data& q,
                                             escript::Data& r) const
{
    if (q.getDataPointSize() != getBlockSize()) {
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

    if (false) {
        // r2=r where q>0, 0 elsewhere
        escript::Data r2(0., q.getDataPointShape(), q.getFunctionSpace());
        r2.copyWithMask(r, q);

        // source -= tp->mass_matrix*r2
        r2.expand();
        source.expand();
        q.expand();
        r2.requireWrite();
        source.requireWrite();
        q.requireWrite();
        double* r2_dp = r2.getSampleDataRW(0);
        double* source_dp = source.getSampleDataRW(0);
        double* q_dp = q.getSampleDataRW(0);
    
        SystemMatrix_MatrixVector(-1., m_transport_problem->mass_matrix,
                                  r2_dp, 1., source_dp);
        checkPasoError();

        // insert 0 rows into transport matrix
        m_transport_problem->transport_matrix->nullifyRows(q_dp, 0.);
        checkPasoError();

        // insert 0 rows and 1 in main diagonal into mass matrix
        m_transport_problem->mass_matrix->nullifyRowsAndCols(q_dp, q_dp, 1.);
        checkPasoError();

        source.copyWithMask(escript::Data(0.,q.getDataPointShape(),q.getFunctionSpace()),q);
    } else {
        r.expand();
        source.expand();
        q.expand();
        r.requireWrite();
        source.requireWrite();
        q.requireWrite();
        double* r_dp = r.getSampleDataRW(0);
        double* source_dp = source.getSampleDataRW(0);
        double* q_dp = q.getSampleDataRW(0);
        m_transport_problem->setUpConstraint(q_dp);
        checkPasoError();
        m_transport_problem->insertConstraint(r_dp, source_dp);
        checkPasoError();
    }
}

double TransportProblemAdapter::getSafeTimeStepSize() const
{
    const double dt = m_transport_problem->getSafeTimeStepSize();
    checkPasoError();
    return dt;
}

double TransportProblemAdapter::getUnlimitedTimeStepSize() const
{
    return LARGE_POSITIVE_FLOAT;
}

int TransportProblemAdapter::getTransportTypeId(int solver, int preconditioner,
                                                int package, bool symmetry,
                                                const esysUtils::JMPI& mpiInfo)
{
    return TransportProblem::getTypeId(
            SystemMatrixAdapter::mapOptionToPaso(solver),
            SystemMatrixAdapter::mapOptionToPaso(preconditioner),
            SystemMatrixAdapter::mapOptionToPaso(package), symmetry, mpiInfo);
}


}  // end of namespace

