
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#ifndef __PASO_TRANSPORT_H__
#define __PASO_TRANSPORT_H__

#define DT_FACTOR_MAX 100000.

#include "Paso.h"
#include "Options.h"
#include "SystemMatrix.h"

#include <escript/AbstractTransportProblem.h>

namespace paso {

class PASO_DLL_API TransportProblem;
typedef boost::shared_ptr<TransportProblem> TransportProblem_ptr;
typedef boost::shared_ptr<const TransportProblem> const_TransportProblem_ptr;

class PASO_DLL_API TransportProblem : public escript::AbstractTransportProblem,
                         public boost::enable_shared_from_this<TransportProblem>
{
public:
    /// Default constructor - throws exception
    TransportProblem();

    TransportProblem(SystemMatrixPattern_ptr pattern, int blocksize,
                     const escript::FunctionSpace& functionspace);

    ~TransportProblem();

    virtual void resetTransport(bool preserveSolverData) const;

    void solve(double* u, double dt, double* u0, double* q, Options* options);

    virtual double getSafeTimeStepSize() const;

    virtual double getUnlimitedTimeStepSize() const;

    void insertConstraint(const double* r,  double* source) const;

    void setUpConstraint(const double* q);

    inline dim_t getBlockSize() const
    {
        return transport_matrix->row_block_size;
    }

    inline SystemMatrix_ptr<double> borrowTransportMatrix() const
    {
        return transport_matrix;
    }

    inline SystemMatrix_ptr<double> borrowMassMatrix() const
    {
        return mass_matrix;
    }

    inline double* borrowLumpedMassMatrix() const
    {
        return lumped_mass_matrix;
    }

    inline dim_t getTotalNumRows() const
    {
        return transport_matrix->getTotalNumRows();
    }

    inline Connector_ptr borrowConnector() const
    {
        return transport_matrix->pattern->col_connector;
    }

    inline index_t* borrowMainDiagonalPointer() const
    {
       return mass_matrix->mainBlock->borrowMainDiagonalPointer();
    }

    inline static int getTypeId(int solver, int preconditioner,
                                int package, bool symmetry,
                                const escript::JMPI& mpi_info)
    {
        return MATRIX_FORMAT_DEFAULT + MATRIX_FORMAT_BLK1;
    }

    SystemMatrix_ptr<double> transport_matrix;
    SystemMatrix_ptr<double> mass_matrix;
    SystemMatrix_ptr<double> iteration_matrix;

    mutable bool valid_matrices;
    /// safe time step size for reactive part
    mutable double dt_max_R;
    /// safe time step size for transport part
    mutable double dt_max_T;
    mutable double* constraint_mask;

    double* main_diagonal_low_order_transport_matrix;
    /// 'relevant' lumped mass matrix is assumed to be positive.
    /// Values with corresponding constraint_mask > 0 value are set to -1
    /// to indicate the value infinity
    double* lumped_mass_matrix;
    double* reactive_matrix;
    double* main_diagonal_mass_matrix;

    escript::JMPI mpi_info;

private:
    virtual void setToSolution(escript::Data& out, escript::Data& u0,
                               escript::Data& source, double dt,
                               boost::python::object& options);

    virtual void copyConstraint(escript::Data& source, escript::Data& q,
                                escript::Data& r);
};

} // namespace paso

#endif // __PASO_TRANSPORT_H__

