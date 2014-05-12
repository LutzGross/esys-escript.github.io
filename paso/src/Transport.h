
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


#ifndef __PASO_TRANSPORT_H__
#define __PASO_TRANSPORT_H__

#define DT_FACTOR_MAX 100000.

#include "SystemMatrix.h"
#include "Options.h"

namespace paso {

struct TransportProblem;
typedef boost::shared_ptr<TransportProblem> TransportProblem_ptr;
typedef boost::shared_ptr<const TransportProblem> const_TransportProblem_ptr;

PASO_DLL_API
struct TransportProblem : boost::enable_shared_from_this<TransportProblem>
{
    TransportProblem(SystemMatrixPattern_ptr pattern, int block_size);
    ~TransportProblem();

    void reset();

    void solve(double* u, double dt, double* u0, double* q, Options* options);

    double getSafeTimeStepSize();

    void insertConstraint(const double* r,  double* source);

    void setUpConstraint(const double* q);

    inline dim_t getBlockSize() const
    {
        return transport_matrix->row_block_size;
    }

    inline SystemMatrix_ptr borrowTransportMatrix() const
    {
        return transport_matrix;
    }

    inline SystemMatrix_ptr borrowMassMatrix() const
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

    inline static index_t getTypeId(index_t solver, index_t preconditioner,
                                    index_t package, bool symmetry,
                                    const esysUtils::JMPI& mpi_info)
    {
        return MATRIX_FORMAT_DEFAULT + MATRIX_FORMAT_BLK1;
    }

    SystemMatrix_ptr transport_matrix;
    SystemMatrix_ptr mass_matrix;
    SystemMatrix_ptr iteration_matrix;

    bool valid_matrices;
    /// safe time step size for reactive part
    double dt_max_R;
    /// safe time step size for transport part
    double dt_max_T;
    double* constraint_mask;

    double* main_diagonal_low_order_transport_matrix;
    /// 'relevant' lumped mass matrix is assumed to be positive.
    /// Values with corresponding constraint_mask > 0 value are set to -1
    /// to indicate the value infinity
    double* lumped_mass_matrix;
    double* reactive_matrix;
    double* main_diagonal_mass_matrix;

    esysUtils::JMPI mpi_info;
};

} // namespace paso

#endif // __PASO_TRANSPORT_H__

