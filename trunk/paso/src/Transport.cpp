
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


/****************************************************************************/

/* Paso: TransportProblem (see TransportSolver::solve)                      */

/****************************************************************************/

/* Author: l.gross@uq.edu.au                                                */

/****************************************************************************/

#include "Transport.h"
#include "PasoUtil.h"
#include "Preconditioner.h"

namespace paso {

TransportProblem::TransportProblem(SystemMatrixPattern_ptr pattern,
                                   int block_size) :
    valid_matrices(false),
    dt_max_R(LARGE_POSITIVE_FLOAT),
    dt_max_T(LARGE_POSITIVE_FLOAT),
    constraint_mask(NULL),
    main_diagonal_low_order_transport_matrix(NULL),
    lumped_mass_matrix(NULL),
    reactive_matrix(NULL),
    main_diagonal_mass_matrix(NULL)
{
    // at the moment only block size 1 is supported
    SystemMatrixType matrix_type = MATRIX_FORMAT_DEFAULT+MATRIX_FORMAT_BLK1;

    transport_matrix.reset(new SystemMatrix(matrix_type, pattern, block_size,
                                            block_size, false));
    mass_matrix.reset(new SystemMatrix(matrix_type, pattern, block_size,
                                       block_size, false));

    mpi_info = Esys_MPIInfo_getReference(pattern->mpi_info);

    if (Esys_noError()) {
        const dim_t n = transport_matrix->getTotalNumRows();
        constraint_mask = new double[n];
        lumped_mass_matrix = new double[n];
        reactive_matrix = new double[n];
        main_diagonal_mass_matrix = new double[n];
        main_diagonal_low_order_transport_matrix = new double[n];

#pragma omp parallel for
        for (dim_t i = 0; i < n; ++i) {
            lumped_mass_matrix[i] = 0.;
            main_diagonal_low_order_transport_matrix[i] = 0.;
            constraint_mask[i] = 0.;
        }
    }
}

TransportProblem::~TransportProblem()
{
    Esys_MPIInfo_free(mpi_info);
    delete[] constraint_mask;
    delete[] reactive_matrix;
    delete[] main_diagonal_mass_matrix;
    delete[] lumped_mass_matrix;
    delete[] main_diagonal_low_order_transport_matrix;
}

void TransportProblem::reset()
{
    const dim_t n = transport_matrix->getTotalNumRows();
    transport_matrix->setValues(0.);
    mass_matrix->setValues(0.);
    Paso_solve_free(iteration_matrix.get());
    util::zeroes(n, constraint_mask);
    valid_matrices = false;
}


void TransportProblem::setUpConstraint(const double* q)
{
    if (valid_matrices) {
        Esys_setError(VALUE_ERROR, "TransportProblem::setUpConstraint: "
                            "Cannot insert a constraint into a valid system.");
        return;
    }

    const dim_t n = transport_matrix->getTotalNumRows();
#pragma omp parallel for
    for (dim_t i=0; i<n; ++i) {
        if (q[i] > 0) {
            constraint_mask[i]=1;
        } else {
            constraint_mask[i]=0;
        }
    }
}

void TransportProblem::insertConstraint(const double* r,  double* source)
{
    const dim_t n = transport_matrix->getTotalNumRows();

#pragma omp parallel for
    for (dim_t i=0; i<n; ++i) {
        if (constraint_mask[i] > 0)
            source[i] = r[i];
    }
}

} // namespace paso

