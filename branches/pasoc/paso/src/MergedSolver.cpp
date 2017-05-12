
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************/

/* Paso: Merged solver for AMG                                              */

/****************************************************************************/

/* Author: lgao@uq.edu.au, l.gross@uq.edu.au                                */

/****************************************************************************/

#include "MergedSolver.h"
#include "Options.h"
#include "PasoUtil.h"
#include "Preconditioner.h"
#include "MKL.h"
#include "UMFPACK.h"

namespace paso {

MergedSolver::MergedSolver(const_SystemMatrix_ptr M, const Options* options)
{
    const index_t rank = M->mpi_info->rank;
    const index_t size = M->mpi_info->size;
    const dim_t global_n = M->getGlobalNumRows();
    const dim_t n_block = M->mainBlock->row_block_size;
    const std::vector<index_t> dist(M->pattern->input_distribution->first_component);

    SparseMatrix_ptr M_temp(M->mergeSystemMatrix());

    mpi_info = M->mpi_info;
    reordering = options->reordering;
    refinements = options->coarse_matrix_refinements;
    //verbose = options->verbose;
    verbose = false;
    sweeps = options->pre_sweeps+options->post_sweeps;

    // First, gather x and b into rank 0
    b = new double[global_n*n_block];
    x = new double[global_n*n_block];
    counts = new int[size];
    offset = new int[size];

#pragma omp parallel for
    for (dim_t i=0; i<size; i++) {
        const dim_t p = dist[i];
        counts[i] = (dist[i+1] - p)*n_block;
        offset[i] = p*n_block;
    }

    if (rank == 0) {
#ifdef ESYS_HAVE_MKL
        A = M_temp->unroll(MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_OFFSET1);
        A->solver_package = PASO_MKL;
#elif defined ESYS_HAVE_UMFPACK
        A = M_temp->unroll(MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_CSC);
        A->solver_package = PASO_UMFPACK;
#else
        A->solver_p = Preconditioner_LocalSmoother_alloc(A,
                            (options->smoother == PASO_JACOBI), verbose);
        A->solver_package = PASO_SMOOTHER;
#endif
    }
}

MergedSolver::~MergedSolver()
{
    delete[] x;
    delete[] b;
    delete[] counts;
    delete[] offset;
}

void MergedSolver::solve(double* local_x, const double* local_b)
{
    const index_t rank = mpi_info->rank;
    const dim_t count = counts[rank];

#ifdef ESYS_MPI
    MPI_Gatherv(const_cast<double*>(local_b), count, MPI_DOUBLE, b, counts,
                offset, MPI_DOUBLE, 0, mpi_info->comm);
#else
#pragma omp parallel for
    for (dim_t i=0; i<count; i++) {
        b[i] = local_b[i];
        x[i] = local_x[i];
    }
#endif
    if (rank == 0) {
        switch (A->solver_package) {
            case PASO_MKL:
                MKL_solve(A, x, b, reordering, refinements, verbose);
                break;
            case PASO_UMFPACK:
                UMFPACK_solve(A, x, b, refinements, verbose);
                break;
            case PASO_SMOOTHER:
                Preconditioner_LocalSmoother_solve(A, reinterpret_cast<Preconditioner_LocalSmoother*>(A->solver_p), x, b, sweeps, false);
                break;
        }
    }
#ifdef ESYS_MPI
    // now we need to distribute the solution to all ranks
    MPI_Scatterv(x, counts, offset, MPI_DOUBLE, local_x, count, MPI_DOUBLE,
                 0, mpi_info->comm);
#else
#pragma omp parallel for
    for (dim_t i=0; i<count; i++)
        local_x[i] = x[i];
#endif
}

} // namespace paso

