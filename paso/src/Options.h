
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

/*   Paso: Options */

/****************************************************************************/

/*   Copyrights by ACcESS Australia 2003,2004,2005 */
/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_OPTIONS_H__
#define __PASO_OPTIONS_H__

#include "esysUtils/types.h"

// valid solver options
#define PASO_DEFAULT 0
#define PASO_DIRECT 1
#define PASO_CHOLEVSKY 2
#define PASO_PCG 3
#define PASO_CR 4
#define PASO_CGS 5
#define PASO_BICGSTAB 6
#define PASO_ILU0 8
#define PASO_ILUT 9
#define PASO_JACOBI 10
#define PASO_GMRES 11
#define PASO_PRES20 12
#define PASO_LUMPING 13
#define PASO_MKL 15
#define PASO_UMFPACK 16
#define PASO_NO_REORDERING 17
#define PASO_MINIMUM_FILL_IN 18
#define PASO_NESTED_DISSECTION 19
#define PASO_ITERATIVE 20
#define PASO_PASO 21
#define PASO_AMG 22
#define PASO_REC_ILU  23
#define PASO_TRILINOS  24
#define PASO_NONLINEAR_GMRES  25
#define PASO_TFQMR  26
#define PASO_MINRES  27
#define PASO_GAUSS_SEIDEL 28
#define PASO_GS PASO_GAUSS_SEIDEL
#define PASO_RILU 29
#define PASO_DEFAULT_REORDERING 30
#define PASO_SUPER_LU 31
#define PASO_PASTIX 32
#define PASO_YAIR_SHAPIRA_COARSENING 33
#define PASO_RUGE_STUEBEN_COARSENING 34
#define PASO_AGGREGATION_COARSENING 35
#define PASO_NO_PRECONDITIONER 36
#define PASO_MIN_COARSE_MATRIX_SIZE 37
#define PASO_AMLI 38
#define PASO_STANDARD_COARSENING 39
#define PASO_CLASSIC_INTERPOLATION_WITH_FF_COUPLING 50
#define PASO_CLASSIC_INTERPOLATION 51
#define PASO_DIRECT_INTERPOLATION 52
#define PASO_BOOMERAMG 60
#define PASO_CIJP_FIXED_RANDOM_COARSENING 61
#define PASO_CIJP_COARSENING 62
#define PASO_FALGOUT_COARSENING 63
#define PASO_PMIS_COARSENING 64
#define PASO_HMIS_COARSENING 65
#define PASO_LINEAR_CRANK_NICOLSON 66
#define PASO_CRANK_NICOLSON 67
#define PASO_BACKWARD_EULER 68

#define PASO_SMOOTHER 99999999

struct Esys_MPIInfo;

namespace paso {

PASO_DLL_API
struct Options
{
    Options() { setDefaults(); }

    /// sets the default values for solver options
    void setDefaults();

    /// prints current option values
    void show() const;

    /// prints diagnostic data
    void showDiagnostics() const;

    static const char* name(index_t key);

    static index_t getPackage(index_t solver, index_t package, bool symmetry,
                              Esys_MPIInfo* mpi_info);

    /// returns the solver to be used with given combination
    static index_t getSolver(index_t solver, index_t package, bool symmetry,
                             Esys_MPIInfo* mpi_info);

    index_t method;
    index_t package;
    bool symmetric;
    double tolerance;
    double absolute_tolerance;
    double inner_tolerance;
    bool adapt_inner_tolerance;
    bool verbose;
    bool reordering;
    index_t preconditioner;
    dim_t iter_max;
    dim_t inner_iter_max;
    double drop_tolerance;
    double drop_storage;
    index_t truncation;
    index_t restart;
    dim_t sweeps;
    dim_t pre_sweeps;
    dim_t post_sweeps;
    dim_t cycle_type;
    dim_t level_max;
    dim_t min_coarse_matrix_size;
    dim_t smoother;
    double coarsening_threshold;
    bool accept_failed_convergence;
    index_t coarsening_method;
    double relaxation_factor;
    bool use_local_preconditioner;
    double min_coarse_sparsity;
    dim_t refinements;
    dim_t coarse_matrix_refinements;
    double diagonal_dominance_threshold;
    bool usePanel;
    index_t interpolation_method;
    index_t ode_solver;
    
    // diagnostic values
    dim_t num_iter;
    dim_t num_level;
    dim_t num_inner_iter;
    double time;
    double set_up_time;
    double coarsening_selection_time;
    double coarsening_matrix_time;
    double net_time;
    double residual_norm;
    bool converged;
    double preconditioner_size; // in Mbytes
    bool time_step_backtracking_used;
    double coarse_level_sparsity;
    dim_t num_coarse_unknowns;
};

} // namespace paso

#endif // __PASO_OPTIONS_H__

