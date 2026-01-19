
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


/****************************************************************************/

/*   Paso: Options */

/****************************************************************************/

/*   Copyrights by ACcESS Australia 2003,2004,2005 */
/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_OPTIONS_H__
#define __PASO_OPTIONS_H__

#include "Paso.h"

#include <boost/python/object.hpp>

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
#define PASO_MKL 15
#define PASO_UMFPACK 16
#define PASO_NO_REORDERING 17
#define PASO_MINIMUM_FILL_IN 18
#define PASO_NESTED_DISSECTION 19
#define PASO_ITERATIVE 20
#define PASO_PASO 21
#define PASO_MUMPS 22
#define PASO_REC_ILU  23
#define PASO_TRILINOS  24
#define PASO_NONLINEAR_GMRES  25
#define PASO_TFQMR  26
#define PASO_MINRES  27
#define PASO_GAUSS_SEIDEL 28
#define PASO_GS PASO_GAUSS_SEIDEL
#define PASO_RILU 29
#define PASO_DEFAULT_REORDERING 30
#define PASO_NO_PRECONDITIONER 36
#define PASO_CLASSIC_INTERPOLATION_WITH_FF_COUPLING 50
#define PASO_CLASSIC_INTERPOLATION 51
#define PASO_DIRECT_INTERPOLATION 52
#define PASO_LINEAR_CRANK_NICOLSON 66
#define PASO_CRANK_NICOLSON 67
#define PASO_BACKWARD_EULER 68

#define PASO_SMOOTHER 99999999

namespace paso {

struct PASO_DLL_API Options
{
    Options() { setDefaults(); }

    /// constructor that fills values from an escript SolverBuddy instance
    Options(const boost::python::object& options);

    /// sets the default values for solver options
    void setDefaults();

    /// prints current option values
    void show() const;

    /// prints diagnostic data
    void showDiagnostics() const;

    /// updates SolverBuddy diagnostics from this
    void updateEscriptDiagnostics(boost::python::object& options) const;

    /// returns the corresponding paso option code for an escript option code
    static int mapEscriptOption(int escriptOption);

    static const char* name(int key);

    static int getPackage(int solver, int package, bool symmetry,
                          const escript::JMPI& mpi_info);

    /// returns the solver to be used with given combination
    static int getSolver(int solver, int package, bool symmetry,
                         const escript::JMPI& mpi_info);

    int method;
    int package;
    bool symmetric;
    bool hermitian;
    double tolerance;
    double absolute_tolerance;
    double inner_tolerance;
    bool adapt_inner_tolerance;
    bool verbose;
    bool reordering;
    int preconditioner;
    dim_t iter_max;
    dim_t inner_iter_max;
    double drop_tolerance;
    double drop_storage;
    index_t truncation;
    index_t restart;
    int sweeps;
    bool accept_failed_convergence;
    double relaxation_factor;
    bool use_local_preconditioner;
    dim_t refinements;
    int ode_solver;

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

