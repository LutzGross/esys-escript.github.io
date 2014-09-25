
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

/*   Paso: solver options */

/****************************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004 */
/*   author: l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "Options.h"

namespace paso {

void Options::setDefaults()
{
    verbose = false;
    method = PASO_DEFAULT;
    package = PASO_DEFAULT;
    symmetric = false;
    reordering = PASO_NO_REORDERING;
    tolerance = 1.e-8;
    absolute_tolerance = 0.;
    inner_tolerance = 0.9;
    adapt_inner_tolerance = true;
    preconditioner = PASO_JACOBI;
    iter_max = 10000;
    inner_iter_max = 10;
    drop_tolerance = 0.01;
    drop_storage = 2.;
    restart = -1;
    truncation = 20;
    sweeps = 2;
    pre_sweeps = 2;
    post_sweeps = 2;
    coarsening_threshold = 0.25;
    min_coarse_matrix_size = 500;
    level_max = 100;
    accept_failed_convergence = false;
    coarsening_method = PASO_DEFAULT;
    relaxation_factor = 0.95;
    smoother = PASO_GS;
    use_local_preconditioner = false;
    min_coarse_sparsity = 0.05;
    refinements = 2;
    coarse_matrix_refinements = 0;
    diagonal_dominance_threshold = 0.5;
    cycle_type = 1;
    usePanel = true;
    interpolation_method = PASO_DIRECT_INTERPOLATION;
    ode_solver = PASO_LINEAR_CRANK_NICOLSON;

    // diagnostic values
    num_iter = -1;
    num_level = -1;
    num_inner_iter = -1;
    time = -1.;
    set_up_time = -1.;
    coarsening_selection_time = -1.;
    coarsening_matrix_time = -1;
    net_time = -1.;
    residual_norm = -1.;
    converged = false;
    preconditioner_size = -1.;
    time_step_backtracking_used = false;
    coarse_level_sparsity = -1.;
    num_coarse_unknowns = -1;
}

void Options::showDiagnostics() const
{
    std::cout << "Paso diagnostics:" << std::endl
        << "\tnum_iter = " << num_iter << std::endl
        << "\tnum_level = " << num_level << std::endl
        << "\tnum_inner_iter = " << num_inner_iter << std::endl
        << "\ttime = " << time << std::endl
        << "\tset_up_time = " << set_up_time << std::endl
        << "\tcoarsening_selection_time = " << coarsening_selection_time << std::endl
        << "\tcoarsening_matrix_time = " << coarsening_matrix_time << std::endl
        << "\tnet_time = " << net_time << std::endl
        << "\tresidual_norm = " << residual_norm << std::endl
        << "\tconverged = " << converged << std::endl
        << "\tpreconditioner_size = " << preconditioner_size << " MBytes" << std::endl
        << "\ttime_step_backtracking_used = " << time_step_backtracking_used << std::endl;
}

void Options::show() const
{
    std::cout << "Paso options settings:" << std::endl
        << "\tverbose = " << verbose << std::endl
        << "\tmethod = " << name(method) << " (" << method << ")" << std::endl
        << "\tpackage = " << name(package) << " (" << package << ")" << std::endl
        << "\tsymmetric = " << symmetric << std::endl
        << "\treordering = " << name(reordering) << " (" << reordering << ")" << std::endl
        << "\ttolerance = " << tolerance << std::endl
        << "\tabsolute_tolerance = " << absolute_tolerance << std::endl
        << "\tinner_tolerance = " << inner_tolerance << std::endl
        << "\tadapt_inner_tolerance = " << adapt_inner_tolerance << std::endl
        << "\tpreconditioner = " << name(preconditioner) << " (" << preconditioner << ")" << std::endl
        << "\titer_max = " << iter_max << std::endl
        << "\tinner_iter_max = " << inner_iter_max << std::endl
        << "\tdrop_tolerance = " << drop_tolerance << std::endl
        << "\tdrop_storage = " << drop_storage << std::endl
        << "\trestart = " << restart << std::endl
        << "\ttruncation = " << truncation << std::endl
        << "\tsweeps = " << sweeps << std::endl
        << "\tpre_sweeps = " << pre_sweeps << std::endl
        << "\tpost_sweeps = " << post_sweeps << std::endl
        << "\tcoarsening_threshold = " << coarsening_threshold << std::endl
        << "\tlevel_max = " << level_max << std::endl
        << "\taccept_failed_convergence = " << accept_failed_convergence << std::endl
        << "\tcoarsening_method = " << name(coarsening_method) << " (" << coarsening_method << ")" << std::endl
        << "\trelaxation_factor = " << relaxation_factor << std::endl
        << "\tuse_local_preconditioner = " << use_local_preconditioner << std::endl
        << "\tmin_coarse_sparsity = " << min_coarse_sparsity << std::endl
        << "\trefinements = " << refinements << std::endl
        << "\tcoarse_matrix_refinements = " << coarse_matrix_refinements << std::endl
        << "\tcycle_type = " << cycle_type << std::endl
        << "\tode_solver = " << ode_solver << std::endl;
}

const char* Options::name(int key)
{
    switch (key) {
       case PASO_DEFAULT:
            return "DEFAULT";
       case PASO_DIRECT:
            return "DIRECT";
       case PASO_CHOLEVSKY:
            return "CHOLEVSKY";
       case PASO_PCG:
            return "PCG";
       case PASO_CR:
            return "CR";
       case PASO_CGS:
            return "CGS";
       case PASO_BICGSTAB:
            return "BICGSTAB";
       case PASO_ILU0:
            return "ILU0";
       case PASO_ILUT:
            return "ILUT";
       case PASO_JACOBI:
            return "JACOBI";
       case PASO_GMRES:
            return "GMRES";
       case PASO_PRES20:
            return "PRES20";
       case PASO_NO_REORDERING:
            return "NO_REORDERING";
       case PASO_MINIMUM_FILL_IN:
            return "MINIMUM_FILL_IN";
       case PASO_NESTED_DISSECTION:
            return "NESTED_DISSECTION";
       case PASO_MKL:
            return "MKL";
       case PASO_UMFPACK:
            return "UMFPACK";
       case PASO_BOOMERAMG:
            return "BOOMERAMG";
       case PASO_ITERATIVE:
            return "ITERATIVE";
       case PASO_PASO:
            return "PASO";
       case PASO_AMG:
            return "AMG";
       case PASO_AMLI:
            return "AMLI";
       case PASO_REC_ILU:
            return "REC_ILU";
       case PASO_TRILINOS:
            return "TRILINOS";
       case PASO_NONLINEAR_GMRES:
            return "NONLINEAR_GMRES";
       case PASO_TFQMR :
            return "TFQMR";
       case PASO_MINRES:
            return "MINRES";
       case PASO_GAUSS_SEIDEL:
            return "GAUSS_SEIDEL";
       case PASO_RILU:
            return "RILU";
       case PASO_DEFAULT_REORDERING:
            return "DEFAULT_REORDERING";
       case PASO_SUPER_LU:
            return "SUPER_LU";
       case PASO_PASTIX:
            return "PASTIX";
       case PASO_YAIR_SHAPIRA_COARSENING:
            return "YAIR_SHAPIRA_COARSENING";
       case PASO_RUGE_STUEBEN_COARSENING:
            return "RUGE_STUEBEN_COARSENING";
       case PASO_AGGREGATION_COARSENING:
            return "AGGREGATION_COARSENING";
       case PASO_STANDARD_COARSENING:
            return "STANDARD_COARSENING";
       case PASO_NO_PRECONDITIONER:
            return "NO_PRECONDITIONER";
       case PASO_CIJP_FIXED_RANDOM_COARSENING:
            return "CIJP_FIXED_RANDOM_COARSENING";
       case PASO_CIJP_COARSENING:
            return "CIJP_COARSENING";
       case PASO_FALGOUT_COARSENING:
            return "FALGOUT_COARSENING";
       case PASO_PMIS_COARSENING:
            return "PMIS_COARSENING";
       case PASO_HMIS_COARSENING:
            return "HMIS_COARSENING";
       case PASO_CRANK_NICOLSON:
            return "PASO_CRANK_NICOLSON";
       case PASO_LINEAR_CRANK_NICOLSON:
            return "PASO_CRANK_NICOLSON";
       case PASO_BACKWARD_EULER:
            return "PASO_BACKWARD_EULER";
       default:
            return "<unknown>";
    }
}

int Options::getSolver(int solver, int pack, bool symmetry,
                       const esysUtils::JMPI& mpi_info)
{
    int out = PASO_DEFAULT;
    // PASO //
    if (pack==PASO_PASO) {
        switch (solver) {
            case PASO_BICGSTAB:
                out=PASO_BICGSTAB;
                break;
            case PASO_PCG:
                out=PASO_PCG;
                break;
            case PASO_PRES20:
                out=PASO_PRES20;
                break;
            case PASO_GMRES:
                out=PASO_GMRES;
                break;
            case PASO_NONLINEAR_GMRES:
                out=PASO_NONLINEAR_GMRES;
                break;
            case PASO_TFQMR:
                out=PASO_TFQMR;
                break;
            case PASO_MINRES:
                out=PASO_MINRES;
                break;
            default:
                if (symmetry) {
                    out=PASO_PCG;
                } else {
                    out=PASO_BICGSTAB;
                }
            break;
        }
    // MKL //
    } else if (pack==PASO_MKL) {
        switch (solver) {
            case PASO_CHOLEVSKY:
                out=PASO_CHOLEVSKY;
            break;
            case PASO_DIRECT:
                out=PASO_DIRECT;
            break;
            default:
                if (symmetry) {
                    out=PASO_CHOLEVSKY;
                } else {
                    out=PASO_DIRECT;
                }
            break;
        }
    // TRILINOS //
    } else if (pack==PASO_TRILINOS) {
        switch (solver) {
            case PASO_BICGSTAB:
                out=PASO_BICGSTAB;
                break;
            case PASO_PCG:
                out=PASO_PCG;
                break;
            case PASO_PRES20:
                out=PASO_PRES20;
                break;
            case PASO_GMRES:
                out=PASO_GMRES;
                break;
            case PASO_TFQMR:
                out=PASO_TFQMR;
                break;
            case PASO_MINRES:
                out=PASO_MINRES;
                break;
            default:
                if (symmetry) {
                    out=PASO_PCG;
                } else {
                    out=PASO_BICGSTAB;
                }
            break;
        }
    // UMFPACK //
    } else if (pack==PASO_UMFPACK) {
        out=PASO_DIRECT;
    } else {
        Esys_setError(VALUE_ERROR, "Options::getSolver: Unidentified package.");
    }
    return out;
}

int Options::getPackage(int solver, int pack, bool symmetry,
                        const esysUtils::JMPI& mpi_info)
{
    int out = PASO_PASO;

    switch (pack) {
        case PASO_DEFAULT:
            if (solver == PASO_DIRECT) {
                // these packages require CSC which is not supported with MPI
                if (mpi_info->size == 1) {
#if defined MKL
                    out = PASO_MKL;
#elif defined UMFPACK
                    out = PASO_UMFPACK;
#elif defined PASTIX
                    out = PASO_PASTIX
#endif
                }
            }
            break;

        case PASO_PASO:
            break;

        case PASO_MKL:
        case PASO_UMFPACK:
        case PASO_PASTIX:
        case PASO_TRILINOS:
            out = pack;
            break;

        default:
            Esys_setError(VALUE_ERROR, "Options::getPackage: Unidentified package.");
    }
    return out;
}

} // namespace paso

