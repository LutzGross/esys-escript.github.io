
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

#include "Paso.h"
#include "Options.h"
#include "PasoException.h"

#include <escript/SolverOptions.h>

#include <boost/python/extract.hpp>
#include <iostream>
#include <sstream>

namespace bp = boost::python;

namespace paso {

Options::Options(const bp::object& options)
{
    escript::SolverBuddy sb = bp::extract<escript::SolverBuddy>(options);

    setDefaults();
    method = mapEscriptOption(sb.getSolverMethod());
    package = mapEscriptOption(sb.getPackage());
    verbose = sb.isVerbose();
    symmetric = sb.isSymmetric();
    tolerance = sb.getTolerance();
    absolute_tolerance = sb.getAbsoluteTolerance();
    inner_tolerance = sb.getInnerTolerance();
    adapt_inner_tolerance = sb.adaptInnerTolerance();
    reordering = mapEscriptOption(sb.getReordering());
    preconditioner = mapEscriptOption(sb.getPreconditioner());
    ode_solver = mapEscriptOption(sb.getODESolver());
    iter_max = sb.getIterMax();
    inner_iter_max = sb.getInnerIterMax();
    drop_tolerance = sb.getDropTolerance();
    drop_storage = sb.getDropStorage();
    truncation = sb.getTruncation();
    restart = sb._getRestartForC();
    sweeps = sb.getNumSweeps();
    pre_sweeps = sb.getNumPreSweeps();
    post_sweeps = sb.getNumPostSweeps();
    level_max = sb.getLevelMax();
    min_coarse_matrix_size = sb.getMinCoarseMatrixSize();
    coarsening_threshold = sb.getCoarseningThreshold();
    accept_failed_convergence = sb.acceptConvergenceFailure();
    coarsening_method = mapEscriptOption(sb.getCoarsening());
    smoother = mapEscriptOption(sb.getSmoother());
    relaxation_factor = sb.getRelaxationFactor();
    use_local_preconditioner = sb.useLocalPreconditioner();
    min_coarse_sparsity = sb.getMinCoarseMatrixSparsity();
    refinements = sb.getNumRefinements();
    coarse_matrix_refinements = sb.getNumCoarseMatrixRefinements();
    usePanel = sb.usePanel();
    interpolation_method = sb.getAMGInterpolation();
    diagonal_dominance_threshold = sb.getDiagonalDominanceThreshold();
}

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
                       const escript::JMPI& mpi_info)
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
        throw PasoException("Options::getSolver: Unidentified package.");
    }
    return out;
}

int Options::getPackage(int solver, int pack, bool symmetry,
                        const escript::JMPI& mpi_info)
{
    int out = PASO_PASO;

    switch (pack) {
        case PASO_DEFAULT:
            if (solver == PASO_DIRECT) {
                // these packages require CSC which is not supported with MPI
                if (mpi_info->size == 1) {
#ifdef ESYS_HAVE_MKL
                    out = PASO_MKL;
#elif defined ESYS_HAVE_UMFPACK
                    out = PASO_UMFPACK;
#endif
                } else{
#ifdef ESYS_HAVE_MKL
                    throw PasoException("MKL does not currently support MPI");
#elif defined ESYS_HAVE_UMFPACK
                    throw PasoException("UMFPACK does not currently support MPI");
#endif
                }
            }
            break;

        case PASO_PASO:
            break;

        case PASO_MKL:
        case PASO_UMFPACK:
        case PASO_TRILINOS:
            out = pack;
            break;

        default:
            throw PasoException("Options::getPackage: Unidentified package.");
    }
    return out;
}

int Options::mapEscriptOption(int escriptOption)
{
    switch (escriptOption) {
        case escript::SO_DEFAULT:
            return PASO_DEFAULT;

        case escript::SO_PACKAGE_MKL:
            return PASO_MKL;
        case escript::SO_PACKAGE_PASO:
            return PASO_PASO;
        case escript::SO_PACKAGE_TRILINOS:
            return PASO_TRILINOS;
        case escript::SO_PACKAGE_UMFPACK:
            return PASO_UMFPACK;

        case escript::SO_METHOD_BICGSTAB:
            return PASO_BICGSTAB;
        case escript::SO_METHOD_CGS:
            return PASO_CGS;
        case escript::SO_METHOD_CHOLEVSKY:
            return PASO_CHOLEVSKY;
        case escript::SO_METHOD_CR:
            return PASO_CR;
        case escript::SO_METHOD_DIRECT:
            return PASO_DIRECT;
        case escript::SO_METHOD_GMRES:
            return PASO_GMRES;
        case escript::SO_METHOD_ITERATIVE:
            return PASO_ITERATIVE;
        case escript::SO_METHOD_MINRES:
            return PASO_MINRES;
        case escript::SO_METHOD_NONLINEAR_GMRES:
            return PASO_NONLINEAR_GMRES;
        case escript::SO_METHOD_PCG:
            return PASO_PCG;
        case escript::SO_METHOD_PRES20:
            return PASO_PRES20;
        case escript::SO_METHOD_TFQMR:
            return PASO_TFQMR;

        case escript::SO_PRECONDITIONER_AMG:
            return PASO_AMG;
        case escript::SO_PRECONDITIONER_AMLI:
            return PASO_AMLI;
        case escript::SO_PRECONDITIONER_BOOMERAMG:
            return PASO_BOOMERAMG;
        case escript::SO_PRECONDITIONER_GAUSS_SEIDEL:
            return PASO_GAUSS_SEIDEL;
        case escript::SO_PRECONDITIONER_ILU0:
            return PASO_ILU0;
        case escript::SO_PRECONDITIONER_ILUT:
            return PASO_ILUT;
        case escript::SO_PRECONDITIONER_JACOBI:
            return PASO_JACOBI;
        case escript::SO_PRECONDITIONER_NONE:
            return PASO_NO_PRECONDITIONER;
        case escript::SO_PRECONDITIONER_REC_ILU:
            return PASO_REC_ILU;
        case escript::SO_PRECONDITIONER_RILU:
            return PASO_RILU;

        case escript::SO_ODESOLVER_BACKWARD_EULER:         
            return PASO_BACKWARD_EULER;
        case escript::SO_ODESOLVER_CRANK_NICOLSON:
            return PASO_CRANK_NICOLSON;
        case escript::SO_ODESOLVER_LINEAR_CRANK_NICOLSON:
            return PASO_LINEAR_CRANK_NICOLSON;

        case escript::SO_INTERPOLATION_CLASSIC:
            return PASO_CLASSIC_INTERPOLATION;
        case escript::SO_INTERPOLATION_CLASSIC_WITH_FF_COUPLING:
            return PASO_CLASSIC_INTERPOLATION_WITH_FF_COUPLING;
        case escript::SO_INTERPOLATION_DIRECT:
            return PASO_DIRECT_INTERPOLATION;

        case escript::SO_COARSENING_AGGREGATION:
            return PASO_AGGREGATION_COARSENING;
        case escript::SO_COARSENING_CIJP:
            return PASO_CIJP_COARSENING;
        case escript::SO_COARSENING_CIJP_FIXED_RANDOM:
            return PASO_CIJP_FIXED_RANDOM_COARSENING;
        case escript::SO_COARSENING_FALGOUT:
            return PASO_FALGOUT_COARSENING;
        case escript::SO_COARSENING_HMIS:
            return PASO_HMIS_COARSENING;
        case escript::SO_COARSENING_PMIS:
            return PASO_PMIS_COARSENING;
        case escript::SO_COARSENING_RUGE_STUEBEN:
            return PASO_RUGE_STUEBEN_COARSENING;
        case escript::SO_COARSENING_STANDARD:
            return PASO_STANDARD_COARSENING;   
        case escript::SO_COARSENING_YAIR_SHAPIRA:
            return PASO_YAIR_SHAPIRA_COARSENING;

        case escript::SO_REORDERING_DEFAULT:
            return PASO_DEFAULT_REORDERING;
        case escript::SO_REORDERING_MINIMUM_FILL_IN:
            return PASO_MINIMUM_FILL_IN;
        case escript::SO_REORDERING_NESTED_DISSECTION:
            return PASO_NESTED_DISSECTION;
        case escript::SO_REORDERING_NONE:
            return PASO_NO_REORDERING;

        default:
            std::stringstream temp;
            temp << "Error - Cannot map option value "<< escriptOption
                 << " onto Paso";
            throw PasoException(temp.str());
    }
}

void Options::updateEscriptDiagnostics(bp::object& options) const
{
#define SET(__key__,__val__,__type__) options.attr("_updateDiagnostics")(__key__,(__type__)__val__)
   SET("num_iter", num_iter, int);
   SET("num_level", num_level, int);
   SET("num_inner_iter", num_inner_iter, int);
   SET("time", time, double);
   SET("set_up_time", set_up_time, double);
   SET("net_time", net_time, double);
   SET("residual_norm", residual_norm, double);
   SET("converged", converged, bool);
   SET("time_step_backtracking_used", time_step_backtracking_used, bool);
   SET("coarse_level_sparsity", coarse_level_sparsity, double);
   SET("num_coarse_unknowns", num_coarse_unknowns, int);
#undef SET
}

} // namespace paso

