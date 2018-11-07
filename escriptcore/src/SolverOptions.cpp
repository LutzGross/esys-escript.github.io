
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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

#include "SolverOptions.h"
#include "EsysException.h"

#include <boost/python.hpp>

namespace bp = boost::python;

template <class R>
bool convert(bp::object bpo, R& result)
{
    bool b = bp::extract<R>(bpo).check();
    if (b)
        result = bp::extract<R>(bpo);
    return b;
}

namespace escript {

SolverBuddy::SolverBuddy() :
    target(SO_TARGET_CPU),
    package(SO_DEFAULT),
    method(SO_DEFAULT),
    preconditioner(SO_PRECONDITIONER_JACOBI),
    ode_solver(SO_ODESOLVER_LINEAR_CRANK_NICOLSON),
    smoother(SO_PRECONDITIONER_GAUSS_SEIDEL),
    reordering(SO_REORDERING_DEFAULT),
    coarsening(SO_DEFAULT),
    amg_interpolation_method(SO_INTERPOLATION_DIRECT),
    level_max(100),
    coarsening_threshold(0.25),
    sweeps(1),
    pre_sweeps(1),
    post_sweeps(1),
    tolerance(1e-8),
    absolute_tolerance(0.),
    inner_tolerance(0.9),
    drop_tolerance(0.0005),
    drop_storage(2.),
    iter_max(100000),
    inner_iter_max(10),
    truncation(20),
    restart(0),
    is_complex(false),
    symmetric(false),
    verbose(false),
    adapt_inner_tolerance(true),
    accept_convergence_failure(false),
    min_coarse_matrix_size(500),
    relaxation(0.3),
    use_local_preconditioner(false),
    min_sparsity(0.05),
    refinements(2),
    coarse_refinements(2),
    use_panel(true),
    diagonal_dominance_threshold(0.5),
    cycle_type(1)
{
    resetDiagnostics(true);
}

SolverBuddy::~SolverBuddy()
{
}


std::string SolverBuddy::getSummary() const
{
    std::stringstream out;
    out << "Solver Package = " << getName(getPackage()) << std::endl
        << "Solver target = " << getName(getSolverTarget()) << std::endl
        << "Verbosity = " << isVerbose() << std::endl
        << "Accept failed convergence = " << acceptConvergenceFailure()
        << std::endl
        << "Relative tolerance = " << getTolerance() << std::endl
        << "Absolute tolerance = " << getAbsoluteTolerance() << std::endl
        << "Symmetric problem = " << isSymmetric() << std::endl
        << "Maximum number of iteration steps = " << getIterMax() << std::endl
        << "Inner tolerance = " << getInnerTolerance() << std::endl
        << "Adapt innner tolerance = " << adaptInnerTolerance() << std::endl;

    if (getPackage() == SO_DEFAULT || getPackage() == SO_PACKAGE_PASO ||
            getPackage() == SO_PACKAGE_CUSP ||
            getPackage() == SO_PACKAGE_TRILINOS) {
        out << "Solver method = " << getName(getSolverMethod()) << std::endl;
        if (getSolverMethod() == SO_METHOD_GMRES) {
            out << "Truncation  = " << getTruncation() << std::endl
                << "Restart  = " << getRestart() << std::endl;
        } else if (getSolverMethod() == SO_PRECONDITIONER_AMG) {
            out << "Number of pre / post sweeps = " << getNumPreSweeps()
                << " / " << getNumPostSweeps() << ", " << getNumSweeps()
                << std::endl
                << "Maximum number of levels = " << getLevelMax() << std::endl
                << "Coarsening threshold = " << getCoarseningThreshold()
                << std::endl
                << "Coarsening method = " << getName(getCoarsening())
                << std::endl;
        }
        out << "Preconditioner = " << getName(getPreconditioner()) << std::endl
            << "Apply preconditioner locally = " << useLocalPreconditioner()
            << std::endl;
        switch (getPreconditioner()) {
            case SO_PRECONDITIONER_AMG:
                out << "Maximum number of levels = " << getLevelMax()
                    << std::endl
                    << "Coarsening threshold = " << getCoarseningThreshold()
                    << std::endl
                    << "Minimal sparsity on coarsest level = "
                    << getMinCoarseMatrixSparsity() << std::endl
                    << "Smoother = " << getName(getSmoother()) << std::endl
                    << "Minimum size of the coarsest level matrix = "
                    << getMinCoarseMatrixSize() << std::endl
                    << "Number of pre / post sweeps = " << getNumPreSweeps()
                    << " / " << getNumPostSweeps() << ", " << getNumSweeps()
                    << std::endl
                    << "Number of refinement steps in coarsest level solver = "
                    << getNumCoarseMatrixRefinements() << std::endl
                    << "Use node panel = " << usePanel() << std::endl
                    << "Interpolation = " << getName(getAMGInterpolation())
                    << std::endl
                    << "Threshold for diagonal dominant rows = "
                    << getDiagonalDominanceThreshold() << std::endl;
                break;
            case SO_PRECONDITIONER_AMLI:
                out << "Maximum number of levels = " << getLevelMax()
                    << std::endl
                    << "Coarsening method = " << getName(getCoarsening())
                    << std::endl
                    << "Coarsening threshold = " << getMinCoarseMatrixSize()
                    << std::endl
                    << "Minimum size of the coarsest level matrix = "
                    << getCoarseningThreshold() << std::endl
                    << "Number of pre / post sweeps = " << getNumPreSweeps()
                    << " / " << getNumPostSweeps() << ", " << getNumSweeps()
                    << std::endl;
                break;
            case SO_PRECONDITIONER_GAUSS_SEIDEL:
                out << "Number of sweeps = " << getNumSweeps() << std::endl;
                break;
            case SO_PRECONDITIONER_ILUT:
                out << "Drop tolerance = " << getDropTolerance() << std::endl;
                out << "Storage increase = " << getDropStorage() << std::endl;
                break;
            case SO_PRECONDITIONER_RILU:
                out << "Relaxation factor = " << getRelaxationFactor()
                    << std::endl;
                break;
            default:
                break;
        } // preconditioner switch
        out << "ODE solver = " << getName(getODESolver()) << std::endl;
    }
    return out.str();
}

const char* SolverBuddy::getName(int key) const
{
    switch (static_cast<SolverOptions>(key)) {
        case SO_DEFAULT: return "DEFAULT";
        case SO_TARGET_CPU: return "CPU";
        case SO_TARGET_GPU: return "GPU";

        case SO_PACKAGE_CUSP: return "CUSP";
        case SO_PACKAGE_MKL: return "MKL";
        case SO_PACKAGE_PASO: return "PASO";
        case SO_PACKAGE_TRILINOS: return "TRILINOS";
        case SO_PACKAGE_UMFPACK: return "UMFPACK";

        case SO_METHOD_BICGSTAB: return "BICGSTAB";
        case SO_METHOD_CGLS: return "CGLS";
        case SO_METHOD_CGS: return "CGS";
        case SO_METHOD_CHOLEVSKY: return "CHOLEVSKY";
        case SO_METHOD_CR: return "CR";
        case SO_METHOD_DIRECT: return "DIRECT";
        case SO_METHOD_DIRECT_MUMPS: return "DIRECT_MUMPS";
        case SO_METHOD_DIRECT_PARDISO: return "DIRECT_PARDISO";
        case SO_METHOD_DIRECT_SUPERLU: return "DIRECT_SUPERLU";
        case SO_METHOD_DIRECT_TRILINOS: return "DIRECT_TRILINOS";
        case SO_METHOD_GMRES: return "GMRES";
        case SO_METHOD_HRZ_LUMPING: return "HRZ_LUMPING";
        case SO_METHOD_ITERATIVE: return "ITERATIVE";
        case SO_METHOD_LSQR: return "LSQR";
        case SO_METHOD_MINRES: return "MINRES";
        case SO_METHOD_NONLINEAR_GMRES: return "NONLINEAR_GMRES";
        case SO_METHOD_PCG: return "PCG";
        case SO_METHOD_PRES20: return "PRES20";
        case SO_METHOD_ROWSUM_LUMPING: return "ROWSUM_LUMPING";
        case SO_METHOD_TFQMR: return "TFQMR";

        case SO_PRECONDITIONER_AMG: return "AMG";
        case SO_PRECONDITIONER_AMLI: return "AMLI";
        case SO_PRECONDITIONER_GAUSS_SEIDEL: return "GAUSS_SEIDEL";
        case SO_PRECONDITIONER_ILU0: return "ILU0";
        case SO_PRECONDITIONER_ILUT: return "ILUT";
        case SO_PRECONDITIONER_JACOBI: return "JACOBI";
        case SO_PRECONDITIONER_NONE: return "NO_PRECONDITIONER";
        case SO_PRECONDITIONER_REC_ILU: return "REC_ILU";
        case SO_PRECONDITIONER_RILU: return "RILU";

        case SO_ODESOLVER_BACKWARD_EULER: return "BACKWARD_EULER";
        case SO_ODESOLVER_CRANK_NICOLSON: return "CRANK_NICOLSON";
        case SO_ODESOLVER_LINEAR_CRANK_NICOLSON: return "LINEAR_CRANK_NICOLSON";

        case SO_INTERPOLATION_CLASSIC: return "CLASSIC_INTERPOLATION";
        case SO_INTERPOLATION_CLASSIC_WITH_FF_COUPLING:
            return "CLASSIC_INTERPOLATION_WITH_FF";
        case SO_INTERPOLATION_DIRECT: return "DIRECT_INTERPOLATION";

        case SO_COARSENING_AGGREGATION: return "AGGREGATION_COARSENING";
        case SO_COARSENING_CIJP: return "CIJP_COARSENING";
        case SO_COARSENING_CIJP_FIXED_RANDOM:
            return "CIJP_FIXED_RANDOM_COARSENING";
        case SO_COARSENING_FALGOUT: return "FALGOUT_COARSENING";
        case SO_COARSENING_HMIS: return "HMIS_COARSENING";
        case SO_COARSENING_PMIS: return "PMIS_COARSENING";
        case SO_COARSENING_RUGE_STUEBEN: return "RUGE_STUEBEN_COARSENING";
        case SO_COARSENING_STANDARD: return "STANDARD_COARSENING";
        case SO_COARSENING_YAIR_SHAPIRA: return "YAIR_SHAPIRA_COARSENING";

        case SO_REORDERING_DEFAULT: return "DEFAULT_REORDERING";
        case SO_REORDERING_MINIMUM_FILL_IN: return "MINIMUM_FILL_IN";
        case SO_REORDERING_NESTED_DISSECTION: return "NESTED_DISSECTION";
        case SO_REORDERING_NONE: return "NO_REORDERING";
        default:
            throw ValueError("getName() invalid option given");
    }
    return "invalid option";
}

void SolverBuddy::resetDiagnostics(bool all)
{
    num_iter = 0;
    num_level = 0;
    num_inner_iter = 0;
    time = 0.;
    set_up_time = 0.;
    net_time = 0.;
    residual_norm = 0.;
    converged = false;
    preconditioner_size = -1;
    time_step_backtracking_used = false;
    coarse_level_sparsity = 0;
    num_coarse_unknowns = 0;
    if (all) {
        cum_num_inner_iter = 0;
        cum_num_iter = 0;
        cum_time = 0.;
        cum_set_up_time = 0.;
        cum_net_time = 0.;
    }
}

void SolverBuddy::updateDiagnostics(const std::string& name, bool value)
{
    if (name == "converged") {
        converged = value;
    } else if (name == "time_step_backtracking_used") {
        time_step_backtracking_used = value;
    } else {
        throw ValueError(std::string("Unknown diagnostic: ") + name);
    }
}

void SolverBuddy::updateDiagnostics(const std::string& name, int value)
{
    if (name == "num_iter") {
        cum_num_iter += num_iter = value;
    } else if (name == "num_level") {
        num_level = value;
    } else if (name == "num_inner_iter") {
        cum_num_inner_iter += num_inner_iter = value;
    } else if (name == "num_coarse_unknowns") {
        num_coarse_unknowns = value;
    } else {
        throw ValueError(std::string("Unknown diagnostic: ") + name);
    }
}

void SolverBuddy::updateDiagnostics(const std::string& name, double value)
{
    if (name == "time") {
        cum_time += time = value;
    } else if (name == "set_up_time") {
        cum_set_up_time += set_up_time = value;
    } else if (name == "net_time") {
        cum_net_time += net_time = value;
    } else if (name == "residual_norm") {
        residual_norm = value;
    } else if (name == "coarse_level_sparsity") {
        coarse_level_sparsity = value;
    } else {
        throw ValueError(std::string("Unknown diagnostic: ") + name);
    }
}

void SolverBuddy::updateDiagnosticsPy(const std::string& name,
                                      const bp::object& value)
{
    int i=0;
    double d=0; // to keep older compilers happy
    bool b=false;
    bool ib = convert<int>(value, i);
    bool db = convert<double>(value, d);
    bool bb = convert<bool>(value, b);

    if (name == "num_iter") {
        if (!ib)
            throw ValueError("setting num_iter to non-int value");
        cum_num_iter += num_iter = i;
    } else if (name == "num_level") {
        if (!ib)
            throw ValueError("setting num_level to non-int value");
        num_level = i;
    } else if (name == "num_inner_iter") {
        if (!ib)
            throw ValueError("setting num_inner_iter to non-int value");
        cum_num_inner_iter += num_inner_iter = i;
    } else if (name == "time") {
        if (!db)
            throw ValueError("setting time to non-double value");
        cum_time += time = d;
    } else if (name == "set_up_time") {
        if (!db)
            throw ValueError("setting set_up_time to non-double value");
        cum_set_up_time += set_up_time = d;
    } else if (name == "net_time") {
        if (!db)
            throw ValueError("setting net_time to non-double value");
        cum_net_time += net_time = d;
    } else if (name == "residual_norm") {
        if (!db)
            throw ValueError("setting residual_norm to non-double value");
        residual_norm = d;
    } else if (name == "converged") {
        if (!bb)
            throw ValueError("setting converged to non-bool value");
        converged = b;
    } else if (name == "time_step_backtracking_used") {
        if (!bb)
            throw ValueError("setting time_step_backtracking_used to non-bool value");
        time_step_backtracking_used = b;
    } else if (name == "coarse_level_sparsity") {
        if (!db)
            throw ValueError("setting coarse_level_sparsity to non-double value");
        coarse_level_sparsity = d;
    } else if (name == "num_coarse_unknowns") {
        if (!ib)
            throw ValueError("setting num_coarse_unknowns to non-int value");
        num_coarse_unknowns = i;
    } else {
        throw ValueError(std::string("Unknown diagnostic: ") + name);
    }
}

double SolverBuddy::getDiagnostics(const std::string name) const
{
    if (name == "num_iter") return num_iter;
    else if (name == "cum_num_iter") return cum_num_iter;
    else if (name == "num_level") return num_level;
    else if (name == "num_inner_iter") return num_inner_iter;
    else if (name == "cum_num_inner_iter") return cum_num_inner_iter;
    else if (name == "time") return time;
    else if (name == "cum_time") return cum_time;
    else if (name == "set_up_time") return set_up_time;
    else if (name == "cum_set_up_time") return cum_set_up_time;
    else if (name == "net_time") return net_time;
    else if (name == "cum_net_time") return cum_net_time;
    else if (name == "residual_norm") return residual_norm;
    else if (name == "converged") return converged;
    else if (name == "preconditioner_size") return preconditioner_size;
    else if (name == "time_step_backtracking_used")
        return  time_step_backtracking_used;
    else if (name == "coarse_level_sparsity") return coarse_level_sparsity;
    else if (name == "num_coarse_unknowns") return num_coarse_unknowns;

    throw ValueError(std::string("unknown diagnostic item: ") + name);
}

bool SolverBuddy::hasConverged() const
{
    return converged;
}

void SolverBuddy::setCoarsening(int method)
{
    SolverOptions meth = static_cast<SolverOptions>(method);
    switch (meth) {
        case SO_DEFAULT:
        case SO_COARSENING_AGGREGATION:
        case SO_COARSENING_CIJP:
        case SO_COARSENING_CIJP_FIXED_RANDOM:
        case SO_COARSENING_FALGOUT:
        case SO_COARSENING_HMIS:
        case SO_COARSENING_PMIS:
        case SO_COARSENING_RUGE_STUEBEN:
        case SO_COARSENING_STANDARD:
        case SO_COARSENING_YAIR_SHAPIRA:
            coarsening = meth;
            break;
        default:
            throw ValueError("unknown coarsening method");
    }
}

SolverOptions SolverBuddy::getCoarsening() const
{
    return coarsening;
}

void SolverBuddy::setMinCoarseMatrixSize(int size)
{
    if (size < 0) {
        throw ValueError("minimum size of the coarsest level "
                                     "matrix must be non-negative.");
    }
    min_coarse_matrix_size = size;
}

int SolverBuddy::getMinCoarseMatrixSize() const
{
    return min_coarse_matrix_size;
}

void SolverBuddy::setPreconditioner(int precon)
{
    SolverOptions preconditioner = static_cast<SolverOptions>(precon);
    switch(preconditioner) {
        case SO_PRECONDITIONER_AMG:
/*
#ifdef ESYS_MPI
            throw ValueError("AMG preconditioner is not supported in MPI builds");
            break;
#endif
*/
        case SO_PRECONDITIONER_AMLI:
        case SO_PRECONDITIONER_GAUSS_SEIDEL:
        case SO_PRECONDITIONER_ILU0:
        case SO_PRECONDITIONER_ILUT:
        case SO_PRECONDITIONER_JACOBI:
        case SO_PRECONDITIONER_NONE:
        case SO_PRECONDITIONER_REC_ILU:
        case SO_PRECONDITIONER_RILU:
            this->preconditioner = preconditioner;
            break;
        default:
            throw ValueError("unknown preconditioner");
    }
}

SolverOptions SolverBuddy::getPreconditioner() const
{
    return preconditioner;
}

void SolverBuddy::setSmoother(int s)
{
    SolverOptions smoother = static_cast<SolverOptions>(s);
    if (smoother != SO_PRECONDITIONER_JACOBI &&
            smoother != SO_PRECONDITIONER_GAUSS_SEIDEL) {
        throw ValueError("unknown smoother");
    }
    this->smoother = smoother;
}

SolverOptions SolverBuddy::getSmoother() const
{
    return smoother;
}

void SolverBuddy::setSolverMethod(int method)
{
    SolverOptions meth = static_cast<SolverOptions>(method);
    switch(meth) {
        case SO_DEFAULT:
        case SO_METHOD_BICGSTAB:
        case SO_METHOD_CGLS:
        case SO_METHOD_CGS:
        case SO_METHOD_CHOLEVSKY:
        case SO_METHOD_CR:
        case SO_METHOD_GMRES:
        case SO_METHOD_HRZ_LUMPING:
        case SO_METHOD_ITERATIVE:
        case SO_METHOD_LSQR:
        case SO_METHOD_MINRES:
        case SO_METHOD_NONLINEAR_GMRES:
        case SO_METHOD_PCG:
        case SO_METHOD_PRES20:
        case SO_METHOD_ROWSUM_LUMPING:
        case SO_METHOD_TFQMR:
            this->method = meth;
            break;
        case SO_METHOD_DIRECT:
        case SO_METHOD_DIRECT_MUMPS:
        case SO_METHOD_DIRECT_PARDISO:
        case SO_METHOD_DIRECT_SUPERLU:
        case SO_METHOD_DIRECT_TRILINOS:
#if defined(ESYS_HAVE_UMFPACK) || defined(ESYS_HAVE_TRILINOS) || defined(ESYS_HAVE_MKL)
#ifndef ESYS_HAVE_TRILINOS
            // translate specific direct solver setting to generic one for PASO
            this->method = SO_METHOD_DIRECT;
#else
            this->method = meth;
#endif
            break;
#else
            throw ValueError("Cannot use DIRECT solver method, the running "
                    "escript was not compiled with a direct solver enabled");
#endif
        default:
            throw ValueError("unknown solver method");
    }
}

SolverOptions SolverBuddy::getSolverMethod() const
{
    return method;
}

void SolverBuddy::setSolverTarget(int target)
{
    SolverOptions targ = static_cast<SolverOptions>(target);
    switch (targ) {
        case SO_TARGET_CPU:
        case SO_TARGET_GPU:
            this->target = targ;
            break;
        default:
            throw ValueError("unknown solver target");
    }
}

SolverOptions SolverBuddy::getSolverTarget() const
{
    return target;
}

void SolverBuddy::setPackage(int package)
{
    SolverOptions pack = static_cast<SolverOptions>(package);
    switch (pack) {
        case SO_DEFAULT:
        case SO_PACKAGE_CUSP:
        case SO_PACKAGE_MKL:
        case SO_PACKAGE_PASO:
        case SO_PACKAGE_TRILINOS:
        case SO_PACKAGE_UMFPACK:
            this->package = pack;
            break;
        default:
            throw ValueError("unknown solver package");
    }
}

SolverOptions SolverBuddy::getPackage() const
{
    return package;
}

void SolverBuddy::setReordering(int ordering)
{
    SolverOptions ord = static_cast<SolverOptions>(ordering);
    switch (ordering) {
        case SO_REORDERING_DEFAULT:
        case SO_REORDERING_MINIMUM_FILL_IN:
        case SO_REORDERING_NESTED_DISSECTION:
        case SO_REORDERING_NONE:
            reordering = ord;
            break;
        default:
            throw ValueError("unknown reordering strategy");
    }
}

SolverOptions SolverBuddy::getReordering() const
{
    return reordering;
}

void SolverBuddy::setRestart(int restart)
{
    if (restart < 0)
        throw ValueError("restart must be non-negative.");

    this->restart = restart;
}

int SolverBuddy::getRestart() const
{
    return restart;
}

int SolverBuddy::_getRestartForC() const
{
    int r = getRestart();
    if (r == 0)
        return -1;
    return r;
}

void SolverBuddy::setDiagonalDominanceThreshold(double value)
{
    if (value < 0. || value > 1.)
        throw ValueError("Diagonal dominance threshold must be between 0 and 1.");
    diagonal_dominance_threshold = value;
}

double SolverBuddy::getDiagonalDominanceThreshold() const
{
    return diagonal_dominance_threshold;
}

void SolverBuddy::setTruncation(int truncation)
{
    if (truncation < 1)
        throw ValueError("truncation must be positive.");
    this->truncation = truncation;
}

int SolverBuddy::getTruncation() const
{
    return truncation;
}

void SolverBuddy::setInnerIterMax(int iter_max)
{
    if (iter_max < 1)
        throw ValueError("maximum number of inner iteration must be positive.");
    inner_iter_max = iter_max;
}

int SolverBuddy::getInnerIterMax() const
{
    return inner_iter_max;
}

void SolverBuddy::setIterMax(int iter_max)
{
    if (iter_max < 1)
        throw ValueError("maximum number of iteration steps must be positive.");
    this->iter_max = iter_max;
}

int SolverBuddy::getIterMax() const
{
    return iter_max;
}

void SolverBuddy::setLevelMax(int level_max)
{
    if (level_max < 0)
        throw ValueError("maximum number of coarsening levels must be non-negative.");
    this->level_max = level_max;
}

int SolverBuddy::getLevelMax() const
{
    return level_max;
}

void SolverBuddy::setCycleType(int cycle_type)
{
    this->cycle_type = cycle_type;
}

int SolverBuddy::getCycleType() const
{
    return cycle_type;
}

void SolverBuddy::setCoarseningThreshold(double theta)
{
    if (theta < 0. || theta > 1.)
        throw ValueError("threshold must be between 0 and 1.");
    coarsening_threshold = theta;
}

double SolverBuddy::getCoarseningThreshold() const
{
    return coarsening_threshold;
}

void SolverBuddy::setNumSweeps(int sweeps)
{
    if (sweeps < 1)
        throw ValueError("number of sweeps must be positive.");
    this->sweeps = sweeps;
}

int SolverBuddy::getNumSweeps() const
{
    return sweeps;
}

void SolverBuddy::setNumPreSweeps(int sweeps)
{
    if (sweeps < 1)
        throw ValueError("number of pre-sweeps must be positive.");
    pre_sweeps = sweeps;
}

int SolverBuddy::getNumPreSweeps() const
{
    return pre_sweeps;
}

void SolverBuddy::setNumPostSweeps(int sweeps)
{
    if (sweeps < 1)
       throw ValueError("number of post-sweeps must be positive.");
    post_sweeps = sweeps;
}

int SolverBuddy::getNumPostSweeps() const
{
    return post_sweeps;
}

void SolverBuddy::setTolerance(double rtol)
{
    if (rtol < 0. || rtol > 1.)
        throw ValueError("tolerance must be between 0 and 1.");
    tolerance = rtol;
}

double SolverBuddy::getTolerance() const
{
    return tolerance;
}

void SolverBuddy::setAbsoluteTolerance(double atol)
{
    if (atol < 0.)
       throw ValueError("absolute tolerance must be non-negative.");
    absolute_tolerance = atol;
}

double SolverBuddy::getAbsoluteTolerance() const
{
    return absolute_tolerance;
}

void SolverBuddy::setInnerTolerance(double rtol)
{
    if (rtol <= 0. || rtol > 1.)
        throw ValueError("tolerance must be positive and less than or equal to 1.");
    inner_tolerance = rtol;
}

double SolverBuddy::getInnerTolerance() const
{
    return inner_tolerance;
}

void SolverBuddy::setDropTolerance(double drop_tol)
{
    if (drop_tol < 0. || drop_tol > 1.)
        throw ValueError("drop tolerance must be between 0 and 1.");
    drop_tolerance = drop_tol;
}

double SolverBuddy::getDropTolerance() const
{
    return drop_tolerance;
}

void SolverBuddy::setDropStorage(double storage)
{
    if (storage < 1.)
        throw ValueError("allowed storage increase must be greater than or equal to 1.");
    drop_storage = storage;
}

double SolverBuddy::getDropStorage() const
{
    return drop_storage;
}

void SolverBuddy::setRelaxationFactor(double factor)
{
    if (factor < 0.)
        throw ValueError("relaxation factor must be non-negative.");
    relaxation = factor;
}

double SolverBuddy::getRelaxationFactor() const
{
    return relaxation;
}

bool SolverBuddy::isComplex() const
{
    return is_complex;
}

void SolverBuddy::setComplex(bool flag)
{
    is_complex = flag;
}

bool SolverBuddy::isSymmetric() const
{
    return symmetric;
}

void SolverBuddy::setSymmetryOn()
{
    symmetric = true;
}

void SolverBuddy::setSymmetryOff()
{
    symmetric = false;
}

void SolverBuddy::setSymmetry(bool flag)
{
    if (flag)
        setSymmetryOn();
    else
        setSymmetryOff();
}

bool SolverBuddy::isVerbose() const
{
    return verbose;
}

void SolverBuddy::setVerbosityOn()
{
    verbose = true;
}

void SolverBuddy::setVerbosityOff()
{
    verbose = false;
}

void SolverBuddy::setVerbosity(bool verbose)
{
    if (verbose)
        setVerbosityOn();
    else
        setVerbosityOff();
}

bool SolverBuddy::adaptInnerTolerance() const
{
    return adapt_inner_tolerance;
}

void SolverBuddy::setInnerToleranceAdaptionOn()
{
    adapt_inner_tolerance = true;
}

void SolverBuddy::setInnerToleranceAdaptionOff()
{
    adapt_inner_tolerance = false;
}

void SolverBuddy::setInnerToleranceAdaption(bool adapt)
{
    if (adapt)
        setInnerToleranceAdaptionOn();
    else
        setInnerToleranceAdaptionOff();
}

bool SolverBuddy::acceptConvergenceFailure() const
{
    return accept_convergence_failure;
}

void SolverBuddy::setAcceptanceConvergenceFailureOn()
{
    accept_convergence_failure = true;
}

void SolverBuddy::setAcceptanceConvergenceFailureOff()
{
    accept_convergence_failure = false;
}

void SolverBuddy::setAcceptanceConvergenceFailure(bool accept)
{
    if (accept)
        setAcceptanceConvergenceFailureOn();
    else
        setAcceptanceConvergenceFailureOff();
}

bool SolverBuddy::useLocalPreconditioner() const
{
    return use_local_preconditioner;
}

void SolverBuddy::setLocalPreconditionerOn()
{
    use_local_preconditioner = true;
}

void SolverBuddy::setLocalPreconditionerOff()
{
    use_local_preconditioner=false;
}

void SolverBuddy::setLocalPreconditioner(bool use)
{
    if (use)
        setLocalPreconditionerOn();
    else
        setLocalPreconditionerOff();
}

void SolverBuddy::setMinCoarseMatrixSparsity(double sparsity)
{
    if (sparsity < 0. || sparsity > 1.)
        throw ValueError("sparsity must be between 0 and 1.");
    min_sparsity = sparsity;
}

double SolverBuddy::getMinCoarseMatrixSparsity() const
{
    return min_sparsity;
}

void SolverBuddy::setNumRefinements(int refinements)
{
    if (refinements < 0)
        throw ValueError("number of refinements must be non-negative.");
    this->refinements = refinements;
}

int SolverBuddy::getNumRefinements() const
{
    return refinements;
}

void SolverBuddy::setNumCoarseMatrixRefinements(int refinements)
{
    if (refinements < 0)
        throw ValueError("number of coarse matrix refinements must be non-negative.");
    coarse_refinements = refinements;
}

int SolverBuddy::getNumCoarseMatrixRefinements() const
{
    return coarse_refinements;
}

bool SolverBuddy::usePanel() const
{
    return use_panel;
}

void SolverBuddy::setUsePanelOn()
{
    use_panel = true;
}

void SolverBuddy::setUsePanelOff()
{
    use_panel = false;
}

void SolverBuddy::setUsePanel(bool use)
{
    if (use)
        setUsePanelOn();
    else
        setUsePanelOff();
}

void SolverBuddy::setAMGInterpolation(int method)
{
    SolverOptions meth = static_cast<SolverOptions>(method);
    switch (meth) {
        case SO_INTERPOLATION_CLASSIC:
        case SO_INTERPOLATION_CLASSIC_WITH_FF_COUPLING:
        case SO_INTERPOLATION_DIRECT:
            amg_interpolation_method = meth;
            break;
        default:
            throw ValueError("unknown AMG interpolation method");
    }
}

SolverOptions SolverBuddy::getAMGInterpolation() const
{
    return amg_interpolation_method;
}

void SolverBuddy::setODESolver(int method)
{
    SolverOptions ode = static_cast<SolverOptions>(method);
    switch (ode) {
        case SO_ODESOLVER_BACKWARD_EULER:
        case SO_ODESOLVER_CRANK_NICOLSON:
        case SO_ODESOLVER_LINEAR_CRANK_NICOLSON:
            ode_solver = ode;
            break;
        default:
            throw ValueError("unknown ODE solver method");
    }
}

SolverOptions SolverBuddy::getODESolver() const
{
    return ode_solver;
}

void SolverBuddy::setTrilinosParameter(const std::string& name,
                                       const bp::object& value)
{
#ifdef ESYS_HAVE_TRILINOS
    trilinosParams[name] = value;
#endif
}

bp::dict SolverBuddy::getTrilinosParameters() const
{
    return trilinosParams;
}

} // namespace escript

