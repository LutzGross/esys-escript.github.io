
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

#include <boost/python.hpp>


#include "SolverOptions.h"
#include "SolverOptionsException.h"


template <class R>
bool convert(boost::python::object bpo, R& result) {
    bool b = boost::python::extract<R>(bpo).check();
    if (b)
        result = boost::python::extract<R>(bpo);
    return b;
}

namespace escript {

SolverBuddy::SolverBuddy()
{
    level_max = 100;
    coarsening_threshold = 0.25;
    smoother = ESCRIPT_GAUSS_SEIDEL;
    sweeps = 1;
    pre_sweeps = 1;
    post_sweeps = 1;
    tolerance = 1e-8;
    absolute_tolerance = 0.;
    inner_tolerance = 0.9;
    drop_tolerance = 0.01;
    drop_storage = 2.;
    iter_max = 100000;
    inner_iter_max = 10;
    truncation = 20;
    restart = 0;
    symmetric = false;
    verbose = false;
    adapt_inner_tolerance = true;
    accept_convergence_failure = false;
    reordering = ESCRIPT_DEFAULT_REORDERING;
    package = ESCRIPT_DEFAULT;
    method = ESCRIPT_DEFAULT;
    preconditioner = ESCRIPT_JACOBI;
    coarsening = ESCRIPT_DEFAULT;
    target = ESCRIPT_TARGET_CPU;
    MinCoarseMatrixSize = 500;
    relaxation = 0.3;
    use_local_preconditioner = false;
    min_sparsity = 0.05;
    refinements = 2;
    coarse_refinements = 2;
    use_panel = true;
    diagonal_dominance_threshold = 0.5;
    amg_interpolation_method = ESCRIPT_DIRECT_INTERPOLATION;
    cycle_type = 1;
    ode_solver = ESCRIPT_LINEAR_CRANK_NICOLSON;

    resetDiagnostics(true);
}

std::string SolverBuddy::getSummary() const {
    std::stringstream out;
    out << "Solver Package: " << getName(getPackage())
        << std::endl
        << "Verbosity = " << isVerbose() << std::endl
        << "Accept failed convergence = " << acceptConvergenceFailure()
        << std::endl
        << "Relative tolerance = " << getTolerance()
        << std::endl
        << "Absolute tolerance = " << getAbsoluteTolerance()
        << std::endl
        << "Symmetric problem = " << isSymmetric()
        << std::endl
        << "Maximum number of iteration steps = " << getIterMax()
        << std::endl
        << "Inner tolerance = " << getInnerTolerance()
        << std::endl
        << "Adapt innner tolerance = " << adaptInnerTolerance()
        << std::endl;

    if (getPackage() == ESCRIPT_DEFAULT || getPackage() == ESCRIPT_PASO) {
        out << "Solver method = " << getName(getSolverMethod())
            << std::endl;
        if (getSolverMethod() == ESCRIPT_GMRES) {
            out << "Truncation  = " << getTruncation()
                << std::endl
                << "Restart  = " << getRestart()
                << std::endl;
        } else if (getSolverMethod() == ESCRIPT_AMG) {
            out << "Number of pre / post sweeps = " << getNumPreSweeps()
                << " / " << getNumPostSweeps() << ", " << getNumSweeps()
                << std::endl
                << "Maximum number of levels = " << getLevelMax()
                << std::endl
                << "Coarsening threshold = " << getCoarseningThreshold()
                << std::endl
                << "Coarsening method = " << getName(getCoarsening())
                << std::endl;
        }
        out << "Preconditioner = " << getName(getPreconditioner())
            << std::endl
            << "Apply preconditioner locally = " << useLocalPreconditioner()
            << std::endl;
        if (getPreconditioner() == ESCRIPT_AMG) {
            out << "Maximum number of levels = " << getLevelMax()
                << std::endl
                << "Coarsening threshold = " << getCoarseningThreshold()
                << std::endl
                << "Minimal sparsity on coarsest level = "
                << getMinCoarseMatrixSparsity()
                << std::endl
                << "Smoother = " << getName(getSmoother())
                << std::endl
                << "Minimum size of the coarsest level matrix = "
                << getMinCoarseMatrixSize()
                << std::endl
                << "Number of pre / post sweeps = " << getNumPreSweeps()
                << " / " << getNumPostSweeps() << ", " << getNumSweeps()
                << std::endl
                << "Number of refinement steps in coarsest level solver = "
                << getNumCoarseMatrixRefinements()
                << std::endl
                << "Use node panel = " << usePanel()
                << std::endl
                << "Interpolation = " << getName(getAMGInterpolation())
                << std::endl
                << "Threshold for diagonal dominant rows = "
                << getDiagonalDominanceThreshold()
                << std::endl;
        }
        if (getPreconditioner() == ESCRIPT_AMLI) {
            out << "Maximum number of levels = " << getLevelMax()
                << std::endl
                << "Coarsening method = " << getName(getCoarsening())
                << std::endl
                << "Coarsening threshold = " << getMinCoarseMatrixSize()
                << std::endl
                << "Minimum size of the coarsest level matrix = "
                << getCoarseningThreshold()
                << std::endl
                << "Number of pre / post sweeps = " << getNumPreSweeps()
                << " / " << getNumPostSweeps() << ", " << getNumSweeps()
                << std::endl;
        } else if (getPreconditioner() == ESCRIPT_BOOMERAMG) {
            out << "Maximum number of levels = " << getLevelMax()
                << std::endl
                << "Coarsening threshold = " << getCoarseningThreshold()
                << std::endl
                << "Threshold for diagonal dominant rows = "
                << getDiagonalDominanceThreshold()
                << std::endl
                << "Coarsening method = " << getName(getCoarsening())
                << std::endl
                << "V-cycle (1) or W-cyle (2) = " << getCycleType()
                << std::endl
                << "Number of pre / post sweeps = " << getNumPreSweeps()
                << " / " << getNumPostSweeps() << ", " << getNumSweeps()
                << std::endl
                << "Smoother = " << getName(getSmoother())
                << std::endl;
        }
        if (getPreconditioner() == ESCRIPT_GAUSS_SEIDEL) {
            out << "Number of sweeps = " << getNumSweeps() << std::endl;
        } else if (getPreconditioner() == ESCRIPT_ILUT) {
            out << "Drop tolerance = " << getDropTolerance() << std::endl;
            out << "Storage increase = " << getDropStorage() << std::endl;
        } else if (getPreconditioner() == ESCRIPT_RILU) {
            out << "Relaxation factor = " << getRelaxationFactor() << std::endl;
        }
        out << "ODE solver = " << getName(getODESolver()) << std::endl;
        out << "Solver target = " << getName(getSolverTarget()) << std::endl;
    }
    return out.str();
}

const char *SolverBuddy::getName(int key) const {
    SolverOptions k = static_cast<SolverOptions>(key);
    switch(k) {
        case ESCRIPT_DEFAULT: return "DEFAULT";
        case ESCRIPT_DIRECT: return "DIRECT";
        case ESCRIPT_CHOLEVSKY: return "CHOLEVSKY";
        case ESCRIPT_PCG: return "PCG";
        case ESCRIPT_CR: return "CR";
        case ESCRIPT_CGS: return "CGS";
        case ESCRIPT_BICGSTAB: return "BICGSTAB";
        case ESCRIPT_ILU0: return "ILU0";
        case ESCRIPT_ILUT: return "ILUT";
        case ESCRIPT_JACOBI: return "JACOBI";
        case ESCRIPT_GMRES: return "GMRES";
        case ESCRIPT_LSQR: return "LSQR";
        case ESCRIPT_PRES20: return "PRES20";
        case ESCRIPT_ROWSUM_LUMPING: return "ROWSUM_LUMPING";
        case ESCRIPT_HRZ_LUMPING: return "HRZ_LUMPING";
        case ESCRIPT_NO_REORDERING: return "NO_REORDERING";
        case ESCRIPT_MINIMUM_FILL_IN: return "MINIMUM_FILL_IN";
        case ESCRIPT_NESTED_DISSECTION: return "NESTED_DISSECTION";
        case ESCRIPT_MKL: return "MKL";
        case ESCRIPT_UMFPACK: return "UMFPACK";
        case ESCRIPT_ITERATIVE: return "ITERATIVE";
        case ESCRIPT_PASO: return "PASO";
        case ESCRIPT_AMG: return "AMG";
        case ESCRIPT_AMLI: return "AMLI";
        case ESCRIPT_REC_ILU: return "REC_ILU";
        case ESCRIPT_TRILINOS: return "TRILINOS";
        case ESCRIPT_NONLINEAR_GMRES: return "NONLINEAR_GMRES";
        case ESCRIPT_TFQMR: return "TFQMR";
        case ESCRIPT_MINRES: return "MINRES";
        case ESCRIPT_GAUSS_SEIDEL: return "GAUSS_SEIDEL";
        case ESCRIPT_RILU: return "RILU";
        case ESCRIPT_DEFAULT_REORDERING: return "DEFAULT_REORDERING";
        case ESCRIPT_SUPER_LU: return "SUPER_LU";
        case ESCRIPT_PASTIX: return "PASTIX";
        case ESCRIPT_YAIR_SHAPIRA_COARSENING: return "YAIR_SHAPIRA_COARSENING";
        case ESCRIPT_RUGE_STUEBEN_COARSENING: return "RUGE_STUEBEN_COARSENING";
        case ESCRIPT_STANDARD_COARSENING: return "STANDARD_COARSENING";
        case ESCRIPT_AGGREGATION_COARSENING: return "AGGREGATION_COARSENING";
        case ESCRIPT_NO_PRECONDITIONER: return "NO_PRECONDITIONER";
        case ESCRIPT_CLASSIC_INTERPOLATION_WITH_FF_COUPLING:
            return "CLASSIC_INTERPOLATION_WITH_FF";
        case ESCRIPT_CLASSIC_INTERPOLATION: return "CLASSIC_INTERPOLATION";
        case ESCRIPT_DIRECT_INTERPOLATION: return "DIRECT_INTERPOLATION";
        case ESCRIPT_BOOMERAMG: return "BOOMERAMG";
        case ESCRIPT_CIJP_FIXED_RANDOM_COARSENING:
            return "CIJP_FIXED_RANDOM_COARSENING";
        case ESCRIPT_CIJP_COARSENING: return "CIJP_COARSENING";
        case ESCRIPT_FALGOUT_COARSENING: return "FALGOUT_COARSENING";
        case ESCRIPT_PMIS_COARSENING: return "PMIS_COARSENING";
        case ESCRIPT_HMIS_COARSENING: return "HMIS_COARSENING";
        case ESCRIPT_LINEAR_CRANK_NICOLSON: return "LINEAR_CRANK_NICOLSON";
        case ESCRIPT_CRANK_NICOLSON: return "CRANK_NICOLSON";
        case ESCRIPT_BACKWARD_EULER: return "BACKWARD_EULER";
        case ESCRIPT_TARGET_CPU: return "TARGET_CPU";
        case ESCRIPT_TARGET_GPU: return "TARGET_GPU";
        default:
            throw SolverOptionsException("getName() invalid option given");
    }
    return "invalid option";
}

void SolverBuddy::resetDiagnostics(bool all) {
    num_iter = 0;
    num_level = 0;
    num_inner_iter = 0;
    time = 0;
    set_up_time = 0;
    net_time = 0;
    residual_norm = 0;
    converged = 0;
    preconditioner_size = -1;
    time_step_backtracking_used = 0;
    if (all) {
        cum_num_inner_iter = 0;
        cum_num_iter = 0;
        cum_time = 0;
        cum_set_up_time = 0;
        cum_net_time = 0;
    }
}

void SolverBuddy::updateDiagnostics(std::string name, boost::python::object value) {

        int i=0;
        double d=0;	// to keep older compilers happy
        bool b=false;
        bool ib = convert<int>(value, i);
        bool db = convert<double>(value, d);
        bool bb = convert<bool>(value, b);

        if (!strcmp(name.c_str(), "num_iter")) {
            if (!ib)
                throw SolverOptionsException("setting num_iter to non-int value");
            cum_num_iter += num_iter = i;
        } else if (!strcmp(name.c_str(), "num_level")) {
            if (!ib)
                throw SolverOptionsException("setting num_level to non-int value");
            num_level = i;
        } else if (!strcmp(name.c_str(), "num_inner_iter")) {
            if (!ib)
                throw SolverOptionsException("setting num_cuM_inner_iter to non-int value");
            cum_num_inner_iter += num_inner_iter = i;
        } else if (!strcmp(name.c_str(), "time")) {
            if (!db)
                throw SolverOptionsException("setting time to non-double value");
            cum_time += time = d;
        } else if (!strcmp(name.c_str(), "set_up_time")) {
            if (!db)
                throw SolverOptionsException("setting set_up_time to non-double value");
            cum_set_up_time += set_up_time = d;
        } else if (!strcmp(name.c_str(), "net_time")) {
            if (!db)
                throw SolverOptionsException("setting net_time to non-double value");
            cum_net_time += net_time = d;
        } else if (!strcmp(name.c_str(), "residual_norm")) {
            if (!db)
                throw SolverOptionsException("setting residual_norm to non-double value");
            residual_norm = d;
        } else if (!strcmp(name.c_str(), "converged")) {
            if (!bb)
                throw SolverOptionsException("setting converged to non-bool value");
            converged = (b == true);
        } else if (!strcmp(name.c_str(), "time_step_backtracking_used")) {
            if (!bb)
                throw SolverOptionsException("setting converged to non-bool value");
            time_step_backtracking_used = (b == true);
        } else if (!strcmp(name.c_str(), "coarse_level_sparsity")) {
            if (!db)
                throw SolverOptionsException("setting coarse_level_sparsity to non-double value");
            coarse_level_sparsity = i;
        } else if (!strcmp(name.c_str(), "num_coarse_unknowns")) {
            if (!ib)
                throw SolverOptionsException("setting num_coarse_unknowns to non-int value");
            num_coarse_unknowns = i;
        } else {
            throw SolverOptionsException(std::string("Unknown diagnostic: ") + name);
        }
}

double SolverBuddy::getDiagnostics(std::string name) const {

    if (!strcmp(name.c_str(), "num_iter")) return num_iter;
    else if (!strcmp(name.c_str(), "cum_num_iter")) return cum_num_iter;
    else if (!strcmp(name.c_str(), "num_level")) return num_level;
    else if (!strcmp(name.c_str(), "num_inner_iter")) return num_inner_iter;
    else if (!strcmp(name.c_str(), "cum_num_inner_iter"))
        return cum_num_inner_iter;
    else if (!strcmp(name.c_str(), "time")) return time;
    else if (!strcmp(name.c_str(), "cum_time")) return cum_time;
    else if (!strcmp(name.c_str(), "set_up_time")) return set_up_time;
    else if (!strcmp(name.c_str(), "cum_set_up_time")) return cum_set_up_time;
    else if (!strcmp(name.c_str(), "net_time")) return net_time;
    else if (!strcmp(name.c_str(), "cum_net_time")) return cum_net_time;
    else if (!strcmp(name.c_str(), "residual_norm")) return residual_norm;
    else if (!strcmp(name.c_str(), "converged")) return converged;
    else if (!strcmp(name.c_str(), "preconditioner_size"))
        return preconditioner_size;
    else if (!strcmp(name.c_str(), "time_step_backtracking_used"))
        return  time_step_backtracking_used;
    else if (!strcmp(name.c_str(), "coarse_level_sparsity"))
        return coarse_level_sparsity;
    else if (!strcmp(name.c_str(), "num_coarse_unknowns"))
        return num_coarse_unknowns;
    else
        throw SolverOptionsException(std::string("unknown diagnostic"
                " item: ") + name);
}

bool SolverBuddy::hasConverged() const {
    return converged;
}

void SolverBuddy::setCoarsening(int method) {
    SolverOptions meth = static_cast<SolverOptions>(method);
        switch (meth) {
            case ESCRIPT_DEFAULT:
            case ESCRIPT_YAIR_SHAPIRA_COARSENING:
            case ESCRIPT_RUGE_STUEBEN_COARSENING:
            case ESCRIPT_AGGREGATION_COARSENING:
            case ESCRIPT_STANDARD_COARSENING:
            case ESCRIPT_CIJP_FIXED_RANDOM_COARSENING:
            case ESCRIPT_CIJP_COARSENING:
            case ESCRIPT_FALGOUT_COARSENING:
            case ESCRIPT_PMIS_COARSENING:
            case ESCRIPT_HMIS_COARSENING:
                coarsening = meth;
                break;
            default:
                throw SolverOptionsException("unknown coarsening method");
        }
}

SolverOptions SolverBuddy::getCoarsening() const {
    return coarsening;
}

void SolverBuddy::setMinCoarseMatrixSize(int size) {
    if (size<0) {
        throw SolverOptionsException("minimum size of the coarsest level matrix must be non-negative.");
    }
    MinCoarseMatrixSize=size;
}

int SolverBuddy::getMinCoarseMatrixSize() const {
    return MinCoarseMatrixSize;
}

void SolverBuddy::setPreconditioner(int precon) {
    SolverOptions preconditioner = static_cast<SolverOptions>(precon);
    switch(preconditioner) {
        case ESCRIPT_ILU0:
        case ESCRIPT_ILUT:
        case ESCRIPT_JACOBI:
        case ESCRIPT_AMLI:
        case ESCRIPT_REC_ILU:
        case ESCRIPT_GAUSS_SEIDEL:
        case ESCRIPT_RILU:
        case ESCRIPT_BOOMERAMG:
        case ESCRIPT_NO_PRECONDITIONER:
#ifndef ESYS_MPI
        case ESCRIPT_AMG:
            this->preconditioner = preconditioner;
            break;
#else
            this->preconditioner = preconditioner;
            break;
        case ESCRIPT_AMG:
            throw SolverOptionsException("AMG preconditioner is not supported in MPI builds");
            break;
#endif
        default:
            throw SolverOptionsException("unknown preconditioner");
    }
}

SolverOptions SolverBuddy::getPreconditioner() const {
    return preconditioner;
}

void SolverBuddy::setSmoother(int s) {
    SolverOptions smoother = static_cast<SolverOptions>(s);
    if (smoother != ESCRIPT_JACOBI && smoother != ESCRIPT_GAUSS_SEIDEL) {
        throw SolverOptionsException("unknown smoother");
    }
    this->smoother = smoother;
}

SolverOptions SolverBuddy::getSmoother() const {
    return smoother;
}

void SolverBuddy::setSolverMethod(int method) {
    SolverOptions meth = static_cast<SolverOptions>(method);
    switch(meth) {
        case ESCRIPT_DEFAULT:
        case ESCRIPT_DIRECT:
        case ESCRIPT_CHOLEVSKY:
        case ESCRIPT_PCG:
        case ESCRIPT_CR:
        case ESCRIPT_CGS:
        case ESCRIPT_BICGSTAB:
        case ESCRIPT_GMRES:
        case ESCRIPT_PRES20:
        case ESCRIPT_ROWSUM_LUMPING:
        case ESCRIPT_HRZ_LUMPING:
        case ESCRIPT_ITERATIVE:
        case ESCRIPT_NONLINEAR_GMRES:
        case ESCRIPT_TFQMR:
        case ESCRIPT_MINRES:
        case ESCRIPT_LSQR:
            this->method = meth;
            break;
        default:
            throw SolverOptionsException("unknown solver method");
    }
}

SolverOptions SolverBuddy::getSolverMethod() const {
    return method;
}

void SolverBuddy::setSolverTarget(int target) {
    SolverOptions targ = static_cast<SolverOptions>(target);
    switch(targ) {
        case ESCRIPT_TARGET_CPU:
        case ESCRIPT_TARGET_GPU:
            this->target = targ;
            break;
        default:
            throw SolverOptionsException("unknown solver target");
    }
}

SolverOptions SolverBuddy::getSolverTarget() const
{
    return target;
}

void SolverBuddy::setPackage(int package) {
    SolverOptions pack = static_cast<SolverOptions>(package);
    switch (pack) {
        case ESCRIPT_DEFAULT:
        case ESCRIPT_PASO:
        case ESCRIPT_SUPER_LU:
        case ESCRIPT_PASTIX:
        case ESCRIPT_MKL:
        case ESCRIPT_UMFPACK:
        case ESCRIPT_TRILINOS:
            this->package = pack;
            break;
        default:
            throw SolverOptionsException("unknown solver package");
    }
}

SolverOptions SolverBuddy::getPackage() const {
    return package;
}

void SolverBuddy::setReordering(int ordering) {
    SolverOptions ord = static_cast<SolverOptions>(ordering);
    switch (ordering) {
        case ESCRIPT_NO_REORDERING:
        case ESCRIPT_MINIMUM_FILL_IN:
        case ESCRIPT_NESTED_DISSECTION:
        case ESCRIPT_DEFAULT_REORDERING:
            reordering = ord;
            break;
        default:
            throw SolverOptionsException("unknown reordering strategy");
    }
}

SolverOptions SolverBuddy::getReordering() const {
    return reordering;
}

void SolverBuddy::setRestart(int restart) {
    if (restart < 0) {
        throw SolverOptionsException("restart must be positive.");
    }
    this->restart = restart;
}

int SolverBuddy::getRestart() const {
    return restart;
}

int SolverBuddy::_getRestartForC() const {
    int r = getRestart();
    if (r == 0)
        return -1;
    return r;
}

void SolverBuddy::setDiagonalDominanceThreshold(double value) {
    if (value < 0. || value > 1.)
        throw SolverOptionsException("Diagonal dominance threshold must be between 0 and 1.");
    diagonal_dominance_threshold = value;
}

double SolverBuddy::getDiagonalDominanceThreshold() const {
    return diagonal_dominance_threshold;
}

void SolverBuddy::setTruncation(int truncation) {
    if (truncation < 1)
        throw SolverOptionsException("truncation must be positive.");
    this->truncation = truncation;
}

int SolverBuddy::getTruncation() const {
    return truncation;
}

void SolverBuddy::setInnerIterMax(int iter_max) {
    if (iter_max < 1)
        throw SolverOptionsException("maximum number of inner iteration must be positive.");
    inner_iter_max = iter_max;
}

int SolverBuddy::getInnerIterMax() const {
    return inner_iter_max;
}

void SolverBuddy::setIterMax(int iter_max) {
    if (iter_max < 1)
        throw SolverOptionsException("maximum number of iteration steps must be positive.");
    this->iter_max = iter_max;
}

int SolverBuddy::getIterMax() const {
    return iter_max;
}

void SolverBuddy::setLevelMax(int level_max) {
    if (level_max < 0)
        throw SolverOptionsException("maximum number of coarsening levels must be non-negative.");
    this->level_max = level_max;
}

int SolverBuddy::getLevelMax() const {
    return level_max;
}

void SolverBuddy::setCycleType(int cycle_type) {
    this->cycle_type = cycle_type;
}

int SolverBuddy::getCycleType() const {
    return cycle_type;
}

void SolverBuddy::setCoarseningThreshold(double theta) {
    if (theta < 0. || theta > 1)
        throw SolverOptionsException("threshold must be non-negative and less or equal 1.");
    coarsening_threshold = theta;
}

double SolverBuddy::getCoarseningThreshold() const {
    return coarsening_threshold;
}

void SolverBuddy::setNumSweeps(int sweeps) {
    if (sweeps < 1)
        throw SolverOptionsException("number of sweeps must be positive.");
    this->sweeps = sweeps;
}

int SolverBuddy::getNumSweeps() const {
    return sweeps;
}

void SolverBuddy::setNumPreSweeps(int sweeps) {
    if (sweeps < 1)
        throw SolverOptionsException("number of sweeps must be positive.");
    pre_sweeps = sweeps;
}

int SolverBuddy::getNumPreSweeps() const {
    return pre_sweeps;
}

void SolverBuddy::setNumPostSweeps(int sweeps) {
    if (sweeps < 1)
       throw SolverOptionsException("number of sweeps must be positive.");
    post_sweeps = sweeps;
}

int SolverBuddy::getNumPostSweeps() const {
    return post_sweeps;
}

void SolverBuddy::setTolerance(double rtol) {
    if (rtol < 0. || rtol > 1.)
        throw SolverOptionsException("tolerance must be non-negative and less or equal 1.");
    tolerance = rtol;
}

double SolverBuddy::getTolerance() const {
    return tolerance;
}

void SolverBuddy::setAbsoluteTolerance(double atol) {
    if (atol < 0.)
       throw SolverOptionsException("tolerance must be non-negative.");
    absolute_tolerance = atol;
}

double SolverBuddy::getAbsoluteTolerance() const {
    return absolute_tolerance;
}

void SolverBuddy::setInnerTolerance(double rtol) {
    if (rtol <= 0. || rtol > 1.)
        throw SolverOptionsException("tolerance must be positive and less than or equal to 1.");
    inner_tolerance = rtol;
}

double SolverBuddy::getInnerTolerance() const {
    return inner_tolerance;
}

void SolverBuddy::setDropTolerance(double drop_tol) {
    if (drop_tol < 0. || drop_tol > 1.)
        throw SolverOptionsException("drop tolerance must be zero or positive and less than or equal to 1.");
    drop_tolerance = drop_tol;
}

double SolverBuddy::getDropTolerance() const {
    return drop_tolerance;
}

void SolverBuddy::setDropStorage(double storage) {
    if (storage < 1.)
        throw SolverOptionsException("allowed storage increase must be greater than or equal to 1.");
    drop_storage = storage;
}

double SolverBuddy::getDropStorage() const {
    return drop_storage;
}

void SolverBuddy::setRelaxationFactor(double factor) {
    if (factor < 0.)
        throw SolverOptionsException("relaxation factor must be non-negative.");
    relaxation = factor;
}

double SolverBuddy::getRelaxationFactor() const {
    return relaxation;
}

bool SolverBuddy::isSymmetric() const {
    return symmetric;
}

void SolverBuddy::setSymmetryOn() {
    symmetric = true;
}

void SolverBuddy::setSymmetryOff() {
    symmetric = false;
}

void SolverBuddy::setSymmetry(bool flag) {
    if (flag)
        setSymmetryOn();
    else
        setSymmetryOff();
}

bool SolverBuddy::isVerbose() const {
    return verbose;
}

void SolverBuddy::setVerbosityOn() {
    verbose = true;
}

void SolverBuddy::setVerbosityOff() {
    verbose = false;
}

void SolverBuddy::setVerbosity(bool verbose) {
    if (verbose)
        setVerbosityOn();
    else
        setVerbosityOff();
}

bool SolverBuddy::adaptInnerTolerance() const {
    return adapt_inner_tolerance;
}

void SolverBuddy::setInnerToleranceAdaptionOn() {
    adapt_inner_tolerance = true;
}

void SolverBuddy::setInnerToleranceAdaptionOff() {
    adapt_inner_tolerance = false;
}

void SolverBuddy::setInnerToleranceAdaption(bool adapt) {
    if (adapt)
        setInnerToleranceAdaptionOn();
    else
        setInnerToleranceAdaptionOff();
}

bool SolverBuddy::acceptConvergenceFailure() const {
    return accept_convergence_failure;
}

void SolverBuddy::setAcceptanceConvergenceFailureOn() {
    accept_convergence_failure = true;
}

void SolverBuddy::setAcceptanceConvergenceFailureOff() {
    accept_convergence_failure = false;
}

void SolverBuddy::setAcceptanceConvergenceFailure(bool accept) {
    if (accept)
        setAcceptanceConvergenceFailureOn();
    else
        setAcceptanceConvergenceFailureOff();
}

bool SolverBuddy::useLocalPreconditioner() const {
    return use_local_preconditioner;
}

void SolverBuddy::setLocalPreconditionerOn() {
    use_local_preconditioner = true;
}

void SolverBuddy::setLocalPreconditionerOff() {
    use_local_preconditioner=false;
}

void SolverBuddy::setLocalPreconditioner(bool use) {
    if (use)
        setLocalPreconditionerOn();
    else
        setLocalPreconditionerOff();
}

void SolverBuddy::setMinCoarseMatrixSparsity(double sparsity) {
    if (sparsity < 0.)
        throw SolverOptionsException("sparsity must be non-negative.");
    else if (sparsity > 1)
        throw SolverOptionsException("sparsity must be less than 1.");
    min_sparsity = sparsity;
}

double SolverBuddy::getMinCoarseMatrixSparsity() const {
    return min_sparsity;
}

void SolverBuddy::setNumRefinements(int refinements) {
    if (refinements < 0)
        throw SolverOptionsException("number of refinements must be non-negative.");
    this->refinements = refinements;
}

int SolverBuddy::getNumRefinements() const {
    return refinements;
}

void SolverBuddy::setNumCoarseMatrixRefinements(int refinements) {
    if (refinements < 0)
        throw SolverOptionsException("number of refinements must be non-negative.");
    coarse_refinements = refinements;
}

int SolverBuddy::getNumCoarseMatrixRefinements() const {
    return coarse_refinements;
}

bool SolverBuddy::usePanel() const {
    return use_panel;
}

void SolverBuddy::setUsePanelOn() {
    use_panel = true;
}

void SolverBuddy::setUsePanelOff() {
    use_panel = false;
}

void SolverBuddy::setUsePanel(bool use) {
    if (use)
        setUsePanelOn();
    else
        setUsePanelOff();
}

void SolverBuddy::setAMGInterpolation(int method) {
    SolverOptions meth = static_cast<SolverOptions>(method);
    switch (meth) {
        case ESCRIPT_CLASSIC_INTERPOLATION_WITH_FF_COUPLING:
        case ESCRIPT_CLASSIC_INTERPOLATION:
        case ESCRIPT_DIRECT_INTERPOLATION:
            amg_interpolation_method = meth;
            break;
        default:
            throw SolverOptionsException("unknown AMG interpolation method");
    }
}

SolverOptions SolverBuddy::getAMGInterpolation() const {
    return amg_interpolation_method;
}

void SolverBuddy::setODESolver(int method) {
    SolverOptions ode = static_cast<SolverOptions>(method);
    switch (ode) {
        case ESCRIPT_CRANK_NICOLSON:
        case ESCRIPT_BACKWARD_EULER:
        case ESCRIPT_LINEAR_CRANK_NICOLSON:
            ode_solver = ode;
            break;
        default:
            throw SolverOptionsException("unknown ODE solver method");
    }
}

SolverOptions SolverBuddy::getODESolver() const {
    return ode_solver;
}

}
