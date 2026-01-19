
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

#include "SolverOptions.h"
#include "EsysException.h"

#include <boost/python.hpp>

#ifdef ESYS_HAVE_TRILINOS
#include <Amesos2.hpp>
#endif

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
    reordering(SO_REORDERING_DEFAULT),
    sweeps(1),
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
    hermitian(false),
    verbose(false),
    adapt_inner_tolerance(true),
    accept_convergence_failure(false),
    relaxation(0.3),
    use_local_preconditioner(false),
    refinements(2),
    dim(2),
    using_default_solver_method(false),
    have_oxley(false)
{
    // setPackage(SO_DEFAULT);
    // setSolverMethod(SO_DEFAULT);
    resetDiagnostics(true);
}

SolverBuddy::~SolverBuddy()
{
}


std::string SolverBuddy::getSummary() const
{
    std::stringstream out;
    out << "Solver Package = " << getName(getPackage()) << std::endl
        << "Verbosity = " << isVerbose() << std::endl
        << "Accept failed convergence = " << acceptConvergenceFailure() << std::endl
        << "Relative tolerance = " << getTolerance() << std::endl
        << "Absolute tolerance = " << getAbsoluteTolerance() << std::endl
        << "Symmetric problem = " << isSymmetric() << std::endl
        << "Hermitian problem = " << isHermitian() << std::endl
        << "Maximum number of iteration steps = " << getIterMax() << std::endl
        << "Inner tolerance = " << getInnerTolerance() << std::endl
        << "Adapt innner tolerance = " << adaptInnerTolerance() << std::endl;

    if (getPackage() == SO_DEFAULT || getPackage() == SO_PACKAGE_PASO ||
            getPackage() == SO_PACKAGE_TRILINOS) {
        out << "Solver method = " << getName(getSolverMethod()) << std::endl;
        if (getSolverMethod() == SO_METHOD_GMRES) {
            out << "Truncation  = " << getTruncation() << std::endl
                << "Restart  = " << getRestart() << std::endl;
        }
        out << "Preconditioner = " << getName(getPreconditioner()) << std::endl
            << "Apply preconditioner locally = " << useLocalPreconditioner()
            << std::endl;
        switch (getPreconditioner()) {
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

        case SO_PACKAGE_MKL: return "MKL";
        case SO_PACKAGE_PASO: return "PASO";
        case SO_PACKAGE_TRILINOS: return "TRILINOS";
        case SO_PACKAGE_UMFPACK: return "UMFPACK";
        case SO_PACKAGE_MUMPS: return "MUMPS";

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
    throw ValueError(std::string("unknown diagnostic item: ") + name);
}

bool SolverBuddy::hasConverged() const
{
    return converged;
}

void SolverBuddy::setPreconditioner(int precon)
{
    SolverOptions preconditioner = static_cast<SolverOptions>(precon);
    switch(preconditioner) {
        case SO_PRECONDITIONER_AMG:
#if !defined(ESYS_HAVE_TRILINOS) && !defined(ESYS_HAVE_MUMPS)
        throw ValueError("escript was not compiled with Trilinos or MUMPS enabled");
#endif
        case SO_PRECONDITIONER_GAUSS_SEIDEL:
        case SO_PRECONDITIONER_JACOBI: // This is the default preconditioner in ifpack2
        case SO_PRECONDITIONER_ILU0:
        case SO_PRECONDITIONER_ILUT:
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

void SolverBuddy::setSolverMethod(int method)
{
    SolverOptions meth = static_cast<SolverOptions>(method);

//     bool havePASODirect = false;
//     using_default_solver_method=false;
// #if defined(ESYS_HAVE_PASO) && (defined(ESYS_HAVE_MKL) || defined(ESYS_HAVE_UMFPACK) || defined(ESYS_HAVE_MUMPS))
//     havePASODirect = true;
// #endif

    switch(meth) {
        case SO_DEFAULT:
        case SO_METHOD_ITERATIVE:
        case SO_METHOD_BICGSTAB:
        case SO_METHOD_CGLS:
        case SO_METHOD_CGS:
        case SO_METHOD_CHOLEVSKY:
        case SO_METHOD_CR:
        case SO_METHOD_GMRES:
        case SO_METHOD_HRZ_LUMPING:
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
#if defined(ESYS_HAVE_UMFPACK) || defined(ESYS_HAVE_TRILINOS) || defined(ESYS_HAVE_MKL) || defined(ESYS_HAVE_MUMPS)
#ifdef ESYS_HAVE_TRILINOS
            // translate specific direct solver setting to generic one for PASO
            this->method = meth;
#else
            this->method = SO_METHOD_DIRECT;
#endif
            break;
#else
            throw ValueError("Cannot use DIRECT solver method, the running "
                    "escript was not compiled with a direct solver enabled");
#endif
        case SO_METHOD_DIRECT_TRILINOS:
#ifdef ESYS_HAVE_TRILINOS
            this->method = meth;
            break;
#else
            throw ValueError("escript was not compiled with Trilinos");
#endif
        case SO_METHOD_DIRECT_MUMPS:
#ifdef ESYS_HAVE_MUMPS
            this->method=meth;
#else
            throw ValueError("escript was not compiled with MUMPS");
#endif
        case SO_METHOD_DIRECT_PARDISO:
#ifdef ESYS_HAVE_TRILINOS
            if(Amesos2::query("pardiso_mkl"))
                this->method = meth;
            else
                throw ValueError("Trilinos was not compiled with MKL Pardiso");
#else
            throw ValueError("escript was not compiled with Trilinos");
#endif
        case SO_METHOD_DIRECT_SUPERLU:
#ifdef ESYS_HAVE_TRILINOS
            if(Amesos2::query("superludist") || Amesos2::query("superlu") || Amesos2::query("superlumt"))
                this->method = meth;
            else
                throw ValueError("Trilinos was not compiled with SuperLU ");
#else
            throw ValueError("escript was not compiled with Trilinos");
#endif
        default:
            throw ValueError("unknown solver method");
    }
}

SolverOptions SolverBuddy::getSolverMethod() const
{
    return method;
}

void SolverBuddy::setPackage(int package)
{
    SolverOptions pack = static_cast<SolverOptions>(package);
    switch (pack) {
        case SO_DEFAULT:
        // Default to trilinos if it is available, else use PASO
        // This should always work because escript cannot be compiled
        // unless at least one of these is available
#ifdef ESYS_HAVE_TRILINOS
            this->package = SO_PACKAGE_TRILINOS;
            setSolverMethod(getSolverMethod());
            break;
#else
            this->package = SO_PACKAGE_PASO;
            setSolverMethod(getSolverMethod());
            break;
#endif
        // Set to PASO iff escript was compiled with PASO
        case SO_PACKAGE_PASO:
#ifdef ESYS_HAVE_PASO
            this->package = SO_PACKAGE_PASO;
            setSolverMethod(getSolverMethod());
            break;
#else
            throw ValueError("escript was not compiled with PASO enabled");
#endif
        // Set to Trilinos iff escript was compiled with Trilinos
        case SO_PACKAGE_TRILINOS:
#ifdef ESYS_HAVE_TRILINOS
            this->package = SO_PACKAGE_TRILINOS;
            setSolverMethod(getSolverMethod());
            break;
#else
            throw ValueError("escript was not compiled with Trilinos enabled");
#endif
        // Set to MKL iff escript was compiled with MKL
        case SO_PACKAGE_MKL:
#ifdef ESYS_HAVE_MKL
            this->package = SO_PACKAGE_MKL;
            setSolverMethod(getSolverMethod());
            break;
#else
            throw ValueError("escript was not compiled with MKL enabled");
#endif
        // Set to Umfpack iff escript was compiled with Umfpack
        case SO_PACKAGE_UMFPACK:
#ifdef ESYS_HAVE_UMFPACK
            this->package = SO_PACKAGE_UMFPACK;
            setSolverMethod(getSolverMethod());
            break;
#else
            throw ValueError("escript was not compiled with UMFPACK enabled");
#endif
        // Set to MUMPS iff escript was compiled with MUMPS
        case SO_PACKAGE_MUMPS:
#ifdef ESYS_HAVE_MUMPS
            this->package = SO_PACKAGE_MUMPS;
            setSolverMethod(getSolverMethod());
            break;
#else
            throw ValueError("escript was not compiled with MUMPS enabled");
#endif
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

bool SolverBuddy::isHermitian() const
{
    return hermitian;
}

void SolverBuddy::setHermitianOn()
{
    hermitian = true;
}

void SolverBuddy::setHermitianOff()
{
    hermitian = false;
}

void SolverBuddy::setHermitian(bool flag)
{
    if (flag)
        setHermitianOn();
    else
        setHermitianOff();
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

void SolverBuddy::setDim(int dim)
{
    if (dim != 2 && dim != 3)
        throw ValueError("Dimension must be either 2 or 3.");
    this->dim = dim;
}

int SolverBuddy::getDim()
{
    return dim;
}

bool SolverBuddy::using_default_method() const
{
    return using_default_solver_method;
}

void SolverBuddy::setOxleyDomain(bool using_oxley)
{
    have_oxley=using_oxley;
}

bool SolverBuddy::getOxleyDomain()
{
    return have_oxley;
}

} // namespace escript
