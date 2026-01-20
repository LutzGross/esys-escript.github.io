
/******************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESCRIPT_SOLVEROPTIONS_H__
#define __ESCRIPT_SOLVEROPTIONS_H__

#include <boost/python/dict.hpp>
#include <boost/python/object.hpp>
#include "system_dep.h"

namespace escript {

/**
This enum defines the options for solving linear/non-linear systems with escript.

SO_DEFAULT: use escript defaults for specific option
SO_TARGET_CPU: use CPUs to solve system
SO_TARGET_GPU: use GPUs to solve system

SO_PACKAGE_MKL: Intel's MKL solver library
SO_PACKAGE_PASO: PASO solver package
SO_PACKAGE_TRILINOS: The TRILINOS parallel solver class library from Sandia National Labs
SO_PACKAGE_UMFPACK: The UMFPACK library
SO_PACKAGE_MUMPS: The MUMPS library

SO_METHOD_BICGSTAB: The stabilized Bi-Conjugate Gradient method
SO_METHOD_CHOLEVSKY: The direct solver based on LDLT factorization (can only be applied for symmetric PDEs)
SO_METHOD_CGS: The conjugate gradient square method
SO_METHOD_CR: The conjugate residual method
SO_METHOD_DIRECT: A direct solver based on LDU factorization
SO_METHOD_DIRECT_MUMPS: MUMPS parallel direct solver
SO_METHOD_DIRECT_PARDISO: MKL Pardiso direct solver
SO_METHOD_DIRECT_SUPERLU: SuperLU direct solver
SO_METHOD_DIRECT_TRILINOS: Trilinos-based direct solver
SO_METHOD_GMRES: The Gram-Schmidt minimum residual method
SO_METHOD_HRZ_LUMPING: Matrix lumping using the HRZ approach
SO_METHOD_ITERATIVE: The default iterative solver
SO_METHOD_LSQR: Least squares with QR factorization
SO_METHOD_MINRES: Minimum residual method
SO_METHOD_PCG: The preconditioned conjugate gradient method (can only be applied for symmetric PDEs)
SO_METHOD_PRES20: Special GMRES with restart after 20 steps and truncation after 5 residuals
SO_METHOD_ROWSUM_LUMPING: Matrix lumping using row sum
SO_METHOD_TFQMR: Transpose Free Quasi Minimal Residual method

SO_PRECONDITIONER_AMG: Algebraic Multi Grid
SO_PRECONDITIONER_GAUSS_SEIDEL: Gauss-Seidel preconditioner
SO_PRECONDITIONER_ILU0: The incomplete LU factorization preconditioner with no fill-in
SO_PRECONDITIONER_ILUT: The incomplete LU factorization preconditioner with fill-in
SO_PRECONDITIONER_JACOBI: The Jacobi preconditioner
SO_PRECONDITIONER_NONE: no preconditioner is applied
SO_PRECONDITIONER_REC_ILU: recursive ILU0
SO_PRECONDITIONER_RILU: relaxed ILU0

SO_ODESOLVER_BACKWARD_EULER: backward Euler scheme
SO_ODESOLVER_CRANK_NICOLSON: Crank-Nicolson scheme
SO_ODESOLVER_LINEAR_CRANK_NICOLSON: linerized Crank-Nicolson scheme

SO_INTERPOLATION_CLASSIC: classical interpolation in AMG
SO_INTERPOLATION_CLASSIC_WITH_FF_COUPLING: classical interpolation in AMG with enforced
SO_INTERPOLATION_DIRECT: direct interpolation in AMG

SO_REORDERING_DEFAULT: the reordering method recommended by the solver
SO_REORDERING_MINIMUM_FILL_IN: Reorder matrix to reduce fill-in during factorization
SO_REORDERING_NESTED_DISSECTION: Reorder matrix to improve load balancing during factorization
SO_REORDERING_NONE: No matrix reordering allowed
*/
enum SolverOptions
{
    SO_DEFAULT,

    // Solver targets
    SO_TARGET_CPU,
    SO_TARGET_GPU,

    // Solver packages
    SO_PACKAGE_MKL,
    SO_PACKAGE_PASO,
    SO_PACKAGE_TRILINOS,
    SO_PACKAGE_UMFPACK,
    SO_PACKAGE_MUMPS,

    // Solver methods
    SO_METHOD_BICGSTAB,
    SO_METHOD_CGLS,
    SO_METHOD_CGS,
    SO_METHOD_CHOLEVSKY,
    SO_METHOD_CR,
    SO_METHOD_DIRECT,
    SO_METHOD_DIRECT_MUMPS,
    SO_METHOD_DIRECT_PARDISO,
    SO_METHOD_DIRECT_SUPERLU,
    SO_METHOD_DIRECT_TRILINOS,
    SO_METHOD_GMRES,
    SO_METHOD_HRZ_LUMPING,
    SO_METHOD_ITERATIVE,
    SO_METHOD_LSQR,
    SO_METHOD_MINRES,
    SO_METHOD_NONLINEAR_GMRES,
    SO_METHOD_PCG,
    SO_METHOD_PRES20,
    SO_METHOD_ROWSUM_LUMPING,
    SO_METHOD_TFQMR,

    // Preconditioners
    SO_PRECONDITIONER_AMG,
    SO_PRECONDITIONER_GAUSS_SEIDEL,
    SO_PRECONDITIONER_ILU0,
    SO_PRECONDITIONER_ILUT,
    SO_PRECONDITIONER_JACOBI,
    SO_PRECONDITIONER_NONE,
    SO_PRECONDITIONER_REC_ILU,
    SO_PRECONDITIONER_RILU,

    // ODE solvers
    SO_ODESOLVER_BACKWARD_EULER,
    SO_ODESOLVER_CRANK_NICOLSON,
    SO_ODESOLVER_LINEAR_CRANK_NICOLSON,

    // Interpolation methods
    SO_INTERPOLATION_CLASSIC,
    SO_INTERPOLATION_CLASSIC_WITH_FF_COUPLING,
    SO_INTERPOLATION_DIRECT,

    SO_REORDERING_DEFAULT,
    SO_REORDERING_MINIMUM_FILL_IN,
    SO_REORDERING_NESTED_DISSECTION,
    SO_REORDERING_NONE
};

/// returns true if the passed solver method refers to a direct solver type
inline bool isDirectSolver(const SolverOptions& method)
{
    switch (method) {
        case SO_METHOD_DIRECT:
        case SO_METHOD_DIRECT_MUMPS:
        case SO_METHOD_DIRECT_PARDISO:
        case SO_METHOD_DIRECT_SUPERLU:
        case SO_METHOD_DIRECT_TRILINOS:
            return true;
        default:
            break;
    }
    return false;
}

class ESCRIPT_DLL_API SolverBuddy
{
public:
    SolverBuddy();
    ~SolverBuddy();

    /**
        Returns a string reporting the current settings
    */
    std::string getSummary() const;

    /**
        Returns the name of a given key

        \param key a valid key from SolverOptions
    */
    const char* getName(int key) const;

    /**
        Resets the diagnostics

        \param all if ``all`` is ``true`` all diagnostics including accumulative
                counters are reset.
    */
    void resetDiagnostics(bool all=false);

    /**
        Updates diagnostic information

        \param key name of diagnostic (a python string in the list "num_iter",
                 "num_level", "num_inner_iter", "time", "set_up_time",
                 "net_time", "residual_norm", "converged").
        \param value new value of the diagnostic information
    */
    void updateDiagnosticsPy(const std::string& key,
                             const boost::python::object& value);

    void updateDiagnostics(const std::string& key, bool value);
    void updateDiagnostics(const std::string& key, int value);
    void updateDiagnostics(const std::string& key, double value);

    /**
        Returns the diagnostic information for the given ``name``.
        Possible values are:

        - "num_iter": the number of iteration steps
        - "cum_num_iter": the cumulative number of iteration steps
        - "num_level": the number of level in multi level solver
        - "num_inner_iter": the number of inner iteration steps
        - "cum_num_inner_iter": the cumulative number of inner iteration steps
        - "time": execution time
        - "cum_time": cumulative execution time
        - "set_up_time": time to set up the solver, typically this includes factorization and reordering
        - "cum_set_up_time": cumulative time to set up the solver
        - "net_time": net execution time, excluding setup time for the solver and execution time for preconditioner
        - "cum_net_time": cumulative net execution time
        - "preconditioner_size": size of preconditioner [Bytes]
        - "converged": true if solution has converged
        - "time_step_backtracking_used": true if time step back tracking has been used
        - "coarse_level_sparsity": sparsity of the matrix on the coarsest level
        - "num_coarse_unknowns": number of unknowns on the coarsest level

        \param name name of diagnostic information to return

        \return requested value or 0 if the value is yet to be defined.

        \note If the solver has thrown an exception diagnostic values have an
              undefined status.
    */
    double getDiagnostics(const std::string name) const;

    /**
        Returns ``true`` if the last solver call has been finalized
        successfully.

        \note if an exception has been thrown by the solver the status of this
                 flag is undefined.
    */
    bool hasConverged() const;

    /**
        Sets the preconditioner to be used.

        \param preconditioner key of the preconditioner to be used, one of
            `SO_PRECONDITIONER_ILU0`, `SO_PRECONDITIONER_ILUT`,
            `SO_PRECONDITIONER_JACOBI`, `SO_PRECONDITIONER_AMG`,
            `SO_PRECONDITIONER_AMLI`, `SO_PRECONDITIONER_REC_ILU`,
            `SO_PRECONDITIONER_GAUSS_SEIDEL`, `SO_PRECONDITIONER_RILU`,
            `SO_PRECONDITIONER_NONE`

        \note Not all packages support all preconditioners. It can be assumed
              that a package makes a reasonable choice if it encounters an
              unknown preconditioner.
    */
    void setPreconditioner(int preconditioner);

    /**
        Returns the key of the preconditioner to be used.
    */
    SolverOptions getPreconditioner() const;

    /**
        Sets the solver method to be used. Use ``method``=``SO_METHOD_DIRECT``
        to indicate that a direct rather than an iterative solver should be
        used and use ``method``=``SO_METHOD_ITERATIVE`` to indicate that an
        iterative rather than a direct solver should be used.

        \param method key of the solver method to be used, should be in
            `SO_DEFAULT`, `SO_METHOD_DIRECT`, `SO_METHOD_DIRECT_MUMPS`,
            `SO_METHOD_DIRECT_PARDISO`, `SO_METHOD_DIRECT_SUPERLU`,
            `SO_METHOD_DIRECT_TRILINOS`, `SO_METHOD_CHOLEVSKY`,
            `SO_METHOD_PCG`, `SO_METHOD_CR`, `SO_METHOD_CGS`,
            `SO_METHOD_BICGSTAB`, `SO_METHOD_GMRES`, `SO_METHOD_PRES20`,
            `SO_METHOD_ROWSUM_LUMPING`, `SO_METHOD_HRZ_LUMPING`,
            `SO_METHOD_ITERATIVE`, `SO_METHOD_LSQR`,
            `SO_METHOD_NONLINEAR_GMRES`, `SO_METHOD_TFQMR`, `SO_METHOD_MINRES`

        \note Not all packages support all solvers. It can be assumed that a
              package makes a reasonable choice if it encounters an unknown
              solver method.
    */
    void setSolverMethod(int method);

    /**
        Returns key of the solver method to be used.
    */
    SolverOptions getSolverMethod() const;

    /**
        Sets the solver package to be used as a solver.

        \param package key of the solver package to be used, should be in
               `SO_DEFAULT`, `SO_PACKAGE_CUSP`, `SO_PACKAGE_PASO`,
               `SO_PACKAGE_MKL`, `SO_PACKAGE_UMFPACK`,
               `SO_PACKAGE_TRILINOS`, `SO_PACKAGE_MUMPS`

        \note Not all packages are supported on all implementation.
              An exception may be thrown on some platforms if the selected
              package is unsupported.
    */
    void setPackage(int package);

    /**
        Returns the solver package key
    */
    SolverOptions getPackage() const;

    /**
        Sets the key of the reordering method to be applied if supported by the
        solver. Some direct solvers support reordering
        to optimize compute time and storage use during elimination.

        \param ordering selects the reordering strategy, should be in
               `SO_REORDERING_NONE`, `SO_REORDERING_MINIMUM_FILL_IN`,
               `SO_REORDERING_NESTED_DISSECTION`, 'SO_REORDERING_DEFAULT`
    */
    void setReordering(int ordering);

    /**
        Returns the key of the reordering method to be applied if supported by
        the solver.
    */
    SolverOptions getReordering() const;

    /**
        Sets the number of iterations steps after which GMRES performs a
        restart.

        \param restart number of iteration steps after which to perform a
               restart. If 0 no restart is performed.
    */
    void setRestart(int restart);

    /**
        Returns the number of iterations steps after which GMRES performs a
        restart. 0 means no restart is performed.
    */
    int getRestart() const;

    /**
        Returns the number of iterations steps after which GMRES performs a
        restart. If -1 is returned no restart is performed.
    */
    int _getRestartForC() const;

    /**
        Sets the number of residuals in GMRES to be stored for
        orthogonalization. The more residuals are stored the faster GMRES
        converges but more memory is required.
    */
    void setTruncation(int truncation);

    /**
        Returns the number of residuals in GMRES to be stored for
        orthogonalization.
    */
    int getTruncation() const;

    /**
        Sets the maximum number of iteration steps for the inner iteration.

        \param iter_max maximum number of inner iterations
    */
    void setInnerIterMax(int iter_max);

    /**
        Returns maximum number of inner iteration steps
    */
    int getInnerIterMax() const;

    /**
        Sets the maximum number of iteration steps

        \param iter_max maximum number of iteration steps
    */
    void setIterMax(int iter_max);

    /**
        Returns maximum number of iteration steps
    */
    int getIterMax() const;

    /**
        Sets the number of sweeps in a Jacobi or Gauss-Seidel/SOR
        preconditioner.

        \param sweeps number of sweeps
    */
    void setNumSweeps(int sweeps);

    /**
        Returns the number of sweeps in a Jacobi or Gauss-Seidel/SOR
        preconditioner.
    */
    int getNumSweeps() const;

    /**
        Sets the relative tolerance for the solver

        \param rtol relative tolerance
    */
    void setTolerance(double rtol);

    /**
        Returns the relative tolerance for the solver
    */
    double getTolerance() const;

    /**
        Sets the absolute tolerance for the solver

        \param atol absolute tolerance
    */
    void setAbsoluteTolerance(double atol);

    /**
        Returns the absolute tolerance for the solver
    */
    double getAbsoluteTolerance() const;

    /**
        Sets the relative tolerance for an inner iteration scheme for instance
        on the coarsest level in a multi-level scheme.

        \param rtol inner relative tolerance
    */
    void setInnerTolerance(double rtol);

    /**
        Returns the relative tolerance for an inner iteration scheme
    */
    double getInnerTolerance() const;

    /**
        Sets the relative drop tolerance in ILUT

        \param drop_tol drop tolerance
    */
    void setDropTolerance(double drop_tol);

    /**
        Returns the relative drop tolerance in ILUT
    */
    double getDropTolerance() const;

    /**
        Sets the maximum allowed increase in storage for ILUT. An increase of
        2 would mean that a doubling of the storage needed for the coefficient
        matrix is allowed during ILUT factorization.

        \param drop allowed storage increase
    */
    void setDropStorage(double drop);

    /**
        Returns the maximum allowed increase in storage for ILUT
    */
    double getDropStorage() const;

    /**
        Sets the relaxation factor used to add dropped elements in RILU to the
        main diagonal.

        \param factor relaxation factor
        \note RILU with a relaxation factor 0 is identical to ILU0
    */
    void setRelaxationFactor(double factor);

    /**
        Returns the relaxation factor used to add dropped elements in RILU to
        the main diagonal.
    */
    double getRelaxationFactor() const;

    /**
        Checks if the coefficient matrix is set to be complex-valued.

        \return true if a complex-valued PDE is indicated, false otherwise
    */
    bool isComplex() const;

    /**
        Sets the complex flag for the coefficient matrix to ``flag``.

        \param complex If true, the complex flag is set otherwise reset.
    */
    void setComplex(bool complex);

    /**
        Checks if symmetry of the coefficient matrix is indicated.

        \return true if a symmetric PDE is indicated, false otherwise
    */
    bool isSymmetric() const;

    /**
        Sets the symmetry flag to indicate that the coefficient matrix is
        symmetric.
    */
    void setSymmetryOn();

    /**
        Clears the symmetry flag for the coefficient matrix.
    */
    void setSymmetryOff();

    /**
        Sets the symmetry flag for the coefficient matrix to ``flag``.

        \param symmetry If true, the symmetry flag is set otherwise reset.
    */
    void setSymmetry(bool symmetry);

    /**
        Checks if the coefficient matrix is indicated to be Hermitian.

        \return true if a Hermitian PDE is indicated, false otherwise
    */
    bool isHermitian() const;

    /**
        Sets the symmetry flag to indicate that the coefficient matrix is
        Hermitian.
    */
    void setHermitianOn();

    /**
        Clears the Hermitian flag for the coefficient matrix.
    */
    void setHermitianOff();

    /**
        Sets the Hermitian flag for the coefficient matrix to ``flag``.

        \param Hermitian If true, the symmetry flag is set otherwise reset.
    */
    void setHermitian(bool hermitian);

    /**
        Returns ``true`` if the solver is expected to be verbose.

        \return true if verbosity is on
    */
    bool isVerbose() const;

    /**
        Switches the verbosity of the solver on.
    */
    void setVerbosityOn();

    /**
        Switches the verbosity of the solver off.
    */
    void setVerbosityOff();

    /**
        Sets the verbosity flag for the solver to ``flag``.

        \param verbose If ``true``, the verbosity of the solver is switched on.
    */
    void setVerbosity(bool verbose);

    /**
        Returns ``true`` if the tolerance of the inner solver is selected
        automatically. Otherwise the inner tolerance set by `setInnerTolerance`
        is used.

        \returns ``true`` if inner tolerance adaption is chosen.
    */
    bool adaptInnerTolerance() const;

    /**
        Switches the automatic selection of inner tolerance on
    */
    void setInnerToleranceAdaptionOn();

    /**
        Switches the automatic selection of inner tolerance off
    */
    void setInnerToleranceAdaptionOff();

    /**
        Sets the flag to indicate automatic selection of the inner tolerance.

        \param adaption If ``true``, the inner tolerance is selected automatically
    */
    void setInnerToleranceAdaption(bool adaption);

    /**
        Returns ``true`` if a failure to meet the stopping criteria within the
        given number of iteration steps does not raise an exception.
        This is useful if a solver is used in a non-linear context where the
        non-linear solver can continue even if the returned solution does not
        necessarily meet the stopping criteria. One can use the `hasConverged`
        method to check if the last call to the solver was successful.

        \returns ``true`` if a failure to achieve convergence is accepted.
    */
    bool acceptConvergenceFailure() const;

    /**
        Switches the acceptance of a failure of convergence on
    */
    void setAcceptanceConvergenceFailureOn();

    /**
        Switches the acceptance of a failure of convergence off
    */
    void setAcceptanceConvergenceFailureOff();

    /**
        Sets the flag to indicate the acceptance of a failure of convergence.

        \param acceptance If ``true``, any failure to achieve convergence is
               accepted.
    */
    void setAcceptanceConvergenceFailure(bool acceptance);

    /**
        Returns ``true`` if the preconditoner is applied locally on each MPI
        rank. This reduces communication costs and speeds up the application of
        the preconditioner but at the cost of more iteration steps.
        This can be an advantage on clusters with slower interconnects.
    */
    bool useLocalPreconditioner() const;

    /**
        Sets the flag to use local preconditioning to on
    */
    void setLocalPreconditionerOn();

    /**
        Sets the flag to use local preconditioning to off
    */
    void setLocalPreconditionerOff();

    /**
        Sets the flag to use local preconditioning

        \param local If ``true``, local proconditioning on each MPI rank is
               applied
    */
    void setLocalPreconditioner(bool local);

    /**
        Sets the number of refinement steps to refine the solution when a
        direct solver is applied.

        \param refinements number of refinements
    */
    void setNumRefinements(int refinements);

    /**
        Returns the number of refinement steps to refine the solution when a
        direct solver is applied.
    */
    int getNumRefinements() const;

    /**
        Sets the solver method for ODEs.

        \param solver key of the ODE solver method to be used, should be in
               `SO_ODESOLVER_CRANK_NICOLSON`, `SO_ODESOLVER_BACKWARD_EULER`,
               `SO_ODESOLVER_LINEAR_CRANK_NICOLSON`
    */
    void setODESolver(int solver);

    /**
        Returns the key of the solver method for ODEs.
    */
    SolverOptions getODESolver() const;

    /**
        Sets a Trilinos preconditioner/solver parameter.
        \note Escript does not check for validity of the parameter name
        (e.g. spelling mistakes). Parameters are passed 1:1 to escript's
        Trilinos wrapper and from there to the relevant Trilinos package.
        See the relevant Trilinos documentation for valid parameter strings
        and values.
        \note This method does nothing in a non-Trilinos build.
    */
    void setTrilinosParameter(const std::string& name,
                              const boost::python::object& value);

    /**
        Returns a boost python dictionary of set Trilinos parameters.
        \note This method returns an empty dictionary in a non-Trilinos build.
    */
    boost::python::dict getTrilinosParameters() const;

    /**
        Sets the dimension of the problem we are solving. 
    */
    void setDim(int dim);

    /**
        Returns the dimension of the problem we are solving. 
    */
    int getDim();

    /**
        True if we are using the default solver settings
    */
    bool using_default_method() const;

    /**
        Sets the parameter Oxley_Domain
    */
    void setOxleyDomain(bool);

    /**
        True if we are using an Oxley domain, False otherwise;
    */
    bool getOxleyDomain();    


protected:
    boost::python::dict trilinosParams;

    SolverOptions target;
    SolverOptions package;
    SolverOptions method;
    SolverOptions preconditioner;
    SolverOptions ode_solver;
    SolverOptions reordering;
    int sweeps;
    double tolerance;
    double absolute_tolerance;
    double inner_tolerance;
    double drop_tolerance;
    double drop_storage;
    int iter_max;
    int inner_iter_max;
    int truncation;
    int restart; //0 will have to be None in python, will get tricky
    bool is_complex;
    bool symmetric;
    bool hermitian;
    bool verbose;
    bool adapt_inner_tolerance;
    bool accept_convergence_failure;
    double relaxation;
    bool use_local_preconditioner;
    int refinements;
    int dim; // Dimension of the problem, either 2 or 3. Used internally

    int num_iter;
    int num_level;
    int num_inner_iter;
    double time;
    double set_up_time;
    double net_time;
    double residual_norm;
    bool converged;
    int preconditioner_size;
    bool time_step_backtracking_used;
    double coarse_level_sparsity;
    int num_coarse_unknowns;
    int cum_num_inner_iter;
    int cum_num_iter;
    double cum_time;
    double cum_set_up_time;
    double cum_net_time;

    bool using_default_solver_method;

    bool have_oxley;
};

typedef boost::shared_ptr<SolverBuddy> SB_ptr;

} // namespace escript

#endif // __ESCRIPT_SOLVEROPTIONS_H__

