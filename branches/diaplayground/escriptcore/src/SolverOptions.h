
/******************************************************************************
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

#ifndef SOLVER_OPTIONS_H
#define SOLVER_OPTIONS_H

#include <boost/python/object.hpp>
#include "system_dep.h"
namespace escript {

/**
This enum defines the solver options for a linear or non-linear solver.

ESCRIPT_DEFAULT: The default method used to solve the system of linear equations
ESCRIPT_DIRECT: The direct solver based on LDU factorization
ESCRIPT_CHOLEVSKY: The direct solver based on LDLT factorization (can only be applied for symmetric PDEs)
ESCRIPT_PCG: The preconditioned conjugate gradient method (can only be applied for symmetric PDEs)
ESCRIPT_CR: The conjugate residual method
ESCRIPT_CGS: The conjugate gradient square method
ESCRIPT_BICGSTAB: The stabilized Bi-Conjugate Gradient method
ESCRIPT_LSQR: Least Squares with QR factorization
ESCRIPT_TFQMR: Transpose Free Quasi Minimal Residual method
ESCRIPT_MINRES: Minimum residual method
ESCRIPT_ILU0: The incomplete LU factorization preconditioner with no fill-in
ESCRIPT_ILUT: The incomplete LU factorization preconditioner with fill-in
ESCRIPT_JACOBI: The Jacobi preconditioner
ESCRIPT_GMRES: The Gram-Schmidt minimum residual method
ESCRIPT_PRES20: Special ESCRIPT_GMRES with restart after 20 steps and truncation after 5 residuals
ESCRIPT_ROWSUM_LUMPING: Matrix lumping using row sum
ESCRIPT_HRZ_LUMPING: Matrix lumping using the HRZ approach
ESCRIPT_NO_REORDERING: No matrix reordering allowed
ESCRIPT_MINIMUM_FILL_IN: Reorder matrix to reduce fill-in during factorization
ESCRIPT_NESTED_DISSECTION: Reorder matrix to improve load balancing during factorization
ESCRIPT_PASO: PASO solver package
ESCRIPT_SCSL: SGI SCSL solver library
ESCRIPT_MKL: Intel's MKL solver library
ESCRIPT_UMFPACK: The UMFPACK library
ESCRIPT_TRILINOS: The TRILINOS parallel solver class library from Sandia National Labs
ESCRIPT_ITERATIVE: The default iterative solver
ESCRIPT_AMG: Algebraic Multi Grid
ESCRIPT_AMLI: Algebraic Multi Level Iteration
ESCRIPT_REC_ILU: recursive ESCRIPT_ILU0
ESCRIPT_RILU: relaxed ESCRIPT_ILU0
ESCRIPT_GAUSS_SEIDEL: Gauss-Seidel preconditioner
ESCRIPT_DEFAULT_REORDERING: the reordering method recommended by the solver
ESCRIPT_SUPER_LU: the Super_LU solver package
ESCRIPT_PASTIX: the Pastix direct solver_package
ESCRIPT_YAIR_SHAPIRA_COARSENING: ESCRIPT_AMG and ESCRIPT_AMLI coarsening method by Yair-Shapira
ESCRIPT_RUGE_STUEBEN_COARSENING: ESCRIPT_AMG and ESCRIPT_AMLI coarsening method by Ruge and Stueben
ESCRIPT_AGGREGATION_COARSENING: ESCRIPT_AMG and ESCRIPT_AMLI coarsening using (symmetric) aggregation
ESCRIPT_STANDARD_COARSENING: ESCRIPT_AMG and ESCRIPT_AMLI standard coarsening using measure of importance of the unknowns
ESCRIPT_MIN_COARSE_MATRIX_SIZE: minimum size of the coarsest level matrix to use direct solver.
ESCRIPT_NO_PRECONDITIONER: no preconditioner is applied.
ESCRIPT_CLASSIC_INTERPOLATION_WITH_FF_COUPLING: classical interpolation in AMG with enforced
ESCRIPT_CLASSIC_INTERPOLATION: classical interpolation in AMG
ESCRIPT_DIRECT_INTERPOLATION: direct interploation in AMG
ESCRIPT_BOOMERAMG: Boomer AMG in hypre library
ESCRIPT_CIJP_FIXED_RANDOM_COARSENING: BoomerAMG parallel coarsening method CIJP by using fixed random vector
ESCRIPT_CIJP_COARSENING: BoomerAMG parallel coarsening method CIJP
ESCRIPT_PAESCRIPT_FALGOUT_COARSENING: BoomerAMG parallel coarsening method falgout
ESCRIPT_PMIS_COARSENING: BoomerAMG parallel coarsening method PMIS
ESCRIPT_HMIS_COARSENING: BoomerAMG parallel coarsening method HMIS
ESCRIPT_BACKWARD_EULER: backward Euler scheme
ESCRIPT_CRANK_NICOLSON: Crank-Nicolson scheme
ESCRIPT_LINEAR_CRANK_NICOLSON: linerized Crank-Nicolson scheme

*/
enum SolverOptions {
    ESCRIPT_DEFAULT,
    ESCRIPT_DIRECT = 1,
    ESCRIPT_CHOLEVSKY = 2,
    ESCRIPT_PCG = 3,
    ESCRIPT_CR = 4,
    ESCRIPT_CGS = 5,
    ESCRIPT_BICGSTAB = 6,
    ESCRIPT_LSQR = 7,
    ESCRIPT_ILU0 = 8,
    ESCRIPT_ILUT = 9,
    ESCRIPT_JACOBI = 10,
    ESCRIPT_GMRES = 11,
    ESCRIPT_PRES20 = 12,
    ESCRIPT_LUMPING = 13,
    ESCRIPT_ROWSUM_LUMPING = 13,
    ESCRIPT_HRZ_LUMPING = 14,
    ESCRIPT_NO_REORDERING = 17,
    ESCRIPT_MINIMUM_FILL_IN = 18,
    ESCRIPT_NESTED_DISSECTION = 19,
    ESCRIPT_MKL = 15,
    ESCRIPT_UMFPACK = 16,
    ESCRIPT_ITERATIVE = 20,
    ESCRIPT_PASO = 21,
    ESCRIPT_AMG = 22,
    ESCRIPT_REC_ILU = 23,
    ESCRIPT_TRILINOS = 24,
    ESCRIPT_NONLINEAR_GMRES = 25,
    ESCRIPT_TFQMR = 26,
    ESCRIPT_MINRES = 27,
    ESCRIPT_GAUSS_SEIDEL = 28,
    ESCRIPT_RILU = 29,
    ESCRIPT_DEFAULT_REORDERING = 30,
    ESCRIPT_SUPER_LU = 31,
    ESCRIPT_PASTIX = 32,
    ESCRIPT_YAIR_SHAPIRA_COARSENING = 33,
    ESCRIPT_RUGE_STUEBEN_COARSENING = 34,
    ESCRIPT_AGGREGATION_COARSENING = 35,
    ESCRIPT_NO_PRECONDITIONER = 36,
    ESCRIPT_MIN_COARSE_MATRIX_SIZE = 37,
    ESCRIPT_AMLI = 38,
    ESCRIPT_STANDARD_COARSENING = 39,
    ESCRIPT_CLASSIC_INTERPOLATION_WITH_FF_COUPLING = 50,
    ESCRIPT_CLASSIC_INTERPOLATION = 51,
    ESCRIPT_DIRECT_INTERPOLATION = 52,
    ESCRIPT_BOOMERAMG = 60,
    ESCRIPT_CIJP_FIXED_RANDOM_COARSENING = 61,
    ESCRIPT_CIJP_COARSENING = 62,
    ESCRIPT_FALGOUT_COARSENING = 63,
    ESCRIPT_PMIS_COARSENING = 64,
    ESCRIPT_HMIS_COARSENING = 65,
    ESCRIPT_LINEAR_CRANK_NICOLSON = 66,
    ESCRIPT_CRANK_NICOLSON = 67,
    ESCRIPT_BACKWARD_EULER = 68,
    ESCRIPT_TARGET_CPU = 69,
    ESCRIPT_TARGET_GPU = 70
};

class SolverBuddy {
public:
    ESCRIPT_DLL_API SolverBuddy();
    /**
        Returns a string reporting the current settings
    */
    ESCRIPT_DLL_API std::string getSummary() const;

    /**
        Returns the name of a given key

        \param key a valid key from SolverOptions
    */
    ESCRIPT_DLL_API const char *getName(int key) const;

    /**
        Resets the diagnostics

        \param all if ``all`` is ``true`` all diagnostics including accumulative
                counters are reset.
    */
    ESCRIPT_DLL_API void resetDiagnostics(bool all=false);

    /**
        Updates diagnostic information

        \param key name of diagnostic (a python string in the list "num_iter",
                 "num_level", "num_inner_iter", "time", "set_up_time",
                 "net_time", "residual_norm", "converged").
        \param value new value of the diagnostic information
    */
    ESCRIPT_DLL_API void updateDiagnostics(std::string key, boost::python::object value);

    /**
    Returns the diagnostic information for the given ``name``. Possible values are:

    - "num_iter": the number of iteration steps
    - "cum_num_iter": the cumulative number of iteration steps
    - "num_level": the number of level in multi level solver
    - "num_inner_iter": the number of inner iteration steps
    - "cum_num_inner_iter": the cumulative number of inner iteration steps
    - "time": execution time
    - "cum_time": cumulative execution time
    - "set_up_time": time to set up of the solver, typically this includes factorization and reordering
    - "cum_set_up_time": cumulative time to set up of the solver
    - "net_time": net execution time, excluding setup time for the solver and execution time for preconditioner
    - "cum_net_time": cumulative net execution time
    - "preconditioner_size": size of preconditioner [Bytes]
    - "converged": return true if solution has converged.
    - "time_step_backtracking_used": returns true if time step back tracking has been used.
    - "coarse_level_sparsity": returns the sparsity of the matrix on the coarsest level
    - "num_coarse_unknowns": returns the number of unknowns on the coarsest level


    \param name name of diagnostic information to return
    
    \returns requested value. 0 is returned if the value is yet to be defined.

    \note If the solver has thrown an exception diagnostic values have an undefined status.
    */
    ESCRIPT_DLL_API double getDiagnostics(std::string name) const;

    /**
        Returns ``true`` if the last solver call has been finalized successfully.

        \note if an exception has been thrown by the solver the status of this
                 flag is undefined.
    */
    ESCRIPT_DLL_API bool hasConverged() const;

    /**
    Sets the key of the coarsening method to be applied in AMG or AMLI or BoomerAMG

    \param coarsening selects the coarsening method, should be in
            `ESCRIPT_DEFAULT`, `ESCRIPT_YAIR_SHAPIRA_COARSENING`,
            `ESCRIPT_RUGE_STUEBEN_COARSENING`, `ESCRIPT_AGGREGATION_COARSENING`,
            `ESCRIPT_CIJP_FIXED_RANDOM_COARSENING`, `ESCRIPT_CIJP_COARSENING`,
            `ESCRIPT_FALGOUT_COARSENING`, `ESCRIPT_PMIS_COARSENING`,
            `ESCRIPT_HMIS_COARSENING`
    */
    ESCRIPT_DLL_API void setCoarsening(int coarsening);

    /**
        Returns the key of the coarsening algorithm to be applied AMG, AMLI
        or BoomerAMG

    */
    ESCRIPT_DLL_API SolverOptions getCoarsening() const;

    /**
        Sets the minimum size of the coarsest level matrix in AMG or AMLI

        \param size minimum size of the coarsest level matrix .
    */
    ESCRIPT_DLL_API void setMinCoarseMatrixSize(int size);

    /**
        Returns the minimum size of the coarsest level matrix in AMG or AMLI
    */
    ESCRIPT_DLL_API int getMinCoarseMatrixSize() const;

    /**
    Sets the preconditioner to be used.

    \param preconditioner key of the preconditioner to be used, should be in
        `ESCRIPT_ILU0`, `ESCRIPT_ILUT`,
        `ESCRIPT_JACOBI`, `ESCRIPT_AMG`, `ESCRIPT_AMLI`,
        `ESCRIPT_REC_ILU`, `ESCRIPT_GAUSS_SEIDEL`,
        `ESCRIPT_RILU`, `ESCRIPT_BOOMERAMG`,
        `ESCRIPT_NO_PRECONDITIONER`
    
    \note Not all packages support all preconditioner. It can be assumed that
            a package makes a reasonable choice if it encounters an unknown
            preconditioner.
    */
    ESCRIPT_DLL_API void setPreconditioner(int preconditioner);

    /**
    Returns the key of the preconditioner to be used.

    */
    ESCRIPT_DLL_API SolverOptions getPreconditioner() const;

    /**
    Sets the smoother to be used.

    \param smoother key of the smoother to be used, should be in 
            `ESCRIPT_JACOBI`, `ESCRIPT_GAUSS_SEIDEL`

    \note Not all packages support all smoothers. It can be assumed that a
            package makes a reasonable choice if it encounters an unknown
            smoother.
    */
    ESCRIPT_DLL_API void setSmoother(int smoother);

    /**
    Returns key of the smoother to be used.

    */
    ESCRIPT_DLL_API SolverOptions getSmoother() const;

    /**
    Sets the solver method to be used. Use ``method``=``ESCRIPT_DIRECT`` to
    indicate that a direct rather than an iterative solver should be used and
    use ``method``=``ESCRIPT_ITERATIVE`` to indicate that an iterative rather
    than a direct solver should be used.

    \param method key of the solver method to be used, should be in
            `ESCRIPT_DEFAULT`, `ESCRIPT_DIRECT`, `ESCRIPT_CHOLEVSKY`,
            `ESCRIPT_PCG`, `ESCRIPT_CR`, `ESCRIPT_CGS`, `ESCRIPT_BICGSTAB`,
            `ESCRIPT_GMRES`, `ESCRIPT_PRES20`, `ROWSUM_ESCRIPT_LUMPING`,
            `HRZ_ESCRIPT_LUMPING`, `ESCRIPT_ITERATIVE`, `ESCRIPT_LSQR`,
            `ESCRIPT_NONLINEAR_GMRES`, `ESCRIPT_TFQMR`, `ESCRIPT_MINRES`
    
    \note Not all packages support all solvers. It can be assumed that a
            package makes a reasonable choice if it encounters an unknown
            solver method.
    */
    ESCRIPT_DLL_API void setSolverMethod(int method);

    /**
    Returns key of the solver method to be used.

    */
    ESCRIPT_DLL_API SolverOptions getSolverMethod() const;

    /**
    Sets the solver target to be used. By default the solver is run on
    the host CPU(s). If escript was compiled with CUDA support then
    `ESCRIPT_TARGET_GPU` is a valid option and the solver will be run
    on GPU(s).

    \param target key of the solver target. Valid settings:
            `ESCRIPT_TARGET_CPU`, `ESCRIPT_TARGET_GPU`
    
    */
    ESCRIPT_DLL_API void setSolverTarget(int target);

    /**
    Returns key of the solver target.

    */
    ESCRIPT_DLL_API SolverOptions getSolverTarget() const;

    /**
    Sets the solver package to be used as a solver.

    \param package key of the solver package to be used, should be in 
            `ESCRIPT_DEFAULT`, `ESCRIPT_PASO`, `ESCRIPT_SUPER_LU`,
            `ESCRIPT_PASTIX`, `ESCRIPT_MKL`, `ESCRIPT_UMFPACK`,
            `ESCRIPT_TRILINOS`

    \note Not all packages are support on all implementation. An exception may
            be thrown on some platforms if a particular is requested.
    */
    ESCRIPT_DLL_API void setPackage(int package);

    /**
    Returns the solver package key

    */
    ESCRIPT_DLL_API SolverOptions getPackage() const;

    /**
    Sets the key of the reordering method to be applied if supported by the
    solver. Some direct solvers support reordering
    to optimize compute time and storage use during elimination.

    \param ordering selects the reordering strategy, should be in 
            'ESCRIPT_NO_REORDERING', 'ESCRIPT_MINIMUM_FILL_IN',
            'ESCRIPT_NESTED_DISSECTION', 'ESCRIPT_DEFAULT_REORDERING'
    */
    ESCRIPT_DLL_API void setReordering(int ordering);

    /**
    Returns the key of the reordering method to be applied if supported by the
    solver.

    */
    ESCRIPT_DLL_API SolverOptions getReordering() const;

    /**
    Sets the number of iterations steps after which GMRES performs a restart.

    \param restart number of iteration steps after which to perform a
        restart. If 0 no restart is performed.
    */
    ESCRIPT_DLL_API void setRestart(int restart);

    /**
    Returns the number of iterations steps after which GMRES performs a
    restart. If 0 is returned no restart is performed.

    */
    ESCRIPT_DLL_API int getRestart() const;

    /**
    Returns the number of iterations steps after which GMRES performs a
    restart. If -1 is returned no restart is performed.

    */
    ESCRIPT_DLL_API int _getRestartForC() const;

    /**
    Sets the threshold for diagonal dominant rows which are eliminated during
    AMG coarsening.

    */
    ESCRIPT_DLL_API void setDiagonalDominanceThreshold(double threshold);

    /**
    Returns the threshold for diagonal dominant rows which are eliminated
    during AMG coarsening.

    */
    ESCRIPT_DLL_API double getDiagonalDominanceThreshold() const;

    /**
    Sets the number of residuals in GMRES to be stored for orthogonalization.
    The more residuals are stored the faster GMRES converged but

    */
    ESCRIPT_DLL_API void setTruncation(int truncation);

    /**
    Returns the number of residuals in GMRES to be stored for orthogonalization

    */
    ESCRIPT_DLL_API int getTruncation() const;

    /**
    Sets the maximum number of iteration steps for the inner iteration.

    \param iter_max maximum number of inner iterations
    */
    ESCRIPT_DLL_API void setInnerIterMax(int iter_max);

    /**
    Returns maximum number of inner iteration steps

    */
    ESCRIPT_DLL_API int getInnerIterMax() const;

    /**
    Sets the maximum number of iteration steps

    \param iter_max maximum number of iteration steps
    */
    ESCRIPT_DLL_API void setIterMax(int iter_max);

    /**
    Returns maximum number of iteration steps

    */
    ESCRIPT_DLL_API int getIterMax() const;

    /**
    Sets the maximum number of coarsening levels to be used in an algebraic
    multi level solver or preconditioner

    \param level_max maximum number of levels
    */
    ESCRIPT_DLL_API void setLevelMax(int level_max);

    /**
    Returns the maximum number of coarsening levels to be used in an algebraic
    multi level solver or preconditioner

    */
    ESCRIPT_DLL_API int getLevelMax() const;

    /**
    Sets the cycle type (V-cycle or W-cycle) to be used in an algebraic multi
    level solver or preconditioner

    \param cycle_type the type of cycle
    */
    ESCRIPT_DLL_API void setCycleType(int cycle_type);

    /**
    Returns the cycle type (V- or W-cycle) to be used in an algebraic multi
    level solver or preconditioner

    */
    ESCRIPT_DLL_API int getCycleType() const;

    /**
    Sets the threshold for coarsening in the algebraic multi level solver or
    preconditioner

    \param theta threshold for coarsening
    */
    ESCRIPT_DLL_API void setCoarseningThreshold(double theta);

    /**
    Returns the threshold for coarsening in the algebraic multi level solver
    or preconditioner

    */
    ESCRIPT_DLL_API double getCoarseningThreshold() const;

    /**
    Sets the number of sweeps in a Jacobi or Gauss-Seidel/SOR preconditioner.

    \param sweeps number of sweeps
    */
    ESCRIPT_DLL_API void setNumSweeps(int sweeps);

    /**
    Returns the number of sweeps in a Jacobi or Gauss-Seidel/SOR preconditioner.

    */
    ESCRIPT_DLL_API int getNumSweeps() const;

    /**
    Sets the number of sweeps in the pre-smoothing step of a multi level
    solver or preconditioner

    \param sweeps number of sweeps
    */
    ESCRIPT_DLL_API void setNumPreSweeps(int sweeps);

    /**
    Returns he number of sweeps in the pre-smoothing step of a multi level
    solver or preconditioner

    */
    ESCRIPT_DLL_API int getNumPreSweeps() const;

    /**
    Sets the number of sweeps in the post-smoothing step of a multi level
    solver or preconditioner

    \param sweeps number of sweeps
    */
    ESCRIPT_DLL_API void setNumPostSweeps(int sweeps);

    /**
    Returns the number of sweeps in the post-smoothing step of a multi level
    solver or preconditioner

    */
    ESCRIPT_DLL_API int getNumPostSweeps() const;

    /**
    Sets the relative tolerance for the solver

    \param rtol relative tolerance
    */
    ESCRIPT_DLL_API void setTolerance(double rtol);

    /**
    Returns the relative tolerance for the solver

    */
    ESCRIPT_DLL_API double getTolerance() const;

    /**
    Sets the absolute tolerance for the solver

    \param atol absolute tolerance
    */
    ESCRIPT_DLL_API void setAbsoluteTolerance(double atol);

    /**
    Returns the absolute tolerance for the solver

    */
    ESCRIPT_DLL_API double getAbsoluteTolerance() const;

    /**
     Sets the relative tolerance for an inner iteration scheme for instance on
     the coarsest level in a multi-level scheme.

    \param rtol inner relative tolerance
    */
    ESCRIPT_DLL_API void setInnerTolerance(double rtol);

    /**
    Returns the relative tolerance for an inner iteration scheme

    */
    ESCRIPT_DLL_API double getInnerTolerance() const;

    /**
    Sets the relative drop tolerance in ILUT

    \param drop_tol drop tolerance
    */
    ESCRIPT_DLL_API void setDropTolerance(double drop_tol);

    /**
    Returns the relative drop tolerance in ILUT

    */
    ESCRIPT_DLL_API double getDropTolerance() const;

    /**
    Sets the maximum allowed increase in storage for ILUT. An increase of 2 would
    mean that a doubling of the storage needed for the coefficient matrix is
    allowed during ILUT factorization.

    \param storage allowed storage increase
    */
    ESCRIPT_DLL_API void setDropStorage(double drop);

    /**
    Returns the maximum allowed increase in storage for ILUT

    */
    ESCRIPT_DLL_API double getDropStorage() const;

    /**
    Sets the relaxation factor used to add dropped elements in RILU to the
    main diagonal.

    \param factor relaxation factor
    \note RILU with a relaxation factor 0 is identical to ILU0
    */
    ESCRIPT_DLL_API void setRelaxationFactor(double factor);


    /**
    Returns the relaxation factor used to add dropped elements in RILU to the
    main diagonal.

    */
    ESCRIPT_DLL_API double getRelaxationFactor() const;

    /**
    Checks if symmetry of the coefficient matrix is indicated.

    \return true if a symmetric PDE is indicated, false otherwise
    */
    ESCRIPT_DLL_API bool isSymmetric() const;

    /**
    Sets the symmetry flag to indicate that the coefficient matrix is symmetric.
    */
    ESCRIPT_DLL_API void setSymmetryOn();

    /**
    Clears the symmetry flag for the coefficient matrix.
    */
    ESCRIPT_DLL_API void setSymmetryOff();

    /**
    Sets the symmetry flag for the coefficient matrix to ``flag``.

    \param flag If true, the symmetry flag is set otherwise reset.
    */ESCRIPT_DLL_API void setSymmetry(bool symmetry);

    /**
    Returns ``true`` if the solver is expected to be verbose.

    \return true if verbosity is switched on.
    */
    ESCRIPT_DLL_API bool isVerbose() const;

    /**
    Switches the verbosity of the solver on.
    */
    ESCRIPT_DLL_API void setVerbosityOn();

    /**
    Switches the verbosity of the solver off.
    */
    ESCRIPT_DLL_API void setVerbosityOff();

    /**
    Sets the verbosity flag for the solver to ``flag``.

    \param verbose If ``true``, the verbosity of the solver is switched on.
    */
    ESCRIPT_DLL_API void setVerbosity(bool verbosity);

    /**
    Returns ``true`` if the tolerance of the inner solver is selected
    automatically. Otherwise the inner tolerance set by `setInnerTolerance` is
    used.

    \returns ``true`` if inner tolerance adaption is chosen.
    */
    ESCRIPT_DLL_API bool adaptInnerTolerance() const;

    /**
    Switches the automatic selection of inner tolerance on
    */
    ESCRIPT_DLL_API void setInnerToleranceAdaptionOn();

    /**
    Switches the automatic selection of inner tolerance off.
    */
    ESCRIPT_DLL_API void setInnerToleranceAdaptionOff();

    /**
    Sets the flag to indicate automatic selection of the inner tolerance.

    \param adapt If ``true``, the inner tolerance is selected automatically.
    */
    ESCRIPT_DLL_API void setInnerToleranceAdaption(bool adaption);

    /**
    Returns ``true`` if a failure to meet the stopping criteria within the
    given number of iteration steps is not raising in exception. This is useful
    if a solver is used in a non-linear context where the non-linear solver can
    continue even if the returned the solution does not necessarily meet the
    stopping criteria. One can use the `hasConverged` method to check if the
    last call to the solver was successful.

    \returns ``true`` if a failure to achieve convergence is accepted.
    */
    ESCRIPT_DLL_API bool acceptConvergenceFailure() const;

    /**
    Switches the acceptance of a failure of convergence on
    */
    ESCRIPT_DLL_API void setAcceptanceConvergenceFailureOn();

    /**
    Switches the acceptance of a failure of convergence off.
    */
    ESCRIPT_DLL_API void setAcceptanceConvergenceFailureOff();

    /**
    Sets the flag to indicate the acceptance of a failure of convergence.

    \param accept If ``true``, any failure to achieve convergence is accepted.
    */
    ESCRIPT_DLL_API void setAcceptanceConvergenceFailure(bool acceptance);

    /**
    Returns ``true`` if the preconditoner is applied locally on each MPI. This
    reduces communication costs and speeds up the application of the
    preconditioner but at the costs of more iteration steps. This can be an
    advantage on clusters with slower interconnects.
    */
    ESCRIPT_DLL_API bool useLocalPreconditioner() const;

    /**
    Sets the flag to use local preconditioning to on
    */
    ESCRIPT_DLL_API void setLocalPreconditionerOn();

    /**
    Sets the flag to use local preconditioning to off
    */
    ESCRIPT_DLL_API void setLocalPreconditionerOff();

    /**
    Sets the flag to use local preconditioning

    \param use If ``true``, local proconditioning on each MPI rank is applied
    */
    ESCRIPT_DLL_API void setLocalPreconditioner(bool local);

    /**
    Sets the minimum sparsity at the coarsest level. Typically
    a direct solver is used when the sparsity becomes bigger than
    the set limit.

    \param sparsity minimal sparsity
    */
    ESCRIPT_DLL_API void setMinCoarseMatrixSparsity(double sparsity);

    /**
    Returns the minimum sparsity on the coarsest level. Typically
    a direct solver is used when the sparsity becomes bigger than
    the set limit.

    \returns minimal sparsity
    */
    ESCRIPT_DLL_API double getMinCoarseMatrixSparsity() const;

    /**
    Sets the number of refinement steps to refine the solution when a direct
    solver is applied.

    \param refinements number of refinements
    */
    ESCRIPT_DLL_API void setNumRefinements(int refinements);

    /**
    Returns the number of refinement steps to refine the solution when a direct
    solver is applied.

    */
    ESCRIPT_DLL_API int getNumRefinements() const;

    /**
    Sets the number of refinement steps to refine the solution on the coarsest
    level when a direct solver is applied.

    \param refinements number of refinements
    */
    ESCRIPT_DLL_API void setNumCoarseMatrixRefinements(int refinements);

    /**
    Returns the number of refinement steps to refine the solution on the
    coarsest level when a direct solver is applied.

    */
    ESCRIPT_DLL_API int getNumCoarseMatrixRefinements() const;

    /**
    Returns ``true`` if a panel is used to search for unknown in the AMG
    coarsening, The panel approach is normally faster
    but can lead to larger coarse level systems.

    */
    ESCRIPT_DLL_API bool usePanel() const;

    /**
    Sets the flag to use a panel to find unknowns in AMG coarsening
    */
    ESCRIPT_DLL_API void setUsePanelOn();

    /**
    Sets the flag to use a panel to find unknowns in AMG coarsening to off
    */
    ESCRIPT_DLL_API void setUsePanelOff();

    /**
    Sets the flag to use  a panel to find unknowns in AMG coarsening

    \param use If ``true``,a panel is used to find unknowns in AMG coarsening
    */
    ESCRIPT_DLL_API void setUsePanel(bool use);

    /**
    Set the interpolation method for the AMG preconditioner.

    \param method key of the interpolation method to be used, should be in 
            `ESCRIPT_CLASSIC_INTERPOLATION_WITH_FF_COUPLING`,
            `ESCRIPT_CLASSIC_INTERPOLATION`, `ESCRIPT_DIRECT_INTERPOLATION`
    */
    ESCRIPT_DLL_API void setAMGInterpolation(int interpolation);

    /**
    Returns key of the interpolation method for the SAMG preconditioner

    */
    ESCRIPT_DLL_API SolverOptions getAMGInterpolation() const;

    /**
    Set the solver method for ODEs.

    \param method key of the ODE solver method to be used, should be in 
            `ESCRIPT_CRANK_NICOLSON`, `ESCRIPT_BACKWARD_EULER`,
            `ESCRIPT_LINEAR_CRANK_NICOLSON`
    */
    ESCRIPT_DLL_API void setODESolver(int solver);

    /**
    Returns key of the solver method for ODEs.

    */
    ESCRIPT_DLL_API SolverOptions getODESolver() const;


protected:
    int level_max;
    double coarsening_threshold;
    SolverOptions smoother;
    int sweeps;
    int pre_sweeps;
    int post_sweeps;
    double tolerance;
    double absolute_tolerance;
    double inner_tolerance;
    double drop_tolerance;
    double drop_storage;
    int iter_max;
    int inner_iter_max;
    int truncation;
    int restart; //0 will have to be None in python, will get tricky
    bool symmetric;
    bool verbose;
    bool adapt_inner_tolerance;
    bool accept_convergence_failure;
    SolverOptions reordering;
    SolverOptions package;
    SolverOptions method;
    SolverOptions preconditioner;
    SolverOptions coarsening;
    SolverOptions target;
    int MinCoarseMatrixSize;
    double relaxation;
    bool use_local_preconditioner;
    double min_sparsity;
    int refinements;
    int coarse_refinements;
    bool use_panel;
    double diagonal_dominance_threshold;
    SolverOptions amg_interpolation_method;
    int cycle_type;
    SolverOptions ode_solver;

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
    int cum_num_inner_iter;
    int cum_num_iter;
    double cum_time;
    double cum_set_up_time;
    double cum_net_time;
    double coarse_level_sparsity;
    int num_coarse_unknowns;
};

typedef boost::shared_ptr<SolverBuddy> SB_ptr;

}

#endif
