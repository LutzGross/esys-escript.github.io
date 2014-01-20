# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2013 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

__copyright__="""Copyright (c) 2003-2013 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
The module provides an interface to define and solve linear partial
differential equations (PDEs) and Transport problems within `escript`.
`linearPDEs` does not provide any solver capabilities in itself but hands the
PDE over to the PDE solver library defined through the `Domain`
of the PDE. The general interface is provided through the `LinearPDE` class.
`TransportProblem` provides an interface to initial value problems dominated
by its advective terms.

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

from . import escriptcpp as escore
from . import util
import math
import numpy

__author__="Lutz Gross, l.gross@uq.edu.au"


class SolverOptions(object):
    """
    this class defines the solver options for a linear or non-linear solver.
    
    The option also supports the handling of diagnostic informations. 
    
    Typical usage is 
    
    ::

      opts=SolverOptions()
      print(opts)
      opts.resetDiagnostics()
      u=solver(opts)
      print("number of iteration steps: =",opts.getDiagnostics("num_iter"))

    :cvar DEFAULT: The default method used to solve the system of linear equations
    :cvar DIRECT: The direct solver based on LDU factorization
    :cvar CHOLEVSKY: The direct solver based on LDLt factorization (can only be applied for symmetric PDEs)
    :cvar PCG: The preconditioned conjugate gradient method (can only be applied for symmetric PDEs)
    :cvar CR: The conjugate residual method
    :cvar CGS: The conjugate gradient square method
    :cvar BICGSTAB: The stabilized Bi-Conjugate Gradient method
    :cvar TFQMR: Transpose Free Quasi Minimal Residual method
    :cvar MINRES: Minimum residual method
    :cvar ILU0: The incomplete LU factorization preconditioner with no fill-in
    :cvar ILUT: The incomplete LU factorization preconditioner with fill-in
    :cvar JACOBI: The Jacobi preconditioner
    :cvar GMRES: The Gram-Schmidt minimum residual method
    :cvar PRES20: Special GMRES with restart after 20 steps and truncation after 5 residuals
    :cvar ROWSUM_LUMPING: Matrix lumping using row sum
    :cvar HRZ_LUMPING: Matrix lumping using the HRZ approach 
    :cvar NO_REORDERING: No matrix reordering allowed
    :cvar MINIMUM_FILL_IN: Reorder matrix to reduce fill-in during factorization
    :cvar NESTED_DISSECTION: Reorder matrix to improve load balancing during factorization
    :cvar PASO: PASO solver package
    :cvar SCSL: SGI SCSL solver library
    :cvar MKL: Intel's MKL solver library
    :cvar UMFPACK: The UMFPACK library
    :cvar TRILINOS: The TRILINOS parallel solver class library from Sandia National Labs
    :cvar ITERATIVE: The default iterative solver
    :cvar AMG: Algebraic Multi Grid
    :cvar AMLI: Algebraic Multi Level Iteration
    :cvar REC_ILU: recursive ILU0
    :cvar RILU: relaxed ILU0
    :cvar GAUSS_SEIDEL: Gauss-Seidel preconditioner
    :cvar DEFAULT_REORDERING: the reordering method recommended by the solver
    :cvar SUPER_LU: the Super_LU solver package
    :cvar PASTIX: the Pastix direct solver_package
    :cvar YAIR_SHAPIRA_COARSENING: AMG and AMLI coarsening method by Yair-Shapira
    :cvar RUGE_STUEBEN_COARSENING: AMG and AMLI coarsening method by Ruge and Stueben
    :cvar AGGREGATION_COARSENING: AMG and AMLI coarsening using (symmetric) aggregation
    :cvar STANDARD_COARSENING: AMG and AMLI standard coarsening using mesure of importance of the unknowns
    :cvar MIN_COARSE_MATRIX_SIZE: minimum size of the coarsest level matrix to use direct solver.
    :cvar NO_PRECONDITIONER: no preconditioner is applied.
    :cvar CLASSIC_INTERPOLATION_WITH_FF_COUPLING: classical interpolation in AMG with enforced 
    :cvar CLASSIC_INTERPOLATION: classical interpolation in AMG
    :cvar DIRECT_INTERPOLATION: direct interploation in AMG
    :cvar BOOMERAMG: Boomer AMG in hypre library
    :cvar CIJP_FIXED_RANDOM_COARSENING: BoomerAMG parallel coarsening method CIJP by using fixed random vector
    :cvar CIJP_COARSENING: BoomerAMG parallel coarsening method CIJP
    :cvar PASO_FALGOUT_COARSENING: BoomerAMG parallel coarsening method falgout
    :cvar PASO_PMIS_COARSENING: BoomerAMG parallel coarsening method PMIS
    :cvar PASO_HMIS_COARSENING: BoomerAMG parallel coarsening method HMIS
    :cvar BACKWARD_EULER: backward Euler scheme
    :cvar CRANK_NICOLSON: Crank-Nicolson scheme
    :cvar LINEAR_CRANK_NICOLSON: linerized Crank-Nicolson scheme

    """
    DEFAULT= 0
    DIRECT= 1
    CHOLEVSKY= 2
    PCG= 3
    CR= 4
    CGS= 5
    BICGSTAB= 6
    ILU0= 8
    ILUT= 9
    JACOBI= 10
    GMRES= 11
    PRES20= 12
    LUMPING= 13
    ROWSUM_LUMPING= 13
    HRZ_LUMPING= 14
    NO_REORDERING= 17
    MINIMUM_FILL_IN= 18
    NESTED_DISSECTION= 19
    MKL= 15
    UMFPACK= 16
    ITERATIVE= 20
    PASO= 21
    AMG= 22
    REC_ILU = 23
    TRILINOS = 24
    NONLINEAR_GMRES = 25
    TFQMR = 26
    MINRES = 27
    GAUSS_SEIDEL=28
    RILU=29
    DEFAULT_REORDERING=30
    SUPER_LU=31
    PASTIX=32
    YAIR_SHAPIRA_COARSENING=33
    RUGE_STUEBEN_COARSENING=34
    AGGREGATION_COARSENING=35
    NO_PRECONDITIONER=36
    MIN_COARSE_MATRIX_SIZE=37
    AMLI=38
    STANDARD_COARSENING=39
    CLASSIC_INTERPOLATION_WITH_FF_COUPLING=50
    CLASSIC_INTERPOLATION=51
    DIRECT_INTERPOLATION=52
    BOOMERAMG=60
    CIJP_FIXED_RANDOM_COARSENING=61
    CIJP_COARSENING=62
    FALGOUT_COARSENING=63
    PMIS_COARSENING=64
    HMIS_COARSENING=65
    LINEAR_CRANK_NICOLSON=66
    CRANK_NICOLSON=67
    BACKWARD_EULER=68
 
    def __init__(self):
        self.setLevelMax()
        self.setCoarseningThreshold()
        self.setSmoother()
        self.setNumSweeps()
        self.setNumPreSweeps()
        self.setNumPostSweeps()
        self.setTolerance()
        self.setAbsoluteTolerance()
        self.setInnerTolerance()
        self.setDropTolerance()
        self.setDropStorage()
        self.setIterMax()
        self.setInnerIterMax()
        self.setTruncation()
        self.setRestart()
        self.setSymmetry()
        self.setVerbosity()
        self.setInnerToleranceAdaption()
        self.setAcceptanceConvergenceFailure()
        self.setReordering()
        self.setPackage()
        self.setSolverMethod()
        self.setPreconditioner()
        self.setCoarsening()
        self.setMinCoarseMatrixSize()
        self.setRelaxationFactor()        
        self.setLocalPreconditionerOff()
        self.resetDiagnostics(all=True)
        self.setMinCoarseMatrixSparsity()
        self.setNumRefinements()
        self.setNumCoarseMatrixRefinements()
        self.setUsePanel()
        self.setDiagonalDominanceThreshold()
        self.setAMGInterpolation()
        self.setCycleType()
        self.setODESolver()
        

    def __str__(self):
        return self.getSummary()
    def getSummary(self):
        """
        Returns a string reporting the current settings
        """
        out="Solver Package: %s"%(self.getName(self.getPackage()))
        out+="\nVerbosity = %s"%self.isVerbose()
        out+="\nAccept failed convergence = %s"%self.acceptConvergenceFailure()
        out+="\nRelative tolerance = %e"%self.getTolerance()
        out+="\nAbsolute tolerance = %e"%self.getAbsoluteTolerance()
        out+="\nSymmetric problem = %s"%self.isSymmetric()
        out+="\nMaximum number of iteration steps = %s"%self.getIterMax()
        # out+="\nInner tolerance = %e"%self.getInnerTolerance()
        # out+="\nAdapt innner tolerance = %s"%self.adaptInnerTolerance()

        if self.getPackage() in (self.DEFAULT, self.PASO):
            out+="\nSolver method = %s"%self.getName(self.getSolverMethod())
            if self.getSolverMethod() == self.GMRES:
                out+="\nTruncation  = %s"%self.getTruncation()
                out+="\nRestart  = %s"%self.getRestart()
            if self.getSolverMethod() == self.AMG:
                out+="\nNumber of pre / post sweeps = %s / %s, %s"%(self.getNumPreSweeps(), self.getNumPostSweeps(), self.getNumSweeps())
                out+="\nMaximum number of levels = %s"%self.getLevelMax()
                out+="\nCoarsening threshold = %e"%self.getCoarseningThreshold()
                out+="\nCoarsening method = %s"%self.getName(self.getCoarsening())
            out+="\nPreconditioner = %s"%self.getName(self.getPreconditioner())
            out+="\nApply preconditioner locally = %s"%self.useLocalPreconditioner()
            if self.getPreconditioner() == self.AMG:
                out+="\nMaximum number of levels = %s"%self.getLevelMax()
                out+="\nCoarsening threshold = %e"%self.getCoarseningThreshold()
                out+="\nMinimal sparsity on coarsest level = %e"%self.getMinCoarseMatrixSparsity()
                out+="\nSmoother = %s"%self.getName(self.getSmoother())
                out+="\nMinimum size of the coarsest level matrix = %e"%self.getMinCoarseMatrixSize()
                out+="\nNumber of pre / post sweeps = %s / %s, %s"%(self.getNumPreSweeps(), self.getNumPostSweeps(), self.getNumSweeps())
                out+="\nNumber of refinement steps in coarsest level solver = %d"%self.getNumCoarseMatrixRefinements()
                out+="\nUse node panel = %s"%self.usePanel()
                out+="\nInterpolation = %s"%(self.getName(self.getAMGInterpolation()))
                out+="\nThreshold for diagonal dominant rows = %s"%self.getDiagonalDominanceThreshold()

            if self.getPreconditioner() == self.AMLI:
                out+="\nMaximum number of levels = %s"%self.getLevelMax()
                out+="\nCoarsening method = %s"%self.getName(self.getCoarsening())
                out+="\nCoarsening threshold = %e"%self.getMinCoarseMatrixSize()
                out+="\nMinimum size of the coarsest level matrix = %e"%self.getCoarseningThreshold()
                out+="\nNumber of pre / post sweeps = %s / %s, %s"%(self.getNumPreSweeps(), self.getNumPostSweeps(), self.getNumSweeps())
            if self.getPreconditioner() == self.BOOMERAMG:
                out+="\nMaximum number of levels = %s"%self.getLevelMax()
                out+="\nCoarsening threshold = %e"%self.getCoarseningThreshold()
                out+="\nThreshold for diagonal dominant rows = %s"%setDiagonalDominanceThreshold()
                out+="\nCoarsening method = %s"%self.getName(self.getCoarsening())
                out+="\nV-cycle (1) or W-cyle (2) = %s"%self.getCycleType()
                out+="\nNumber of pre / post sweeps = %s / %s, %s"%(self.getNumPreSweeps(), self.getNumPostSweeps(), self.getNumSweeps())
                out+="\nSmoother = %s"%self.getName(self.getSmoother())
            if self.getPreconditioner() == self.GAUSS_SEIDEL:
                out+="\nNumber of sweeps = %s"%self.getNumSweeps()
            if self.getPreconditioner() == self.ILUT:
                out+="\nDrop tolerance = %e"%self.getDropTolerance()
                out+="\nStorage increase = %e"%self.getDropStorage()
            if self.getPreconditioner() == self.RILU:
                out+="\nRelaxation factor = %e"%self.getRelaxationFactor()
            out+="\nODE solver = %s"%self.getName(self.getODESolver())
        return out
        
    def getName(self,key):
        """
        returns the name of a given key
        
        :param key: a valid key
        """
        if key == self.DEFAULT: return "DEFAULT"
        if key == self.DIRECT: return "DIRECT"
        if key == self.CHOLEVSKY: return "CHOLEVSKY"
        if key == self.PCG: return "PCG"
        if key == self.CR: return "CR"
        if key == self.CGS: return "CGS"
        if key == self.BICGSTAB: return "BICGSTAB"
        if key == self.ILU0: return "ILU0:"
        if key == self.ILUT: return "ILUT"
        if key == self.JACOBI: return "JACOBI"
        if key == self.GMRES: return "GMRES"
        if key == self.PRES20: return "PRES20"
        if key == self.ROWSUM_LUMPING: return "ROWSUM_LUMPING"
        if key == self.HRZ_LUMPING: return "HRZ_LUMPING"
        if key == self.NO_REORDERING: return "NO_REORDERING"
        if key == self.MINIMUM_FILL_IN: return "MINIMUM_FILL_IN"
        if key == self.NESTED_DISSECTION: return "NESTED_DISSECTION"
        if key == self.MKL: return "MKL"
        if key == self.UMFPACK: return "UMFPACK"
        if key == self.ITERATIVE: return "ITERATIVE"
        if key == self.PASO: return "PASO"
        if key == self.AMG: return "AMG"
        if key == self.AMLI: return "AMLI"
        if key == self.REC_ILU: return "REC_ILU"
        if key == self.TRILINOS: return "TRILINOS"
        if key == self.NONLINEAR_GMRES: return "NONLINEAR_GMRES"
        if key == self.TFQMR: return "TFQMR"
        if key == self.MINRES: return "MINRES"
        if key == self.GAUSS_SEIDEL: return "GAUSS_SEIDEL"
        if key == self.RILU: return "RILU"
        if key == self.DEFAULT_REORDERING: return "DEFAULT_REORDERING"
        if key == self.SUPER_LU: return "SUPER_LU"
        if key == self.PASTIX: return "PASTIX"
        if key == self.YAIR_SHAPIRA_COARSENING: return "YAIR_SHAPIRA_COARSENING"
        if key == self.RUGE_STUEBEN_COARSENING: return "RUGE_STUEBEN_COARSENING"
        if key == self.STANDARD_COARSENING: return "STANDARD_COARSENING"
        if key == self.AGGREGATION_COARSENING: return "AGGREGATION_COARSENING"
        if key == self.NO_PRECONDITIONER: return "NO_PRECONDITIONER"
        if key == self.CLASSIC_INTERPOLATION_WITH_FF_COUPLING: return "CLASSIC_INTERPOLATION_WITH_FF"
        if key == self.CLASSIC_INTERPOLATION: return "CLASSIC_INTERPOLATION"
        if key == self.DIRECT_INTERPOLATION: return "DIRECT_INTERPOLATION"
        if key == self.BOOMERAMG: return "BOOMERAMG"
        if key == self.CIJP_FIXED_RANDOM_COARSENING: return "CIJP_FIXED_RANDOM_COARSENING"
        if key == self.CIJP_COARSENING: return "CIJP_COARSENING"
        if key == self.FALGOUT_COARSENING: return "FALGOUT_COARSENING"
        if key == self.PMIS_COARSENING: return "PMIS_COARSENING"
        if key == self.HMIS_COARSENING: return "HMIS_COARSENING"
        if key == self.LINEAR_CRANK_NICOLSON: return "LINEAR_CRANK_NICOLSON"
        if key == self.CRANK_NICOLSON: return "CRANK_NICOLSON"
        if key == self.BACKWARD_EULER: return "BACKWARD_EULER"


     
    def resetDiagnostics(self,all=False):
        """
        resets the diagnostics
        
        :param all: if ``all`` is ``True`` all diagnostics including accumulative counters are reset.
        :type all: ``bool``
        """
        self.__num_iter=None
        self.__num_level=None
        self.__num_inner_iter=None
        self.__time=None
        self.__set_up_time=None
        self.__net_time=None
        self.__residual_norm=None
        self.__converged=None
        self.__preconditioner_size=-1
        self.__time_step_backtracking_used=None
        if all: 
            self.__cum_num_inner_iter=0
            self.__cum_num_iter=0
            self.__cum_time=0
            self.__cum_set_up_time=0
            self.__cum_net_time=0

    def _updateDiagnostics(self, name, value):
        """
        Updates diagnostic information 
        
        :param name: name of  diagnostic information
        :type name: ``str`` in the list "num_iter", "num_level", "num_inner_iter", "time", "set_up_time", "net_time", "residual_norm", "converged".
        :param value: new value of the diagnostic information
        :note: this function is used by a solver to report diagnostics informations.
        """
        if name == "num_iter": 
            self.__num_iter=int(value)
            self.__cum_num_iter+=self.__num_iter
        if name == "num_level":
            self.__num_level=int(value)
        if name == "num_inner_iter":
            self.__num_inner_iter=int(value)
            self.__cum_num_inner_iter+=self.__num_inner_iter
        if name == "time":
            self.__time=float(value)
            self.__cum_time+=self.__time
        if name == "set_up_time":
            self.__set_up_time=float(value)
            self.__cum_set_up_time+=self.__set_up_time
        if name == "net_time":
            self.__net_time=float(value)
            self.__cum_net_time+=self.__net_time
        if name == "residual_norm":
            self.__residual_norm=float(value)
        if name == "converged":
            self.__converged = (value == True)
        if name == "time_step_backtracking_used":
            self.__time_step_backtracking_used = (value == True)
        if name == "coarse_level_sparsity":
            self.__coarse_level_sparsity = value
        if name == "num_coarse_unknowns": 
            self.__num_coarse_unknowns= value
            
    def getDiagnostics(self, name):
        """
        Returns the diagnostic information ``name``. Possible values are:
            
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
        - "converged": return True if solution has converged.
        - "time_step_backtracking_used": returns True if time step back tracking has been used.
        - "coarse_level_sparsity": returns the sparsity of the matrix on the coarsest level
        - "num_coarse_unknowns": returns the number of unknowns on the coarsest level
        
        
        :param name: name of diagnostic information to return
        :type name: ``str`` in the list above.
        :return: requested value. ``None`` is returned if the value is undefined.
        :note: If the solver has thrown an exception diagnostic values have an undefined status.
        """
        if name == "num_iter": return self.__num_iter
        elif name == "cum_num_iter": return self.__cum_num_iter
        elif name == "num_level": return self.__num_level
        elif name == "num_inner_iter": return self.__num_inner_iter
        elif name == "cum_num_inner_iter": return self.__cum_num_inner_iter
        elif name == "time": return self.__time
        elif name == "cum_time": return self.__cum_time
        elif name == "set_up_time": return self.__set_up_time
        elif name == "cum_set_up_time": return self.__cum_set_up_time
        elif name == "net_time": return self.__net_time
        elif name == "cum_net_time": return self.__cum_net_time
        elif name == "residual_norm": return self.__residual_norm
        elif name == "converged": return self.__converged      
        elif name == "preconditioner_size": return  self.__preconditioner_size
        elif name == "time_step_backtracking_used": return  self.__time_step_backtracking_used
        elif name == "coarse_level_sparsity": return self.__coarse_level_sparsity
        elif name == "num_coarse_unknowns": return self.__num_coarse_unknowns
        else:
            raise ValueError("unknown diagnostic item %s"%name)
    def hasConverged(self):
        """
        Returns ``True`` if the last solver call has been finalized successfully.
        :note: if an exception has been thrown by the solver the status of this flag is undefined.
        """
        return self.getDiagnostics("converged")
    def setCoarsening(self,method=0):
        """
        Sets the key of the coarsening method to be applied in AMG or AMLI or BoomerAMG

        :param method: selects the coarsening method .
        :type method: in {SolverOptions.DEFAULT}, `SolverOptions.YAIR_SHAPIRA_COARSENING`,  `SolverOptions.RUGE_STUEBEN_COARSENING`, `SolverOptions.AGGREGATION_COARSENING`, `SolverOptions.CIJP_FIXED_RANDOM_COARSENING`, `SolverOptions.CIJP_COARSENING`, `SolverOptions.FALGOUT_COARSENING`, `SolverOptions.PMIS_COARSENING`, `SolverOptions.HMIS_COARSENING`
        """
        if method==None: method=0
        if not method in [self.DEFAULT, self.YAIR_SHAPIRA_COARSENING, self.RUGE_STUEBEN_COARSENING, self.AGGREGATION_COARSENING, self.STANDARD_COARSENING, self.CIJP_FIXED_RANDOM_COARSENING, self.CIJP_COARSENING, self.FALGOUT_COARSENING, self.PMIS_COARSENING, self.HMIS_COARSENING,]:
             raise ValueError("unknown coarsening method %s"%method)
        self.__coarsening=method
    
    def getCoarsening(self):
        """
        Returns the key of the coarsening algorithm to be applied AMG or AMLI or BoomerAMG

        :rtype: in the list `SolverOptions.DEFAULT`, `SolverOptions.YAIR_SHAPIRA_COARSENING`, `SolverOptions.RUGE_STUEBEN_COARSENING`, `SolverOptions.AGGREGATION_COARSENING`, `SolverOptions.CIJP_FIXED_RANDOM_COARSENING`, `SolverOptions.CIJP_COARSENING`, `SolverOptions.FALGOUT_COARSENING`, `SolverOptions.PMIS_COARSENING`, `SolverOptions.HMIS_COARSENING`
        """
        return self.__coarsening
      
    def setMinCoarseMatrixSize(self,size=None):
        """
        Sets the minumum size of the coarsest level matrix in AMG or AMLI

        :param size: minumum size of the coarsest level matrix .
        :type size: positive ``int`` or ``None``
        """
        if size==None: size=500
        size=int(size)
        if size<0:
           raise ValueError("minimum size of the coarsest level matrix must be non-negative.")
        self.__MinCoarseMatrixSize=size
        
    def getMinCoarseMatrixSize(self):
        """
        Returns the minimum size of the coarsest level matrix in AMG or AMLI

        :rtype: ``int``
        """
        return self.__MinCoarseMatrixSize
      
    def setPreconditioner(self, preconditioner=10):
        """
        Sets the preconditioner to be used. 

        :param preconditioner: key of the preconditioner to be used.
        :type preconditioner: in `SolverOptions.ILU0`, `SolverOptions.ILUT`,
            `SolverOptions.JACOBI`, `SolverOptions.AMG`, `SolverOptions.AMLI`,
            `SolverOptions.REC_ILU`, `SolverOptions.GAUSS_SEIDEL`,
            `SolverOptions.RILU`, `SolverOptions.BOOMERAMG`,
            `SolverOptions.NO_PRECONDITIONER`
        :note: Not all packages support all preconditioner. It can be assumed that a package makes a reasonable choice if it encounters an unknown preconditioner. 
        """
        if preconditioner==None: preconditioner=10
        if not preconditioner in [  SolverOptions.ILU0, SolverOptions.ILUT, SolverOptions.JACOBI, 
                                    SolverOptions.AMG, SolverOptions.AMLI, SolverOptions.REC_ILU, SolverOptions.GAUSS_SEIDEL, SolverOptions.RILU, SolverOptions.BOOMERAMG,
                                    SolverOptions.NO_PRECONDITIONER] :
             raise ValueError("unknown preconditioner %s"%preconditioner)
        if preconditioner==SolverOptions.AMG and escore.getEscriptParamInt('DISABLE_AMG',0):
             raise ValueError("AMG preconditioner is not supported in MPI builds")
        self.__preconditioner=preconditioner    

    def getPreconditioner(self):
        """
        Returns the key of the preconditioner to be used. 

        :rtype: in the list `SolverOptions.ILU0`, `SolverOptions.ILUT`,
            `SolverOptions.JACOBI`, `SolverOptions.AMLI`, `SolverOptions.AMG`,
            `SolverOptions.REC_ILU`, `SolverOptions.GAUSS_SEIDEL`,
            `SolverOptions.RILU`, `SolverOptions.BOOMERAMG`,
            `SolverOptions.NO_PRECONDITIONER`
        """
        return self.__preconditioner

    def setSmoother(self, smoother=None):
        """
        Sets the smoother to be used. 

        :param smoother: key of the smoother to be used.
        :type smoother: in `SolverOptions.JACOBI`, `SolverOptions.GAUSS_SEIDEL`
        :note: Not all packages support all smoothers. It can be assumed that a package makes a reasonable choice if it encounters an unknown smoother. 
        """
        if smoother==None: smoother=28
        if not smoother in [ SolverOptions.JACOBI, SolverOptions.GAUSS_SEIDEL ] :
             raise ValueError("unknown smoother %s"%smoother)
        self.__smoother=smoother

    def getSmoother(self):
        """
        Returns key of the smoother to be used. 

        :rtype: in the list `SolverOptions.JACOBI`, `SolverOptions.GAUSS_SEIDEL`
        """
        return self.__smoother  

    def setSolverMethod(self, method=0):
        """
        Sets the solver method to be used. Use ``method``=``SolverOptions.DIRECT`` to indicate that a direct rather than an iterative
        solver should be used and Use ``method``=``SolverOptions.ITERATIVE`` to indicate that an iterative rather than a direct
        solver should be used. 

        :param method: key of the solver method to be used.
        :type method: in `SolverOptions.DEFAULT`, `SolverOptions.DIRECT`, `SolverOptions.CHOLEVSKY`, `SolverOptions.PCG`, 
                        `SolverOptions.CR`, `SolverOptions.CGS`, `SolverOptions.BICGSTAB`, 
                        `SolverOptions.GMRES`, `SolverOptions.PRES20`, `SolverOptions.ROWSUM_LUMPING`, `SolverOptions.HRZ_LUMPING`, `SolverOptions.ITERATIVE`, 
                        `SolverOptions.NONLINEAR_GMRES`, `SolverOptions.TFQMR`, `SolverOptions.MINRES`
        :note: Not all packages support all solvers. It can be assumed that a package makes a reasonable choice if it encounters an unknown solver method. 
        """
        if method==None: method=0
        if not method in [ SolverOptions.DEFAULT, SolverOptions.DIRECT, SolverOptions.CHOLEVSKY, SolverOptions.PCG, 
                           SolverOptions.CR, SolverOptions.CGS, SolverOptions.BICGSTAB, 
                           SolverOptions.GMRES, SolverOptions.PRES20, SolverOptions.ROWSUM_LUMPING, SolverOptions.HRZ_LUMPING,
                           SolverOptions.ITERATIVE,
                           SolverOptions.NONLINEAR_GMRES, SolverOptions.TFQMR, SolverOptions.MINRES ]:
             raise ValueError("unknown solver method %s"%method)
        self.__method=method
    def getSolverMethod(self):
        """
        Returns key of the solver method to be used. 

        :rtype: in the list `SolverOptions.DEFAULT`, `SolverOptions.DIRECT`, `SolverOptions.CHOLEVSKY`, `SolverOptions.PCG`, 
                        `SolverOptions.CR`, `SolverOptions.CGS`, `SolverOptions.BICGSTAB`, 
                        `SolverOptions.GMRES`, `SolverOptions.PRES20`, `SolverOptions.ROWSUM_LUMPING`, `SolverOptions.HRZ_LUMPING`, `SolverOptions.ITERATIVE`, 
                        `SolverOptions.NONLINEAR_GMRES`, `SolverOptions.TFQMR`, `SolverOptions.MINRES`
        """
        return self.__method
        
    def setPackage(self, package=0):
        """
        Sets the solver package to be used as a solver.  

        :param package: key of the solver package to be used.
        :type package: in `SolverOptions.DEFAULT`, `SolverOptions.PASO`, `SolverOptions.SUPER_LU`, `SolverOptions.PASTIX`, `SolverOptions.MKL`, `SolverOptions.UMFPACK`, `SolverOptions.TRILINOS`
        :note: Not all packages are support on all implementation. An exception may be thrown on some platforms if a particular is requested. 
        """
        if package==None: package=0
        if not package in [SolverOptions.DEFAULT, SolverOptions.PASO, SolverOptions.SUPER_LU, SolverOptions.PASTIX, SolverOptions.MKL, SolverOptions.UMFPACK, SolverOptions.TRILINOS]:
             raise ValueError("unknown solver package %s"%package)
        self.__package=package
    def getPackage(self):
        """
        Returns the solver package key

        :rtype: in the list `SolverOptions.DEFAULT`, `SolverOptions.PASO`, `SolverOptions.SUPER_LU`, `SolverOptions.PASTIX`, `SolverOptions.MKL`, `SolverOptions.UMFPACK`, `SolverOptions.TRILINOS`
        """
        return self.__package
    def setReordering(self,ordering=30):
        """
        Sets the key of the reordering method to be applied if supported by the solver. Some direct solvers support reordering 
        to optimize compute time and storage use during elimination. 

        :param ordering: selects the reordering strategy.
        :type ordering: in 'SolverOptions.NO_REORDERING', 'SolverOptions.MINIMUM_FILL_IN', 'SolverOptions.NESTED_DISSECTION', 'SolverOptions.DEFAULT_REORDERING'
        """
        if not ordering in [self.NO_REORDERING, self.MINIMUM_FILL_IN, self.NESTED_DISSECTION, self.DEFAULT_REORDERING]:
             raise ValueError("unknown reordering strategy %s"%ordering)
        self.__reordering=ordering

    def getReordering(self):
        """
        Returns the key of the reordering method to be applied if supported by the solver.

        :rtype: in `SolverOptions.NO_REORDERING`,
            `SolverOptions.MINIMUM_FILL_IN`, `SolverOptions.NESTED_DISSECTION`,
            `SolverOptions.DEFAULT_REORDERING`
        """
        return self.__reordering

    def setRestart(self,restart=None):
        """
        Sets the number of iterations steps after which GMRES perfors a restart.

        :param restart: number of iteration steps after which to perform a
            restart. If ``None`` no restart is performed.
        :type restart: ``int`` or ``None``
        """
        if restart == None:
            self.__restart=restart
        else:
            restart=int(restart)
            if restart<1:
                raise ValueError("restart must be positive.")
            self.__restart=restart
            
    def getRestart(self):
        """
        Returns the number of iterations steps after which GMRES is performing a restart.
        If ``None`` is returned no restart is performed.

        :rtype: ``int`` or ``None``
        """
        if (self.__restart is None) or (self.__restart < 0):
            return None
        else:
            return self.__restart
    def _getRestartForC(self):
            r=self.getRestart()
            if r == None: 
               return -1
            else: 
               return r
   
    def setDiagonalDominanceThreshold(self,value=0.5):
        """
        Sets the threshold for diagonal dominant rows which are eliminated during AMG coarsening. 

        :param value: threshold
        :type value: ``float``
        """
        value=float(value)
        if value<0 or value>1.:
           raise ValueError("Diagonal dominance threshold must be between 0 and 1.")
        self.__diagonal_dominance_threshold=value
        
    def getDiagonalDominanceThreshold(self):
        """
        Returns the threshold for diagonal dominant rows which are eliminated during AMG coarsening. 

        :rtype: ``float``
        """
        return self.__diagonal_dominance_threshold
        
    def setTruncation(self,truncation=20):
        """
        Sets the number of residuals in GMRES to be stored for orthogonalization.  The more residuals are stored
        the faster GMRES converged but

        :param truncation: truncation
        :type truncation: ``int``
        """
        truncation=int(truncation)
        if truncation<1:
           raise ValueError("truncation must be positive.")
        self.__truncation=truncation
    def getTruncation(self):
        """
        Returns the number of residuals in GMRES to be stored for orthogonalization

        :rtype: ``int``
        """
        return self.__truncation
    def setInnerIterMax(self,iter_max=10):
        """
        Sets the maximum number of iteration steps for the inner iteration.

        :param iter_max: maximum number of inner iterations
        :type iter_max: ``int``
        """
        iter_max=int(iter_max)
        if iter_max<1:
           raise ValueError("maximum number of inner iteration must be positive.")
        self.__inner_iter_max=iter_max
    def getInnerIterMax(self):
        """
        Returns maximum number of inner iteration steps

        :rtype: ``int``
        """
        return self.__inner_iter_max
    def setIterMax(self,iter_max=100000):
        """
        Sets the maximum number of iteration steps

        :param iter_max: maximum number of iteration steps
        :type iter_max: ``int``
        """
        iter_max=int(iter_max)
        if iter_max<1:
           raise ValueError("maximum number of iteration steps must be positive.")
        self.__iter_max=iter_max
    def getIterMax(self):
        """
        Returns maximum number of iteration steps

        :rtype: ``int``
        """
        return self.__iter_max
    def setLevelMax(self,level_max=100):
        """
        Sets the maximum number of coarsening levels to be used in an algebraic multi level solver or preconditioner

        :param level_max: maximum number of levels
        :type level_max: ``int``
        """
        level_max=int(level_max)
        if level_max<0:
           raise ValueError("maximum number of coarsening levels must be non-negative.")
        self.__level_max=level_max
    def getLevelMax(self):
        """
        Returns the maximum number of coarsening levels to be used in an algebraic multi level solver or preconditioner

        :rtype: ``int``
        """
        return self.__level_max
    def setCycleType(self, cycle_type=1): 
        """
        Sets the cycle type (V-cycle or W-cycle) to be used in an algebraic multi level solver or preconditioner

        :param cycle_type: the type of cycle
        :type cycle_type: ``int``
        """
        cycle_type=int(cycle_type)
        self.__cycle_type=cycle_type
    def getCycleType(self):
        """
        Returns the cyle type (V- or W-cycle) to be used in an algebraic multi level solver or preconditioner

        :rtype: ``int``
        """
        return self.__cycle_type
    def setCoarseningThreshold(self,theta=0.25):
        """
        Sets the threshold for coarsening in the algebraic multi level solver or preconditioner

        :param theta: threshold for coarsening
        :type theta: positive ``float``
        """
        theta=float(theta)
        if theta<0 or theta>1:
           raise ValueError("threshold must be non-negative and less or equal 1.")
        self.__coarsening_threshold=theta
    def getCoarseningThreshold(self):
        """
        Returns the threshold for coarsening in the algebraic multi level solver or preconditioner

        :rtype: ``float``
        """
        return self.__coarsening_threshold
    def setNumSweeps(self,sweeps=1):
        """
        Sets the number of sweeps in a Jacobi or Gauss-Seidel/SOR preconditioner.

        :param sweeps: number of sweeps
        :type sweeps: positive ``int``
        """
        sweeps=int(sweeps)
        if sweeps<1:
           raise ValueError("number of sweeps must be positive.")
        self.__sweeps=sweeps
    def getNumSweeps(self):
        """
        Returns the number of sweeps in a Jacobi or Gauss-Seidel/SOR preconditioner.

        :rtype: ``int``
        """
        return self.__sweeps
    def setNumPreSweeps(self,sweeps=1):
        """
        Sets the number of sweeps in the pre-smoothing step of a multi level solver or preconditioner

        :param sweeps: number of sweeps
        :type sweeps: positive ``int``
        """
        sweeps=int(sweeps)
        if sweeps<1:
           raise ValueError("number of sweeps must be positive.")
        self.__pre_sweeps=sweeps
    def getNumPreSweeps(self):
        """
        Returns he number of sweeps in the pre-smoothing step of a multi level solver or preconditioner

        :rtype: ``int``
        """
        return self.__pre_sweeps
    def setNumPostSweeps(self,sweeps=1):
        """
        Sets the number of sweeps in the post-smoothing step of a multi level solver or preconditioner

        :param sweeps: number of sweeps
        :type sweeps: positive ``int``
        """
        sweeps=int(sweeps)
        if sweeps<1:
           raise ValueError("number of sweeps must be positive.")
        self.__post_sweeps=sweeps
    def getNumPostSweeps(self):
        """
        Returns he number of sweeps in the post-smoothing step of a multi level solver or preconditioner

        :rtype: ``int``
        """
        return self.__post_sweeps

    def setTolerance(self,rtol=1.e-8):
        """
        Sets the relative tolerance for the solver

        :param rtol: relative tolerance
        :type rtol: non-negative ``float``
        """
        rtol=float(rtol)
        if rtol<0 or rtol>1:
           raise ValueError("tolerance must be non-negative and less or equal 1.")
        self.__tolerance=rtol
    def getTolerance(self):
        """
        Returns the relative tolerance for the solver

        :rtype: ``float``
        """
        return self.__tolerance
    def setAbsoluteTolerance(self,atol=0.):
        """
        Sets the absolute tolerance for the solver

        :param atol:  absolute tolerance
        :type atol: non-negative ``float``
        """
        atol=float(atol)
        if atol<0:
           raise ValueError("tolerance must be non-negative.")
        self.__absolute_tolerance=atol
    def getAbsoluteTolerance(self):
        """
        Returns the absolute tolerance for the solver

        :rtype: ``float``
        """
        return self.__absolute_tolerance

    def setInnerTolerance(self,rtol=0.9):
        """
         Sets the relative tolerance for an inner iteration scheme for instance on the coarsest level in a multi-level scheme.

        :param rtol: inner relative tolerance
        :type rtol: positive ``float``
        """
        rtol=float(rtol)
        if rtol<=0 or rtol>1:
            raise ValueError("tolerance must be positive and less or equal 1.")
        self.__inner_tolerance=rtol
    def getInnerTolerance(self):
        """
        Returns the relative tolerance for an inner iteration scheme

        :rtype: ``float``
        """
        return self.__inner_tolerance
    def setDropTolerance(self,drop_tol=0.01):
        """
        Sets the relative drop tolerance in ILUT

        :param drop_tol: drop tolerance
        :type drop_tol: positive ``float``
        """
        drop_tol=float(drop_tol)
        if drop_tol<=0 or drop_tol>1:
            raise ValueError("drop tolerance must be positive and less or equal 1.")
        self.__drop_tolerance=drop_tol
    def getDropTolerance(self):
        """
        Returns the relative drop tolerance in ILUT

        :rtype: ``float``
        """
        return self.__drop_tolerance
    def setDropStorage(self,storage=2.):
        """
        Sets the maximum allowed increase in storage for ILUT. ``storage`` =2 would mean that
        a doubling of the storage needed for the coefficient matrix is allowed in the ILUT factorization.

        :param storage: allowed storage increase
        :type storage: ``float``
        """
        storage=float(storage)
        if storage<1:
            raise ValueError("allowed storage increase must be greater or equal to 1.")
        self.__drop_storage=storage
    def getDropStorage(self):
    
        """
        Returns the maximum allowed increase in storage for ILUT

        :rtype: ``float``
        """
        return self.__drop_storage
    def setRelaxationFactor(self,factor=0.3):
        """
        Sets the relaxation factor used to add dropped elements in RILU to the main diagonal.

        :param factor: relaxation factor
        :type factor: ``float``
        :note: RILU with a relaxation factor 0 is identical to ILU0
        """
        factor=float(factor)
        if factor<0: 
            raise ValueError("relaxation factor must be non-negative.")
        self.__relaxation=factor
    def getRelaxationFactor(self):
    
        """
        Returns the relaxation factor used to add dropped elements in RILU to the main diagonal.

        :rtype: ``float``
        """
        return self.__relaxation
    def isSymmetric(self):
        """
        Checks if symmetry of the  coefficient matrix is indicated.

        :return: True if a symmetric PDE is indicated, False otherwise
        :rtype: ``bool``
        """
        return self.__symmetric
    def setSymmetryOn(self):
        """
        Sets the symmetry flag to indicate that the coefficient matrix is symmetric.
        """
        self.__symmetric=True
    def setSymmetryOff(self):
        """
        Clears the symmetry flag for the coefficient matrix.
        """
        self.__symmetric=False
    def setSymmetry(self,flag=False):
        """
        Sets the symmetry flag for the coefficient matrix to ``flag``.

        :param flag: If True, the symmetry flag is set otherwise reset.
        :type flag: ``bool``
        """
        if flag:
            self.setSymmetryOn()
        else:
            self.setSymmetryOff()
    def isVerbose(self):
        """
        Returns ``True`` if the solver is expected to be verbose.

        :return: True if verbosity of switched on.
        :rtype: ``bool``
        """
        return self.__verbose

    def setVerbosityOn(self):
        """
        Switches the verbosity of the solver on.
        """
        self.__verbose=True
    def setVerbosityOff(self):
        """
        Switches the verbosity of the solver off.
        """
        self.__verbose=False
    def setVerbosity(self,verbose=False):
        """
        Sets the verbosity flag for the solver to ``flag``.

        :param verbose: If ``True``, the verbosity of the solver is switched on.
        :type verbose: ``bool``
        """
        if verbose:
            self.setVerbosityOn()
        else:
            self.setVerbosityOff()
            
    def adaptInnerTolerance(self):
        """
        Returns ``True`` if the tolerance of the inner solver is selected automatically. 
        Otherwise the inner tolerance set by `setInnerTolerance` is used.

        :return: ``True`` if inner tolerance adaption is chosen.
        :rtype: ``bool``
        """
        return self.__adapt_inner_tolerance

    def setInnerToleranceAdaptionOn(self):
        """
        Switches the automatic selection of inner tolerance on 
        """
        self.__adapt_inner_tolerance=True
    def setInnerToleranceAdaptionOff(self):
        """
        Switches the automatic selection of inner tolerance off.
        """
        self.__adapt_inner_tolerance=False
    def setInnerToleranceAdaption(self,adapt=True):
        """
        Sets the flag to indicate automatic selection of the inner tolerance. 

        :param adapt: If ``True``, the inner tolerance is selected automatically.
        :type adapt: ``bool``
        """
        if adapt:
            self.setInnerToleranceAdaptionOn()
        else:
            self.setInnerToleranceAdaptionOff()

    def acceptConvergenceFailure(self):
        """
        Returns ``True`` if a failure to meet the stopping criteria within the
        given number of iteration steps is not raising in exception. This is useful 
        if a solver is used in a non-linear context where the non-linear solver can 
        continue even if the returned the solution does not necessarily meet the
        stopping criteria. One can use the `hasConverged` method to check if the
        last call to the solver was successful.

        :return: ``True`` if a failure to achieve convergence is accepted.
        :rtype: ``bool``
        """
        return self.__accept_convergence_failure

    def setAcceptanceConvergenceFailureOn(self):
        """
        Switches the acceptance of a failure of convergence on 
        """
        self.__accept_convergence_failure=True
    def setAcceptanceConvergenceFailureOff(self):
        """
        Switches the acceptance of a failure of convergence off.
        """
        self.__accept_convergence_failure=False
    def setAcceptanceConvergenceFailure(self,accept=False):
        """
        Sets the flag to indicate the acceptance of a failure of convergence. 

        :param accept: If ``True``, any failure to achieve convergence is accepted.
        :type accept: ``bool``
        """
        if accept:
            self.setAcceptanceConvergenceFailureOn()
        else:
            self.setAcceptanceConvergenceFailureOff()

    def useLocalPreconditioner(self):
        """
        Returns ``True`` if the preconditoner is applied locally on each MPI. This reducess communication costs
        and speeds up the application of the preconditioner but at the costs of more iteration steps. This can be an 
        advantage on clusters with slower interconnects.

        :return: ``True`` if local preconditioning is applied
        :rtype: ``bool``
        """
        return self.__use_local_preconditioner

    def setLocalPreconditionerOn(self):
        """
        Sets the flag to use  local preconditioning to on 
        """
        self.__use_local_preconditioner=True
    def setLocalPreconditionerOff(self):
        """
        Sets the flag to use  local preconditioning to off
        """
        self.__use_local_preconditioner=False
        
    def setLocalPreconditioner(self,use=False):
        """
        Sets the flag to use  local preconditioning

        :param use: If ``True``, local proconditioning on each MPI rank is applied
        :type use: ``bool``
        """
        if use:
            self.setUseLocalPreconditionerOn()
        else:
            self.setUseLocalPreconditionerOff()
            
    def setMinCoarseMatrixSparsity(self,sparsity=0.05):
       """
       Sets the minimum sparsity on the coarsest level. Typically 
       a direct solver is used when the sparsity becomes bigger than 
       the set limit.
       
       :param sparsity: minimual sparsity
       :type sparsity: ``float``
       """
       sparsity=float(sparsity)
       if sparsity<0: 
          raise ValueError("sparsity must be non-negative.")
       if sparsity>1: 
          raise ValueError("sparsity must be less than 1.")
       self.__min_sparsity=sparsity

    def getMinCoarseMatrixSparsity(self):
      """
      Returns the minimum sparsity on the coarsest level. Typically 
      a direct solver is used when the sparsity becomes bigger than 
      the set limit.
   
      :return: minimual sparsity
      :rtype: ``float``
      """
      return self.__min_sparsity

    def setNumRefinements(self,refinements=2):
      """
      Sets the number of refinement steps to refine the solution when a direct solver is applied.
   
      :param refinements: number of refinements
      :type refinements: non-negative ``int``
      """
      refinements=int(refinements)
      if refinements<0:
        raise ValueError("number of refinements must be non-negative.")
      self.__refinements=refinements
   
    def getNumRefinements(self):
       """
       Returns the number of refinement steps to refine the solution when a direct solver is applied.
       
       :rtype: non-negative ``int``
       """
       return self.__refinements

    def setNumCoarseMatrixRefinements(self,refinements=2):
      """
      Sets the number of refinement steps to refine the solution on the coarset level when a direct solver is applied.
   
      :param refinements: number of refinements
      :type refinements: non-negative ``int``
      """
      refinements=int(refinements)
      if refinements<0:
         raise ValueError("number of refinements must be non-negative.")
      self.__coarse_refinements=refinements
   
    def getNumCoarseMatrixRefinements(self):
      """
      Returns the number of resfinement steps to refine the solution on the coarset level when a direct solver is applied.
      
      :rtype: non-negative ``int``
      """
      return self.__coarse_refinements

    def usePanel(self):
        """
        Returns ``True`` if a panel is used to serach for unknown in the AMG coarsening, The panel approach is normally faster
        but can lead to larger coarse level systems.
        
        :return: ``True`` if a panel is used to find unknowns in AMG coarsening
        :rtype: ``bool``
        """
        return self.__use_panel

    def setUsePanelOn(self):
        """
        Sets the flag to use a panel to find unknowns in AMG coarsening 
        """
        self.__use_panel=True
        
    def setUsePanelOff(self):
        """
        Sets the flag to use a panel to find unknowns in AMG coarsening to off
        """
        self.__use_panel=False
        
    def setUsePanel(self,use=True):
        """
        Sets the flag to use  a panel to find unknowns in AMG coarsening 

        :param use: If ``True``,a panel is used to find unknowns in AMG coarsening
        :type use: ``bool``
        """
        if use:
            self.setUsePanelOn()
        else:
            self.setUsePanelOff()
            
    def setAMGInterpolation(self, method=None):
        """
        Set the interpolation method for the AMG preconditioner. 

        :param method: key of the interpolation method to be used.
        :type method: in `SolverOptions.CLASSIC_INTERPOLATION_WITH_FF_COUPLING`, `SolverOptions.CLASSIC_INTERPOLATION`, `SolverOptions.DIRECT_INTERPOLATION`
        """
        if method==None: method=self.DIRECT_INTERPOLATION
        if not method in [ SolverOptions.CLASSIC_INTERPOLATION_WITH_FF_COUPLING, SolverOptions.CLASSIC_INTERPOLATION, SolverOptions.DIRECT_INTERPOLATION ]:
             raise ValueError("unknown AMG interpolation method %s"%method)
        self.__amg_interpolation_method=method

    def getAMGInterpolation(self):
        """
        Returns key of the interpolation method for the AMG preconditioner

        :rtype: in the list `SolverOptions.CLASSIC_INTERPOLATION_WITH_FF_COUPLING`, `SolverOptions.CLASSIC_INTERPOLATION`, `SolverOptions.DIRECT_INTERPOLATION`
        """
        return self.__amg_interpolation_method
    
    def setODESolver(self, method=None):
        """
        Set the solver method for ODEs. 

        :param method: key of the ODE solver method to be used.
        :type method: in `SolverOptions.CRANK_NICOLSON`, `SolverOptions.BACKWARD_EULER`, `SolverOptions.LINEAR_CRANK_NICOLSON`
        """
        if method is None: method=self.LINEAR_CRANK_NICOLSON
        if not method in [ SolverOptions.CRANK_NICOLSON, SolverOptions.BACKWARD_EULER, SolverOptions.LINEAR_CRANK_NICOLSON ]:
             raise ValueError("unknown ODE solver method %s"%method)
        self.__ode_solver=method

    def getODESolver(self, method=None):
        """
        Returns key of the solver method for ODEs. 

        :param method: key of the ODE solver method to be used.
        :type method: in `SolverOptions.CRANK_NICOLSON`, `SolverOptions.BACKWARD_EULER`, `SolverOptions.LINEAR_CRANK_NICOLSON`
        """
        return self.__ode_solver
        
        
class IllegalCoefficient(ValueError):
  
   """
   Exception that is raised if an illegal coefficient of the general or
   particular PDE is requested.
   """
   pass

class IllegalCoefficientValue(ValueError):
   """
   Exception that is raised if an incorrect value for a coefficient is used.
   """
   pass

class IllegalCoefficientFunctionSpace(ValueError):
   """
   Exception that is raised if an incorrect function space for a coefficient
   is used.
   """

class UndefinedPDEError(ValueError):
   """
   Exception that is raised if a PDE is not fully defined yet.
   """
   pass

class PDECoef(object):
    """
    A class for describing a PDE coefficient.

    :cvar INTERIOR: indicator that coefficient is defined on the interior of
                    the PDE domain
    :cvar BOUNDARY: indicator that coefficient is defined on the boundary of
                    the PDE domain
    :cvar CONTACT: indicator that coefficient is defined on the contact region
                   within the PDE domain
    :cvar INTERIOR_REDUCED: indicator that coefficient is defined on the
                            interior of the PDE domain using a reduced
                            integration order
    :cvar BOUNDARY_REDUCED: indicator that coefficient is defined on the
                            boundary of the PDE domain using a reduced
                            integration order
    :cvar CONTACT_REDUCED: indicator that coefficient is defined on the contact
                           region within the PDE domain using a reduced
                           integration order
    :cvar SOLUTION: indicator that coefficient is defined through a solution of
                    the PDE
    :cvar REDUCED: indicator that coefficient is defined through a reduced
                   solution of the PDE
    :cvar DIRACDELTA: indicator that coefficient is defined as Dirac delta functions
    :cvar BY_EQUATION: indicator that the dimension of the coefficient shape is
                       defined by the number of PDE equations
    :cvar BY_SOLUTION: indicator that the dimension of the coefficient shape is
                       defined by the number of PDE solutions
    :cvar BY_DIM: indicator that the dimension of the coefficient shape is
                  defined by the spatial dimension
    :cvar OPERATOR: indicator that the coefficient alters the operator of
                    the PDE
    :cvar RIGHTHANDSIDE: indicator that the coefficient alters the right
                         hand side of the PDE
    :cvar BOTH: indicator that the coefficient alters the operator as well
                as the right hand side of the PDE

    """
    INTERIOR=0
    BOUNDARY=1
    CONTACT=2
    SOLUTION=3
    REDUCED=4
    BY_EQUATION=5
    BY_SOLUTION=6
    BY_DIM=7
    OPERATOR=10
    RIGHTHANDSIDE=11
    BOTH=12
    INTERIOR_REDUCED=13
    BOUNDARY_REDUCED=14
    CONTACT_REDUCED=15
    DIRACDELTA=16

    def __init__(self, where, pattern, altering):
       """
       Initialises a PDE coefficient type.

       :param where: describes where the coefficient lives
       :type where: one of `INTERIOR`, `BOUNDARY`, `CONTACT`, `SOLUTION`,
                    `REDUCED`, `INTERIOR_REDUCED`, `BOUNDARY_REDUCED`,
                    `CONTACT_REDUCED`, 'DIRACDELTA'
       :param pattern: describes the shape of the coefficient and how the shape
                       is built for a given spatial dimension and numbers of
                       equations and solutions in then PDE. For instance,
                       (`BY_EQUATION`,`BY_SOLUTION`,`BY_DIM`) describes a
                       rank 3 coefficient which is instantiated as shape (3,2,2)
                       in case of three equations and two solution components
                       on a 2-dimensional domain. In the case of single equation
                       and a single solution component the shape components
                       marked by `BY_EQUATION` or `BY_SOLUTION` are dropped.
                       In this case the example would be read as (2,).
       :type pattern: ``tuple`` of `BY_EQUATION`, `BY_SOLUTION`, `BY_DIM`
       :param altering: indicates what part of the PDE is altered if the
                        coefficient is altered
       :type altering: one of `OPERATOR`, `RIGHTHANDSIDE`, `BOTH`
       """
       super(PDECoef, self).__init__()
       self.what=where
       self.pattern=pattern
       self.altering=altering
       self.resetValue()

    def resetValue(self):
       """
       Resets the coefficient value to the default.
       """
       self.value=escore.Data()

    def getFunctionSpace(self,domain,reducedEquationOrder=False,reducedSolutionOrder=False):
       """
       Returns the `FunctionSpace` of the coefficient.

       :param domain: domain on which the PDE uses the coefficient
       :type domain: `Domain`
       :param reducedEquationOrder: True to indicate that reduced order is used
                                    to represent the equation
       :type reducedEquationOrder: ``bool``
       :param reducedSolutionOrder: True to indicate that reduced order is used
                                    to represent the solution
       :type reducedSolutionOrder: ``bool``
       :return: `FunctionSpace` of the coefficient
       :rtype: `FunctionSpace`
       """
       if self.what==self.INTERIOR:
            return escore.Function(domain)
       elif self.what==self.INTERIOR_REDUCED:
            return escore.ReducedFunction(domain)
       elif self.what==self.BOUNDARY:
            return escore.FunctionOnBoundary(domain)
       elif self.what==self.BOUNDARY_REDUCED:
            return escore.ReducedFunctionOnBoundary(domain)
       elif self.what==self.CONTACT:
            return escore.FunctionOnContactZero(domain)
       elif self.what==self.CONTACT_REDUCED:
            return escore.ReducedFunctionOnContactZero(domain)
       elif self.what==self.DIRACDELTA:
            return escore.DiracDeltaFunctions(domain)
       elif self.what==self.SOLUTION:
            if reducedEquationOrder and reducedSolutionOrder:
                return escore.ReducedSolution(domain)
            else:
                return escore.Solution(domain)
       elif self.what==self.REDUCED:
            return escore.ReducedSolution(domain)

    def getValue(self):
       """
       Returns the value of the coefficient.

       :return: value of the coefficient
       :rtype: `Data`
       """
       return self.value

    def setValue(self,domain,numEquations=1,numSolutions=1,reducedEquationOrder=False,reducedSolutionOrder=False,newValue=None):
       """
       Sets the value of the coefficient to a new value.

       :param domain: domain on which the PDE uses the coefficient
       :type domain: `Domain`
       :param numEquations: number of equations of the PDE
       :type numEquations: ``int``
       :param numSolutions: number of components of the PDE solution
       :type numSolutions: ``int``
       :param reducedEquationOrder: True to indicate that reduced order is used
                                    to represent the equation
       :type reducedEquationOrder: ``bool``
       :param reducedSolutionOrder: True to indicate that reduced order is used
                                    to represent the solution
       :type reducedSolutionOrder: ``bool``
       :param newValue: number of components of the PDE solution
       :type newValue: any object that can be converted into a
                       `Data` object with the appropriate shape
                       and `FunctionSpace`
       :raise IllegalCoefficientValue: if the shape of the assigned value does
                                       not match the shape of the coefficient
       :raise IllegalCoefficientFunctionSpace: if unable to interpolate value
                                               to appropriate function space
       """
       if newValue==None:
           newValue=escore.Data()
       elif isinstance(newValue,escore.Data):
           if not newValue.isEmpty():
              if not newValue.getFunctionSpace() == self.getFunctionSpace(domain,reducedEquationOrder,reducedSolutionOrder):
                try:
                  newValue=escore.Data(newValue,self.getFunctionSpace(domain,reducedEquationOrder,reducedSolutionOrder))
                except:
                  raise IllegalCoefficientFunctionSpace("Unable to interpolate coefficient to function space %s"%self.getFunctionSpace(domain))
       else:
           newValue=escore.Data(newValue,self.getFunctionSpace(domain,reducedEquationOrder,reducedSolutionOrder))
       if not newValue.isEmpty():
           if not self.getShape(domain,numEquations,numSolutions)==newValue.getShape():
               raise IllegalCoefficientValue("Expected shape of coefficient is %s but actual shape is %s."%(self.getShape(domain,numEquations,numSolutions),newValue.getShape()))
       self.value=newValue

    def isAlteringOperator(self):
        """
        Checks if the coefficient alters the operator of the PDE.

        :return: True if the operator of the PDE is changed when the
                 coefficient is changed
        :rtype: ``bool``
        """
        if self.altering==self.OPERATOR or self.altering==self.BOTH:
            return not None
        else:
            return None

    def isAlteringRightHandSide(self):
        """
        Checks if the coefficient alters the right hand side of the PDE.

        :rtype: ``bool``
        :return: True if the right hand side of the PDE is changed when the
                 coefficient is changed, ``None`` otherwise.
        """
        if self.altering==self.RIGHTHANDSIDE or self.altering==self.BOTH:
            return not None
        else:
            return None

    def estimateNumEquationsAndNumSolutions(self,domain,shape=()):
       """
       Tries to estimate the number of equations and number of solutions if
       the coefficient has the given shape.

       :param domain: domain on which the PDE uses the coefficient
       :type domain: `Domain`
       :param shape: suggested shape of the coefficient
       :type shape: ``tuple`` of ``int`` values
       :return: the number of equations and number of solutions of the PDE if
                the coefficient has given shape. If no appropriate numbers
                could be identified, ``None`` is returned
       :rtype: ``tuple`` of two ``int`` values or ``None``
       """
       dim=domain.getDim()
       if len(shape)>0:
           num=max(shape)+1
       else:
           num=1
       search=[]
       if self.definesNumEquation() and self.definesNumSolutions():
          for u in range(num):
             for e in range(num):
                search.append((e,u))
          search.sort(key=lambda x: -(x[0]+x[1]))
          for item in search:
             s=self.getShape(domain,item[0],item[1])
             if len(s)==0 and len(shape)==0:
                 return (1,1)
             else:
                 if s==shape: return item
       elif self.definesNumEquation():
          for e in range(num,0,-1):
             s=self.getShape(domain,e,0)
             if len(s)==0 and len(shape)==0:
                 return (1,None)
             else:
                 if s==shape: return (e,None)

       elif self.definesNumSolutions():
          for u in range(num,0,-1):
             s=self.getShape(domain,0,u)
             if len(s)==0 and len(shape)==0:
                 return (None,1)
             else:
                 if s==shape: return (None,u)
       return None

    def definesNumSolutions(self):
       """
       Checks if the coefficient allows to estimate the number of solution
       components.

       :return: True if the coefficient allows an estimate of the number of
                solution components, False otherwise
       :rtype: ``bool``
       """
       for i in self.pattern:
             if i==self.BY_SOLUTION: return True
       return False

    def definesNumEquation(self):
       """
       Checks if the coefficient allows to estimate the number of equations.

       :return: True if the coefficient allows an estimate of the number of
                equations, False otherwise
       :rtype: ``bool``
       """
       for i in self.pattern:
             if i==self.BY_EQUATION: return True
       return False

    def __CompTuple2(self,t1,t2):
      """
      Compares two tuples of possible number of equations and number of
      solutions.

      :param t1: the first tuple
      :param t2: the second tuple
      :return: 0, 1, or -1
      """

      dif=t1[0]+t1[1]-(t2[0]+t2[1])
      if dif<0: return 1
      elif dif>0: return -1
      else: return 0

    def getShape(self,domain,numEquations=1,numSolutions=1):
       """
       Builds the required shape of the coefficient.

       :param domain: domain on which the PDE uses the coefficient
       :type domain: `Domain`
       :param numEquations: number of equations of the PDE
       :type numEquations: ``int``
       :param numSolutions: number of components of the PDE solution
       :type numSolutions: ``int``
       :return: shape of the coefficient
       :rtype: ``tuple`` of ``int`` values
       """
       dim=domain.getDim()
       s=()
       for i in self.pattern:
             if i==self.BY_EQUATION:
                if numEquations>1: s=s+(numEquations,)
             elif i==self.BY_SOLUTION:
                if numSolutions>1: s=s+(numSolutions,)
             else:
                s=s+(dim,)
       return s

#====================================================================================================================

class LinearProblem(object):
   """
   This is the base class to define a general linear PDE-type problem for
   for an unknown function *u* on a given domain defined through a
   `Domain` object. The problem can be given as a single
   equation or as a system of equations.

   The class assumes that some sort of assembling process is required to form
   a problem of the form

   *L u=f*

   where *L* is an operator and *f* is the right hand side. This operator
   problem will be solved to get the unknown *u*.

   """
   def __init__(self,domain,numEquations=None,numSolutions=None,debug=False):
     """
     Initializes a linear problem.

     :param domain: domain of the PDE
     :type domain: `Domain`
     :param numEquations: number of equations. If ``None`` the number of
                          equations is extracted from the coefficients.
     :param numSolutions: number of solution components. If ``None`` the number
                          of solution components is extracted from the
                          coefficients.
     :param debug: if True debug information is printed

     """
     super(LinearProblem, self).__init__()

     self.__debug=debug
     self.__domain=domain
     self.__numEquations=numEquations
     self.__numSolutions=numSolutions
     self.__altered_coefficients=False
     self.__reduce_equation_order=False
     self.__reduce_solution_order=False
     self.__sym=False
     self.setSolverOptions()
     self.__is_RHS_valid=False
     self.__is_operator_valid=False
     self.__COEFFICIENTS={}
     self.__solution_rtol=1.e99
     self.__solution_atol=1.e99
     self.setSymmetryOff()
     # initialize things:
     self.resetAllCoefficients()
     self.initializeSystem()
   # ==========================================================================
   #    general stuff:
   # ==========================================================================
   def __str__(self):
     """
     Returns a string representation of the PDE.

     :return: a simple representation of the PDE
     :rtype: ``str``
     """
     return "<LinearProblem %d>"%id(self)
   # ==========================================================================
   #    debug :
   # ==========================================================================
   def setDebugOn(self):
     """
     Switches debug output on.
     """
     self.__debug=not None

   def setDebugOff(self):
     """
     Switches debug output off.
     """
     self.__debug=None

   def setDebug(self, flag):
     """
     Switches debug output on if ``flag`` is True otherwise it is switched off.

     :param flag: desired debug status
     :type flag: ``bool``
     """
     if flag:
         self.setDebugOn()
     else:
         self.setDebugOff()

   def trace(self,text):
     """
     Prints the text message if debug mode is switched on.

     :param text: message to be printed
     :type text: ``string``
     """
     if self.__debug: print(("%s: %s"%(str(self),text)))

   # ==========================================================================
   # some service functions:
   # ==========================================================================
   def introduceCoefficients(self,**coeff):
       """
       Introduces new coefficients into the problem.

       Use:

       p.introduceCoefficients(A=PDECoef(...), B=PDECoef(...))

       to introduce the coefficients *A* and *B*.
       """
       for name, type in list(coeff.items()):
           if not isinstance(type,PDECoef):
              raise ValueError("coefficient %s has no type."%name)
           self.__COEFFICIENTS[name]=type
           self.__COEFFICIENTS[name].resetValue()
           self.trace("coefficient %s has been introduced."%name)
   def resetRightHandSideCoefficients(self):
       """
       Resets all coefficients defining the right hand side
       """
       for name in self.__COEFFICIENTS:
         if self.__COEFFICIENTS[name].altering == PDECoef.RIGHTHANDSIDE  :
              self.__COEFFICIENTS[name].resetValue()
              self.trace("coefficient %s has been reset."%name)

   def getDomain(self):
     """
     Returns the domain of the PDE.

     :return: the domain of the PDE
     :rtype: `Domain`
     """
     return self.__domain
   def getDomainStatus(self):
     """
     Return the status indicator of the domain
     """
     return self.getDomain().getStatus()

   def getSystemStatus(self):
     """
     Return the domain status used to build the current system
     """
     return self.__system_status
   def setSystemStatus(self,status=None):
     """
     Sets the system status to ``status`` if ``status`` is not present the 
     current status of the domain is used.
     """
     if status == None:
         self.__system_status=self.getDomainStatus()
     else: 
         self.__system_status=status

   def getDim(self):
     """
     Returns the spatial dimension of the PDE.

     :return: the spatial dimension of the PDE domain
     :rtype: ``int``
     """
     return self.getDomain().getDim()

   def getNumEquations(self):
     """
     Returns the number of equations.

     :return: the number of equations
     :rtype: ``int``
     :raise UndefinedPDEError: if the number of equations is not specified yet
     """
     if self.__numEquations==None:
         if self.__numSolutions==None:
            raise UndefinedPDEError("Number of equations is undefined. Please specify argument numEquations.")
         else:
            self.__numEquations=self.__numSolutions
     return self.__numEquations

   def getNumSolutions(self):
     """
     Returns the number of unknowns.

     :return: the number of unknowns
     :rtype: ``int``
     :raise UndefinedPDEError: if the number of unknowns is not specified yet
     """
     if self.__numSolutions==None:
        if self.__numEquations==None:
            raise UndefinedPDEError("Number of solution is undefined. Please specify argument numSolutions.")
        else:
            self.__numSolutions=self.__numEquations
     return self.__numSolutions

   def reduceEquationOrder(self):
     """
     Returns the status of order reduction for the equation.

     :return: True if reduced interpolation order is used for the
              representation of the equation, False otherwise
     :rtype: `bool`
     """
     return self.__reduce_equation_order

   def reduceSolutionOrder(self):
     """
     Returns the status of order reduction for the solution.

     :return: True if reduced interpolation order is used for the
              representation of the solution, False otherwise
     :rtype: `bool`
     """
     return self.__reduce_solution_order

   def getFunctionSpaceForEquation(self):
     """
     Returns the `FunctionSpace` used to discretize
     the equation.

     :return: representation space of equation
     :rtype: `FunctionSpace`
     """
     if self.reduceEquationOrder():
         return escore.ReducedSolution(self.getDomain())
     else:
         return escore.Solution(self.getDomain())

   def getFunctionSpaceForSolution(self):
     """
     Returns the `FunctionSpace` used to represent
     the solution.

     :return: representation space of solution
     :rtype: `FunctionSpace`
     """
     if self.reduceSolutionOrder():
         return escore.ReducedSolution(self.getDomain())
     else:
         return escore.Solution(self.getDomain())

   # ==========================================================================
   #   solver settings:
   # ==========================================================================
   def setSolverOptions(self,options=None):
       """
       Sets the solver options.

       :param options: the new solver options. If equal ``None``, the solver options are set to the default.
       :type options: `SolverOptions` or ``None``
       :note: The symmetry flag of options is overwritten by the symmetry flag of the `LinearProblem`.
       """
       if options==None:
          self.__solver_options=SolverOptions()
       elif isinstance(options, SolverOptions):
          self.__solver_options=options
       else:
          raise ValueError("options must be a SolverOptions object.")
       self.__solver_options.setSymmetry(self.__sym)
     
   def getSolverOptions(self):
       """
       Returns the solver options
   
       :rtype: `SolverOptions`
       """
       self.__solver_options.setSymmetry(self.__sym)
       return self.__solver_options
       
   def isUsingLumping(self):
      """
      Checks if matrix lumping is the current solver method.

      :return: True if the current solver method is lumping
      :rtype: ``bool``
      """
      return self.getSolverOptions().getSolverMethod() in [ SolverOptions.ROWSUM_LUMPING, SolverOptions.HRZ_LUMPING ]
   # ==========================================================================
   #    symmetry  flag:
   # ==========================================================================
   def isSymmetric(self):
      """
      Checks if symmetry is indicated.

      :return: True if a symmetric PDE is indicated, False otherwise
      :rtype: ``bool``
      :note: the method is equivalent to use getSolverOptions().isSymmetric()
      """
      self.getSolverOptions().isSymmetric()

   def setSymmetryOn(self):
      """
      Sets the symmetry flag. 
      :note: The method overwrites the symmetry flag set by the solver options
      """
      self.__sym=True
      self.getSolverOptions().setSymmetryOn()

   def setSymmetryOff(self):
      """
      Clears the symmetry flag.
      :note: The method overwrites the symmetry flag set by the solver options
      """
      self.__sym=False
      self.getSolverOptions().setSymmetryOff()

   def setSymmetry(self,flag=False):
      """
      Sets the symmetry flag to ``flag``.

      :param flag: If True, the symmetry flag is set otherwise reset.
      :type flag: ``bool``
      :note: The method overwrites the symmetry flag set by the solver options
      """
      self.getSolverOptions().setSymmetry(flag)
   # ==========================================================================
   # function space handling for the equation as well as the solution
   # ==========================================================================
   def setReducedOrderOn(self):
     """
     Switches reduced order on for solution and equation representation.

     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     self.setReducedOrderForSolutionOn()
     self.setReducedOrderForEquationOn()

   def setReducedOrderOff(self):
     """
     Switches reduced order off for solution and equation representation

     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     self.setReducedOrderForSolutionOff()
     self.setReducedOrderForEquationOff()

   def setReducedOrderTo(self,flag=False):
     """
     Sets order reduction state for both solution and equation representation
     according to flag.

     :param flag: if True, the order reduction is switched on for both solution
                  and equation representation, otherwise or if flag is not
                  present order reduction is switched off
     :type flag: ``bool``
     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     self.setReducedOrderForSolutionTo(flag)
     self.setReducedOrderForEquationTo(flag)


   def setReducedOrderForSolutionOn(self):
     """
     Switches reduced order on for solution representation.

     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     if not self.__reduce_solution_order:
         if self.__altered_coefficients:
              raise RuntimeError("order cannot be altered after coefficients have been defined.")
         self.trace("Reduced order is used for solution representation.")
         self.__reduce_solution_order=True
         self.initializeSystem()

   def setReducedOrderForSolutionOff(self):
     """
     Switches reduced order off for solution representation

     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set.
     """
     if self.__reduce_solution_order:
         if self.__altered_coefficients:
              raise RuntimeError("order cannot be altered after coefficients have been defined.")
         self.trace("Full order is used to interpolate solution.")
         self.__reduce_solution_order=False
         self.initializeSystem()

   def setReducedOrderForSolutionTo(self,flag=False):
     """
     Sets order reduction state for solution representation according to flag.

     :param flag: if flag is True, the order reduction is switched on for
                  solution representation, otherwise or if flag is not present
                  order reduction is switched off
     :type flag: ``bool``
     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     if flag:
        self.setReducedOrderForSolutionOn()
     else:
        self.setReducedOrderForSolutionOff()

   def setReducedOrderForEquationOn(self):
     """
     Switches reduced order on for equation representation.

     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     if not self.__reduce_equation_order:
         if self.__altered_coefficients:
              raise RuntimeError("order cannot be altered after coefficients have been defined.")
         self.trace("Reduced order is used for test functions.")
         self.__reduce_equation_order=True
         self.initializeSystem()

   def setReducedOrderForEquationOff(self):
     """
     Switches reduced order off for equation representation.

     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     if self.__reduce_equation_order:
         if self.__altered_coefficients:
              raise RuntimeError("order cannot be altered after coefficients have been defined.")
         self.trace("Full order is used for test functions.")
         self.__reduce_equation_order=False
         self.initializeSystem()

   def setReducedOrderForEquationTo(self,flag=False):
     """
     Sets order reduction state for equation representation according to flag.

     :param flag: if flag is True, the order reduction is switched on for
                  equation representation, otherwise or if flag is not present
                  order reduction is switched off
     :type flag: ``bool``
     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     if flag:
        self.setReducedOrderForEquationOn()
     else:
        self.setReducedOrderForEquationOff()
   def getOperatorType(self):
      """
      Returns the current system type.
      """
      return self.__operator_type

   def checkSymmetricTensor(self,name,verbose=True):
      """
      Tests a coefficient for symmetry.

      :param name: name of the coefficient
      :type name: ``str``
      :param verbose: if set to True or not present a report on coefficients
                      which break the symmetry is printed.
      :type verbose: ``bool``
      :return: True if coefficient ``name`` is symmetric
      :rtype: ``bool``
      """
      SMALL_TOLERANCE=util.EPSILON*10.
      A=self.getCoefficient(name)
      verbose=verbose or self.__debug
      out=True
      if not A.isEmpty():
         tol=util.Lsup(A)*SMALL_TOLERANCE
         s=A.getShape()
         if A.getRank() == 4:
            if s[0]==s[2] and s[1] == s[3]:
               for i in range(s[0]):
                  for j in range(s[1]):
                     for k in range(s[2]):
                        for l in range(s[3]):
                            if util.Lsup(A[i,j,k,l]-A[k,l,i,j])>tol:
                               if verbose: print(("non-symmetric problem as %s[%d,%d,%d,%d]!=%s[%d,%d,%d,%d]"%(name,i,j,k,l,name,k,l,i,j)))
                               out=False
            else:
                 if verbose: print(("non-symmetric problem because of inappropriate shape %s of coefficient %s."%(s,name)))
                 out=False
         elif A.getRank() == 2:
            if s[0]==s[1]:
               for j in range(s[0]):
                  for l in range(s[1]):
                     if util.Lsup(A[j,l]-A[l,j])>tol:
                        if verbose: print(("non-symmetric problem because %s[%d,%d]!=%s[%d,%d]"%(name,j,l,name,l,j)))
                        out=False
            else:
                 if verbose: print(("non-symmetric problem because of inappropriate shape %s of coefficient %s."%(s,name)))
                 out=False
         elif A.getRank() == 0:
            pass
         else:
             raise ValueError("Cannot check rank %s of %s."%(A.getRank(),name))
      return out

   def checkReciprocalSymmetry(self,name0,name1,verbose=True):
      """
      Tests two coefficients for reciprocal symmetry.

      :param name0: name of the first coefficient
      :type name0: ``str``
      :param name1: name of the second coefficient
      :type name1: ``str``
      :param verbose: if set to True or not present a report on coefficients
                      which break the symmetry is printed
      :type verbose: ``bool``
      :return: True if coefficients ``name0`` and ``name1`` are reciprocally
               symmetric.
      :rtype: ``bool``
      """
      SMALL_TOLERANCE=util.EPSILON*10.
      B=self.getCoefficient(name0)
      C=self.getCoefficient(name1)
      verbose=verbose or self.__debug
      out=True
      if B.isEmpty() and not C.isEmpty():
         if verbose: print(("non-symmetric problem because %s is not present but %s is"%(name0,name1)))
         out=False
      elif not B.isEmpty() and C.isEmpty():
         if verbose: print(("non-symmetric problem because %s is not present but %s is"%(name0,name1)))
         out=False
      elif not B.isEmpty() and not C.isEmpty():
         sB=B.getShape()
         sC=C.getShape()
         tol=(util.Lsup(B)+util.Lsup(C))*SMALL_TOLERANCE/2.
         if len(sB) != len(sC):
             if verbose: print(("non-symmetric problem because ranks of %s (=%s) and %s (=%s) are different."%(name0,len(sB),name1,len(sC))))
             out=False
         else:
             if len(sB)==0:
               if util.Lsup(B-C)>tol:
                  if verbose: print(("non-symmetric problem because %s!=%s"%(name0,name1)))
                  out=False
             elif len(sB)==1:
               if sB[0]==sC[0]:
                  for j in range(sB[0]):
                     if util.Lsup(B[j]-C[j])>tol:
                        if verbose: print(("non-symmetric PDE because %s[%d]!=%s[%d]"%(name0,j,name1,j)))
                        out=False
               else:
                 if verbose: print(("non-symmetric problem because of inappropriate shapes %s and %s of coefficients %s and %s, respectively."%(sB,sC,name0,name1)))
             elif len(sB)==3:
               if sB[0]==sC[1] and sB[1]==sC[2] and sB[2]==sC[0]:
                   for i in range(sB[0]):
                      for j in range(sB[1]):
                         for k in range(sB[2]):
                            if util.Lsup(B[i,j,k]-C[k,i,j])>tol:
                                 if verbose: print(("non-symmetric problem because %s[%d,%d,%d]!=%s[%d,%d,%d]"%(name0,i,j,k,name1,k,i,j)))
                                 out=False
               else:
                 if verbose: print(("non-symmetric problem because of inappropriate shapes %s and %s of coefficients %s and %s, respectively."%(sB,sC,name0,name1)))
             else:
                 raise ValueError("Cannot check rank %s of %s and %s."%(len(sB),name0,name1))
      return out

   def getCoefficient(self,name):
     """
     Returns the value of the coefficient ``name``.

     :param name: name of the coefficient requested
     :type name: ``string``
     :return: the value of the coefficient
     :rtype: `Data`
     :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
     """
     if self.hasCoefficient(name):
         return self.__COEFFICIENTS[name].getValue()
     else:
        raise IllegalCoefficient("illegal coefficient %s requested for general PDE."%name)

   def hasCoefficient(self,name):
     """
     Returns True if ``name`` is the name of a coefficient.

     :param name: name of the coefficient enquired
     :type name: ``string``
     :return: True if ``name`` is the name of a coefficient of the general PDE,
              False otherwise
     :rtype: ``bool``
     """
     return name in self.__COEFFICIENTS

   def createCoefficient(self, name):
     """
     Creates a `Data` object corresponding to coefficient
     ``name``.

     :return: the coefficient ``name`` initialized to 0
     :rtype: `Data`
     :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
     """
     if self.hasCoefficient(name):
        return escore.Data(0.,self.getShapeOfCoefficient(name),self.getFunctionSpaceForCoefficient(name))
     else:
        raise IllegalCoefficient("illegal coefficient %s requested for general PDE."%name)

   def getFunctionSpaceForCoefficient(self,name):
     """
     Returns the `FunctionSpace` to be used for
     coefficient ``name``.

     :param name: name of the coefficient enquired
     :type name: ``string``
     :return: the function space to be used for coefficient ``name``
     :rtype: `FunctionSpace`
     :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
     """
     if self.hasCoefficient(name):
        return self.__COEFFICIENTS[name].getFunctionSpace(self.getDomain())
     else:
        raise ValueError("unknown coefficient %s requested"%name)

   def getShapeOfCoefficient(self,name):
     """
     Returns the shape of the coefficient ``name``.

     :param name: name of the coefficient enquired
     :type name: ``string``
     :return: the shape of the coefficient ``name``
     :rtype: ``tuple`` of ``int``
     :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
     """
     if self.hasCoefficient(name):
        return self.__COEFFICIENTS[name].getShape(self.getDomain(),self.getNumEquations(),self.getNumSolutions())
     else:
        raise IllegalCoefficient("illegal coefficient %s requested for general PDE."%name)

   def resetAllCoefficients(self):
     """
     Resets all coefficients to their default values.
     """
     for i in list(self.__COEFFICIENTS.keys()):
         self.__COEFFICIENTS[i].resetValue()

   def alteredCoefficient(self,name):
     """
     Announces that coefficient ``name`` has been changed.

     :param name: name of the coefficient affected
     :type name: ``string``
     :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
     :note: if ``name`` is q or r, the method will not trigger a rebuild of the
            system as constraints are applied to the solved system.
     """
     if self.hasCoefficient(name):
        self.trace("Coefficient %s has been altered."%name)
        if not ((name=="q" or name=="r") and self.isUsingLumping()):
           if self.__COEFFICIENTS[name].isAlteringOperator(): self.invalidateOperator()
           if self.__COEFFICIENTS[name].isAlteringRightHandSide(): self.invalidateRightHandSide()
     else:
        raise IllegalCoefficient("illegal coefficient %s requested for general PDE."%name)

   def validSolution(self):
       """
       Marks the solution as valid.
       """
       self.__is_solution_valid=True

   def invalidateSolution(self):
       """
       Indicates the PDE has to be resolved if the solution is requested.
       """
       self.trace("System will be resolved.")
       self.__is_solution_valid=False

   def isSolutionValid(self):
       """
       Returns True if the solution is still valid.
       """
       if not self.getDomainStatus()==self.getSystemStatus(): self.invalidateSolution()
       if self.__solution_rtol>self.getSolverOptions().getTolerance() or \
          self.__solution_atol>self.getSolverOptions().getAbsoluteTolerance():
            self.invalidateSolution()  
       return self.__is_solution_valid

   def validOperator(self):
       """
       Marks the operator as valid.
       """
       self.__is_operator_valid=True

   def invalidateOperator(self):
       """
       Indicates the operator has to be rebuilt next time it is used.
       """
       self.trace("Operator will be rebuilt.")
       self.invalidateSolution()
       self.__is_operator_valid=False

   def isOperatorValid(self):
       """
       Returns True if the operator is still valid.
       """
       if not self.getDomainStatus()==self.getSystemStatus(): self.invalidateOperator()
       if not self.getRequiredOperatorType()==self.getOperatorType(): self.invalidateOperator()
       return self.__is_operator_valid

   def validRightHandSide(self):
       """
       Marks the right hand side as valid.
       """
       self.__is_RHS_valid=True

   def invalidateRightHandSide(self):
       """
       Indicates the right hand side has to be rebuilt next time it is used.
       """
       self.trace("Right hand side has to be rebuilt.")
       self.invalidateSolution()
       self.__is_RHS_valid=False

   def isRightHandSideValid(self):
       """
       Returns True if the operator is still valid.
       """
       if not self.getDomainStatus()==self.getSystemStatus(): self.invalidateRightHandSide()
       return self.__is_RHS_valid

   def invalidateSystem(self):
       """
       Announces that everything has to be rebuilt.
       """
       self.invalidateSolution()
       self.invalidateOperator()
       self.invalidateRightHandSide()

   def isSystemValid(self):
       """
       Returns True if the system (including solution) is still vaild.
       """
       return self.isSolutionValid() and self.isOperatorValid() and self.isRightHandSideValid()

   def initializeSystem(self):
       """
       Resets the system clearing the operator, right hand side and solution.
       """
       self.trace("New System has been created.")
       self.__operator_type=None
       self.setSystemStatus()
       self.__operator=escore.Operator()
       self.__righthandside=escore.Data()
       self.__solution=escore.Data()
       self.invalidateSystem()

   def getOperator(self):
     """
     Returns the operator of the linear problem.

     :return: the operator of the problem
     """
     return self.getSystem()[0]

   def getRightHandSide(self):
     """
     Returns the right hand side of the linear problem.

     :return: the right hand side of the problem
     :rtype: `Data`
     """
     return self.getSystem()[1]

   def createRightHandSide(self):
       """
       Returns an instance of a new right hand side.
       """
       self.trace("New right hand side is allocated.")
       if self.getNumEquations()>1:
           return escore.Data(0.,(self.getNumEquations(),),self.getFunctionSpaceForEquation(),True)
       else:
           return escore.Data(0.,(),self.getFunctionSpaceForEquation(),True)

   def createSolution(self):
       """
       Returns an instance of a new solution.
       """
       self.trace("New solution is allocated.")
       if self.getNumSolutions()>1:
           return escore.Data(0.,(self.getNumSolutions(),),self.getFunctionSpaceForSolution(),True)
       else:
           return escore.Data(0.,(),self.getFunctionSpaceForSolution(),True)

   def resetSolution(self):
       """
       Sets the solution to zero.
       """
       if self.__solution.isEmpty():
           self.__solution=self.createSolution()
       else:
           self.__solution.setToZero()
           self.trace("Solution is reset to zero.")

   def setSolution(self,u, validate=True):
       """
       Sets the solution assuming that makes the system valid with the tolrance
       defined by the solver options
       """
       if validate:
          self.__solution_rtol=self.getSolverOptions().getTolerance()
          self.__solution_atol=self.getSolverOptions().getAbsoluteTolerance()
          self.validSolution()
       self.__solution=u
   def getCurrentSolution(self):
       """
       Returns the solution in its current state.
       """
       if self.__solution.isEmpty(): self.__solution=self.createSolution()
       return self.__solution

   def resetRightHandSide(self):
       """
       Sets the right hand side to zero.
       """
       if self.__righthandside.isEmpty():
           self.__righthandside=self.createRightHandSide()
       else:
           self.__righthandside.setToZero()
           self.trace("Right hand side is reset to zero.")

   def getCurrentRightHandSide(self):
       """
       Returns the right hand side in its current state.
       """
       if self.__righthandside.isEmpty(): self.__righthandside=self.createRightHandSide()
       return self.__righthandside

   def resetOperator(self):
       """
       Makes sure that the operator is instantiated and returns it initialized
       with zeros.
       """
       if self.getOperatorType() == None:
           if self.isUsingLumping():
               self.__operator=self.createSolution()
           else:
               self.__operator=self.createOperator()
           self.__operator_type=self.getRequiredOperatorType()
       else:
           if self.isUsingLumping():
               self.__operator.setToZero()
           else:
               if self.getOperatorType() == self.getRequiredOperatorType():
                   self.__operator.resetValues()
               else:
                   self.__operator=self.createOperator()
                   self.__operator_type=self.getRequiredOperatorType()
           self.trace("Operator reset to zero")

   def getCurrentOperator(self):
       """
       Returns the operator in its current state.
       """
       return self.__operator

   def setValue(self,**coefficients):
      """
      Sets new values to coefficients.

      :raise IllegalCoefficient: if an unknown coefficient keyword is used
      """
      # check if the coefficients are  legal:
      for i in list(coefficients.keys()):
         if not self.hasCoefficient(i):
            raise IllegalCoefficient("Attempt to set unknown coefficient %s"%i)
      # if the number of unknowns or equations is still unknown we try to estimate them:
      if self.__numEquations==None or self.__numSolutions==None:
         for i,d in list(coefficients.items()):
            if hasattr(d,"shape"):
                s=d.shape
            elif hasattr(d,"getShape"):
                s=d.getShape()
            else:
                s=numpy.array(d).shape
            if s!=None:
                # get number of equations and number of unknowns:
                res=self.__COEFFICIENTS[i].estimateNumEquationsAndNumSolutions(self.getDomain(),s)
                if res==None:
                    raise IllegalCoefficientValue("Illegal shape %s of coefficient %s"%(s,i))
                else:
                    if self.__numEquations==None: self.__numEquations=res[0]
                    if self.__numSolutions==None: self.__numSolutions=res[1]
      if self.__numEquations==None: raise UndefinedPDEError("unidentified number of equations")
      if self.__numSolutions==None: raise UndefinedPDEError("unidentified number of solutions")
      # now we check the shape of the coefficient if numEquations and numSolutions are set:
      for i,d in list(coefficients.items()):
        try:
           self.__COEFFICIENTS[i].setValue(self.getDomain(),
                     self.getNumEquations(),self.getNumSolutions(),
                     self.reduceEquationOrder(),self.reduceSolutionOrder(),d)
           self.alteredCoefficient(i)
        except IllegalCoefficientFunctionSpace as m:
            # if the function space is wrong then we try the reduced version:
            i_red=i+"_reduced"
            if (not i_red in list(coefficients.keys())) and i_red in list(self.__COEFFICIENTS.keys()):
                try:
                    self.__COEFFICIENTS[i_red].setValue(self.getDomain(),
                                                      self.getNumEquations(),self.getNumSolutions(),
                                                      self.reduceEquationOrder(),self.reduceSolutionOrder(),d)
                    self.alteredCoefficient(i_red)
                except IllegalCoefficientValue as m:
                    raise IllegalCoefficientValue("Coefficient %s:%s"%(i,m))
                except IllegalCoefficientFunctionSpace as m:
                    raise IllegalCoefficientFunctionSpace("Coefficient %s:%s"%(i,m))
            else:
                raise IllegalCoefficientFunctionSpace("Coefficient %s:%s"%(i,m))
        except IllegalCoefficientValue as m:
           raise IllegalCoefficientValue("Coefficient %s:%s"%(i,m))
      self.__altered_coefficients=True

   # ==========================================================================
   # methods that are typically overwritten when implementing a particular
   # linear problem
   # ==========================================================================
   def getRequiredOperatorType(self):
      """
      Returns the system type which needs to be used by the current set up.

      :note: Typically this method is overwritten when implementing a
             particular linear problem.
      """
      return None

   def createOperator(self):
       """
       Returns an instance of a new operator.

       :note: This method is overwritten when implementing a particular
              linear problem.
       """
       return escore.Operator()

   def checkSymmetry(self,verbose=True):
      """
      Tests the PDE for symmetry.

      :param verbose: if set to True or not present a report on coefficients
                      which break the symmetry is printed
      :type verbose: ``bool``
      :return: True if the problem is symmetric
      :rtype: ``bool``
      :note: Typically this method is overwritten when implementing a
             particular linear problem.
      """
      out=True
      return out

   def getSolution(self,**options):
       """
       Returns the solution of the problem.

       :return: the solution
       :rtype: `Data`

       :note: This method is overwritten when implementing a particular
              linear problem.
       """
       return self.getCurrentSolution()

   def getSystem(self):
       """
       Returns the operator and right hand side of the PDE.

       :return: the discrete version of the PDE
       :rtype: ``tuple`` of `Operator` and `Data`.

       :note: This method is overwritten when implementing a particular
              linear problem.
       """
       return (self.getCurrentOperator(), self.getCurrentRightHandSide())

class LinearPDE(LinearProblem):
   """
   This class is used to define a general linear, steady, second order PDE
   for an unknown function *u* on a given domain defined through a
   `Domain` object.

   For a single PDE having a solution with a single component the linear PDE
   is defined in the following form:

   *-(grad(A[j,l]+A_reduced[j,l])*grad(u)[l]+(B[j]+B_reduced[j])u)[j]+(C[l]+C_reduced[l])*grad(u)[l]+(D+D_reduced)=-grad(X+X_reduced)[j,j]+(Y+Y_reduced)*

   where *grad(F)* denotes the spatial derivative of *F*. Einstein's
   summation convention, ie. summation over indexes appearing twice in a term
   of a sum performed, is used.
   The coefficients *A*, *B*, *C*, *D*, *X* and *Y* have to be specified
   through `Data` objects in `Function` and
   the coefficients *A_reduced*, *B_reduced*, *C_reduced*, *D_reduced*,
   *X_reduced* and *Y_reduced* have to be specified through
   `Data` objects in `ReducedFunction`.
   It is also allowed to use objects that can be converted into such
   `Data` objects. *A* and *A_reduced* are rank two, *B*,
   *C*, *X*, *B_reduced*, *C_reduced* and *X_reduced* are rank one and
   *D*, *D_reduced*, *Y* and *Y_reduced* are scalar.

   The following natural boundary conditions are considered:

   *n[j]*((A[i,j]+A_reduced[i,j])*grad(u)[l]+(B+B_reduced)[j]*u)+(d+d_reduced)*u=n[j]*(X[j]+X_reduced[j])+y*

   where *n* is the outer normal field. Notice that the coefficients *A*,
   *A_reduced*, *B*, *B_reduced*, *X* and *X_reduced* are defined in the
   PDE. The coefficients *d* and *y* are each a scalar in
   `FunctionOnBoundary` and the coefficients
   *d_reduced* and *y_reduced* are each a scalar in
   `ReducedFunctionOnBoundary`.

   Constraints for the solution prescribe the value of the solution at certain
   locations in the domain. They have the form

   *u=r* where *q>0*

   *r* and *q* are each scalar where *q* is the characteristic function
   defining where the constraint is applied. The constraints override any
   other condition set by the PDE or the boundary condition.

   The PDE is symmetrical if

   *A[i,j]=A[j,i]*  and *B[j]=C[j]* and *A_reduced[i,j]=A_reduced[j,i]*
   and *B_reduced[j]=C_reduced[j]*

   For a system of PDEs and a solution with several components the PDE has the
   form

   *-grad((A[i,j,k,l]+A_reduced[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k])[j]+(C[i,k,l]+C_reduced[i,k,l])*grad(u[k])[l]+(D[i,k]+D_reduced[i,k]*u[k] =-grad(X[i,j]+X_reduced[i,j])[j]+Y[i]+Y_reduced[i]*

   *A* and *A_reduced* are of rank four, *B*, *B_reduced*, *C* and
   *C_reduced* are each of rank three, *D*, *D_reduced*, *X_reduced* and
   *X* are each of rank two and *Y* and *Y_reduced* are of rank one.
   The natural boundary conditions take the form:

   *n[j]*((A[i,j,k,l]+A_reduced[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k])+(d[i,k]+d_reduced[i,k])*u[k]=n[j]*(X[i,j]+X_reduced[i,j])+y[i]+y_reduced[i]*

   The coefficient *d* is of rank two and *y* is of rank one both in
   `FunctionOnBoundary`. The coefficients
   *d_reduced* is of rank two and *y_reduced* is of rank one both in
   `ReducedFunctionOnBoundary`.

   Constraints take the form

   *u[i]=r[i]*  where  *q[i]>0*

   *r* and *q* are each rank one. Notice that at some locations not
   necessarily all components must have a constraint.

   The system of PDEs is symmetrical if

      - *A[i,j,k,l]=A[k,l,i,j]*
      - *A_reduced[i,j,k,l]=A_reduced[k,l,i,j]*
      - *B[i,j,k]=C[k,i,j]*
      - *B_reduced[i,j,k]=C_reduced[k,i,j]*
      - *D[i,k]=D[i,k]*
      - *D_reduced[i,k]=D_reduced[i,k]*
      - *d[i,k]=d[k,i]*
      - *d_reduced[i,k]=d_reduced[k,i]*

   `LinearPDE` also supports solution discontinuities over a contact region
   in the domain. To specify the conditions across the discontinuity we are
   using the generalised flux *J* which, in the case of a system of PDEs
   and several components of the solution, is defined as

   *J[i,j]=(A[i,j,k,l]+A_reduced[[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k]-X[i,j]-X_reduced[i,j]*

   For the case of single solution component and single PDE *J* is defined as

   *J[j]=(A[i,j]+A_reduced[i,j])*grad(u)[j]+(B[i]+B_reduced[i])*u-X[i]-X_reduced[i]*

   In the context of discontinuities *n* denotes the normal on the
   discontinuity pointing from side 0 towards side 1 calculated from
   `FunctionSpace.getNormal` of `FunctionOnContactZero`.
   For a system of PDEs the contact condition takes the form

   *n[j]*J0[i,j]=n[j]*J1[i,j]=(y_contact[i]+y_contact_reduced[i])- (d_contact[i,k]+d_contact_reduced[i,k])*jump(u)[k]*

   where *J0* and *J1* are the fluxes on side 0 and side 1 of the
   discontinuity, respectively. *jump(u)*, which is the difference of the
   solution at side 1 and at side 0, denotes the jump of *u* across
   discontinuity along the normal calculated by `jump`.
   The coefficient *d_contact* is of rank two and *y_contact* is of rank one
   both in `FunctionOnContactZero` or
   `FunctionOnContactOne`.
   The coefficient *d_contact_reduced* is of rank two and *y_contact_reduced*
   is of rank one both in `ReducedFunctionOnContactZero`
   or `ReducedFunctionOnContactOne`.
   In case of a single PDE and a single component solution the contact
   condition takes the form

   *n[j]*J0_{j}=n[j]*J1_{j}=(y_contact+y_contact_reduced)-(d_contact+y_contact_reduced)*jump(u)*

   In this case the coefficient *d_contact* and *y_contact* are each scalar
   both in `FunctionOnContactZero` or
   `FunctionOnContactOne` and the coefficient
   *d_contact_reduced* and *y_contact_reduced* are each scalar both in
   `ReducedFunctionOnContactZero` or
   `ReducedFunctionOnContactOne`.

   Typical usage::

       p = LinearPDE(dom)
       p.setValue(A=kronecker(dom), D=1, Y=0.5)
       u = p.getSolution()

   """

   def __init__(self,domain,numEquations=None,numSolutions=None,debug=False):
     """
     Initializes a new linear PDE.

     :param domain: domain of the PDE
     :type domain: `Domain`
     :param numEquations: number of equations. If ``None`` the number of
                          equations is extracted from the PDE coefficients.
     :param numSolutions: number of solution components. If ``None`` the number
                          of solution components is extracted from the PDE
                          coefficients.
     :param debug: if True debug information is printed

     """
     super(LinearPDE, self).__init__(domain,numEquations,numSolutions,debug)
     #
     #   the coefficients of the PDE:
     #
     self.introduceCoefficients(
       A=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       B=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       C=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       D=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       X=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM),PDECoef.RIGHTHANDSIDE),
       Y=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       d=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       d_contact=PDECoef(PDECoef.CONTACT,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y_contact=PDECoef(PDECoef.CONTACT,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       A_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       B_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       C_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       D_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       X_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_DIM),PDECoef.RIGHTHANDSIDE),
       Y_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       d_reduced=PDECoef(PDECoef.BOUNDARY_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y_reduced=PDECoef(PDECoef.BOUNDARY_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       d_contact_reduced=PDECoef(PDECoef.CONTACT_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y_contact_reduced=PDECoef(PDECoef.CONTACT_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       d_dirac=PDECoef(PDECoef.DIRACDELTA,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y_dirac=PDECoef(PDECoef.DIRACDELTA,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       r=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.RIGHTHANDSIDE),
       q=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.BOTH) )

   def __str__(self):
     """
     Returns the string representation of the PDE.

     :return: a simple representation of the PDE
     :rtype: ``str``
     """
     return "<LinearPDE %d>"%id(self)

   def getRequiredOperatorType(self):
      """
      Returns the system type which needs to be used by the current set up.
      """
      if self.isUsingLumping():
         return "__ESCRIPT_DATA"
      else:
         solver_options=self.getSolverOptions()
         return self.getDomain().getSystemMatrixTypeId(solver_options.getSolverMethod(), solver_options.getPreconditioner(),solver_options.getPackage(), solver_options.isSymmetric())

   def checkSymmetry(self,verbose=True):
      """
      Tests the PDE for symmetry.

      :param verbose: if set to True or not present a report on coefficients
                      which break the symmetry is printed.
      :type verbose: ``bool``
      :return: True if the PDE is symmetric
      :rtype: `bool`
      :note: This is a very expensive operation. It should be used for
             degugging only! The symmetry flag is not altered.
      """
      out=True
      out=out and self.checkSymmetricTensor("A", verbose)
      out=out and self.checkSymmetricTensor("A_reduced", verbose)
      out=out and self.checkReciprocalSymmetry("B","C", verbose)
      out=out and self.checkReciprocalSymmetry("B_reduced","C_reduced", verbose)
      out=out and self.checkSymmetricTensor("D", verbose)
      out=out and self.checkSymmetricTensor("D_reduced", verbose)
      out=out and self.checkSymmetricTensor("d", verbose)
      out=out and self.checkSymmetricTensor("d_reduced", verbose)
      out=out and self.checkSymmetricTensor("d_contact", verbose)
      out=out and self.checkSymmetricTensor("d_contact_reduced", verbose)
      out=out and self.checkSymmetricTensor("d_dirac", verbose)
      return out

   def createOperator(self):
       """
       Returns an instance of a new operator.
       """
       optype=self.getRequiredOperatorType()
       self.trace("New operator of type %s is allocated."%optype)
       return self.getDomain().newOperator( \
                           self.getNumEquations(), \
                           self.getFunctionSpaceForEquation(), \
                           self.getNumSolutions(), \
                           self.getFunctionSpaceForSolution(), \
                           optype)

   def getSolution(self):
       """
       Returns the solution of the PDE.

       :return: the solution
       :rtype: `Data`
       """
       option_class=self.getSolverOptions()
       if not self.isSolutionValid():
          mat,f=self.getSystem()
          if self.isUsingLumping():
             if not util.inf(abs(mat)) > 0.:
                 raise ZeroDivisionError("Lumped mass matrix has zero entry (try order 1 elements or HRZ lumping).")
             self.setSolution(f*1/mat)
          else:
             self.trace("PDE is resolved.")
             self.trace("solver options: %s"%str(option_class))
             self.setSolution(mat.solve(f,option_class))
       return self.getCurrentSolution()

   def getSystem(self):
       """
       Returns the operator and right hand side of the PDE.

       :return: the discrete version of the PDE
       :rtype: ``tuple`` of `Operator` and
               `Data`
       """
       if not self.isOperatorValid() or not self.isRightHandSideValid():
          if self.isUsingLumping():
              if not self.isOperatorValid():
                 if not self.getFunctionSpaceForEquation()==self.getFunctionSpaceForSolution():
                      raise TypeError("Lumped matrix requires same order for equations and unknowns")
                 if not self.getCoefficient("A").isEmpty():
                      raise ValueError("coefficient A in lumped matrix may not be present.")
                 if not self.getCoefficient("B").isEmpty():
                      raise ValueError("coefficient B in lumped matrix may not be present.")
                 if not self.getCoefficient("C").isEmpty():
                      raise ValueError("coefficient C in lumped matrix may not be present.")
                 if not self.getCoefficient("d_contact").isEmpty():
                      raise ValueError("coefficient d_contact in lumped matrix may not be present.")
                 if not self.getCoefficient("A_reduced").isEmpty():
                      raise ValueError("coefficient A_reduced in lumped matrix may not be present.")
                 if not self.getCoefficient("B_reduced").isEmpty():
                      raise ValueError("coefficient B_reduced in lumped matrix may not be present.")
                 if not self.getCoefficient("C_reduced").isEmpty():
                      raise ValueError("coefficient C_reduced in lumped matrix may not be present.")
                 if not self.getCoefficient("d_contact_reduced").isEmpty():
                      raise ValueError("coefficient d_contact_reduced in lumped matrix may not be present.")
                 D=self.getCoefficient("D")
                 d=self.getCoefficient("d")
                 D_reduced=self.getCoefficient("D_reduced")
                 d_reduced=self.getCoefficient("d_reduced")
                 d_dirac=self.getCoefficient("d_dirac")

                 if not D.isEmpty():
                     if self.getNumSolutions()>1:
                        D_times_e=util.matrix_mult(D,numpy.ones((self.getNumSolutions(),)))
                     else:
                        D_times_e=D
                 else:
                    D_times_e=escore.Data()
                 if not d.isEmpty():
                     if self.getNumSolutions()>1:
                        d_times_e=util.matrix_mult(d,numpy.ones((self.getNumSolutions(),)))
                     else:
                        d_times_e=d
                 else:
                    d_times_e=escore.Data()

                 if not D_reduced.isEmpty():
                     if self.getNumSolutions()>1:
                        D_reduced_times_e=util.matrix_mult(D_reduced,numpy.ones((self.getNumSolutions(),)))
                     else:
                        D_reduced_times_e=D_reduced
                 else:
                    D_reduced_times_e=escore.Data()
                    
                 if not d_reduced.isEmpty():
                     if self.getNumSolutions()>1:
                        d_reduced_times_e=util.matrix_mult(d_reduced,numpy.ones((self.getNumSolutions(),)))
                     else:
                        d_reduced_times_e=d_reduced
                 else:
                    d_reduced_times_e=escore.Data()
                    
                 if not d_dirac.isEmpty():
                     if self.getNumSolutions()>1:
                        d_dirac_times_e=util.matrix_mult(d_dirac,numpy.ones((self.getNumSolutions(),)))
                     else:
                        d_reduced_dirac_e=d_dirac
                 else:
                    d_dirac_times_e=escore.Data()

                 self.resetOperator()
                 operator=self.getCurrentOperator()
                 if hasattr(self.getDomain(), "addPDEToLumpedSystem") :
                    hrz_lumping=( self.getSolverOptions().getSolverMethod() ==  SolverOptions.HRZ_LUMPING )
                    self.getDomain().addPDEToLumpedSystem(operator, D_times_e, d_times_e, d_dirac_times_e,  hrz_lumping )
                    self.getDomain().addPDEToLumpedSystem(operator, D_reduced_times_e, d_reduced_times_e, escore.Data(), hrz_lumping)
                 else:
                    self.getDomain().addPDEToRHS(operator, \
                                                 escore.Data(), \
                                                 D_times_e, \
                                                 d_times_e,\
                                                 escore.Data(),\
                                                 d_dirac_times_e)
                    self.getDomain().addPDEToRHS(operator, \
                                                 escore.Data(), \
                                                 D_reduced_times_e, \
                                                 d_reduced_times_e,\
                                                 escore.Data(), \
                                                 escore.Data())
                 self.trace("New lumped operator has been built.")
              if not self.isRightHandSideValid():
                 self.resetRightHandSide()
                 righthandside=self.getCurrentRightHandSide()
                 self.getDomain().addPDEToRHS(righthandside, \
                               self.getCoefficient("X"), \
                               self.getCoefficient("Y"),\
                               self.getCoefficient("y"),\
                               self.getCoefficient("y_contact"), \
                               self.getCoefficient("y_dirac"))
                 self.getDomain().addPDEToRHS(righthandside, \
                               self.getCoefficient("X_reduced"), \
                               self.getCoefficient("Y_reduced"),\
                               self.getCoefficient("y_reduced"),\
                               self.getCoefficient("y_contact_reduced"), \
                               escore.Data())
                 self.trace("New right hand side has been built.")
                 self.validRightHandSide()
              self.insertConstraint(rhs_only=False)
              self.validOperator()
          else:
             if not self.isOperatorValid() and not self.isRightHandSideValid():
                 self.resetRightHandSide()
                 righthandside=self.getCurrentRightHandSide()
                 self.resetOperator()
                 operator=self.getCurrentOperator()
                 self.getDomain().addPDEToSystem(operator,righthandside, \
                               self.getCoefficient("A"), \
                               self.getCoefficient("B"), \
                               self.getCoefficient("C"), \
                               self.getCoefficient("D"), \
                               self.getCoefficient("X"), \
                               self.getCoefficient("Y"), \
                               self.getCoefficient("d"), \
                               self.getCoefficient("y"), \
                               self.getCoefficient("d_contact"), \
                               self.getCoefficient("y_contact"), \
                               self.getCoefficient("d_dirac"), \
                               self.getCoefficient("y_dirac"))
                 self.getDomain().addPDEToSystem(operator,righthandside, \
                               self.getCoefficient("A_reduced"), \
                               self.getCoefficient("B_reduced"), \
                               self.getCoefficient("C_reduced"), \
                               self.getCoefficient("D_reduced"), \
                               self.getCoefficient("X_reduced"), \
                               self.getCoefficient("Y_reduced"), \
                               self.getCoefficient("d_reduced"), \
                               self.getCoefficient("y_reduced"), \
                               self.getCoefficient("d_contact_reduced"), \
                               self.getCoefficient("y_contact_reduced"), \
                               escore.Data(), \
                               escore.Data())
                 self.insertConstraint(rhs_only=False)
                 self.trace("New system has been built.")
                 self.validOperator()
                 self.validRightHandSide()
             elif not self.isRightHandSideValid():
                 self.resetRightHandSide()
                 righthandside=self.getCurrentRightHandSide()
                 self.getDomain().addPDEToRHS(righthandside,
                               self.getCoefficient("X"), \
                               self.getCoefficient("Y"),\
                               self.getCoefficient("y"),\
                               self.getCoefficient("y_contact"), \
                               self.getCoefficient("y_dirac") )
                 self.getDomain().addPDEToRHS(righthandside,
                               self.getCoefficient("X_reduced"), \
                               self.getCoefficient("Y_reduced"),\
                               self.getCoefficient("y_reduced"),\
                               self.getCoefficient("y_contact_reduced"), \
                               escore.Data())
                 self.insertConstraint(rhs_only=True)
                 self.trace("New right hand side has been built.")
                 self.validRightHandSide()
             elif not self.isOperatorValid():
                 self.resetOperator()
                 operator=self.getCurrentOperator()
                 self.getDomain().addPDEToSystem(operator,escore.Data(), \
                            self.getCoefficient("A"), \
                            self.getCoefficient("B"), \
                            self.getCoefficient("C"), \
                            self.getCoefficient("D"), \
                            escore.Data(), \
                            escore.Data(), \
                            self.getCoefficient("d"), \
                            escore.Data(),\
                            self.getCoefficient("d_contact"), \
                            escore.Data(),                   \
                            self.getCoefficient("d_dirac"),   \
                            escore.Data())
                 self.getDomain().addPDEToSystem(operator,escore.Data(), \
                            self.getCoefficient("A_reduced"), \
                            self.getCoefficient("B_reduced"), \
                            self.getCoefficient("C_reduced"), \
                            self.getCoefficient("D_reduced"), \
                            escore.Data(), \
                            escore.Data(), \
                            self.getCoefficient("d_reduced"), \
                            escore.Data(),\
                            self.getCoefficient("d_contact_reduced"), \
                            escore.Data(),  \
                            escore.Data(),  \
                            escore.Data())
                 self.insertConstraint(rhs_only=False)
                 self.trace("New operator has been built.")
                 self.validOperator()
       self.setSystemStatus()
       self.trace("System status is %s."%self.getSystemStatus())
       return (self.getCurrentOperator(), self.getCurrentRightHandSide())

   def insertConstraint(self, rhs_only=False):
      """
      Applies the constraints defined by q and r to the PDE.

      :param rhs_only: if True only the right hand side is altered by the
                       constraint
      :type rhs_only: ``bool``
      """
      q=self.getCoefficient("q")
      r=self.getCoefficient("r")
      righthandside=self.getCurrentRightHandSide()
      operator=self.getCurrentOperator()

      if not q.isEmpty():
         if r.isEmpty():
            r_s=self.createSolution()
         else:
            r_s=r
         if not rhs_only and not operator.isEmpty():
             if self.isUsingLumping():
                 operator.copyWithMask(escore.Data(1.,q.getShape(),q.getFunctionSpace()),q)
             else:
                 row_q=escore.Data(q,self.getFunctionSpaceForEquation())
                 col_q=escore.Data(q,self.getFunctionSpaceForSolution())
                 u=self.createSolution()
                 u.copyWithMask(r_s,col_q)
                 righthandside-=operator*u
                 operator.nullifyRowsAndCols(row_q,col_q,1.)
         righthandside.copyWithMask(r_s,q)

   def setValue(self,**coefficients):
      """
      Sets new values to coefficients.

      :param coefficients: new values assigned to coefficients
      :keyword A: value for coefficient ``A``
      :type A: any type that can be cast to a `Data` object on
               `Function`
      :keyword A_reduced: value for coefficient ``A_reduced``
      :type A_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword B: value for coefficient ``B``
      :type B: any type that can be cast to a `Data` object on
               `Function`
      :keyword B_reduced: value for coefficient ``B_reduced``
      :type B_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword C: value for coefficient ``C``
      :type C: any type that can be cast to a `Data` object on
               `Function`
      :keyword C_reduced: value for coefficient ``C_reduced``
      :type C_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword D: value for coefficient ``D``
      :type D: any type that can be cast to a `Data` object on
               `Function`
      :keyword D_reduced: value for coefficient ``D_reduced``
      :type D_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword X: value for coefficient ``X``
      :type X: any type that can be cast to a `Data` object on
               `Function`
      :keyword X_reduced: value for coefficient ``X_reduced``
      :type X_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword Y: value for coefficient ``Y``
      :type Y: any type that can be cast to a `Data` object on
               `Function`
      :keyword Y_reduced: value for coefficient ``Y_reduced``
      :type Y_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword d: value for coefficient ``d``
      :type d: any type that can be cast to a `Data` object on
               `FunctionOnBoundary`
      :keyword d_reduced: value for coefficient ``d_reduced``
      :type d_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunctionOnBoundary`
      :keyword y: value for coefficient ``y``
      :type y: any type that can be cast to a `Data` object on
               `FunctionOnBoundary`
      :keyword d_contact: value for coefficient ``d_contact``
      :type d_contact: any type that can be cast to a `Data`
                       object on `FunctionOnContactOne`
                       or `FunctionOnContactZero`
      :keyword d_contact_reduced: value for coefficient ``d_contact_reduced``
      :type d_contact_reduced: any type that can be cast to a `Data`
                               object on `ReducedFunctionOnContactOne`
                               or `ReducedFunctionOnContactZero`
      :keyword y_contact: value for coefficient ``y_contact``
      :type y_contact: any type that can be cast to a `Data`
                       object on `FunctionOnContactOne`
                       or `FunctionOnContactZero`
      :keyword y_contact_reduced: value for coefficient ``y_contact_reduced``
      :type y_contact_reduced: any type that can be cast to a `Data`
                               object on `ReducedFunctionOnContactOne`
                               or `ReducedFunctionOnContactZero`
      :keyword d_dirac: value for coefficient ``d_dirac``
      :type d_dirac: any type that can be cast to a `Data` object on `DiracDeltaFunctions`
      :keyword y_dirac: value for coefficient ``y_dirac``
      :type y_dirac: any type that can be cast to a `Data` object on `DiracDeltaFunctions`
      :keyword r: values prescribed to the solution at the locations of
                  constraints
      :type r: any type that can be cast to a `Data` object on
               `Solution` or `ReducedSolution`
               depending on whether reduced order is used for the solution
      :keyword q: mask for location of constraints
      :type q: any type that can be cast to a `Data` object on
               `Solution` or `ReducedSolution`
               depending on whether reduced order is used for the
               representation of the equation
      :raise IllegalCoefficient: if an unknown coefficient keyword is used
      """
      super(LinearPDE,self).setValue(**coefficients)
      # check if the systrem is inhomogeneous:
      if len(coefficients)>0 and not self.isUsingLumping():
         q=self.getCoefficient("q")
         r=self.getCoefficient("r")
         if not q.isEmpty() and not r.isEmpty():
             if util.Lsup(q*r)>0.:
               self.trace("Inhomogeneous constraint detected.")
               self.invalidateSystem()


   def getResidual(self,u=None):
     """
     Returns the residual of u or the current solution if u is not present.

     :param u: argument in the residual calculation. It must be representable
               in `self.getFunctionSpaceForSolution()`. If u is not present
               or equals ``None`` the current solution is used.
     :type u: `Data` or None
     :return: residual of u
     :rtype: `Data`
     """
     if u==None:
        return self.getOperator()*self.getSolution()-self.getRightHandSide()
     else:
        return self.getOperator()*escore.Data(u,self.getFunctionSpaceForSolution())-self.getRightHandSide()

   def getFlux(self,u=None):
     """
     Returns the flux *J* for a given *u*.

     *J[i,j]=(A[i,j,k,l]+A_reduced[A[i,j,k,l]]*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])u[k]-X[i,j]-X_reduced[i,j]*

     or

     *J[j]=(A[i,j]+A_reduced[i,j])*grad(u)[l]+(B[j]+B_reduced[j])u-X[j]-X_reduced[j]*

     :param u: argument in the flux. If u is not present or equals `None` the
               current solution is used.
     :type u: `Data` or None
     :return: flux
     :rtype: `Data`
     """
     if u==None: u=self.getSolution()
     if self.getNumEquations()>1:
       out = escore.Data(0.,(self.getNumEquations(),self.getDim()),self.getFunctionSpaceForCoefficient("X"))
     else:
       out = escore.Data(0.,(self.getDim(), ),self.getFunctionSpaceForCoefficient("X"))

     A=self.getCoefficient("A")
     if not A.isEmpty():
           out+=util.tensormult(A,util.grad(u,self.getFunctionSpaceForCoefficient("A")))
      
     B=self.getCoefficient("B")
     if not B.isEmpty():
           if B.getRank() == 1:
               out+=B * u
           else:
               out+=util.generalTensorProduct(B,u,axis_offset=1)

     X=self.getCoefficient("X") 
     if not X.isEmpty():
           out-=X

     A_reduced=self.getCoefficient("A_reduced")
     if not A_reduced.isEmpty():
           out+=util.tensormult(A_reduced, util.grad(u,self.getFunctionSpaceForCoefficient("A_reduced"))) \

     B_reduced=self.getCoefficient("B_reduced")
     if not B_reduced.isEmpty():
           if B_reduced.getRank() == 1:
                out+=B_reduced*u
           else:
                out+=util.generalTensorProduct(B_reduced,u,axis_offset=1)

     X_reduced=self.getCoefficient("X_reduced")
     if not X_reduced.isEmpty():
           out-=X_reduced
     return out

class Poisson(LinearPDE):
   """
   Class to define a Poisson equation problem. This is generally a
   `LinearPDE` of the form

   *-grad(grad(u)[j])[j] = f*

   with natural boundary conditions

   *n[j]*grad(u)[j] = 0*

   and constraints:

   *u=0* where *q>0*

   """

   def __init__(self,domain,debug=False):
     """
     Initializes a new Poisson equation.

     :param domain: domain of the PDE
     :type domain: `Domain`
     :param debug: if True debug information is printed

     """
     super(Poisson, self).__init__(domain,1,1,debug)
     self.introduceCoefficients(
                        f=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
                        f_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE))
     self.setSymmetryOn()

   def setValue(self,**coefficients):
     """
     Sets new values to coefficients.

     :param coefficients: new values assigned to coefficients
     :keyword f: value for right hand side *f*
     :type f: any type that can be cast to a `Scalar` object
              on `Function`
     :keyword q: mask for location of constraints
     :type q: any type that can be cast to a rank zero `Data`
              object on `Solution` or
              `ReducedSolution` depending on whether
              reduced order is used for the representation of the equation
     :raise IllegalCoefficient: if an unknown coefficient keyword is used
     """
     super(Poisson, self).setValue(**coefficients)


   def getCoefficient(self,name):
     """
     Returns the value of the coefficient ``name`` of the general PDE.

     :param name: name of the coefficient requested
     :type name: ``string``
     :return: the value of the coefficient ``name``
     :rtype: `Data`
     :raise IllegalCoefficient: invalid coefficient name
     :note: This method is called by the assembling routine to map the Poisson
            equation onto the general PDE.
     """
     if name == "A" :
         return escore.Data(util.kronecker(self.getDim()),escore.Function(self.getDomain()))
     elif name == "Y" :
         return self.getCoefficient("f")
     elif name == "Y_reduced" :
         return self.getCoefficient("f_reduced")
     else:
         return super(Poisson, self).getCoefficient(name)

class Helmholtz(LinearPDE):
   """
   Class to define a Helmholtz equation problem. This is generally a
   `LinearPDE` of the form

   *omega*u - grad(k*grad(u)[j])[j] = f*

   with natural boundary conditions

   *k*n[j]*grad(u)[j] = g- alphau*

   and constraints:

   *u=r* where *q>0*

   """

   def __init__(self,domain,debug=False):
     """
     Initializes a new Helmholtz equation.

     :param domain: domain of the PDE
     :type domain: `Domain`
     :param debug: if True debug information is printed

     """
     super(Helmholtz, self).__init__(domain,1,1,debug)
     self.introduceCoefficients(
                        omega=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.OPERATOR),
                        k=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.OPERATOR),
                        f=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
                        f_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
                        alpha=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,),PDECoef.OPERATOR),
                        g=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
                        g_reduced=PDECoef(PDECoef.BOUNDARY_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE))
     self.setSymmetryOn()

   def setValue(self,**coefficients):
     """
     Sets new values to coefficients.

     :param coefficients: new values assigned to coefficients
     :keyword omega: value for coefficient *omega*
     :type omega: any type that can be cast to a `Scalar`
                  object on `Function`
     :keyword k: value for coefficient *k*
     :type k: any type that can be cast to a `Scalar` object
              on `Function`
     :keyword f: value for right hand side *f*
     :type f: any type that can be cast to a `Scalar` object
              on `Function`
     :keyword alpha: value for right hand side *alpha*
     :type alpha: any type that can be cast to a `Scalar`
                  object on `FunctionOnBoundary`
     :keyword g: value for right hand side *g*
     :type g: any type that can be cast to a `Scalar` object
              on `FunctionOnBoundary`
     :keyword r: prescribed values *r* for the solution in constraints
     :type r: any type that can be cast to a `Scalar` object
              on `Solution` or
              `ReducedSolution` depending on whether
              reduced order is used for the representation of the equation
     :keyword q: mask for the location of constraints
     :type q: any type that can be cast to a `Scalar` object
              on `Solution` or
              `ReducedSolution` depending on whether
              reduced order is used for the representation of the equation
     :raise IllegalCoefficient: if an unknown coefficient keyword is used
     """
     super(Helmholtz, self).setValue(**coefficients)

   def getCoefficient(self,name):
     """
     Returns the value of the coefficient ``name`` of the general PDE.

     :param name: name of the coefficient requested
     :type name: ``string``
     :return: the value of the coefficient ``name``
     :rtype: `Data`
     :raise IllegalCoefficient: invalid name
     """
     if name == "A" :
         if self.getCoefficient("k").isEmpty():
              return escore.Data(numpy.identity(self.getDim()),escore.Function(self.getDomain()))
         else:
              return escore.Data(numpy.identity(self.getDim()),escore.Function(self.getDomain()))*self.getCoefficient("k")
     elif name == "D" :
         return self.getCoefficient("omega")
     elif name == "Y" :
         return self.getCoefficient("f")
     elif name == "d" :
         return self.getCoefficient("alpha")
     elif name == "y" :
         return self.getCoefficient("g")
     elif name == "Y_reduced" :
         return self.getCoefficient("f_reduced")
     elif name == "y_reduced" :
        return self.getCoefficient("g_reduced")
     else:
        return super(Helmholtz, self).getCoefficient(name)

class VTIWavePDE(LinearPDE):
    """
    A class specifically for waves, passes along values to native implementation
    to save computational time.
    """
    def __init__(self,domain,c,numEquations=None,numSolutions=None,debug=False):
        """
        Initializes a new linear PDE.

        :param domain: domain of the PDE
        :type domain: `Domain`
        :param numEquations: number of equations. If ``None`` the number of
                          equations is extracted from the PDE coefficients.
        :param numSolutions: number of solution components. If ``None`` the number
                          of solution components is extracted from the PDE
                          coefficients.
        :param debug: if True debug information is printed

        """
        super(VTIWavePDE, self).__init__(domain,numEquations,numSolutions,debug)
        #
        #   the coefficients of the PDE:
        #
        self.introduceCoefficients(
           A=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
           B=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
           C=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
           D=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
           du=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM),PDECoef.RIGHTHANDSIDE),
           Y=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
           d=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
           y=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
           d_dirac=PDECoef(PDECoef.DIRACDELTA,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
           y_dirac=PDECoef(PDECoef.DIRACDELTA,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
           r=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.RIGHTHANDSIDE),

#X=PDECoef(PDECoef.INTERIOR, (PDECoef.BY_EQUATION,PDECoef.BY_DIM),PDECoef.RIGHTHANDSIDE),

           q=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.BOTH))

        domain.setAssembler("WaveAssembler", c)

    
    def getSystem(self):
        """
        Returns the operator and right hand side of the PDE.

        :return: the discrete version of the PDE
        :rtype: ``tuple`` of `Operator` and
                `Data`
        """
        if not self.isOperatorValid() or not self.isRightHandSideValid():
          if self.isUsingLumping():
              if not self.isOperatorValid():
                 if not self.getFunctionSpaceForEquation()==self.getFunctionSpaceForSolution():
                      raise TypeError("Lumped matrix requires same order for equations and unknowns")
                 if not self.getCoefficient("A").isEmpty():
                      raise ValueError("coefficient A in lumped matrix may not be present.")
                 if not self.getCoefficient("B").isEmpty():
                      raise ValueError("coefficient B in lumped matrix may not be present.")
                 if not self.getCoefficient("C").isEmpty():
                      raise ValueError("coefficient C in lumped matrix may not be present.")

                 D=self.getCoefficient("D")
                 d=self.getCoefficient("d")
                 d_dirac=self.getCoefficient("d_dirac")
                 
                 if not D.isEmpty():
                     if self.getNumSolutions()>1:
                        D_times_e=util.matrix_mult(D,numpy.ones((self.getNumSolutions(),)))
                     else:
                        D_times_e=D
                 else:
                    D_times_e=escore.Data()
                 if not d.isEmpty():
                     if self.getNumSolutions()>1:
                        d_times_e=util.matrix_mult(d,numpy.ones((self.getNumSolutions(),)))
                     else:
                        d_times_e=d
                 else:
                    d_times_e=escore.Data()
                    
                 if not d_dirac.isEmpty():
                     if self.getNumSolutions()>1:
                        d_dirac_times_e=util.matrix_mult(d_dirac,numpy.ones((self.getNumSolutions(),)))
                 else:
                    d_dirac_times_e=escore.Data()

                 self.resetOperator()
                 operator=self.getCurrentOperator()
                 if hasattr(self.getDomain(), "addPDEToLumpedSystem") :
                    hrz_lumping=( self.getSolverOptions().getSolverMethod() ==  SolverOptions.HRZ_LUMPING )
                    self.getDomain().addPDEToLumpedSystem(operator, D_times_e, d_times_e, d_dirac_times_e,  hrz_lumping )
                 else:
                    self.getDomain().addToRHS(operator, 
                        [("Y", D_times_e), ("y", d_times_e),
                         ("y_dirac", d_dirac_times_e)])
                 self.trace("New lumped operator has been built.")
              if not self.isRightHandSideValid():
                 self.resetRightHandSide()
                 righthandside=self.getCurrentRightHandSide()
                 self.getDomain().addToRHS(righthandside,
                                [(i, self.getCoefficient(i)) for i in 
                                    ["du", "Y", "y", "y_dirac"]
                                ])
                 self.trace("New right hand side has been built.")
                 self.validRightHandSide()
              self.insertConstraint(rhs_only=False)
              self.validOperator()
          else:
             if not self.isOperatorValid() and not self.isRightHandSideValid():
                 self.resetRightHandSide()
                 righthandside=self.getCurrentRightHandSide()
                 self.resetOperator()
                 operator=self.getCurrentOperator()
                 data = [(i, self.getCoefficient(i)) for i in ["A", "B", "C",
                                "D", "Y", "d", "y", "d_contact",
                                "y_contact", "d_dirac", "y_dirac", "du"]
                            ]
                 self.getDomain().addToSystem(operator, righthandside, data)
                 self.insertConstraint(rhs_only=False)
                 self.trace("New system has been built.")
                 self.validOperator()
                 self.validRightHandSide()
             elif not self.isRightHandSideValid():
                 self.resetRightHandSide()
                 righthandside=self.getCurrentRightHandSide()
                 self.getDomain().addToRHS(righthandside,
                                [(i, self.getCoefficient(i)) for i in 
                                    ["du", "Y", "y", "y_contact", "y_dirac"]
                                ])
                 self.insertConstraint(rhs_only=True)
                 self.trace("New right hand side has been built.")
                 self.validRightHandSide()
             elif not self.isOperatorValid():
                 self.resetOperator()
                 operator=self.getCurrentOperator()
                 data = [(i, self.getCoefficient(i)) for i in ["A", "B", "C",
                        "D", "d", "d_contact", "d_dirac", "du"]]
                 self.getDomain().addToSystem(operator, escore.Data(), data)
                 self.insertConstraint(rhs_only=False)
                 self.trace("New operator has been built.")
                 self.validOperator()
        self.setSystemStatus()
        self.trace("System status is %s."%self.getSystemStatus())
        return (self.getCurrentOperator(), self.getCurrentRightHandSide())


class LameEquation(LinearPDE):
   """
   Class to define a Lame equation problem. This problem is defined as:

   *-grad(mu*(grad(u[i])[j]+grad(u[j])[i]))[j] - grad(lambda*grad(u[k])[k])[j] = F_i -grad(sigma[ij])[j]*

   with natural boundary conditions:

   *n[j]*(mu*(grad(u[i])[j]+grad(u[j])[i]) + lambda*grad(u[k])[k]) = f_i +n[j]*sigma[ij]*

   and constraints:

   *u[i]=r[i]* where *q[i]>0*

   """

   def __init__(self,domain,debug=False):
      """
      Initializes a new Lame equation.

      :param domain: domain of the PDE
      :type domain: `Domain`
      :param debug: if True debug information is printed

      """
      super(LameEquation, self).__init__(domain,\
                                         domain.getDim(),domain.getDim(),debug)
      self.introduceCoefficients(lame_lambda=PDECoef(PDECoef.INTERIOR,(),PDECoef.OPERATOR),
                                 lame_mu=PDECoef(PDECoef.INTERIOR,(),PDECoef.OPERATOR),
                                 F=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
                                 sigma=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM),PDECoef.RIGHTHANDSIDE),
                                 f=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE))
      self.setSymmetryOn()

   def setValues(self,**coefficients):
     """
     Sets new values to coefficients.

     :param coefficients: new values assigned to coefficients
     :keyword lame_mu: value for coefficient *mu*
     :type lame_mu: any type that can be cast to a `Scalar`
                    object on `Function`
     :keyword lame_lambda: value for coefficient *lambda*
     :type lame_lambda: any type that can be cast to a `Scalar`
                        object on `Function`
     :keyword F: value for internal force *F*
     :type F: any type that can be cast to a `Vector` object
              on `Function`
     :keyword sigma: value for initial stress *sigma*
     :type sigma: any type that can be cast to a `Tensor`
                  object on `Function`
     :keyword f: value for external force *f*
     :type f: any type that can be cast to a `Vector` object
              on `FunctionOnBoundary`
     :keyword r: prescribed values *r* for the solution in constraints
     :type r: any type that can be cast to a `Vector` object
              on `Solution` or
              `ReducedSolution` depending on whether
              reduced order is used for the representation of the equation
     :keyword q: mask for the location of constraints
     :type q: any type that can be cast to a `Vector` object
              on `Solution` or
              `ReducedSolution` depending on whether
              reduced order is used for the representation of the equation
     :raise IllegalCoefficient: if an unknown coefficient keyword is used
     """
     super(LameEquation, self).setValues(**coefficients)

   def getCoefficient(self,name):
     """
     Returns the value of the coefficient ``name`` of the general PDE.

     :param name: name of the coefficient requested
     :type name: ``string``
     :return: the value of the coefficient ``name``
     :rtype: `Data`
     :raise IllegalCoefficient: invalid coefficient name
     """
     out =self.createCoefficient("A")
     if name == "A" :
         if self.getCoefficient("lame_lambda").isEmpty(): 
            if self.getCoefficient("lame_mu").isEmpty():
                pass
            else:
                for i in range(self.getDim()):
                  for j in range(self.getDim()):
                    out[i,j,j,i] += self.getCoefficient("lame_mu")
                    out[i,j,i,j] += self.getCoefficient("lame_mu")
         else: 
            if self.getCoefficient("lame_mu").isEmpty():
                for i in range(self.getDim()):
                  for j in range(self.getDim()):
                    out[i,i,j,j] += self.getCoefficient("lame_lambda")
            else:
                for i in range(self.getDim()):
                  for j in range(self.getDim()):
                    out[i,i,j,j] += self.getCoefficient("lame_lambda")
                    out[i,j,j,i] += self.getCoefficient("lame_mu")
                    out[i,j,i,j] += self.getCoefficient("lame_mu")
         return out
     elif name == "X" :
         return self.getCoefficient("sigma")
     elif name == "Y" :
         return self.getCoefficient("F")
     elif name == "y" :
         return self.getCoefficient("f")
     else:
        return super(LameEquation, self).getCoefficient(name)


def LinearSinglePDE(domain,debug=False):
   """
   Defines a single linear PDE.

   :param domain: domain of the PDE
   :type domain: `Domain`
   :param debug: if True debug information is printed
   :rtype: `LinearPDE`
   """
   return LinearPDE(domain,numEquations=1,numSolutions=1,debug=debug)

def LinearPDESystem(domain,debug=False):
   """
   Defines a system of linear PDEs.

   :param domain: domain of the PDEs
   :type domain: `Domain`
   :param debug: if True debug information is printed
   :rtype: `LinearPDE`
   """
   return LinearPDE(domain,numEquations=domain.getDim(),numSolutions=domain.getDim(),debug=debug)


class TransportPDE(LinearProblem):
   """
   This class is used to define a transport problem given by a general linear,
   time dependent, second order PDE for an unknown, non-negative,
   time-dependent function *u* on a given domain defined through a
   `Domain` object.

   For a single equation with a solution with a single component the transport
   problem is defined in the following form:

   *(M+M_reduced)*u_t=-(grad(A[j,l]+A_reduced[j,l]) * grad(u)[l]+(B[j]+B_reduced[j])u)[j]+(C[l]+C_reduced[l])*grad(u)[l]+(D+D_reduced)-grad(X+X_reduced)[j,j]+(Y+Y_reduced)*

   where *u_t* denotes the time derivative of *u* and *grad(F)* denotes the
   spatial derivative of *F*. Einstein's summation convention,  ie. summation
   over indexes appearing twice in a term of a sum performed, is used.
   The coefficients *M*, *A*, *B*, *C*, *D*, *X* and *Y* have to be
   specified through `Data` objects in `Function`
   and the coefficients *M_reduced*, *A_reduced*, *B_reduced*, *C_reduced*,
   *D_reduced*, *X_reduced* and *Y_reduced* have to be specified through
   `Data` objects in `ReducedFunction`.
   It is also allowed to use objects that can be converted into such
   `Data` objects. *M* and *M_reduced* are scalar, *A* and
   *A_reduced* are rank two, *B*, *C*, *X*, *B_reduced*, *C_reduced* and
   *X_reduced* are rank one and *D*, *D_reduced*, *Y* and *Y_reduced* are
   scalar.

   The following natural boundary conditions are considered:

   *n[j]*((A[i,j]+A_reduced[i,j])*grad(u)[l]+(B+B_reduced)[j]*u+X[j]+X_reduced[j])+(d+d_reduced)*u+y+y_reduced=(m+m_reduced)*u_t*

   where *n* is the outer normal field. Notice that the coefficients *A*,
   *A_reduced*, *B*, *B_reduced*, *X* and *X_reduced* are defined in the
   transport problem. The coefficients *m*, *d* and *y* are each a scalar in
   `FunctionOnBoundary` and the coefficients
   *m_reduced*, *d_reduced* and *y_reduced* are each a scalar in
   `ReducedFunctionOnBoundary`.

   Constraints for the solution prescribing the value of the solution at
   certain locations in the domain have the form

   *u_t=r* where *q>0*

   *r* and *q* are each scalar where *q* is the characteristic function
   defining where the constraint is applied. The constraints override any other
   condition set by the transport problem or the boundary condition.

   The transport problem is symmetrical if

   *A[i,j]=A[j,i]* and *B[j]=C[j]* and *A_reduced[i,j]=A_reduced[j,i]* and
   *B_reduced[j]=C_reduced[j]*

   For a system and a solution with several components the transport problem
   has the form

   *(M[i,k]+M_reduced[i,k]) * u[k]_t=-grad((A[i,j,k,l]+A_reduced[i,j,k,l]) * grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k]) * u[k])[j]+(C[i,k,l]+C_reduced[i,k,l]) * grad(u[k])[l]+(D[i,k]+D_reduced[i,k] * u[k]-grad(X[i,j]+X_reduced[i,j])[j]+Y[i]+Y_reduced[i]*

   *A* and *A_reduced* are of rank four, *B*, *B_reduced*, *C* and
   *C_reduced* are each of rank three, *M*, *M_reduced*, *D*, *D_reduced*,
   *X_reduced* and *X* are each of rank two and *Y* and *Y_reduced* are of
   rank one. The natural boundary conditions take the form:

   *n[j]*((A[i,j,k,l]+A_reduced[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k]+X[i,j]+X_reduced[i,j])+(d[i,k]+d_reduced[i,k])*u[k]+y[i]+y_reduced[i]= (m[i,k]+m_reduced[i,k])*u[k]_t*

   The coefficient *d* and *m* are of rank two and *y* is of rank one with
   all in `FunctionOnBoundary`. The coefficients
   *d_reduced* and *m_reduced* are of rank two and *y_reduced* is of rank
   one all in `ReducedFunctionOnBoundary`.

   Constraints take the form

   *u[i]_t=r[i]* where *q[i]>0*

   *r* and *q* are each rank one. Notice that at some locations not
   necessarily all components must have a constraint.

   The transport problem is symmetrical if

      - *M[i,k]=M[i,k]*
      - *M_reduced[i,k]=M_reduced[i,k]*
      - *A[i,j,k,l]=A[k,l,i,j]*
      - *A_reduced[i,j,k,l]=A_reduced[k,l,i,j]*
      - *B[i,j,k]=C[k,i,j]*
      - *B_reduced[i,j,k]=C_reduced[k,i,j]*
      - *D[i,k]=D[i,k]*
      - *D_reduced[i,k]=D_reduced[i,k]*
      - *m[i,k]=m[k,i]*
      - *m_reduced[i,k]=m_reduced[k,i]*
      - *d[i,k]=d[k,i]*
      - *d_reduced[i,k]=d_reduced[k,i]*
      - *d_dirac[i,k]=d_dirac[k,i]*

   `TransportPDE` also supports solution discontinuities over a contact region
   in the domain. To specify the conditions across the discontinuity we are
   using the generalised flux *J* which, in the case of a system of PDEs and
   several components of the solution, is defined as

   *J[i,j]=(A[i,j,k,l]+A_reduced[[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k]+X[i,j]+X_reduced[i,j]*

   For the case of single solution component and single PDE *J* is defined as

   *J[j]=(A[i,j]+A_reduced[i,j])*grad(u)[j]+(B[i]+B_reduced[i])*u+X[i]+X_reduced[i]*

   In the context of discontinuities *n* denotes the normal on the
   discontinuity pointing from side 0 towards side 1 calculated from
   `FunctionSpace.getNormal` of `FunctionOnContactZero`.
   For a system of transport problems the contact condition takes the form

   *n[j]*J0[i,j]=n[j]*J1[i,j]=(y_contact[i]+y_contact_reduced[i])- (d_contact[i,k]+d_contact_reduced[i,k])*jump(u)[k]*

   where *J0* and *J1* are the fluxes on side 0 and side 1 of the
   discontinuity, respectively. *jump(u)*, which is the difference of the
   solution at side 1 and at side 0, denotes the jump of *u* across
   discontinuity along the normal calculated by `jump`.
   The coefficient *d_contact* is of rank two and *y_contact* is of rank one
   both in `FunctionOnContactZero` or `FunctionOnContactOne`.
   The coefficient *d_contact_reduced* is of rank two and *y_contact_reduced*
   is of rank one both in `ReducedFunctionOnContactZero` or `ReducedFunctionOnContactOne`.
   In case of a single PDE and a single component solution the contact
   condition takes the form

   *n[j]*J0_{j}=n[j]*J1_{j}=(y_contact+y_contact_reduced)-(d_contact+y_contact_reduced)*jump(u)*

   In this case the coefficient *d_contact* and *y_contact* are each scalar
   both in `FunctionOnContactZero` or
   `FunctionOnContactOne` and the coefficient
   *d_contact_reduced* and *y_contact_reduced* are each scalar both in
   `ReducedFunctionOnContactZero` or
   `ReducedFunctionOnContactOne`.

   Typical usage::

       p = TransportPDE(dom)
       p.setValue(M=1., C=[-1.,0.])
       p.setInitialSolution(u=exp(-length(dom.getX()-[0.1,0.1])**2)
       t = 0
       dt = 0.1
       while (t < 1.):
           u = p.solve(dt)

   """
   def __init__(self,domain,numEquations=None,numSolutions=None, useBackwardEuler=None, debug=False):
     """
     Initializes a transport problem.

     :param domain: domain of the PDE
     :type domain: `Domain`
     :param numEquations: number of equations. If ``None`` the number of
                          equations is extracted from the coefficients.
     :param numSolutions: number of solution components. If ``None`` the number
                          of solution components is extracted from the
                          coefficients.
     :param debug: if True debug information is printed
     """
     super(TransportPDE, self).__init__(domain,numEquations,numSolutions,debug)
     #
     #   the coefficients of the transport problem
     #
     self.introduceCoefficients(
       M=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       A=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       B=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       C=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       D=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       X=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM),PDECoef.RIGHTHANDSIDE),
       Y=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       m=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       d=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       d_contact=PDECoef(PDECoef.CONTACT,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y_contact=PDECoef(PDECoef.CONTACT,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       M_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       A_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       B_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       C_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       D_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       X_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_DIM),PDECoef.RIGHTHANDSIDE),
       Y_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       m_reduced=PDECoef(PDECoef.BOUNDARY_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       d_reduced=PDECoef(PDECoef.BOUNDARY_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y_reduced=PDECoef(PDECoef.BOUNDARY_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       d_contact_reduced=PDECoef(PDECoef.CONTACT_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y_contact_reduced=PDECoef(PDECoef.CONTACT_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       d_dirac=PDECoef(PDECoef.DIRACDELTA,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y_dirac=PDECoef(PDECoef.DIRACDELTA,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       r=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.RIGHTHANDSIDE),
       q=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.BOTH) )
     if not useBackwardEuler == None:
        import warnings
        warnings.warn("Argument useBackwardEuler has expired and will be removed in a later release. Please use SolverOptions.setODESolver() instead.", PendingDeprecationWarning, stacklevel=2)
        if useBackwardEuler: self.getSolverOptions().setODESolver(SolverOptions.BACKWARD_EULER)

   def __str__(self):
     """
     Returns the string representation of the transport problem.

     :return: a simple representation of the transport problem
     :rtype: ``str``
     """
     return "<TransportPDE %d>"%id(self)

   def checkSymmetry(self,verbose=True):
      """
      Tests the transport problem for symmetry.

      :param verbose: if set to True or not present a report on coefficients
                      which break the symmetry is printed.
      :type verbose: ``bool``
      :return:  True if the PDE is symmetric
      :rtype: ``bool``
      :note: This is a very expensive operation. It should be used for
             degugging only! The symmetry flag is not altered.
      """
      out=True
      out=out and self.checkSymmetricTensor("M", verbose)
      out=out and self.checkSymmetricTensor("M_reduced", verbose)
      out=out and self.checkSymmetricTensor("A", verbose)
      out=out and self.checkSymmetricTensor("A_reduced", verbose)
      out=out and self.checkReciprocalSymmetry("B","C", verbose)
      out=out and self.checkReciprocalSymmetry("B_reduced","C_reduced", verbose)
      out=out and self.checkSymmetricTensor("D", verbose)
      out=out and self.checkSymmetricTensor("D_reduced", verbose)
      out=out and self.checkSymmetricTensor("m", verbose)
      out=out and self.checkSymmetricTensor("m_reduced", verbose)
      out=out and self.checkSymmetricTensor("d", verbose)
      out=out and self.checkSymmetricTensor("d_reduced", verbose)
      out=out and self.checkSymmetricTensor("d_contact", verbose)
      out=out and self.checkSymmetricTensor("d_contact_reduced", verbose)
      out=out and self.checkSymmetricTensor("d_dirac", verbose)
      return out

   def setValue(self,**coefficients):
      """
      Sets new values to coefficients.

      :param coefficients: new values assigned to coefficients
      :keyword M: value for coefficient ``M``
      :type M: any type that can be cast to a `Data` object on
               `Function`
      :keyword M_reduced: value for coefficient ``M_reduced``
      :type M_reduced: any type that can be cast to a `Data`
                       object on `Function`
      :keyword A: value for coefficient ``A``
      :type A: any type that can be cast to a `Data` object on
               `Function`
      :keyword A_reduced: value for coefficient ``A_reduced``
      :type A_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword B: value for coefficient ``B``
      :type B: any type that can be cast to a `Data` object on
               `Function`
      :keyword B_reduced: value for coefficient ``B_reduced``
      :type B_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword C: value for coefficient ``C``
      :type C: any type that can be cast to a `Data` object on
               `Function`
      :keyword C_reduced: value for coefficient ``C_reduced``
      :type C_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword D: value for coefficient ``D``
      :type D: any type that can be cast to a `Data` object on
               `Function`
      :keyword D_reduced: value for coefficient ``D_reduced``
      :type D_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword X: value for coefficient ``X``
      :type X: any type that can be cast to a `Data` object on
               `Function`
      :keyword X_reduced: value for coefficient ``X_reduced``
      :type X_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword Y: value for coefficient ``Y``
      :type Y: any type that can be cast to a `Data` object on
               `Function`
      :keyword Y_reduced: value for coefficient ``Y_reduced``
      :type Y_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword m: value for coefficient ``m``
      :type m: any type that can be cast to a `Data` object on
               `FunctionOnBoundary`
      :keyword m_reduced: value for coefficient ``m_reduced``
      :type m_reduced: any type that can be cast to a `Data`
                       object on `FunctionOnBoundary`
      :keyword d: value for coefficient ``d``
      :type d: any type that can be cast to a `Data` object on
               `FunctionOnBoundary`
      :keyword d_reduced: value for coefficient ``d_reduced``
      :type d_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunctionOnBoundary`
      :keyword y: value for coefficient ``y``
      :type y: any type that can be cast to a `Data` object on
               `FunctionOnBoundary`
      :keyword d_contact: value for coefficient ``d_contact``
      :type d_contact: any type that can be cast to a `Data`
                       object on `FunctionOnContactOne` or `FunctionOnContactZero`
      :keyword d_contact_reduced: value for coefficient ``d_contact_reduced``
      :type d_contact_reduced: any type that can be cast to a `Data` object on `ReducedFunctionOnContactOne` or `ReducedFunctionOnContactZero`
      :keyword y_contact: value for coefficient ``y_contact``
      :type y_contact: any type that can be cast to a `Data`
                       object on `FunctionOnContactOne` or `FunctionOnContactZero`
      :keyword y_contact_reduced: value for coefficient ``y_contact_reduced``
      :type y_contact_reduced: any type that can be cast to a `Data` object on `ReducedFunctionOnContactOne` or `ReducedFunctionOnContactZero`
      
      :keyword d_dirac: value for coefficient ``d_dirac``
      :type d_dirac: any type that can be cast to a `Data` object on `DiracDeltaFunctions`
      :keyword y_dirac: value for coefficient ``y_dirac``
      :type y_dirac: any type that can be cast to a `Data` object on `DiracDeltaFunctions`

      :keyword r: values prescribed to the solution at the locations of constraints
      :type r: any type that can be cast to a `Data` object on
               `Solution` or `ReducedSolution`
               depending on whether reduced order is used for the solution
      :keyword q: mask for the location of constraints
      :type q: any type that can be cast to a `Data` object on
               `Solution` or
               `ReducedSolution` depending on whether
               reduced order is used for the representation of the equation
      :raise IllegalCoefficient: if an unknown coefficient keyword is used
      """
      super(TransportPDE,self).setValue(**coefficients)

   def createOperator(self):
       """
       Returns an instance of a new transport operator.
       """
       optype=self.getRequiredOperatorType()
       self.trace("New Transport problem pf type %s is allocated."%optype)
       return self.getDomain().newTransportProblem( \
                               self.getNumEquations(), \
                               self.getFunctionSpaceForSolution(), \
                               optype)


   def getRequiredOperatorType(self):
      """
      Returns the system type which needs to be used by the current set up.

      :return: a code to indicate the type of transport problem scheme used
      :rtype: ``float``
      """
      solver_options=self.getSolverOptions()
      return self.getDomain().getTransportTypeId(solver_options.getSolverMethod(), solver_options.getPreconditioner(),solver_options.getPackage(), solver_options.isSymmetric())

   def getUnlimitedTimeStepSize(self):
      """
      Returns the value returned by the ``getSafeTimeStepSize`` method to
      indicate no limit on the safe time step size.

       :return: the value used to indicate that no limit is set to the time
                step size
       :rtype: ``float``
       :note: Typically the biggest positive float is returned
      """
      return self.getOperator().getUnlimitedTimeStepSize()

   def getSafeTimeStepSize(self):
       """
       Returns a safe time step size to do the next time step.

       :return: safe time step size
       :rtype: ``float``
       :note: If not ``getSafeTimeStepSize()`` < ``getUnlimitedTimeStepSize()``
              any time step size can be used.
       """
       return self.getOperator().getSafeTimeStepSize()

   #====================================================================
   def getSolution(self, dt=None, u0=None):
       """
       Returns the solution by marching forward by time step dt. if ''u0'' is present,
       ''u0'' is used as the initial value otherwise the solution from the last call is used.

       :param dt: time step size. If ``None`` the last solution is returned.
       :type dt: positive ``float`` or ``None``
       :param u0: new initial solution or ``None``
       :type u0: any object that can be interpolated to a `Data`
                object on `Solution` or `ReducedSolution`
       :return: the solution
       :rtype: `Data`
       """
       if not dt == None:
          option_class=self.getSolverOptions()
          if dt<=0:
              raise ValueError("step size needs to be positive.")
          if u0 == None:
              u0=self.getCurrentSolution()
          else:
              u0=util.interpolate(u0,self.getFunctionSpaceForSolution())
              if self.getNumSolutions() == 1:
                if u0.getShape()!=():
                  raise ValueError("Illegal shape %s of initial solution."%(u0.getShape(),))
              else:
                 if u0.getShape()!=(self.getNumSolutions(),):
                   raise ValueError("Illegal shape %s of initial solution."%(u0.getShape(),))
          self.setSolution(self.getOperator().solve(u0, self.getRightHandSide(),dt,option_class))
          self.validSolution()
       return self.getCurrentSolution()

   def setInitialSolution(self,u):
       """
       Sets the initial solution.

       :param u: initial solution
       :type u: any object that can be interpolated to a `Data`
                object on `Solution` or `ReducedSolution`
       """
       u2=util.interpolate(u,self.getFunctionSpaceForSolution())
       if self.getNumSolutions() == 1:
          if u2.getShape()!=():
              raise ValueError("Illegal shape %s of initial solution."%(u2.getShape(),))
       else:
          if u2.getShape()!=(self.getNumSolutions(),):
              raise ValueError("Illegal shape %s of initial solution."%(u2.getShape(),))
       self.setSolution(u2,validate=False)


   def getSystem(self):
       """
       Returns the operator and right hand side of the PDE.

       :return: the discrete version of the PDE
       :rtype: ``tuple`` of `Operator` and
               `Data`

       """
       if not self.isOperatorValid() or not self.isRightHandSideValid():
          self.resetRightHandSide()
          righthandside=self.getCurrentRightHandSide()
          self.resetOperator()
          operator=self.getCurrentOperator()
          self.getDomain().addPDEToTransportProblem(
                            operator,
                            righthandside,
                            self.getCoefficient("M"),
                            self.getCoefficient("A"),
                            self.getCoefficient("B"),
                            self.getCoefficient("C"),
                            self.getCoefficient("D"),
                            self.getCoefficient("X"),
                            self.getCoefficient("Y"),
                            self.getCoefficient("d"),
                            self.getCoefficient("y"),
                            self.getCoefficient("d_contact"),
                            self.getCoefficient("y_contact"),
                            self.getCoefficient("d_dirac"),
                            self.getCoefficient("y_dirac") )
          self.getDomain().addPDEToTransportProblem(
                            operator,
                            righthandside,
                            self.getCoefficient("M_reduced"),
                            self.getCoefficient("A_reduced"),
                            self.getCoefficient("B_reduced"),
                            self.getCoefficient("C_reduced"),
                            self.getCoefficient("D_reduced"),
                            self.getCoefficient("X_reduced"),
                            self.getCoefficient("Y_reduced"),
                            self.getCoefficient("d_reduced"),
                            self.getCoefficient("y_reduced"),
                            self.getCoefficient("d_contact_reduced"),
                            self.getCoefficient("y_contact_reduced"),
                            escore.Data(),
                            escore.Data() )
          operator.insertConstraint(righthandside,self.getCoefficient("q"),self.getCoefficient("r"))
          self.trace("New system has been built.")
          self.validOperator()
          self.validRightHandSide()
       self.setSystemStatus()
       self.trace("System status is %s."%self.getSystemStatus())
       return (self.getCurrentOperator(), self.getCurrentRightHandSide())

   def setDebug(self, flag):
     """
     Switches debug output on if ``flag`` is True,
     otherwise it is switched off.

     :param flag: desired debug status
     :type flag: ``bool``
     """
     if flag:
         self.setDebugOn()
     else:
         self.setDebugOff()

   def setDebugOn(self):
     """
     Switches debug output on.
     """
     super(TransportPDE,self).setDebugOn()
     
   def setDebugOff(self):
     """
     Switches debug output off.
     """
     super(TransportPDE,self).setDebugOff()
     
def SingleTransportPDE(domain, debug=False):
   """
   Defines a single transport problem

   :param domain: domain of the PDE
   :type domain: `Domain`
   :param debug: if True debug information is printed
   :rtype: `TransportPDE`
   """
   return TransportPDE(domain,numEquations=1,numSolutions=1, debug=debug)
