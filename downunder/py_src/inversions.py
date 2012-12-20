
##############################################################################
#
# Copyright (c) 2003-2012 by University of Queensland
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

"""Higher-level classes that allow running inversions with minimal set up"""

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['InversionBase', 'GravityInversion','MagneticInversion']

import logging
from esys.escript import *
import esys.escript.unitsSI as U
from esys.weipa import createDataset

from .inversioncostfunctions import InversionCostFunction
from .forwardmodels import GravityModel, MagneticModel
from .mappings import *
from .minimizers import *
from .regularizations import Regularization
from .datasources import DataSource

class InversionBase(object):
    """
    Base class for running an inversion
    """
    def __init__(self):
        self.logger=logging.getLogger('inv.%s'%self.__class__.__name__)
        self.__costfunction = None
        self.setSolverMaxIterations()
        self.setSolverTolerance()
        self.setSolverCallback()
        self.setSolverClass()
        self._solver_opts = {}
        self.initial_value = None
        # set default solver


    def setCostFunction(self, costfunction):
        """
        sets the cost function of the inversion. This function needs to be called 
        before the inversion iteration can be started.

        :param costfunction: domain of the inversion
        :type costfunction: 'InversionCostFunction'
        """
        self.__costfunction=costfunction

    def getCostFunction(self):
        """
        returns the domain of the inversion

        :rtype: 'InversionCostFunction'
        """
        if self.isSetUp():
           return  self.__costfunction
        else:
	   raise RuntimeError("Inversion is not set up.")
    
    def setSolverClass(self, solverclass=None):
        """
        The solver to be used in the inversion process. See the minimizers
        module for available solvers.
        By default, the L-BFGS minimizer is used.
        """
        if solverclass == None:
            self.solverclass=MinimizerLBFGS
        else:
            self.solverclass=solverclass
            
    def isSetUp(self):
        """
        returns True if the inversion is set up and is ready to run.
        
        :rtype: `bool`
        """
        if  self.__costfunction:
	    return True
	else:
	    return False
	    
    def getDomain(self):
        """
        returns the domain of the inversion

        :rtype: `Domain`
        """
        self.getCostFunction().getDomain()
        
    def setSolverCallback(self, callback=None):
        """
        Sets the callback function which is called after every solver
        iteration.
        """
        self._solver_callback=callback
        
    def setSolverMaxIterations(self, maxiter=None):
        """
        Sets the maximum number of solver iterations to run.
        """
        if maxiter == None: maxiter = 200
        if maxiter>0:
            self._solver_maxiter=maxiter
        else:
            raise ValueError("maxiter must be positive.")

    def setSolverOptions(self, **opts):
        """
        Sets additional solver options. The valid options depend on the
        solver being used.
        """
        self._solver_opts.update(**opts)

    def setSolverTolerance(self, tol=None, atol=None):
	"""
	Sets the error tolerance for the solver. An acceptable solution is
	considered to be found once the tolerance is reached.
	
	:param tol: tolerance for changes to level set function. If `None` changes to the 
		    level set function are not checked for convergence during iteration.
	:param atol: tolerance for changes to cost function. If `None` changes to the 
		    cost function are not checked for convergence during iteration.
	:type tol: `float` or `None`
	:type atol: `float` or `None` 
	:note: if both arguments are equal to `None` the default setting tol=1e-4, atol=None is used.

	"""
	if tol == None and atol==None:
	    tol=1e-4

	self._solver_xtol=tol
	self._solver_atol=atol

 
    def setInitialGuess(self, *args):
        """
        set the initial guess for the inversion iteration. By default zero is used.
        
        """
        self.initial_value=self.getCostFunction().createLevelSetFunction(*args)
        
    def run(self):
        """
        this function runs the inversion. 
        
        
        """
        if not self.isSetUp():
            raise RuntimeError("Inversion is not setup.")           

        if self.initial_value == None: self.setInitialGuess()
        solver=self.solverclass(self.getCostFunction())
        solver.setCallback(self._solver_callback)
        solver.setMaxIterations(self._solver_maxiter)
        solver.setOptions(**self._solver_opts)
        solver.setTolerance(x_tol=self._solver_xtol)
        self.logger.info("Starting solver...")
        try:
            solver.run(self.initial_value)
            self.m=solver.getResult() 
            self.p=self.getCostFunction().getProperties(self.m)
        except MinimizerException as e:
            self.m=solver.getResult()
            self.p=self.getCostFunction().getProperties(self.m)
            self.logger.info("iteration failed.")
            raise e
        self.logger.info("iteration completed.")
        solver.logSummary()
        return self.p

    def setup(self, *args, **k_args): 
        """
        returns True if the inversion is set up and ready to run.
        """
        pass

class GravityInversion(InversionBase):
    """
    Inversion of Gravity (Bouguer) anomaly data.
    """
    def setup(self, domainbuilder,
                    rho0=None, drho=None, z0=None, beta=None,
                    w0=None, w1=None):
        """
        Sets up the inversion parameters from a `DomainBuilder`.

        :param domainbuilder: Domain builder object with gravity source(s)
        :type domainbuilder: `DomainBuilder`
        :param rho0: reference density. If not specified, zero is used.
        :type rho0: ``float`` or `Scalar`
        """
        self.logger.info('Retrieving domain...')
        dom=domainbuilder.getDomain()
        DIM=dom.getDim()
        #========================
        self.logger.info('Creating mapping...')
        rho_mapping=DensityMapping(dom, rho0=rho0, drho=drho, z0=z0, beta=beta)
        scale_mapping=rho_mapping.getTypicalDerivative()
        print " scale_mapping = ",scale_mapping
        #========================
        self.logger.info("Setting up regularization...")
        if w1 is None:
            w1=[1.]*DIM
        rho_mask = domainbuilder.getSetDensityMask()
        regularization=Regularization(dom, numLevelSets=1,\
                               w0=w0, w1=w1, location_of_set_m=rho_mask)
        #====================================================================
        self.logger.info("Retrieving gravity surveys...")
        surveys=domainbuilder.getGravitySurveys()
        g=[]
        w=[]
        for g_i,sigma_i in surveys:
            w_i=safeDiv(1., sigma_i)
            if g_i.getRank()==0:
                g_i=g_i*kronecker(DIM)[DIM-1]
            if w_i.getRank()==0:
                w_i=w_i*kronecker(DIM)[DIM-1]
            g.append(g_i)
            w.append(w_i)
            self.logger.debug("Added gravity survey:")
            self.logger.debug("g = %s"%g_i)
            self.logger.debug("sigma = %s"%sigma_i)
            self.logger.debug("w = %s"%w_i)
        #====================================================================

        self.logger.info("Setting up model...")
        forward_model=GravityModel(dom, w, g)
        forward_model.rescaleWeights(rho_scale=scale_mapping)

 
        #====================================================================
        self.logger.info("Setting cost function...")
        self.setCostFunction(InversionCostFunction(regularization, rho_mapping, forward_model))
        
    def setInitialGuess(self, rho=None):
        """
        set the initial guess *rho* for density the inversion iteration. If no *rho* present
        then an appropriate initial guess is chosen.
        
        :param rho: initial value for the density anomaly.
        :type rho: `Scalar`
        """
        if rho:
	    super(GravityInversion,self).setInitialGuess(rho)
	else:
	    super(GravityInversion,self).setInitialGuess()
        
    def siloWriterCallback(self, k, m, Jm, g_Jm):
        """
        callback function that can be used to track the solution

        :param k: iteration count
        :param m: current m approximation
        :param Jm: value of cost function
        :param g_Jm: gradient of f at x
        """
        fn='inv.%d'%k
        ds=createDataset(rho=self.getCostFunction().mappings[0].getValue(m))
        ds.setCycleAndTime(k,k)
        ds.saveSilo(fn)
        self.logger.debug("J(m) = %e"%Jm)

class MagneticInversion(InversionBase):
    """
    Inversion of magnetic data.
    """
    def setup(self, domainbuilder,
                    k0=None, dk=None, z0=None, beta=None,
                    w0=None, w1=None):
        """
        Sets up the inversion parameters from a `DomainBuilder`.

        :param domainbuilder: Domain builder object with magnetic source(s)
        :type domainbuilder: `DomainBuilder`

        """
        self.logger.info('Retrieving domain...')
        dom=domainbuilder.getDomain()
        DIM=dom.getDim()

        #========================
        self.logger.info('Creating mapping ...')
        susc_mapping=SusceptibilityMapping(dom, k0=k0, dk=dk, z0=z0, beta=beta)
        scale_mapping=susc_mapping.getTypicalDerivative()
        print " scale_mapping = ",scale_mapping
        #========================
        self.logger.info("Setting up regularization...")
        if w1 is None:
            w1=[1.]*DIM
        k_mask = domainbuilder.getSetSusceptibilityMask()
        regularization=Regularization(dom, numLevelSets=1,w0=w0, w1=w1, location_of_set_m=k_mask)

        #====================================================================
        self.logger.info("Retrieving magnetic field surveys...")
        surveys=domainbuilder.getMagneticSurveys()
        B=[]
        w=[]
        for B_i,sigma_i in surveys:
            w_i=safeDiv(1., sigma_i)
            if B_i.getRank()==0:
                B_i=B_i*kronecker(DIM)[DIM-1]
            if w_i.getRank()==0:
                w_i=w_i*kronecker(DIM)[DIM-1]
            B.append(B_i)
            w.append(w_i)
            self.logger.debug("Added magnetic survey:")
            self.logger.debug("B = %s"%B_i)
            self.logger.debug("sigma = %s"%sigma_i)
            self.logger.debug("w = %s"%w_i)
        #====================================================================
        self.logger.info("Setting up model...")
        forward_model=MagneticModel(dom, w, B, domainbuilder.getBackgroundMagneticFluxDensity())
        forward_model.rescaleWeights(k_scale=scale_mapping)

        #====================================================================
        self.logger.info("Setting cost function...")
        self.setCostFunction(InversionCostFunction(regularization, susc_mapping, forward_model))
        
    def setInitialGuess(self, k=None):
        """
        set the initial guess *k* for susceptibility for the inversion iteration. If no *k* present
        then an appropriate initial guess is chosen.
        
        :param k: initial value for the susceptibility anomaly.
        :type k: `Scalar`
        """
        if k:
	    super(MagneticInversion,self).setInitialGuess(k)
	else:
	    super(MagneticInversion,self).setInitialGuess()
	    
    def siloWriterCallback(self, k, m, Jm, g_Jm):
        """
        callback function that can be used to track the solution

        :param k: iteration count
        :param m: current m approximation
        :param Jm: value of cost function
        :param g_Jm: gradient of f at x
        """
        fn='inv.%d'%k
        ds=createDataset(susceptibility=self.getCostFunction().mappings[0].getValue(m))
        ds.setCycleAndTime(k,k)
        ds.saveSilo(fn)
        self.logger.debug("J(m) = %e"%Jm)
        