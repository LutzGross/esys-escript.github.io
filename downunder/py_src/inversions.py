
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

import logging

from esys.escript import *
from esys.weipa import createDataset

from costfunctions import SimpleCostFunction
from forwardmodels import GravityModel, ForwardModel
from mappings import *
from minimizers import *
from regularizations import Regularization
from datasources import DataSource

class InversionBase(object):
    """
    Base class for inversions.
    """
    def __init__(self):
        self.logger=logging.getLogger('inv.%s'%self.__class__.__name__)
        self._solver_callback = None
        self._solver_opts = {}
        self._solver_tol = 1e-9
        self._solver_maxiter = 200
        # use identity mapping by default
        self.setMapping()
        self.setSolverClass()        
        self.__domain=None
        self.__regularization=None
        self.__forwardmodel=None
        
    def setDomain(self, domain):
        """
        sets the domain of the inversion
        
        :param domain: domain of the inversion
        :type domain: `Domain`
        """
        self.__domain=domain
        
    def getDomain(self):
        """
        returns the domain of the inversion
        
        :rtype domain: domain of the inversion
        """
        return  self.__domain
        
    def setMapping(self, mapping=None):
        """
        Sets the mapping oobject to map between model parameters and the data.
        If no mapping is provided, an identity mapping is used (`ScalingMapping`
        with constant 1).
        
        :param mapping: parameter mapping
        :type mapping: `Mapping`
        """
        if mapping == None:
	     self.__mapping=ScalingMapping(1)
	else:
             self.__mapping=mapping
             
    def getMapping(self):
        """
        return the mapping(s) used in the inversion
        
        :rtype: `Mapping`
        """
        return self.__mapping
        
        
    def setRegularization(self, regularization):
        """
        sets the regularization for the inversion. 
        
        :param regularization: regularization
        :type regularization: `Regularization`
        """
        self.__regularization=regularization
        
    def getRegularization(self):
        """
        returns the regularization method(s)
        :rtype: `Regularization`
        """
        return self.__regularization
        
    def setForwardModel(self, forwardmodel):
        """
        sets the forward model(s) for the inversion. 
        
        :param forwardmodel: forward model
        :type forwardmodel: `ForwardModel`
        """
        self.__forwardmodel=forwardmodel
        
    def getForwardModel(self):
        """
        returns the forward model
        :rtype: `ForwardModel`
        """
        return self.__forwardmodel
        
    def setSolverClass(self, solverclass=None):
        """
        The solver to be used in the inversion process. See the minimizers
        module for available solvers. By default, the L-BFGS minimizer is used.
        """
        if solverclass == None:
	   self.solverclass=MinimizerLBFGS
	else:
           self.solverclass=solverclass
            
    def setSolverCallback(self, callback):
        """
        Sets the callback function which is called after every solver iteration
        """
        self._solver_callback=callback

    def setSolverMaxIterations(self, maxiter):
        """
        Sets the maximum number of solver iterations to run
        """
        if maxiter>0:
           self._solver_maxiter=maxiter 
        else:
	   raise ValueError("maxiter must be positive.")

    def setSolverOptions(self, **opts):
        """
        Sets additional solver options. The valid options depend on the solver
        being used.
        """
        self._solver_opts.update(**opts)

    def setSolverTolerance(self, tol):
        """
        Sets the error tolerance for the solver. An acceptable solution is
        considered to be found once the tolerance is reached.
        """
        if tol>0:
           self._solver_tol=tol
        else:
	   raise ValueError("tolerance must be positive.")
        

    def isSetup(self):
        """
        return True if the inversion is set up.
        """
        raise NotImplementedError
        
    def run(self, *args):
        """
        This method starts the inversion process and must be overwritten.
        Preferably, users should be able to set a starting point so
        inversions can be continued after stopping.
        """
        raise NotImplementedError


class GravityInversion(InversionBase):
    """
    """
    def __init__(self):
        super(GravityInversion,self).__init__()
        self.__is_setup=False
        self.setWeights()
        
    def setWeights(self, mu_reg=1., mu_model=1.):
        """
        sets the weighting factors for the cost funtion
        """
        self._mu_reg=mu_reg
        self._mu_model=mu_model

    def isSetUp(self):
        """
        return True if the inversion is set up.
        """
        if self.getRegularization() and self.getMapping() and self.getForwardModel() and self.getDomain():
	    return True
	else:
	    return False

        
    def siloWriterCallback(self, k, x, fx, gfx):
        """
        callback function for solver to track solution
         
        :param k: iteration count
        :param x: current m approximation
        :param fx: value of cost function
        :param gfx: gradient of f at x
        """
        fn='inv.%d'%k
        ds=createDataset(rho=self.getMapping().getValue(x))
        ds.setCycleAndTime(k,k)
        ds.saveSilo(fn)
        self.logger.debug("Jreg(m) = %e"%self.regularization.getValue(x))
        self.logger.debug("f(m) = %e"%fx)

    def setup(self, source, rho_ref=None, w0=None, w=None):
        """
        sets up the inversion parameters from a downunder.DataSource.
        
        :param source: suitable data source
        :type source: `DataSource'
        :param rho_ref: reference density. If not specified, zero is used.
        :type rho_ref: ``float`` or `Scalar`
        :param w0: weighting factor for L2 term in the regularization
        :type rho_ref: ``float`` or `Scalar`
         :param w: weighting factor for H1 term in the regularization
        :type rho_ref: list of float or `Vector`
        """
        self.logger.info('Retrieving domain...')
        self.setDomain(source.getDomain())
        DIM=self.getDomain().getDim()

        if rho_ref == None: 
            rho_ref=0.
        if w==None: 
            w=[1.]*DIM

        self.logger.info("Retrieving density mask...")
        rho_mask = source.getSetDensityMask()
        m_ref=self.getMapping().getInverse(rho_ref)
        
        self.logger.info("Retrieving gravity and standard deviation data...")
        g, sigma=source.getGravityAndStdDev()
        chi=safeDiv(1., sigma*sigma)        
        self.logger.debug("g = %s"%g)
        self.logger.debug("sigma = %s"%sigma)
        self.logger.debug("chi = %s"%chi)
        chi=chi*kronecker(DIM)[DIM-1]
        
        self.logger.info("Set up regularization and forward model...")
        self.setRegularization(Regularization(self.getDomain(), m_ref=m_ref, w0=w0[:DIM], w=w, location_of_set_m=rho_mask))
        self.setForwardModel(GravityModel(self.getDomain(), chi, g))
            
        x=self.getDomain().getX()
        l=0
        for i in range(DIM-1):
              l=max(l, sup(x[i])-inf(x[i]))
        G=6.6742e-11
        mu_reg=0.5*(l*l*G)**2
        mu_model=1.
        #CIHAN: reg is not used!
        self.logger.debug("Set weighting factors...")    
        self.logger.debug("mu_reg = %s"%mu_reg)
        self.logger.debug("mu_model = %s"%mu_model)
        self.setWeights(mu_reg=mu_reg, mu_model=mu_model)

    def run(self, rho_init=0.):
      
        if not self.isSetUp():
	    raise ValueError("Inversion is no setup properly.")
        
        f=SimpleCostFunction(self.getRegularization(), self.getMapping(), self.getForwardModel())
        f.setWeights(mu_reg=self._mu_reg, mu_model=self._mu_model)
             
        solver=self.solverclass(f)
        solver.setCallback(self._solver_callback)
        solver.setMaxIterations(self._solver_maxiter)
        solver.setOptions(**self._solver_opts)
        solver.setTolerance(self._solver_tol)
        if not isinstance(rho_init, Data):
            rho_init=Scalar(rho_init, ContinuousFunction(self.getDomain()))
        m_init=self.getMapping().getInverse(rho_init)

        self.logger.info("Starting solver...")
        solver.run(m_init)
        m_star=solver.getResult()
        self.logger.info("m* = %s"%m_star)
        rho_star=self.getMapping().getValue(m_star)
        self.logger.info("rho* = %s"%rho_star)
        solver.logSummary()
        return rho_star

