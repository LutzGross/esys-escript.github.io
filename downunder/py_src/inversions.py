
########################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import logging

from esys.escript import *
from esys.weipa import createDataset

from costfunctions import SimpleCostFunction
from forwardmodels import GravityModel
from mappings import *
from minimizers import *
from regularizations import Regularization

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
        self.mapping=ScalingMapping(1)
        self.source=None
        self.solverclass=MinimizerLBFGS

    def setSolverCallback(self, callback):
        """
        Sets the callback function which is called after every solver iteration
        """
        self._solver_callback=callback

    def setSolverMaxIterations(self, maxiter):
        """
        Sets the maximum number of solver iterations to run
        """
        self._solver_maxiter=maxiter

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
        self._solver_tol=tol

    def setSolverClass(self, solverclass):
        """
        The solver to be used in the inversion process. See the minimizers
        module for available solvers. By default, the L-BFGS minimizer is used.
        """
        self.solverclass=solverclass

    def setDataSource(self, source):
        """
        Sets the data source which is used to get the survey data to be
        inverted.
        """
        self.source=source

    def setMapping(self, mapping):
        """
        Sets the mapping class to map between model parameters and the data.
        If no mapping is provided, an identity mapping is used (ScalingMapping
        with constant 1).
        """
        self.mapping=mapping

    def getRegularization(self):
        """
        Returns the regularization object used in this inversion.
        The regularization is only valid once setup() has been called.
        """
        raise NotImplementedError

    def getForwardModel(self):
        """
        Returns the forward model object used in this inversion.
        The forward model is only valid once setup() has been called.
        """
        raise NotImplementedError

    def getCostFunction(self):
        """
        Returns the cost function object used in this inversion.
        The cost function is only valid once setup() has been called.
        """
        raise NotImplementedError

    def setup(self):
        """
        This method must be overwritten to perform any setup needed by the
        solver to run. The relevant objects for the inversion (e.g. forward
        model, regularization etc.) should be created and ready to use.
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

    def setWeights(self, mu_reg=None, mu_model=1.):
        self._mu_reg=mu_reg
        self._mu_model=mu_model

    def siloWriterCallback(self, k, x, fx, gfx):
        fn='inv.%d'%k
        ds=createDataset(rho=self.mapping.getValue(x))
        ds.setCycleAndTime(k,k)
        ds.saveSilo(fn)
        self.logger.debug("Jreg(m) = %e"%self.regularization.getValue(x))
        self.logger.debug("f(m) = %e"%fx)

    def getRegularization(self):
        if self.__is_setup:
            return self.regularization
        return None

    def getForwardModel(self):
        if self.__is_setup:
            return self.forwardmodel
        return None

    def getCostFunction(self):
        if self.__is_setup:
            return self.f
        return None

    def setup(self):
        if self.source is None:
            raise ValueError("No data source set!")

        self.logger.info('Retrieving domain...')
        domain=self.source.getDomain()
        DIM=domain.getDim()
        self.logger.info("Retrieving density mask...")
        rho_mask = self.source.getDensityMask()
        self.logger.info("Retrieving gravity and standard deviation data...")
        g, sigma=self.source.getGravityAndStdDev()
        chi=safeDiv(1., sigma*sigma)
        chi=interpolate(chi, Function(domain))
        g=interpolate(g, Function(domain))
        self.logger.debug("g = %s"%g)
        self.logger.debug("sigma = %s"%sigma)
        self.logger.debug("chi = %s"%chi)
        chi=chi*kronecker(DIM)[DIM-1]
        m_ref=self.mapping.getInverse(0.)
        self.regularization=Regularization(domain, m_ref=m_ref, w0=None, w=[1]*DIM, location_of_set_m=rho_mask)
        self.forwardmodel=GravityModel(domain, chi, g)
        self.f=SimpleCostFunction(self.regularization, self.mapping, self.forwardmodel)
        if self._mu_reg is None:
            x=domain.getX()
            l=0
            for i in range(DIM-1):
                l=max(l, sup(x[i])-inf(x[i]))
            G=6.6742e-11
            self._mu_reg=0.5*(l*l*G)**2
        self.logger.debug("mu_reg = %s"%self._mu_reg)
        self.logger.debug("mu_model = %s"%self._mu_model)
        self.f.setWeights(mu_reg=self._mu_reg, mu_model=self._mu_model)
        self.__is_setup=True
        #args={'rho_mask':rho_mask,'g':g[DIM-1],'chi':chi[DIM-1],'sigma':sigma}
        #try:
        #    args['rho_ref']=self.source.getReferenceDensity()
        #except:
        #    pass
        #saveSilo('init', **args)


    def run(self, rho_init=0.):
        if not self.__is_setup:
            self.setup()
        solver=self.solverclass(self.f)
        solver.setCallback(self._solver_callback)
        solver.setMaxIterations(self._solver_maxiter)
        solver.setOptions(**self._solver_opts)
        solver.setTolerance(self._solver_tol)
        if not isinstance(rho_init, Data):
            rho_init=Scalar(rho_init, ContinuousFunction(self.source.getDomain()))
        m_init=self.mapping.getInverse(rho_init)

        self.logger.info("Starting solver...")
        solver.run(m_init)
        m_star=solver.getResult()
        self.logger.info("m* = %s"%m_star)
        rho_star=self.mapping.getValue(m_star)
        self.logger.info("rho* = %s"%rho_star)
        solver.logSummary()
        return rho_star

