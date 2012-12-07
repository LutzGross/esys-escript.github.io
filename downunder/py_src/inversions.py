
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

__all__ = ['InversionBase','SingleParameterInversionBase','GravityInversion','MagneticInversion']

import logging
from esys.escript import *
import esys.escript.unitsSI as U
from esys.weipa import createDataset

from .inversioncostfunctions import SimpleInversionCostFunction
from .forwardmodels import GravityModel, MagneticModel
from .mappings import *
from .minimizers import *
from .regularizations import Regularization
from .datasources import DataSource

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
        # set default solver
        self.setSolverClass()
        self.__domain=None
        self.__regularization=None
        self.__forwardmodel=None
        self.__mapping=None

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

        :rtype: `Domain`
        """
        return  self.__domain

    def setMapping(self, mapping):
        """
        Sets the mapping object to map between model parameters and the data.

        :param mapping: Parameter mapping object
        :type mapping: `Mapping`
        """
        self.__mapping=mapping

    def getMapping(self):
        """
        return the mapping(s) used in the inversion

        :rtype: `Mapping`
        """
        return self.__mapping


    def setRegularization(self, regularization):
        """
        Sets the regularization for the inversion.

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
        Sets the forward model(s) for the inversion.

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
        module for available solvers.
        By default, the L-BFGS minimizer is used.
        """
        if solverclass == None:
            self.solverclass=MinimizerLBFGS
        else:
            self.solverclass=solverclass

    def setSolverCallback(self, callback):
        """
        Sets the callback function which is called after every solver
        iteration.
        """
        self._solver_callback=callback

    def setSolverMaxIterations(self, maxiter):
        """
        Sets the maximum number of solver iterations to run.
        """
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
        returns True if the inversion is set up and ready to run.
        """
        raise NotImplementedError

    def run(self, *args):
        """
        This method starts the inversion process and must be overwritten.
        Preferably, users should be able to set a starting point so
        inversions can be continued after stopping.
        """
        raise NotImplementedError


class SingleParameterInversionBase(InversionBase):
    """
    Base class for inversions with a single parameter to be found.
    """
    def __init__(self):
        super(SingleParameterInversionBase,self).__init__()
        self.setWeights()

    def setWeights(self, mu_reg=None, mu_model=1.):
        """
        Sets the weighting factors for the cost function.
        """
        self.logger.debug("Setting weighting factors...")
        self.logger.debug("mu_reg = %s"%mu_reg)
        self.logger.debug("mu_model = %s"%mu_model)
        self._mu_reg=mu_reg
        self._mu_model=mu_model

    def siloWriterCallback(self, k, x, fx, gfx):
        """
        callback function that can be used to track the solution

        :param k: iteration count
        :param x: current m approximation
        :param fx: value of cost function
        :param gfx: gradient of f at x
        """
        fn='inv.%d'%k
        ds=createDataset(rho=self.getMapping().getValue(x))
        ds.setCycleAndTime(k,k)
        ds.saveSilo(fn)
        self.logger.debug("f(m) = %e"%fx)


    def isSetup(self):
        if self.getRegularization() and self.getMapping() \
                and self.getForwardModel() and self.getDomain():
            return True
        else:
            return False

    def run(self, initial_value=0.):
        if not self.isSetup():
            raise RuntimeError("Inversion is not setup properly.")

        f=SimpleInversionCostFunction(self.getRegularization(), self.getMapping(), self.getForwardModel())
        f.setWeights(mu_reg_1=self._mu_reg, mu_model=self._mu_model)

        solver=self.solverclass(f)
        solver.setCallback(self._solver_callback)
        solver.setMaxIterations(self._solver_maxiter)
        solver.setOptions(**self._solver_opts)
        solver.setTolerance(self._solver_tol)
        if not isinstance(initial_value, Data):
            initial_value=Scalar(initial_value, ContinuousFunction(self.getDomain()))
        m_init=self.getMapping().getInverse(initial_value)

        self.logger.info("Starting solver...")
        solver.run(m_init)
        m_star=solver.getResult()
        self.logger.info("m* = %s"%m_star)
        value_star=self.getMapping().getValue(m_star)
        self.logger.info("result* = %s"%value_star)
        solver.logSummary()
        return value_star


class GravityInversion(SingleParameterInversionBase):
    """
    Inversion of Gravity (Bouguer) anomaly data.
    """
    def setup(self, domainbuilder,
                    rho0=None, drho=None, z0=None, beta=None,
                    s0=None, s1=None):
        """
        Sets up the inversion parameters from a `DomainBuilder`.

        :param domainbuilder: Domain builder object with gravity source(s)
        :type domainbuilder: `DomainBuilder`
        :param rho0: reference density. If not specified, zero is used.
        :type rho0: ``float`` or `Scalar`
        """
        self.logger.info('Retrieving domain...')
        self.setDomain(domainbuilder.getDomain())
        DIM=self.getDomain().getDim()
        #========================
        self.logger.info('Creating mapping...')
        self.setMapping(DensityMapping(self.getDomain(), rho0=rho0, drho=drho, z0=z0, beta=beta))
        #========================
        self.logger.info("Setting up regularization...")
        if s1 is None:
            s1=[1.]*DIM
        rho_mask = domainbuilder.getSetDensityMask()
        self.setRegularization(Regularization(self.getDomain(), numLevelSets=1,\
                               s0=s0, s1=s1, location_of_set_m=rho_mask))
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
        self.setForwardModel(GravityModel(self.getDomain(), w, g))

        # this is switched off for now:
        if self._mu_reg is None and False:
            x=self.getDomain().getX()
            l=0.
            for i in range(DIM-1):
                l=max(l, sup(x[i])-inf(x[i]))
            G=U.Gravitational_Constant
            mu_reg=0.5*(l*l*G)**2
            self.setWeights(mu_reg=mu_reg)


class MagneticInversion(SingleParameterInversionBase):
    """
    Inversion of magnetic data.
    """
    def setup(self, domainbuilder,
                    k0=None, dk=None, z0=None, beta=None,
                    s0=None, s1=None):
        """
        Sets up the inversion parameters from a `DomainBuilder`.

        :param domainbuilder: Domain builder object with magnetic source(s)
        :type domainbuilder: `DomainBuilder`

        """
        self.logger.info('Retrieving domain...')
        self.setDomain(domainbuilder.getDomain())
        DIM=self.getDomain().getDim()

        #========================
        self.logger.info('Creating mapping ...')
        self.setMapping(SusceptibilityMapping(self.getDomain(), k0=k0, dk=dk, z0=z0, beta=beta))

        #========================
        self.logger.info("Setting up regularization...")
        if s1 is None:
            s1=[1.]*DIM
        k_mask = domainbuilder.getSetSusceptibilityMask()
        self.setRegularization(Regularization(self.getDomain(), numLevelSets=1,\
                               s0=s0, s1=s1, location_of_set_m=k_mask))

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
        self.setForwardModel(MagneticModel(self.getDomain(), w, B, domainbuilder.getBackgroundMagneticField()))

        # this is switched off for now:
        if self._mu_reg is None and False:
            x=self.getDomain().getX()
            l=0.
            for i in range(DIM-1):
                l=max(l, sup(x[i])-inf(x[i]))
            mu_reg=0.5*l**2
            self.setWeights(mu_reg=mu_reg)

