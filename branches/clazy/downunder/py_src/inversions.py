
##############################################################################
#
# Copyright (c) 2003-2016 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

"""Higher-level classes that allow running inversions with minimal set up"""

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2016 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['InversionDriver', 'GravityInversion','MagneticInversion', 'JointGravityMagneticInversion', 'StrongJointGravityMagneticInversion']

import logging
import numpy as np

import esys.escript as es
from esys.escript import unitsSI as U
from esys.weipa import createDataset

from .inversioncostfunctions import InversionCostFunction
from .forwardmodels import GravityModel, MagneticModel, MagneticIntensityModel, SelfDemagnetizationModel
from .mappings import *
from .minimizers import *
from .regularizations import Regularization
from .datasources import DataSource
from .coordinates import makeTransformation

class InversionDriver(object):
    """
    Base class for running an inversion
    """
    def __init__(self, solverclass=None, debug=False, self_demagnetization=False, magnetic_intensity_data=False):
        """
        creates an instance of an inversion driver.

        :param solverclass: class of the solver to be used.
        :param self_demagnetization: if True self-demagnitization is applied.
        :param magnetic_intensity_data: if True magnetic intensity is used in the cost function. 
        :type solverclass: 'AbstractMinimizer'.
        """
        # This configures basic logging to get some meaningful output
        if debug:
            level=logging.DEBUG
        else:
            level=logging.INFO
        logging.getLogger('inv').setLevel(level)
        self.logger=logging.getLogger('inv.%s'%self.__class__.__name__)
        if solverclass is None: solverclass=MinimizerLBFGS
        self.__solver=solverclass()
        self.initial_value = None
        self.self_demagnetization=self_demagnetization
        self.magnetic_intensity_data=magnetic_intensity_data
        self.m=None
        self.fixGravityPotentialAtBottom()
        self.fixMagneticPotentialAtBottom()

    def fixGravityPotentialAtBottom(self, status=False):
        """
        indicates to fix the gravity potential at the bottom to zero
        (in addition to the top)

        :param status: if True gravity potential at the bottom is set to zero
        :type status: ``bool``
        """
        self._fixGravityPotentialAtBottom=status

    def fixMagneticPotentialAtBottom(self, status=True):
        """
        indicates to fix the magnetic potential at the bottom to zero
        (in addition to the top)

        :param status: if True magnetic potential at the bottom is set to zero
        :type status: ``bool``
        """
        self._fixMagneticPotentialAtBottom=status

    def setCostFunction(self, costfunction):
        """
        sets the cost function of the inversion. This function needs to be
        called before the inversion iteration can be started.

        :param costfunction: domain of the inversion
        :type costfunction: 'InversionCostFunction'
        """
        self.getSolver().setCostFunction(costfunction)

    def getCostFunction(self):
        """
        returns the cost function of the inversion.

        :rtype: 'InversionCostFunction'
        """
        if self.isSetUp():
            return self.getSolver().getCostFunction()
        else:
           raise RuntimeError("Inversion is not set up.")

    def getSolver(self):
        """
        The solver to be used in the inversion process. See the minimizers
        module for available solvers. By default, the L-BFGS minimizer is used.

        :rtype: 'AbstractMinimizer'.
        """
        return self.__solver

    def isSetUp(self):
        """
        returns True if the inversion is set up and is ready to run.

        :rtype: `bool`
        """
        if self.getSolver().getCostFunction():
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
        self.getSolver().setCallback(callback)

    def setSolverMaxIterations(self, maxiter=None):
        """
        Sets the maximum number of solver iterations to run.
        If `maxiter` is reached the iteration is terminated and
        `MinimizerMaxIterReached` is thrown.

        :param maxiter: maximum number of iteration steps.
        :type maxiter: positive ``int``
        """
        if maxiter == None: maxiter = 200
        if maxiter > 0:
            self.getSolver().setMaxIterations(maxiter)
        else:
            raise ValueError("maxiter must be positive.")

    def setSolverTolerance(self, m_tol=None, J_tol=None):
        """
        Sets the error tolerance for the solver. An acceptable solution is
        considered to be found once the tolerance is reached.

        :param m_tol: tolerance for changes to level set function. If `None`
                      changes to the level set function are not checked for
                      convergence during iteration.
        :param J_tol: tolerance for changes to cost function. If `None` changes
                      to the cost function are not checked for convergence
                      during iteration.
        :type m_tol: `float` or `None`
        :type J_tol: `float` or `None`
        :note: if both arguments are `None` the default setting m_tol=1e-4,
               J_tol=None is used.
        """
        if m_tol is None and J_tol is None:
            m_tol=1e-4
        self.getSolver().setTolerance(m_tol=m_tol, J_tol=J_tol)

    def setInitialGuess(self, *props):
        """
        Sets the initial guess for the inversion iteration.
        By default zero is used.
        """
        self.initial_value=self.getCostFunction().createLevelSetFunction(*props)

    def run(self):
        """
        Executes the inversion.

        :return: physical parameters as result of the inversion
        :rtype: ``list`` of physical parameters or a physical parameter
        """
        if not self.isSetUp():
            raise RuntimeError("Inversion is not setup.")

        if self.initial_value is None: self.setInitialGuess()
        self.logger.info("Starting solver...")
        try:
            self.getSolver().run(self.initial_value)
            self.m=self.getSolver().getResult()
            self.p=self.getCostFunction().getProperties(self.m)
        except MinimizerException as e:
            self.m=self.getSolver().getResult()
            self.p=self.getCostFunction().getProperties(self.m)
            self.logger.info("Solver iteration failed!")
            raise e
        self.logger.info("Solver completed successfully!")
        self.getSolver().logSummary()
        return self.p

    def getLevelSetFunction(self):
        """
        returns the level set function as solution of the optimization problem

        :rtype: `Data`
        """
        return self.m

    def setup(self, *args, **k_args):
        """
        sets up the inversion. The default implementation does nothing but
        it is advised to call this method before calling `run()`.
        """
        pass

class GravityInversion(InversionDriver):
    """
     Driver class to perform an inversion of Gravity (Bouguer) anomaly data.
     The class uses the standard `Regularization` class for a single level set
     function, `DensityMapping` mapping, and the gravity forward model
     `GravityModel`.
    """
    def setup(self, domainbuilder,
                    rho0=None, drho=None, z0=None, beta=None,
                    w0=None, w1=None, rho_at_depth=None):
        """
        Sets up the inversion parameters from a `DomainBuilder`.

        :param domainbuilder: Domain builder object with gravity source(s)
        :type domainbuilder: `DomainBuilder`
        :param rho0: reference density, see `DensityMapping`. If not specified, zero is used.
        :type rho0: ``float`` or ``Scalar``
        :param drho: density scale, see `DensityMapping`. If not specified, 2750kg/m^3 is used.
        :type drho: ``float`` or ``Scalar``
        :param z0: reference depth for depth weighting, see `DensityMapping`. If not specified, zero is used.
        :type z0: ``float`` or ``Scalar``
        :param beta: exponent for  depth weighting, see `DensityMapping`. If not specified, zero is used.
        :type beta: ``float`` or ``Scalar``
        :param w0: weighting factor for level set term regularization. If not set zero is assumed.
        :type w0: ``Scalar`` or ``float``
        :param w1: weighting factor for the gradient term in the regularization. If not set zero is assumed
        :type w1: ``Vector`` or list of ``float``
        :param rho_at_depth: value for density at depth, see `DomainBuilder`.
        :type rho_at_depth: ``float`` or ``None``
        """
        self.logger.info('Retrieving domain...')
        dom=domainbuilder.getDomain()
        trafo=makeTransformation(dom, domainbuilder.getReferenceSystem())
        DIM=dom.getDim()
        rho_mask = domainbuilder.getSetDensityMask()
        #========================
        self.logger.info('Creating mapping...')
        if rho_at_depth:
            rho2 = rho_mask * rho_at_depth + (1-rho_mask) * rho0
        elif rho0:
            rho2 = (1-rho_mask) * rho0
        else:
            rho2 = 0.

        rho_mapping=DensityMapping(dom, rho0=rho2, drho=drho, z0=z0, beta=beta)
        scale_mapping=rho_mapping.getTypicalDerivative()
        #========================
        self.logger.info("Setting up regularization...")
        if w1 is None:
            w1=[1.]*DIM

        regularization=Regularization(dom, numLevelSets=1,\
                               w0=w0, w1=w1, location_of_set_m=rho_mask,
                               coordinates=trafo)
        #====================================================================
        self.logger.info("Retrieving gravity surveys...")
        surveys=domainbuilder.getGravitySurveys()
        g=[]
        w=[]
        for g_i,sigma_i in surveys:
            w_i=es.safeDiv(1., sigma_i)
            if g_i.getRank()==0:
                g_i=g_i*es.kronecker(DIM)[DIM-1]
            if w_i.getRank()==0:
                w_i=w_i*es.kronecker(DIM)[DIM-1]
            g.append(g_i)
            w.append(w_i)
            self.logger.debug("Added gravity survey:")
            self.logger.debug("g = %s"%g_i)
            self.logger.debug("sigma = %s"%sigma_i)
            self.logger.debug("w = %s"%w_i)
        #====================================================================

        self.logger.info("Setting up model...")
        forward_model=GravityModel(dom, w, g, fixPotentialAtBottom=self._fixGravityPotentialAtBottom, coordinates=trafo)
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

    def siloWriterCallback(self, k, x, Jx, g_Jx, norm_dJ=None, norm_dx=None):
        """
        callback function that can be used to track the solution

        :param k: iteration count
        :param x: current approximation
        :param Jx: value of cost function
        :param g_Jx: gradient of f at x
        """
        fn='inv.%d'%k
        ds=createDataset(density=self.getCostFunction().getProperties(x))
        ds.setCycleAndTime(k,k)
        ds.saveSilo(fn)
        self.logger.debug("J(m) = %e"%Jx)

class MagneticInversion(InversionDriver):
    """
    Driver class to perform an inversion of magnetic anomaly data. The class
    uses the standard `Regularization` class for a single level set function,
    `SusceptibilityMapping` mapping and the linear magnetic forward model
    `MagneticModel`.
    """
    def setup(self, domainbuilder,
                    k0=None, dk=None, z0=None, beta=None,
                    w0=None, w1=None,k_at_depth=None ):
        """
        Sets up the inversion from a `DomainBuilder`.
        If magnetic data are given as scalar it is assumed that values are
        collected in direction of the background magnetic field.

        :param domainbuilder: Domain builder object with gravity source(s)
        :type domainbuilder: `DomainBuilder`
        :param k0: reference susceptibility, see `SusceptibilityMapping`. If not specified, zero is used.
        :type k0: ``float`` or ``Scalar``
        :param dk: susceptibility scale, see `SusceptibilityMapping`. If not specified, 1. is used.
        :type dk: ``float`` or ``Scalar``
        :param z0: reference depth for depth weighting, see `SusceptibilityMapping`. If not specified, zero is used.
        :type z0: ``float`` or ``Scalar``
        :param beta: exponent for  depth weighting, see `SusceptibilityMapping`. If not specified, zero is used.
        :type beta: ``float`` or ``Scalar``
        :param w0: weighting factor for level set term regularization. If not set zero is assumed.
        :type w0: ``Scalar`` or ``float``
        :param w1: weighting factor for the gradient term in the regularization. If not set zero is assumed
        :type w1: ``Vector`` or list of ``float``
        :param k_at_depth: value for susceptibility at depth, see `DomainBuilder`.
        :type k_at_depth: ``float`` or ``None``
        """
        self.logger.info('Retrieving domain...')
        dom=domainbuilder.getDomain()
        DIM=dom.getDim()
        trafo=makeTransformation(dom, domainbuilder.getReferenceSystem())
        #========================
        self.logger.info('Creating mapping ...')
        k_mask = domainbuilder.getSetSusceptibilityMask()
        if k_at_depth:
             k2= k_mask *  k_at_depth + (1-k_mask) * k0
        elif k0:
             k2= (1-k_mask) * k0
        else:
             k2=0

        k_mapping=SusceptibilityMapping(dom, k0=k2, dk=dk, z0=z0, beta=beta)
        scale_mapping=k_mapping.getTypicalDerivative()
        #========================
        self.logger.info("Setting up regularization...")
        if w1 is None:
            w1=[1.]*DIM

        regularization=Regularization(dom, numLevelSets=1,w0=w0, w1=w1, location_of_set_m=k_mask, coordinates=trafo)

        #====================================================================
        self.logger.info("Retrieving magnetic field surveys...")
        d_b=es.normalize(domainbuilder.getBackgroundMagneticFluxDensity())
        surveys=domainbuilder.getMagneticSurveys()
        B=[]
        w=[]
        for B_i,sigma_i in surveys:
            w_i=es.safeDiv(1., sigma_i)
            if self.magnetic_intensity_data:
                if not B_i.getRank()==0:
                    B_i=es.length(B_i)
                if not w_i.getRank()==0:
                    w_i=length(w_i)
            else:
                if B_i.getRank()==0:
                    B_i=B_i*d_b
                if w_i.getRank()==0:
                    w_i=w_i*d_b
            B.append(B_i)
            w.append(w_i)
            self.logger.debug("Added magnetic survey:")
            self.logger.debug("B = %s"%B_i)
            self.logger.debug("sigma = %s"%sigma_i)
            self.logger.debug("w = %s"%w_i)
        #====================================================================
        self.logger.info("Setting up model...")
        if self.self_demagnetization:
            forward_model=SelfDemagnetizationModel(dom, w, B, domainbuilder.getBackgroundMagneticFluxDensity(), fixPotentialAtBottom=self._fixMagneticPotentialAtBottom, coordinates=trafo)
        else:
            if self.magnetic_intensity_data:
              forward_model=MagneticIntensityModel(dom, w, B, domainbuilder.getBackgroundMagneticFluxDensity(), fixPotentialAtBottom=self._fixMagneticPotentialAtBottom, coordinates=trafo)
            else:
              forward_model=MagneticModel(dom, w, B, domainbuilder.getBackgroundMagneticFluxDensity(), fixPotentialAtBottom=self._fixMagneticPotentialAtBottom, coordinates=trafo)
        forward_model.rescaleWeights(k_scale=scale_mapping)

        #====================================================================
        self.logger.info("Setting cost function...")
        self.setCostFunction(InversionCostFunction(regularization, k_mapping, forward_model))

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

    def siloWriterCallback(self, k, x, Jx, g_Jx, norm_dJ=None, norm_dx=None):
        """
        callback function that can be used to track the solution

        :param k: iteration count
        :param x: current approximation
        :param Jx: value of cost function
        :param g_Jx: gradient of f at x
        """
        fn='inv.%d'%k
        ds=createDataset(susceptibility=self.getCostFunction().getProperties(x))
        ds.setCycleAndTime(k,k)
        ds.saveSilo(fn)
        self.logger.debug("J(m) = %e"%Jx)

class JointGravityMagneticInversion(InversionDriver):
    """
    Driver class to perform a joint inversion of Gravity (Bouguer) and
    magnetic anomaly data. The class uses the standard `Regularization` class
    for two level set functions with cross-gradient correlation,
    `DensityMapping` and `SusceptibilityMapping` mappings, the gravity forward
    model `GravityModel` and the linear magnetic forward model `MagneticModel`.
    """
    DENSITY=0
    SUSCEPTIBILITY=1

    def setup(self, domainbuilder,
                    rho0=None, drho=None, rho_z0=None, rho_beta=None,
                    k0=None, dk=None, k_z0=None, k_beta=None, w0=None, w1=None,
                    w_gc=None,rho_at_depth=None, k_at_depth=None):
        """
        Sets up the inversion from an instance ``domainbuilder`` of a
        `DomainBuilder`. Gravity and magnetic data attached to the
        ``domainbuilder`` are considered in the inversion.
        If magnetic data are given as scalar it is assumed that values are
        collected in direction of the background magnetic field.

        :param domainbuilder: Domain builder object with gravity source(s)
        :type domainbuilder: `DomainBuilder`
        :param rho0: reference density, see `DensityMapping`. If not specified, zero is used.
        :type rho0: ``float`` or `Scalar`
        :param drho: density scale, see `DensityMapping`. If not specified, 2750kg/m^3 is used.
        :type drho: ``float`` or `Scalar`
        :param rho_z0: reference depth for depth weighting for density, see `DensityMapping`. If not specified, zero is used.
        :type rho_z0: ``float`` or `Scalar`
        :param rho_beta: exponent for  depth weighting  for density, see `DensityMapping`. If not specified, zero is used.
        :type rho_beta: ``float`` or `Scalar`
        :param k0: reference susceptibility, see `SusceptibilityMapping`. If not specified, zero is used.
        :type k0: ``float`` or `Scalar`
        :param dk: susceptibility scale, see `SusceptibilityMapping`. If not specified, 2750kg/m^3 is used.
        :type dk: ``float`` or `Scalar`
        :param k_z0: reference depth for depth weighting for susceptibility, see `SusceptibilityMapping`. If not specified, zero is used.
        :type k_z0: ``float`` or `Scalar`
        :param k_beta: exponent for  depth weighting for susceptibility, see `SusceptibilityMapping`. If not specified, zero is used.
        :type k_beta: ``float`` or `Scalar`
        :param w0: weighting factors for level set term regularization, see `Regularization`. If not set zero is assumed.
        :type w0: `es.Data` or ``ndarray`` of shape (2,)
        :param w1: weighting factor for the gradient term in the regularization see `Regularization`. If not set zero is assumed
        :type w1: `es.Data` or ``ndarray`` of shape (2,DIM)
        :param w_gc: weighting factor for the cross gradient term in the regularization, see `Regularization`. If not set one is assumed
        :type w_gc: `Scalar` or `float`
        :param k_at_depth: value for susceptibility at depth, see `DomainBuilder`.
        :type k_at_depth: ``float`` or ``None``
        :param rho_at_depth: value for density at depth, see `DomainBuilder`.
        :type rho_at_depth: ``float`` or ``None``
        """
        self.logger.info('Retrieving domain...')
        dom=domainbuilder.getDomain()
        DIM=dom.getDim()
        trafo=makeTransformation(dom, domainbuilder.getReferenceSystem())
        #========================
        self.logger.info('Creating mappings ...')
        rho_mask=domainbuilder.getSetDensityMask()
        if rho_at_depth:
             rho2= rho_mask *  rho_at_depth + (1-rho_mask) * rho0
        elif rho0:
             rho2= (1-rho_mask) * rho0
        else:
             rho2=0

        k_mask = domainbuilder.getSetSusceptibilityMask()
        if k_at_depth:
             k2= k_mask *  k_at_depth + (1-k_mask) * k0
        elif k0:
             k2= (1-k_mask) * k0
        else:
             k2=0

        rho_mapping=DensityMapping(dom, rho0=rho2, drho=drho, z0=rho_z0, beta=rho_beta)
        rho_scale_mapping=rho_mapping.getTypicalDerivative()
        self.logger.debug("rho_scale_mapping = %s"%rho_scale_mapping)
        k_mapping=SusceptibilityMapping(dom, k0=k2, dk=dk, z0=k_z0, beta=k_beta)
        k_scale_mapping=k_mapping.getTypicalDerivative()
        self.logger.debug("k_scale_mapping = %s"%k_scale_mapping)
        #========================
        self.logger.info("Setting up regularization...")

        if w1 is None:
            w1=np.ones((2,DIM))

        wc=es.Data(0.,(2,2), es.Function(dom))
        if w_gc is  None:
            wc[0,1]=1
        else:
            wc[0,1]=w_gc

        reg_mask=es.Data(0.,(2,), es.Solution(dom))
        reg_mask[self.DENSITY] = rho_mask
        reg_mask[self.SUSCEPTIBILITY] = k_mask
        regularization=Regularization(dom, numLevelSets=2,\
                               w0=w0, w1=w1,wc=wc, location_of_set_m=reg_mask, coordinates=trafo)
        #====================================================================
        self.logger.info("Retrieving gravity surveys...")
        surveys=domainbuilder.getGravitySurveys()
        g=[]
        w=[]
        for g_i,sigma_i in surveys:
            w_i=es.safeDiv(1., sigma_i)
            if g_i.getRank()==0:
                g_i=g_i*es.kronecker(DIM)[DIM-1]
            if w_i.getRank()==0:
                w_i=w_i*es.kronecker(DIM)[DIM-1]
            g.append(g_i)
            w.append(w_i)
            self.logger.debug("Added gravity survey:")
            self.logger.debug("g = %s"%g_i)
            self.logger.debug("sigma = %s"%sigma_i)
            self.logger.debug("w = %s"%w_i)

        self.logger.info("Setting up gravity model...")
        gravity_model=GravityModel(dom, w, g, fixPotentialAtBottom=self._fixGravityPotentialAtBottom, coordinates=trafo)
        gravity_model.rescaleWeights(rho_scale=rho_scale_mapping)
        #====================================================================
        self.logger.info("Retrieving magnetic field surveys...")
        d_b=es.normalize(domainbuilder.getBackgroundMagneticFluxDensity())
        surveys=domainbuilder.getMagneticSurveys()
        B=[]
        w=[]
        for B_i,sigma_i in surveys:
            w_i=es.safeDiv(1., sigma_i)
            if self.magnetic_intensity_data:
                if not B_i.getRank()==0:
                    B_i=es.length(B_i)
                if not w_i.getRank()==0:
                    w_i=length(w_i)
            else:
              if B_i.getRank()==0:
                  B_i=B_i*d_b
              if w_i.getRank()==0:
                  w_i=w_i*d_b
            B.append(B_i)
            w.append(w_i)
            self.logger.debug("Added magnetic survey:")
            self.logger.debug("B = %s"%B_i)
            self.logger.debug("sigma = %s"%sigma_i)
            self.logger.debug("w = %s"%w_i)

        self.logger.info("Setting up magnetic model...")
        if self.self_demagnetization:
            magnetic_model=SelfDemagnetizationModel(dom, w, B, domainbuilder.getBackgroundMagneticFluxDensity(), fixPotentialAtBottom=self._fixMagneticPotentialAtBottom, coordinates=trafo)
        else:
            if self.magnetic_intensity_data:
              magnetic_model=MagneticIntensityModel(dom, w, B, domainbuilder.getBackgroundMagneticFluxDensity(), fixPotentialAtBottom=self._fixMagneticPotentialAtBottom, coordinates=trafo)
            else:
              magnetic_model=MagneticModel(dom, w, B, domainbuilder.getBackgroundMagneticFluxDensity(), fixPotentialAtBottom=self._fixMagneticPotentialAtBottom, coordinates=trafo)
        magnetic_model.rescaleWeights(k_scale=k_scale_mapping)
        #====================================================================
        self.logger.info("Setting cost function...")
        self.setCostFunction(InversionCostFunction(regularization,
             ((rho_mapping,self.DENSITY), (k_mapping, self.SUSCEPTIBILITY)),
               ((gravity_model,0), (magnetic_model,1)) ))

    def setInitialGuess(self, rho=None, k=None):
        """
        set the initial guess *rho* for density and *k* for susceptibility for the inversion iteration.

        :param rho: initial value for the density anomaly.
        :type rho: `Scalar`
        :param k: initial value for the susceptibility anomaly.
        :type k: `Scalar`
        """
        super(JointGravityMagneticInversion,self).setInitialGuess(rho, k)

    def siloWriterCallback(self, k, x, Jx, g_Jx, norm_dJ=None, norm_dx=None):
        """
        callback function that can be used to track the solution

        :param k: iteration count
        :param x: current approximation
        :param Jx: value of cost function
        :param g_Jx: gradient of f at x
        """
        fn='inv.%d'%k
        p=self.getCostFunction().getProperties(x)
        ds=createDataset(density=p[self.DENSITY], susceptibility=p[self.SUSCEPTIBILITY])
        ds.setCycleAndTime(k,k)
        ds.saveSilo(fn)
        self.logger.debug("J(m) = %e"%Jx)

class StrongJointGravityMagneticInversion(InversionDriver):
    """
    Driver class to perform a joint inversion of Gravity (Bouguer) and magnetic
    anomaly data with the assumption that there is a functional relationship
    between density and susceptibility.

    The class uses the standard `Regularization` class for a single level set
    function, `DensityMapping` and `SusceptibilityMapping` mappings, the
    gravity forward model `GravityModel` and the linear magnetic forward model
    `MagneticModel`.
    """
    DENSITY=0
    SUSCEPTIBILITY=1

    def setup(self, domainbuilder,
                    rho0=None, drho=None, rho_z0=None, rho_beta=None,
                    k0=None, dk=None, k_z0=None, k_beta=None, w0=None, w1=None,
                    w_gc=None,rho_at_depth=None, k_at_depth=None):
        """
        Sets up the inversion from an instance ``domainbuilder`` of a
        `DomainBuilder`. Gravity and magnetic data attached to the
        ``domainbuilder`` are considered in the inversion.
        If magnetic data are given as scalar it is assumed that values are
        collected in direction of the background magnetic field.

        :param domainbuilder: Domain builder object with gravity source(s)
        :type domainbuilder: `DomainBuilder`
        :param rho0: reference density, see `DensityMapping`. If not specified,
                     zero is used.
        :type rho0: ``float`` or `Scalar`
        :param drho: density scale, see `DensityMapping`. If not specified,
                     2750 kg/m^3 is used.
        :type drho: ``float`` or `Scalar`
        :param rho_z0: reference depth for depth weighting for density, see
                       `DensityMapping`. If not specified, zero is used.
        :type rho_z0: ``float`` or `Scalar`
        :param rho_beta: exponent for density depth weighting, see
                         `DensityMapping`. If not specified, zero is used.
        :type rho_beta: ``float`` or `Scalar`
        :param k0: reference susceptibility, see `SusceptibilityMapping`.
                   If not specified, zero is used.
        :type k0: ``float`` or `Scalar`
        :param dk: susceptibility scale, see `SusceptibilityMapping`. If not
                   specified, 1. is used.
        :type dk: ``float`` or `Scalar`
        :param k_z0: reference depth for susceptibility depth weighting, see
                     `SusceptibilityMapping`. If not specified, zero is used.
        :type k_z0: ``float`` or `Scalar`
        :param k_beta: exponent for susceptibility depth weighting, see
                       `SusceptibilityMapping`. If not specified, zero is used.
        :type k_beta: ``float`` or `Scalar`
        :param w0: weighting factor for level set term in the regularization.
                   If not set zero is assumed.
        :type w0: ``Scalar`` or ``float``
        :param w1: weighting factor for the gradient term in the regularization
                   see `Regularization`.  If not set zero is assumed.
        :type w1: `es.Data` or ``ndarray`` of shape (DIM,)
        :param w_gc: weighting factor for the cross gradient term in the
                     regularization, see `Regularization`. If not set one is
                     assumed.
        :type w_gc: `Scalar` or `float`
        :param k_at_depth: value for susceptibility at depth, see `DomainBuilder`.
        :type k_at_depth: ``float`` or ``None``
        :param rho_at_depth: value for density at depth, see `DomainBuilder`.
        :type rho_at_depth: ``float`` or ``None``
        """
        self.logger.info('Retrieving domain...')
        dom=domainbuilder.getDomain()
        DIM=dom.getDim()
        trafo=makeTransformation(dom, domainbuilder.getReferenceSystem())

        rock_mask=wherePositive(domainbuilder.getSetDensityMask() + domainbuilder.getSetSusceptibilityMask())
        #========================
        self.logger.info('Creating mappings ...')
        if rho_at_depth:
             rho2= rock_mask *  rho_at_depth + (1-rock_mask) * rho0
        elif rho0:
             rho2= (1-rock_mask) * rho0
        else:
             rho2=0

        if k_at_depth:
             k2= rock_mask *  k_at_depth + (1-rock_mask) * k0
        elif k0:
             k2= (1-rock_mask) * k0
        else:
             k2=0

        rho_mapping=DensityMapping(dom, rho0=rho2, drho=drho, z0=rho_z0, beta=rho_beta)
        rho_scale_mapping=rho_mapping.getTypicalDerivative()
        self.logger.debug("rho_scale_mapping = %s"%rho_scale_mapping)
        k_mapping=SusceptibilityMapping(dom, k0=k2, dk=dk, z0=k_z0, beta=k_beta)
        k_scale_mapping=k_mapping.getTypicalDerivative()
        self.logger.debug("k_scale_mapping = %s"%k_scale_mapping)
        #========================
        self.logger.info("Setting up regularization...")
        if w1 is None:
            w1=[1.]*DIM

        regularization=Regularization(dom, numLevelSets=1,w0=w0, w1=w1, location_of_set_m=rock_mask, coordinates=trafo)
        #====================================================================
        self.logger.info("Retrieving gravity surveys...")
        surveys=domainbuilder.getGravitySurveys()
        g=[]
        w=[]
        for g_i,sigma_i in surveys:
            w_i=es.safeDiv(1., sigma_i)
            if g_i.getRank()==0:
                g_i=g_i*es.kronecker(DIM)[DIM-1]
            if w_i.getRank()==0:
                w_i=w_i*es.kronecker(DIM)[DIM-1]
            g.append(g_i)
            w.append(w_i)
            self.logger.debug("Added gravity survey:")
            self.logger.debug("g = %s"%g_i)
            self.logger.debug("sigma = %s"%sigma_i)
            self.logger.debug("w = %s"%w_i)

        self.logger.info("Setting up gravity model...")
        gravity_model=GravityModel(dom, w, g, fixPotentialAtBottom=self._fixGravityPotentialAtBottom, coordinates=trafo)
        gravity_model.rescaleWeights(rho_scale=rho_scale_mapping)
        #====================================================================
        self.logger.info("Retrieving magnetic field surveys...")
        d_b=es.normalize(domainbuilder.getBackgroundMagneticFluxDensity())
        surveys=domainbuilder.getMagneticSurveys()
        B=[]
        w=[]
        for B_i,sigma_i in surveys:
            w_i=es.safeDiv(1., sigma_i)
            if self.magnetic_intensity_data:
              if not B_i.getRank()==0:
                B_i=es.length(B_i)
              if not w_i.getRank()==0:
                w_i=length(w_i)
            else:
              if B_i.getRank()==0:
                B_i=B_i*d_b
              if w_i.getRank()==0:
                w_i=w_i*d_b
            B.append(B_i)
            w.append(w_i)
            self.logger.debug("Added magnetic survey:")
            self.logger.debug("B = %s"%B_i)
            self.logger.debug("sigma = %s"%sigma_i)
            self.logger.debug("w = %s"%w_i)

        self.logger.info("Setting up magnetic model...")
        if self.self_demagnetization:
            magnetic_model=SelfDemagnetizationModel(dom, w, B, domainbuilder.getBackgroundMagneticFluxDensity(), fixPotentialAtBottom=self._fixMagneticPotentialAtBottom, coordinates=trafo)
        else:
            if self.magnetic_intensity_data:
              magnetic_model=MagneticIntensityModel(dom, w, B, domainbuilder.getBackgroundMagneticFluxDensity(), fixPotentialAtBottom=self._fixMagneticPotentialAtBottom, coordinates=trafo)
            else:
              magnetic_model=MagneticModel(dom, w, B, domainbuilder.getBackgroundMagneticFluxDensity(), fixPotentialAtBottom=self._fixMagneticPotentialAtBottom, coordinates=trafo)
        magnetic_model.rescaleWeights(k_scale=k_scale_mapping)
        #====================================================================
        self.logger.info("Setting cost function...")

        self.setCostFunction(InversionCostFunction(regularization,
               (rho_mapping, k_mapping),
               ((gravity_model,0), (magnetic_model,1)) ))

    def setInitialGuess(self, rho=None, k=None):
        """
        set the initial guess *rho* for density and *k* for susceptibility for
        the inversion iteration.

        :param rho: initial value for the density anomaly.
        :type rho: `Scalar`
        :param k: initial value for the susceptibility anomaly.
        :type k: `Scalar`
        """
        super(StrongJointGravityMagneticInversion,self).setInitialGuess(rho, k)

    def siloWriterCallback(self, k, x, Jx, g_Jx, norm_dJ=None, norm_dx=None):
        """
        callback function that can be used to track the solution

        :param k: iteration count
        :param x: current m approximation
        :param Jx: value of cost function
        :param g_Jx: gradient of f at x
        """
        fn='inv.%d'%k
        p=self.getCostFunction().getProperties(x)
        ds=createDataset(density=p[self.DENSITY], susceptibility=p[self.SUSCEPTIBILITY])
        ds.setCycleAndTime(k,k)
        ds.saveSilo(fn)
        self.logger.debug("J(m) = %e"%Jx)

