
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

"""Collection of cost functions for inversion"""

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['SimpleInversionCostFunction']

from costfunctions import MeteredCostFunction
from esys.escript.pdetools import ArithmeticTuple
from esys.escript import Data


class SimpleInversionCostFunction(MeteredCostFunction):
    """
    This is a simple cost function with a single continuous (mapped) variable.
    It is the sum of two weighted terms, a single forward model and a single
    regularization term. This cost function is used in the gravity inversion.
    """
    provides_inverse_Hessian_approximation=True

    def __init__(self, regularization, mapping, forwardmodel):
        """
        constructor stores the supplied object references and sets default
        weights.

        :param regularization: The regularization part of the cost function
        :param mapping: Parametrization object
        :param forwardmodel: The forward model part of the cost function
        """
        super(SimpleInversionCostFunction, self).__init__()
        self.forwardmodel=forwardmodel
        self.regularization=regularization
        self.mapping=mapping
        self.setWeights()

    def setWeights(self, mu_model=1., mu_reg_0=None,mu_reg_1=None):
        """
        sets the weighting factors for the forward model and regularization
        terms.
        
        :param mu_model: Weighting factor for the forward model (default=1.)
        :type mu_model: non-negative `float`
        :param mu_reg_0: Weighting factor for the regularization (default=1.)
        :type mu_reg_0: non-negative `float`
        """
        if mu_model<0:
            raise ValueError("weighting factors must be non-negative.")
        self.mu_model=mu_model
        self.regularization.setWeightsForS0(mu_reg_0)
        self.regularization.setWeightsForS1(mu_reg_1)
        
    def _getDualProduct(self, x, r):
        """
        returns ``regularization.getDualProduct(x, r)``

        :rtype: `float`
        """
        return self.regularization.getDualProduct(x, r)

    def _getArguments(self, m):
        """
        returns precalculated values that are shared in the calculation of
        *f(x)* and *grad f(x)*. In this implementation returns a tuple with the
        mapped value of ``m``, the arguments from the forward model and the
        arguments from the regularization.

        :rtype: `tuple`
        """
        p=self.mapping(m)
        return p, self.forwardmodel.getArguments(p), self.regularization.getArguments(m)

    def _getValue(self, m, *args):
        """
        returns the function value at m.
        If the precalculated values are not supplied `getArguments()` is called.

        :rtype: `float`
        """
        # if there is more than one forward_model and/or regularization their
        # contributions need to be added up. But this implementation allows
        # only one of each...
        if len(args)==0:
            args=self.getArguments(m)

        return  self.mu_model * self.forwardmodel.getValue(args[0],*args[1]) \
               +  self.regularization.getValue(m, *args[2])

    def _getGradient(self, m, *args):
        """
        returns the gradient of *f* at *m*.
        If the precalculated values are not supplied `getArguments()` is called.

        :rtype: `esys.escript.Data`
        """
        dpdm = self.mapping.getDerivative(m)
        if len(args)==0:
            args = self.getArguments(m)

        Y = self.forwardmodel.getGradient(args[0],*args[1]) * dpdm
        g_reg = self.regularization.getGradient(m, *args[2])
        print "grad forward = ", Y
        print "grad regularization Y  = ", g_reg[0]
        print "grad regularization X = ", g_reg[1]
        
         
        return self.mu_model * ArithmeticTuple(Y, Data()) + g_reg


    def _getInverseHessianApproximation(self, m, r, *args):
        """
        returns an approximative evaluation *p* of the inverse of the Hessian operator of the cost function
        for a given gradient type *r* at a given location *m*: *H(m) p = r*
        
        :param m: level set approximation where to calculate Hessian inverse
        :type m: ``Data``
        :param r: a given gradient 
        :type r: ``ArithmeticTuple``
        :param args: pre-calculated values for ``m`` from ``getArguments()``
        :rtype: ``Data``
        :note: in the current implementation only the regularization term is considered in the 
          inverse Hessian approximation. 
          
        """
        print "inverse Hessian approximation:"
        print "Y  = ",r[0]
        print "X  = ",r[1]
        m=self.regularization.getInverseHessianApproximation(m, r, *args[2])
        print "m  = ",m
        return m
        
    def updateHessian(self):
        """
        notifies the class that the Hessian operator needs to be updated.
        """
        self.regularization.updateHessian()

