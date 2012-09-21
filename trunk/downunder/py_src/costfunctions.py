
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

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

try:
    # only needed for getDirectionalDerivative(), so ignore error
    from esys.escript import grad
except:
    pass

class CostFunction(object):
    """
    A function f(x) that can be minimized (base class).

    Example of usage:
        cf=CostFunction()
        ... calculate x ...
        args=cf.getArguments(x) # this could be potentially expensive!
        f=cf.getValue(x, *args)
        ... it could be required to update x without using the gradient...
        ... but then ...
        gf=cf.getGradient(x, *args)

    The function calls update statistical information.
    The actual work is done by the methods with corresponding name and a
    leading underscore. These functions need to be overwritten for a particular
    cost function.
    """

    def __init__(self):
        self.resetCounters()

    def resetCounters(self):
        """
        resets all statistical counters
        """
        self.Inner_calls=0
        self.Value_calls=0
        self.Gradient_calls=0
        self.DirectionalDerivative_calls=0
        self.Arguments_calls=0

    def getInner(self, f0, f1):
        """
        returns the inner product of two gradients.
        """
        self.Inner_calls+=1
        return self._getInner(f0, f1)

    def getValue(self, x, *args):
        """
        returns the value f(x) using the precalculated values for x.
        """
        self.Value_calls+=1
        return self._getValue(x, *args)

    def __call__(self, x, *args):
        """
        short for getValue(x, *args).
        """
        return self.getValue(x, *args)

    def getGradient(self, x, *args):
        """
        returns the gradient of f at x using the precalculated values for x.
        """
        self.Gradient_calls+=1
        return self._getGradient(x, *args)

    def getDirectionalDerivative(self, x, d, *args):
        """
        returns inner(grad f(x), d) using the precalculated values for x.
        """
        self.DirectionalDerivative_calls+=1
        return self._getDirectionalDerivative(x, d, *args)

    def getArguments(self, x):
        """
        returns precalculated values that are shared in the calculation of
        f(x) and grad f(x).
        """
        self.Arguments_calls+=1
        return self._getArguments(x)

    def _getInner(self, f0, f1):
        """
        Empty implementation of getInner().

        Needs to be overwritten.
        """
        raise NotImplementedError

    def _getValue(self, x, *args):
        """
        Empty implementation of getValue().

        Needs to be overwritten.
        """
        raise NotImplementedError

    def _getGradient(self, x, *args):
        """
        Empty implementation of getGradient().

        Needs to be overwritten.
        """
        raise NotImplementedError

    def _getDirectionalDerivative(self, x, d, *args):
        """
        returns getInner(grad f(x), d) using the precalculated values for x.

        This function may be overwritten as there might be more efficient ways
        of calculating the return value rather than using a self.getGradient()
        call.
        """
        return self.getInner(self.getGradient(x, *args), d)

    def _getArguments(self, x):
        """
        returns precalculated values that are shared in the calculation of f(x)
        and grad f(x).
        """
        return ()


class SimpleCostFunction(CostFunction):
    """
    This is a simple cost function with a single continuous (mapped) variable.
    It is the sum of two weighted terms, a single forward model and a single
    regularization term. This cost function is used in the gravity inversion.
    """
    def __init__(self, regularization, mapping, forwardmodel):
        super(SimpleCostFunction, self).__init__()
        self.forwardmodel=forwardmodel
        self.regularization=regularization
        self.mapping=mapping
        self.setWeights()

    def setWeights(self, mu_model=1., mu_reg=1.):
        """
        sets the weighting factors for the forward model and regularization
        terms.
        """
        if mu_model<0. or mu_reg<0.:
            raise ValueError("weighting factors must be non-negative.")
        self.mu_model=mu_model
        self.mu_reg=mu_reg
 
    def _getInner(self, f0, f1):
        """
        returns regularization.getInner(f0,f1)
        """
        # if there is more than one regularization involved their contributions
        # need to be added up.
        return self.regularization.getInner(f0, f1)

    def _getArguments(self, m):
        """
        returns precalculated values that are shared in the calculation of f(x)
        and grad f(x).
        """
        rho=self.mapping(m)
        return rho, self.forwardmodel.getArguments(rho), self.regularization.getArguments(m)

    def _getValue(self, m, *args):
        """
        returns the function value at m.
        If the precalculated values are not supplied getArguments() is called.

        """
        # if there is more than one forward_model and/or regularization their
        # contributions need to be added up.
        if len(args)==0:
            args=self.getArguments(m)
        return self.mu_model * self.forwardmodel.getValue(args[0],*args[1]) \
               + self.mu_reg * self.regularization.getValue(m)

    def _getGradient(self, m, *args):
        """
        returns the gradient of f at m.
        If the precalculated values are not supplied getArguments() is called.
        """
        drhodm = self.mapping.getDerivative(m)
        if len(args)==0:
            args = self.getArguments(m)
        Y0 = self.forwardmodel.getGradient(args[0],*args[1])
        Y1, X1 = self.regularization.getGradient(m)
        return self.regularization.project(Y=self.mu_reg*Y1 + self.mu_model*Y0*drhodm, X=self.mu_reg*X1)

    def _getDirectionalDerivative(self, m, d, *args):
        """
        returns the directional derivative at m in direction d.
        """
        drhodm = self.mapping.getDerivative(m)
        Y0 = self.forwardmodel.getGradient(args[0],*args[1])
        Y1, X1 = self.regularization.getGradient(m)
        return self.regularization.getInner(d, self.mu_reg*Y1 + self.mu_model*Y0*drhodm) \
                + self.mu_reg*self.regularization.getInner(grad(d), X1)

