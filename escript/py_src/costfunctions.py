##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

"""General cost functions for minimization"""

from __future__ import print_function, division

__copyright__ = """Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__ = "https://launchpad.net/escript-finley"

__all__ = ['CostFunction']

class CostFunction(object):
    """
    A cost function *F(m)* that can be minimized (base class).
    Near a solution *m* in the search space the cost function is approximated 
    in direction *p* as       

    *F(m+p) ~ F(m) + <grad F(m), p> + < H p, p>*
    
    where *F(m)* is the value of the function at *m*, *grad F* is the gradient, *H* is the Hessian and
    *<.,.>* is dual product. These four procedures are defined in this class.

    The class distinguishes between the representation of the solution
    m (m-type) and the cost function gradient (g-type).

    Example of usage:
        #
        class DerivedCostFunction(CostFunction)
            # overwrite: getArguments, etc.
        cf=DerivedCostFunction()
        # ... calculate m ...
        args=cf.getArguments(m) # this could be potentially expensive!
        f=cf.getValue(m, *args)
        # ... it could be required to update m without using the gradient...
        # ... but then ...
        gf=cf.getGradient(m, *args)

    Use the "AndCount" versions (e.g. getValueAndCount) if you want to count calls.

    """
    def __init__(self):
        """
        Constructor. Initializes logger.
        """
        #self.logger = logging.getLogger('esys.%s' % self.__class__.__name__)
        self.resetStatistics()

    def resetStatistics(self):
        """
        resets all counters
        """
        self.DualProduct_calls = 0
        self.Value_calls = 0
        self.Gradient_calls = 0
        self.Arguments_calls = 0
        self.InverseHessianApproximation_calls = 0
        self.Norm_calls = 0
        
    def getStatistics(self):
        """
        return the call statistics as a string:
        """
        out="Number of cost function evaluations: %d\n" % self.Value_calls
        out+="Number of gradient evaluations: %d\n" % self.Gradient_calls
        out+="Number of inverse Hessian evaluations: %d\n" % self.InverseHessianApproximation_calls
        out+="Number of gradient evaluations: %d\n" % self.Gradient_calls
        out+="Number of inner product evaluations: %d\n" % self.DualProduct_calls
        out+="Number of argument evaluations: %d\n" % self.Arguments_calls
        out+="Number of norm evaluations: %d" % self.Norm_calls
        return out

    def __call__(self, m, *args):
        """
        short for `sum(getValueAndCount(m, *args))`.
        """
        return self.getValueAndCount(m, *args)
    
    def getDualProductAndCount(self, m, r):
        """
        returns the dual product *<.,.>* of ``r`` and ``m``.

        When calling this method the calling statistics is updated.

        :type m: m-type
        :type r: g-type
        :rtype: ``float``
        """
        self.DualProduct_calls += 1
        return self.getDualProduct(m, r)

    def getNormAndCount(self, m):
        """
        returns the norm of ``m``.

        When calling this method the calling statistics is updated.

        :type m: m-type
        :rtype: ``float``
        """
        self.Norm_calls += 1
        return self.getNorm(m)

    def getArgumentsAndCount(self, m):
        """
        returns precalculated values that are shared in the calculation of
        *F(m)* and *grad F(m)* and the Hessian operator. The default
        implementation returns an empty tuple.

        When calling this method the calling statistics is updated.

        :param m: location of derivative
        :type m: m-type
        :rtype: ``tuple``

        :note: The tuple returned by this call will be passed back to this `CostFunction` in other
           calls (eg: `getGradient`). Its contents are not specified at this level because no code,
           other than the `CostFunction` which created it, will be interacting with it.
           That is, the implementor can put whatever information they find useful in it.

        """
        return ()

    def getValueAndCount(self, m, *args):
        """
        returns the value *F(m)* using the precalculated values for *m*.

        If the cost function is a composition the values of these components can be
        returned. The `sum` function is automatically applied by any solver if a
        single value is required.

        When calling this method the calling statistics is updated.

        :param m: a solution approximation
        :type m: m-type
        :param args: pre-calculated values for ``m`` from `getArgumentsAndCount()`
        :rtype: ``float`` or valid argument for ``sum``
        """
        self.Value_calls += 1
        if not args:
            args = self.getArgumentsAndCount(m)
        return self.getValue(m, *args)

    def getGradientAndCount(self, m, *args):
        """
        returns the gradient of *F* at *m* using the precalculated values for *m*.

        When calling this method the calling statistics is updated.

        :param m: location of derivative
        :type m: m-type
        :param args: pre-calculated values for ``m`` from `getArgumentsAndCount()`
        :rtype: g-type
        """
        self.Gradient_calls += 1
        if not args:
            args = self.getArgumentsAndCount(m)
        return self.getGradient(m, *args)

    def getInverseHessianApproximationAndCount(self, r, m,  *args, initializeHessian=True):
        """
        returns an evaluation *p* of the inverse of the Hessian
        operator of the cost function for a given gradient *r*: *H p = r*

        When calling this method the calling statistics is updated.

        :param r: a given gradient
        :type r: g-type
        :poram initializeHessian: indicates if the Hessian operator should be initialized using `m`.
                                    If updating the Hessian is expensive it should only be done
                                    when initializeHessian is True. If this method provides an approximation
                                    only and building the new Hessian approximation is expensive
                                    it is typically more efficient to update the Hessian operator
                                    occasionally.
        :type initializeHessian: bool
        :param args: pre-calculated values for ``m`` from `getArgumentsAndCount()`
        :returns: new search direction p.
        :rtype: m-type
        """
        if not args:
            args = self.getArgumentsAndCount(m)
        self.InverseHessianApproximation_calls += 1
        return self.getInverseHessianApproximation(r, m,  *args, initializeHessian=initializeHessian)

    def getDualProduct(self, m, r):
        """
        returns the dual product of ``m`` and ``r``

        :type m: m-type
        :type r: g-type
        :rtype: ``float``
        :note: Overwrite this method to implement a cost function.
        """
        raise NotImplementedError

    def getNorm(self, m):
        """
        returns the norm of ``m``.
        Typically, this is the 'Lsup' function.

        :type m: m-type
        :rtype: ``float``
        :note: Overwrite this method to implement a cost function.
        """
        raise NotImplementedError

    def getArguments(self, m):
        """
        returns precalculated values that are shared in the calculation of
        *F(m)* and *grad F(m)* and the Hessian operator.

        :note: Overwrite this method to implement a cost function.
        :note: The tuple returned by this call will be passed back to this `CostFunction` in other
           calls(eg: ``getGradientAndCount``). Its contents are not specified at this level because no code,
           other than the `CostFunction` which created it, will be interacting with it.
           That is, the implementor can put whatever information they find useful in it.

        :param m: location of derivative
        :type m: m-type
        :returns: pre-calculated arguments
        :rtype: ``tuple``

        :note: Overwrite this method to implement a cost function.
        """
        return ()

    def getValue(self, m, *args):
        """
        returns the value *F(m)* using the precalculated values for *m*.

        If the cost function is a composition the values of these components can be
        returned. The `sum` function is automatically applied by any solver if a
        single value is required.

        :param m: a solution approximation
        :type m: m-type
        :param args: pre-calculated values for ``m`` from `getArgumentsAndCount()`
        :rtype: ``float`` or valid argument for ``sum``

        :note: Overwrite this method to implement a cost function.
        """
        raise NotImplementedError

    def getGradient(self, m, *args):
        """
        returns the gradient of *F* at *m* using the precalculated values for
        *m*.

        :param m: location of derivative
        :type m: m-type
        :param args: pre-calculated values for ``m`` from `getArgumentsAndCount()`
        :rtype: g-type
        :note: Overwrite this method to implement a cost function.
        """
        raise NotImplementedError

    def updateHessian(self, m, *args):
        """
        this function is called to update the (potentially approximate) Hessian H
        at a given location *m*.

        :param m: location of Hessian operator to be set
        :type m: m-type
        :param args: pre-calculated values for ``m`` from `getArgumentsAndCount()`
        :returns: None
        :note: Overwrite this method to implement a cost function.

        :note: Typically this method is called by the solver when the iteration is restarted.
        """
        pass

    def getInverseHessianApproximation(self, r, m, *args, initializeHessian = True):
        """
        returns an approximate evaluation *p* of the inverse of the Hessian
        operator of the cost function for a given gradient *r*: *H p = r*
        The Hessian is set and updated by the solver calling ``updateHessian``
        :note: by default this method is returning *r*. In this case it is assumed that m-type==g-type
        :note: Overwrite this method to implement a cost function.

        :param r: a given gradient
        :type r: g-type
        :poram initializeHessian: indicates if the Hessian operator should be initialized using `m`.
                                    If this method provides an approximation
                                    only and building the new Hessian approximation is expensive
                                    it is typically more efficient to update the Hessian operator
                                    occasionally only when on input initializeHessian is `True`.
                                    If the Hessian should be updated in each step ignore the initializeHessian
                                    value otherwise update the Hessian only if `initializeHessian == True`.
        :type initializeHessian: bool
        :returns: new search direction p.
        :rtype: m-type

        """
        return r