##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

"""General cost functions for minimization"""


__copyright__ = """Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__ = "https://github.com/LutzGross/esys-escript.github.io"

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
        Returns the call statistics as a string.

        :return: formatted string with counts of all function calls
        :rtype: ``str``
        """
        out="Number of cost function evaluations: %d\n" % self.Value_calls
        out+="Number of gradient evaluations: %d\n" % self.Gradient_calls
        out+="Number of inverse Hessian evaluations: %d\n" % self.InverseHessianApproximation_calls
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
        self.Arguments_calls += 1
        return self.getArguments(m)

    def getValueAndCount(self, m, *args):
        """
        returns the value *F(m)* using the precalculated values for *m*.

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
        :param initializeHessian: indicates if the Hessian operator should be initialized using `m`.
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

    def getSqueezeFactor(self, m, p, *args):
        """
        The new solution is calculated as m+a*p with a>0. This function allows to provide an upper bound
        for a to make sure that m+a*p is valid typically to avoid overflow when the cost function is evaluated.
        the solver will take action to make sure that the value of a is not too small.

        :param m: a solution approximation
        :type m: m-type
        :param p: an increment to the solution
        :type m: m-type
        :param args: pre-calculated values for ``m`` from `getArgumentsAndCount()`
        :rtype: positive ``float`` or None

        :note: Overwrite this method to implement a cost function.
        """
        return None

    def getInverseHessianApproximation(self, r, m, *args, initializeHessian = True):
        """
        returns an approximate evaluation *p* of the inverse of the Hessian
        operator of the cost function for a given gradient *r*: *H p = r*
        :note: by default this method is returning *r*. In this case it is assumed that m-type==g-type
        :note: Overwrite this method to implement a cost function.

        :param r: a given gradient
        :type r: g-type
        :param initializeHessian: indicates if the Hessian operator should be initialized using `m`.
                                    If this method provides an approximation
                                    only and building the new Hessian approximation is expensive
                                    it is typically more efficient to update the Hessian operator
                                    occasionally only when on input initializeHessian is `True`.
                                    If the Hessian should be updated in each step ignore the
                                    value of `initializeHessian`.
        :type initializeHessian: bool
        :returns: new search direction p.
        :rtype: m-type

        """
        return r