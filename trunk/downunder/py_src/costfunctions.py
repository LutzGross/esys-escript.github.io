
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

"""General cost functions for minimization"""

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = [ 'CostFunction', 'MeteredCostFunction' ]


class CostFunction(object):
    """
    A function *f(x)* that can be minimized (base class).

    Example of usage::

        cf=DerivedCostFunction()
        # ... calculate x ...
        args=cf.getArguments(x) # this could be potentially expensive!
        f=cf.getValue(x, *args)
        # ... it could be required to update x without using the gradient...
        # ... but then ...
        gf=cf.getGradient(x, *args) 

    The class distinguishes between the representation of the solution
    x (x-type) and the gradients (r-type). 

    :cvar provides_inverse_Hessian_approximation: This member should be set
            to ``True`` in subclasses that provide a valid implementation of
            `getInverseHessianApproximation()`
    """

    provides_inverse_Hessian_approximation=False

    def __call__(self, x, *args):
        """
        short for `getValue(x, *args)`.
        """
        return self.getValue(x, *args)

    def getArguments(self, x):
        """
        returns precalculated values that are shared in the calculation of
        *f(x)* and *grad f(x)* and the Hessian operator. The default
        implementation returns an empty tuple.
        
        :param x: location of derivative
        :type x: x-type 
        :rtype: ``tuple``
        """
        return ()
        
    def getDualProduct(self, x, r):
        """
        returns the dual product of ``x`` and ``r``
        
        :type x: x-type
        :type r: r-type
        :rtype: ``float``
        """
        raise NotImplementedError

    def getGradient(self, x, *args):
        """
        returns the gradient of *f* at *x* using the precalculated values for
        *x*.
        
        :param x: location of derivative
        :type x: x-type
        :param args: pre-calculated values for ``x`` from `getArguments()`
        :rtype: r-type
        """
        raise NotImplementedError
        
    def getInverseHessianApproximation(self, x, r, *args):
        """
        returns an approximative evaluation *p* of the inverse of the Hessian
        operator of the cost function for a given gradient *r* at a given
        location *x*: *H(x) p = r*
        
        :param x: location of Hessian operator to be evaluated
        :type x: x-type
        :param r: a given gradient 
        :type r: r-type
        :param args: pre-calculated values for ``x`` from `getArguments()`
        :rtype: x-type
        :note: In general it is assumed that the Hessian *H(x)* needs to be
               calculated in each call for a new location *x*. However, the
               solver may suggest that this is not required, typically when
               the iteration is close to completeness.
        :note: Subclasses that implement this method should set the class
               variable `provides_inverse_Hessian_approximation` to ``True`` to
               enable the solver to call this method.
        """
        raise NotImplementedError

    def getValue(self, x, *args):
        """
        returns the value *f(x)* using the precalculated values for *x*.
        
        :param x: a solution approximation
        :type x: x-type
        :rtype: ``float``
        """
        raise NotImplementedError

    def updateHessian(self):
        """
        notifies the class that the Hessian operator needs to be updated.
        This method is called by the solver class.
        """
        pass


class MeteredCostFunction(CostFunction):
    """
    This an intrumented version of the `CostFunction` class. The function
    calls update statistical information.
    The actual work is done by the methods with corresponding name and a
    leading underscore. These functions need to be overwritten for a particular
    cost function implementation.
    """

    def __init__(self):
        """
        the base constructor initializes the counters so subclasses should
        ensure the super class constructor is called.
        """
        super(MeteredCostFunction, self).__init__()
        self.resetCounters()

    def resetCounters(self):
        """
        resets all statistical counters
        """
        self.DualProduct_calls=0
        self.Value_calls=0
        self.Gradient_calls=0
        self.Arguments_calls=0
        self.InverseHessianApproximation_calls=0
        
    def getDualProduct(self, x, r):
        """
        returns the dual product of ``x`` and ``r``
        
        :type x: x-type
        :type r: r-type
        :rtype: ``float``
        """
        self.DualProduct_calls+=1
        return self._getDualProduct(x, r)
        
    def _getDualProduct(self, x, r):
        """
        returns the dual product of ``x`` and ``r``
        
        :type x: x-type
        :type r: r-type
        :rtype: ``float``
        :note: This is the worker for `getDualProduct()`, needs to be overwritten.
        """
        raise NotImplementedError
        
    def getValue(self, x, *args):
        """
        returns the value *f(x)* using the precalculated values for *x*.
        
        :param x: a solution approximation
        :type x: x-type
        :rtype: ``float``
        """
        self.Value_calls+=1
        return self._getValue(x, *args)
    
    def _getValue(self, x, *args):
        """
        returns the value *f(x)* using the precalculated values for *x*.
        
        :param x: a solution approximation
        :type x: x-type
        :rtype: ``float``
        :note: This is the worker for ``getValue()``, needs to be overwritten.
        """
        raise NotImplementedError
        
    def getGradient(self, x, *args):
        """
        returns the gradient of *f* at *x* using the precalculated values for
        *x*.
        
        :param x: location of derivative
        :type x: x-type
        :param args: pre-calculated values for ``x`` from `getArguments()`
        :rtype: r-type        
        """
        self.Gradient_calls+=1
        return self._getGradient(x, *args)

    def _getGradient(self, x, *args):
        """
        returns the gradient of *f* at *x* using the precalculated values for
        *x*.
        
        :param x: location of derivative
        :type x: x-type
        :param args: pre-calculated values for ``x`` from `getArguments()`
        :rtype: r-type  
        :note: This is the worker for `getGradient()`, needs to be overwritten.
        """      
        raise NotImplementedError
        

    def getArguments(self, x):
        """
        returns precalculated values that are shared in the calculation of
        *f(x)* and *grad f(x)* and the Hessian operator
        
        :param x: location of derivative
        :type x: x-type 
        """
        self.Arguments_calls+=1
        return self._getArguments(x)
    
    def _getArguments(self, x):
        """
        returns precalculated values that are shared in the calculation of
        *f(x)* and *grad f(x)* and the Hessian operator
        
        :param x: location of derivative
        :type x: x-type 
        :note: Overwrite this function to implement a specific cost function
        """
        return ()    
        
    def getInverseHessianApproximation(self, x, r,*args):
        """
        returns an approximative evaluation *p* of the inverse of the Hessian
        operator of the cost function for a given gradient *r* at a given
        location *x*: *H(x) p = r*
        
        :param x: location of Hessian operator to be evaluated.
        :type x: x-type
        :param r: a given gradient 
        :type r: r-type
        :param args: pre-calculated values for ``x`` from `getArguments()`
        :rtype: x-type
        :note: In general it is assumed that the Hessian *H(x)* needs to be calculate in each call for a new
          location *x*. However, the solver may suggest that this is not required, typically when the iteration 
          ius close to completness. 
        
        """
        self.InverseHessianApproximation_calls+=1
        return self._getInverseHessianApproximation(x, r, *args)
   
    def _getInverseHessianApproximation(self, x, r, *args):
        """
        returns an approximative evaluation *p* of the inverse of the Hessian operator of the costfunction
        for a given gradient type *r* at a given location *x*: *H(x) p = r*
        
        :param x: location of Hessian operator to be evaluated.
        :type x: x-type
        :param r: a given gradient 
        :type r: r-type
        :param args: pre-calculated values for ``x`` from ``getArguments()``
        :rtype: x-type
        :note: In general it is assumed that the Hessian *H(x)* needs to be
               calculate in each call for a new location *x*. However, the
               solver may suggest that this is not required, typically when
               the iteration is close to completness. 
        :note: This is the worker for getInverseHessianApproximation()`, needs
               to be overwritten.
        :note: Subclasses that implement this method should set the class
               variable `provides_inverse_Hessian_approximation` to ``True`` to
               enable the solver to call this method.
        """
        raise NotImplementedError

