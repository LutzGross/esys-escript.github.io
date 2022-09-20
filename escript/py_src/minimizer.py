#############################################################################
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

"""Generic minimization algorithms"""

from __future__ import print_function, division

__copyright__ = """Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__ = "https://launchpad.net/escript-finley"

__all__ = ['MinimizerException', 'MinimizerIterationIncurableBreakDown', 'LineSearchTerminated',
           'MinimizerMaxIterReached', 'AbstractMinimizer', 'MinimizerLBFGS', 'MinimizerNLCG', 'CostFunction']

import numpy as np
from .costfunctions import *
import logging

lslogger = logging.getLogger('esys.minimizer.linesearch')
zoomlogger = logging.getLogger('esys.minimizer.linesearch.zoom')


class MinimizerException(Exception):
    """
    This is a generic exception thrown by a minimizer.
    """
    pass


class MinimizerMaxIterReached(MinimizerException):
    """
    Exception thrown if the maximum number of iteration steps is reached.
    """
    pass


class MinimizerIterationIncurableBreakDown(MinimizerException):
    """
    Exception thrown if the iteration scheme encountered an incurable
    breakdown.
    """
    pass


class LineSearchTerminated(MinimizerException):
    """
    Exception thrown if the line serach was unsucessful
    """
    pass


def _zoom(Phi, gradPhi, getPhiArgs, alpha_lo, alpha_hi, phi_lo, phi_hi,
          phi_0, grad_phi_0,  changeMin=1.e-6, reduceMin=1.e-4, interpolationOrder=4, c1=1e-4, c2=0.9, iterMax=25, alphaMin=0.):
    """
    Helper function for `line_search` below which tries to tighten the range
    alpha_lo...alpha_hi. See Chapter 3 of 'Numerical Optimization' by
    J. Nocedal for an explanation.
    interpolation options are linear, quadratic, cubic or Newton interpolation.
    """

    def linearinterpolate(alpha_lo, alpha_hi, old_alpha):
        zoomlogger.debug("Linear interpolation is applied.")
        alpha = 0.5 * (alpha_lo + alpha_hi)
        return alpha, 1

    def quadinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi):
        if old_alpha is None:
            return linearinterpolate(alpha_lo, alpha_hi, old_alpha)
        zoomlogger.debug("Quadratic interpolation is applied.")
        denom = 2.0 * (old_phi - phi_0 - grad_phi_0 * old_alpha)
        if denom == 0:
            zoomlogger.debug("Root calculation failed (denom == 0). Falling back to linear interpolation.")
            return linearinterpolate(alpha_lo, alpha_hi, old_alpha)
        alpha = -grad_phi_0 * old_alpha * old_alpha / denom
        if np.abs(alpha - old_alpha) < changeMin or np.abs(alpha) < reduceMin * np.abs(old_alpha):
            zoomlogger.debug("Insufficient change. Falling back to linear interpolation.")
            return linearinterpolate(alpha_lo, alpha_hi, old_alpha)
        if alpha < alphaMin:
            zoomlogger.debug("alpha = %g is too small. Falling back to linear interpolation." % alpha)
            return linearinterpolate(alpha_lo, alpha_hi, old_alpha)
        return alpha, 2

    def cubicinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi, very_old_alpha, very_old_phi):
        if very_old_alpha is None:
            return quadinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi)
        zoomlogger.debug("Cubic interpolation is applied.")
        a0s, a1s = very_old_alpha * very_old_alpha, old_alpha * old_alpha
        denom = a0s * a1s * (old_alpha - very_old_alpha)
        if denom == 0:
            zoomlogger.debug("Root calculation failed (denom == 0). Falling back to quadratic interpolation.")
            return quadinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi)
        a0c, a1c = a0s * very_old_alpha, a1s * old_alpha
        tmpA = old_phi - phi_0 - grad_phi_0 * old_alpha
        tmpB = very_old_phi - phi_0 - grad_phi_0 * very_old_alpha
        a = (a0s * tmpA - a1s * tmpB) / denom
        b = (-a0c * tmpA + a1c * tmpB) / denom
        deter = b * b - 3.0 * a * grad_phi_0
        if deter < 0:
            zoomlogger.debug("Root calculation failed (deter < 0). Falling back to quadratic interpolation.")
            return quadinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi)
        alpha = (-b + np.sqrt(deter)) / (3.0 * a)
        if np.abs(alpha - old_alpha) < changeMin or np.abs(alpha) < reduceMin * np.abs(old_alpha):
            zoomlogger.debug("Insufficient change. Falling back to quadratic interpolation.")
            return quadinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi)
        if alpha < alphaMin:
            zoomlogger.debug("alpha = %g is too small. Falling back to quadratic interpolation." % alpha)
            return quadinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi)
        return alpha, 3

    def newtoninterpolate(alpha_data, phi_data, old_alpha, old_phi):
        # Interpolates using a polynomial of increasing order
        # The coefficients of the interpolated polynomial are found using the Newton method
        if very_old_alpha is None:
            return quadinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi)
        if old_alpha in alpha_data:
            return quadinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi)
        alpha_data.append(old_alpha)
        phi_data.append(old_phi)
        coefs = newton_poly(alpha_data, phi_data)
        alpha, ii = newtonroot(coefs, alpha_data, old_alpha)
        if alpha < alpha_lo or alpha > alpha_hi:  # Root is outside the domain
            zoomlogger.debug("Newton interpolation converged on a root outside of [alpha_lo,alpha_hi]. Falling back on cubic interpolation")
            return cubicinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi, very_old_alpha, very_old_phi)
        if np.abs(alpha - old_alpha) < changeMin or np.abs(alpha) < reduceMin * np.abs(old_alpha):
            alpha = 0.5 * old_alpha
            zoomlogger.debug("Insufficient change.  Falling back on cubic interpolation.")
            return cubicinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi, very_old_alpha, very_old_phi)
        if abs(alpha) <= alphaMin:
            zoomlogger.debug("Newton interpolation returned alpha <= alphaMin = %g. Falling back on cubic interpolation."%(alphaMin, ))
            return cubicinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi, very_old_alpha, very_old_phi)
        return alpha, len(alpha_data)

    def newton_poly(alpha_data, phi_data):
        # Returns the coefficients of the newton form polynomial of Phi(alpha)
        # for the points alpha_data and function values phi_data
        m = len(alpha_data)
        x = np.copy(alpha_data)
        a = np.copy(phi_data)
        for k in range(1, m):
            a[k:m] = (a[k:m] - a[k - 1]) / (x[k:m] - x[k - 1])
        return a

    def newtonroot(coefs, alpha_data, startingGuess):
        # Solves for the root of a polynomial using the newton method
        dfcoefs = list(range(1, len(coefs)))
        for k in range(0, len(coefs) - 1):
            dfcoefs[k] *= coefs[k + 1]
        tol = 1e-6
        maxiterations = 100
        point = startingGuess
        counter = 0
        while counter < maxiterations:
            counter += 1
            product = 1
            numer = coefs[0]
            denom = dfcoefs[0]
            for k in range(1, len(coefs)):
                product *= (point - alpha_data[k - 1])
                numer += product * coefs[k]
                if k < len(coefs) - 1:
                    denom += product * dfcoefs[k]
            if denom == 0:
                zoomlogger.debug("Higher order interpolation failed (denom==0). Falling back on cubic interpolation")
                return cubicinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi, very_old_alpha, very_old_phi)
            newpoint = float(point - numer / denom)
            error = abs(newpoint - point)
            if error < tol:
                zoomlogger.debug("Higher order interpolation converged in %d iterations. Got alpha = %g" % (counter, newpoint))
                return newpoint, counter
            else:
                point = newpoint
        zoomlogger.debug("Higher order interpolation failed (exceeded max iterations). Falling back on cubic interpolation")
        return cubicinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi, very_old_alpha, very_old_phi)

    iterCount = 0
    # Setup
    old_alpha = None
    old_phi = None
    very_old_alpha = None
    very_old_phi = None
    alpha_data = [0.0]
    phi_data = [phi_0]

    while iterCount <= iterMax:
        if interpolationOrder == 1:
            alpha, o = linearinterpolate(alpha_lo, alpha_hi, old_alpha)
        elif interpolationOrder == 2:
            alpha, o = quadinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi)
        elif interpolationOrder == 3:
            alpha, o = cubicinterpolate(alpha_lo, alpha_hi, old_alpha, old_phi, very_old_alpha, very_old_phi)
        else:
            alpha, o = newtoninterpolate(alpha_data, phi_data, old_alpha, old_phi)

        getPhiArgs(alpha)
        phi_a = Phi(alpha)
        zoomlogger.debug("Iteration %d, alpha=%g, phi(alpha)=%g (interpolation order = %d)" % (iterCount+1, alpha, phi_a, o))
        if phi_a > phi_0 + c1 * alpha * grad_phi_0 or phi_a >= phi_lo:
            very_old_alpha, very_old_phi = old_alpha, old_phi
            old_alpha, old_phi = alpha_hi, phi_hi
            alpha_hi, phi_hi = alpha, phi_a
        else:
            grad_phi_a = gradPhi(alpha)
            zoomlogger.debug("phi'(alpha)=%g" % (grad_phi_a, ))
            if np.abs(grad_phi_a) <= -c2 * grad_phi_0:
                zoomlogger.info("Zoom completed after %d steps: alpha=%g, phi(alpha)=%g." % (iterCount+1, alpha, phi_a))
                return alpha, phi_a
            if grad_phi_a * (alpha_hi - alpha_lo) >= 0:
                alpha_hi = alpha_lo
            very_old_alpha, very_old_phi = old_alpha, old_phi
            old_alpha, old_phi = alpha_lo, phi_lo
            alpha_lo, phi_lo = alpha, phi_a
        zoomlogger.debug("alpha range =[%g, %g]" % (alpha_lo, alpha_hi))
        if not alpha_hi > alpha_lo:
            raise LineSearchTerminated("Zoom is terminated (void search interval).")
        iterCount += 1
    raise LineSearchTerminated("Zoom is terminated (exceeded max iterations).")
    return alpha, phi_a


def line_search(F, m, p, grad_Fm, Fm, args_m, alpha=1.0, alphaMax=50.0, alphaMin=0.,
                c1=1e-4, c2=0.9, iterMax=15, zoom_interpolationOrder=1,
                zoom_relChangeMin=1e-6, zoom_reduceMin=1e-5, zoom_iterMax=25):
    """
    Line search method to minimize F(m+alpha*p)
    The iteration is terminated when the strong Wolfe conditions is fulfilled.
    See Chapter 3 of 'Numerical Optimization' by J. Nocedal for an explanation.

    :param F: callable objective function F(m) of type ``Costfunction``
    :param m: offset for the line search
    :param p: search direction
    :param grad_Fm: value for the gradient of F at m
    :param Fm: value of F(m)
    :param args_m: arguments for m
    :param alpha: initial step length. If grad_Fm is properly scaled alpha=1 is a
                  reasonable starting value.
    :param alphaMax: maximum value for alpha. Iteration is terminated if alpha > alphaMax
    :param alphaMin: minimum value for alpha. Iteration is terminated if alpha < alphaMin
    :param c1: value for Armijo condition (see reference)
    :param c2: value for curvature condition (see reference)
    :param iterMax: maximum number of line search iterations to perform
    :param zoom_interpolationOrder: the order of the interpolation used is Zoominb (1, 2, 3, 4)
    :param zoom_relChangeMin: the interpolated value of alpha, alpha_i, must differ from the previous
    :             value by at least this much i.e. abs(alpha_i-alpha_{i-1}) < zoom_relChangeMin
    :param zoom_reduceMin: minimal reduction alpha allowed
    :             If abs(alpha_i) < reduceMin*abs(alpha_{i-1}) the zooming reduces the interpolation order.
    :param zoom_iterMax: maximum number of zoom steps to perform
    :return: alpha, F, grad F, args for m+alpha*p (grad can be None)
    """
    # Fm_info this stores the latest grad F(m+a*p) which is returned
    if args_m is None:
        args_m = F.getArgumentsAndCount(m)
        Fm = None
        grad_Fm = None
    if Fm is None:
        Fm = F(m, *args_m)
        lslogger.debug("initial F(m) calculated= %s" % str(Fm))
    if grad_Fm is None:
        grad_Fm = F.getGradientAndCount(m, *args_m)
        lslogger.debug("initial grad F(m) calculated.")
    # this stores the latest gradf(x+a*p) which is returned
    Fm_info = [args_m, Fm, grad_Fm]

    def Phi(a):
        """ Phi(a):=F(x+a*p) """
        if Fm_info[0] is None:
            getPhiArgs(a)
        Fm_info[1] = F(m + a * p, *Fm_info[0])
        return Fm_info[1]

    def gradPhi(a):
        if Fm_info[0] is None:
            getPhiArgs(a)
        Fm_info[2] = F.getGradientAndCount(m + a * p, *Fm_info[0])
        return F.getDualProductAndCount(p, Fm_info[2])

    def getPhiArgs(a):
        args = F.getArgumentsAndCount(m + a * p)
        Fm_info[0] = args
        Fm_info[1] = None
        Fm_info[2] = None
        return args

    old_alpha = 0.
    phi_0 = Fm
    grad_phi_0 = F.getDualProductAndCount(p, grad_Fm)  # gradPhi(0., *args0)
    lslogger.info("Initial values: phi(0)=%g, phi'(0)=%g" % (phi_0, grad_phi_0))

    old_phi_a = phi_0
    phi_a = phi_0
    iterCount = 1
    if zoom_interpolationOrder not in [1, 2, 3, 4]:
        lslogger.debug("Setting interpolationOrder = 1")
        interpolationOrder = 1
    else:
        interpolationOrder = zoom_interpolationOrder
    sucess = False
    while iterCount < max(iterMax, 2) and alpha > alphaMin and not sucess:
        getPhiArgs(alpha)
        phi_a = Phi(alpha)
        lslogger.info("Iteration %d, alpha=%g, phi(alpha)=%g" % (iterCount, alpha, phi_a) )
        if (phi_a > phi_0 + c1 * alpha * grad_phi_0) or ((phi_a >= old_phi_a) and (iterCount > 1)):
            alpha, phi_a = _zoom(Phi, gradPhi, getPhiArgs, old_alpha, alpha, old_phi_a, phi_a, phi_0, grad_phi_0,
                                 interpolationOrder=interpolationOrder, changeMin=zoom_relChangeMin*alphaMax,
                                 reduceMin=zoom_reduceMin, c1=c1, c2=c2, iterMax=zoom_iterMax, alphaMin=alphaMin)
            sucess = True
            lslogger.info("XXXX")
        else:
            lslogger.info("ZZZZ")
            grad_phi_a = gradPhi(alpha)
            lslogger.info("\tphi'(alpha)=%g" % (grad_phi_a, ))
            if np.abs(grad_phi_a) <= -c2 * grad_phi_0:
                sucess = True
                lslogger.info("Strong Wolfe condition fulfilled.")
            elif grad_phi_a >= 0:
                alpha, phi_a = _zoom(Phi, gradPhi, getPhiArgs, alpha, old_alpha, phi_a, old_phi_a, phi_0, grad_phi_0,
                                 interpolationOrder=zoom_interpolationOrder, changeMin=zoom_relChangeMin*alphaMax,
                                 reduceMin=zoom_reduceMin, c1=c1, c2=c2, iterMax=zoom_iterMax, alphaMin=alphaMin)
                sucess = True
                lslogger.info("YYY")
            else:
                old_alpha = alpha
                old_phi_a = phi_a
                # the factor is arbitrary as long as there is sufficient increase
                alpha = 2. * alpha
                if alpha > alphaMax:
                    raise LineSearchTerminated("Line search terminated (Maximum alpha reached). ")
        iterCount += 1
    if not sucess:
        raise LineSearchTerminated("Line search terminated (Minimal alpha reached). ")
    lslogger.info("Line serach completed after %d steps (alpha=%g, phi(alpha)=%g)." % (iterCount, alpha, phi_a))
    return alpha, phi_a, Fm_info[2], Fm_info[0]  # returns F, grad F, args for m+alpha*p (grad can be None)


class AbstractMinimizer(object):
    """
    Base class for function minimization methods.
    """
    __F = None
    _m_tol = 1e-3
    _F_tol = None
    # maximum number of iteration steps
    _iterMax = 300

    # History size
    _truncation = 30

    # rescale search direction
    _scaleSearchDirection = True

    # Initial Hessian multiplier (if applicable)
    _initial_H = 1.

    # Restart after this many iteration steps
    _restart = 60

    # maximum number of line search steps
    _linesearch_iterMax = 20

    # interpolation order for search direction scaling alpha when zooming
    _zoom_interpolationOrder = 1

    # maximum number of zoom steps in line search
    _zoom_iterMax = 20

    # minimal reduction alpha allowed
    # If abs(alpha_i) < reduceMin * abs(alpha_{i - 1}) the zooming reduces
    # the interpolation order.
    _zoom_reduceMin = 1e-5

    # search direction scaling in line search  must differ from the previous
    # value by at least by _zoom_relChangeMin much i.e. abs(alpha_i-alpha_{i-1}) < zoom_relChangeMin
    _zoom_relChangeMin = 1e-6

    # range for alpha on line search
    _alphaMin = 0
    _alphaMax = 50
    
    # Wolfe condition paramters
    _c1 = 1e-4
    _c2 = 0.9

    def __init__(self, F=None, m_tol=1e-4, F_tol=None, iterMax=300):
        """
        Initializes a new minimizer for a given cost function.

        :param F: the cost function to be minimized
        :type F: `CostFunction`
        :param m_tol: terminate interations when relative change of the level set
                      function is less than or equal m_tol
        :type m_tol: float

        """
        self.setCostFunction(F)
        self._result = None
        self._callback = None
        self.logger = logging.getLogger('esys.minimizer.%s' % self.__class__.__name__)
        self.setOptions(m_tol=m_tol, F_tol=F_tol, iterMax=iterMax)

    def setCostFunction(self, F):
        """
        set the cost function to be minimized

        :param F: the cost function to be minimized
        :type F: `CostFunction`
        """
        self.__F = F

    def getCostFunction(self):
        """
        return the cost function to be minimized

        :rtype: `CostFunction`
        """
        return self.__F

    def setTolerance(self, m_tol=1e-4, F_tol=None):
        """
        Sets the tolerance for the stopping criterion. The minimizer stops
        when an appropriate norm is less than `m_tol`.
        """
        self.setOptions(m_tol=m_tol, F_tol=F_tol)

    def setMaxIterations(self, iterMax):
        """
        Sets the maximum number of iterations before the minimizer terminates.
        """
        self.setOptions(iterMax=iterMax)

    def getOptions(self):
        """
        returns a dictionary of LBFGS options
        rtype: dictionary with keys 'truncation', 'initialHessian', 'restart', 'linesearch_iterMax', 'zoom_iterMax'
        'interpolationOrder', 'zoom_relChangeMin', 'zoom_reduceMin'
        """
        return {'m_tol': self._m_tol, 'F_tol': self._F_tol, 'iterMax': self._iterMax, 'truncation': self._truncation, 
                'scaleSearchDirection': self._scaleSearchDirection, 'initialHessian': self._initial_H,
                'restart': self._restart,
                'linesearch_iterMax': self._linesearch_iterMax, 'zoom_iterMax': self._zoom_iterMax,
                'interpolationOrder': self._zoom_interpolationOrder, 'zoom_relChangeMin': self._zoom_relChangeMin, 'zoom_reduceMin': self._zoom_reduceMin,
                'alphaMin': self._alphaMin, 'alphaMax': self._alphaMax, 'c1': self._c2, 'c2': self._c2}

    def setOptions(self, **opts):
        """
        setOptions for LBFGS.  use       solver.setOptions( key = value)
        :key truncation: sets the number of previous LBFGS iterations to keep
        :type truncation : int
        :default truncation: 30
        :key scaleSearchDirection: if set the search direction is rescaled using an estimation of the norm of the Hessian
        :type scaleSearchDirection: ``bool``
        :default scaleSearchDirection: True
        :key initialHessian: first value used to rescale search direction.
                            This can be useful to a accelerate the line search in the first iteration step.
        :type initialHessian: positive ``float``
        :default initialHessian: 1
        :key restart: restart after this many iteration steps
        :type restart: int
        :default restart: 60
        :key linesearch_iterMax: maximum number of line search iterations
        :type linesearch_iterMax: int
        :default linesearch_iterMax: 25
        :key zoom_iterMax: maximum number of zoom iterations
        :type zoom_iterMax: int
        :default zoom_iterMax: 60
        :key interpolationOrder: order of the interpolation used for line search
        :type interpolationOrder: 1,2,3
        :default interpolationOrder: 1
        :key alphaMin: minimum search length 
        :type alphaMin: float
        :default alphaMin: 0.
        :key alphaMax: maximum search length 
        :type alphaMax: float
        :default alphaMax: 50.
        :key c1 : Wolfe condition parameter c1
        :type c1: float
        :default c1: 1e-4
        :key c2 : Wolfe condition parameter c2
        :type c2: float
        :default c2: 0.9
        :key zoom_relChangeMin: interpolated value of alpha, alpha_i, must differ
                   : from the previous value by at least this much
                   : abs(alpha_i-alpha_{i-1}) > zoom_relChangeMin * alphaMax  (line search)
        :type zoom_relChangeMin: float
        :default zoom_relChangeMin: 1e-6
        :key zoom_reduceMin: interpolated value of alpha must not be less than zoom_reduceMin
                   : abs(alpha_i) < zoom_reduceMin*abs(alpha_{i-1})
        :type zoom_reduceMin: float
        :default zoom_reduceMin: 1e-5

            Example of usage::
              cf=DerivedCostFunction()
              solver=MinimizerLBFGS(J=cf, m_tol = 1e-5, F_tol = 1e-5, iterMax=300)
              solver.setOptions(truncation=20, zoom_relChangeMin =1e-7)
              solver.run(initial_m)
              result=solver.getResult()
        """
        self.logger.debug("Setting options: %s" % (str(opts)))
        for o in opts:
            if o == "iterMax":
                self._iterMax = max(1, int(opts[o]))
            elif o == "m_tol":
                mtol = opts[o]
                if mtol is None:
                    self._m_tol = None
                else:
                    assert float(mtol) > 0, "m_tol must be positive or None"
                    self._m_tol = float(mtol)
            elif o == "F_tol":
                Ftol = opts[o]
                if Ftol is None:
                    self._F_tol = None
                else:
                    assert float(Ftol) > 0, "m_tol must be positive or None"
                    self._F_tol = float(Ftol)
            elif o == 'historySize' or o == 'truncation':
                assert opts[o] > 2, "Trancation must be greater than 2."
                self._truncation = max(0, int(opts[o]))
            elif o == 'scaleSearchDirection':
                self._scaleSearchDirection = opts[o]
            elif o == 'initialHessian':
                self._initial_H = opts[o]
            elif o == 'restart':
                self._restart = max(int(opts[o]), 1)
            elif o == 'linesearch_iterMax':
                self._linesearch_iterMax = max(int(opts[o]), 1)
            elif o == 'zoom_iterMax':
                self._zoom_iterMax = max(int(opts[o]), 1)
            elif o == 'alphaMax':
                self._alphaMax = max(float(opts[o]), 1.)
            elif o == 'alphaMin':
                self._alphaMin = max(float(opts[o]), 0.)
            elif o == 'c1':
                self._c1 = max(float(opts[o]), 0.)
            elif o == 'c2':
                self._c2 = max(float(opts[o]), 0.)
            elif o == 'interpolationOrder':
                self._zoom_interpolationOrder = max(min(int(opts[o]), 4), 1)
            elif o == 'zoom_relChangeMin':
                self._zoom_relChangeMin = max(float(opts[o]), 0.)
            elif o == 'zoom_reduceMin':
                self._zoom_reduceMin = max(float(opts[o]), 0.)
            else:
                raise KeyError("esysalid option '%s'" % o)

    def setCallback(self, callback):
        """
        Sets a callback function to be called after every iteration.
        It is up to the specific implementation what arguments are passed
        to the callback. Subclasses should at least pass the current
        iteration number k, the current estimate x, and possibly F(m),
        grad F(m), and the current error.
        """
        if callback is not None and not callable(callback):
            raise TypeError("Callback function not callable.")
        self._callback = callback

    def _doCallback(self, **args):
        if self._callback is not None:
            self._callback(**args)

    def getResult(self):
        """
        Returns the result of the minimization.
        """
        return self._result

    def run(self, m0):
        """
        Executes the minimization algorithm for *f* starting with the initial
        guess ``m0``.

        :return: the result of the minimization
        """
        raise NotImplementedError

    def logSummary(self):
        """
        Outputs a summary of the completed minimization process to the logger.
        """
        for l in self.getCostFunction().getStatistics().split("\n"):
            self.logger.info(l)


class MinimizerLBFGS(AbstractMinimizer):
    """
    Minimizer that uses the limited-memory Broyden-Fletcher-Goldfarb-Shanno
    method.
    See Chapter 6 of 'Numerical Optimization' by J. Nocedal for an explanation.
    """
    def run(self, m):
        """
        The callback function is called with the following arguments:
            k       - iteration number
            x       - current estimate
            Fm      - value of cost function at x
            grad_Fm    - gradient of cost function at x
            norm_dJ - ||Fm_k - Fm_{k-1}|| (only if F_tol is set)
            norm_dx - ||x_k - x_{k-1}|| (only if m_tol is set)

        :param m: initial guess
        :type m: m-type
        :return: solution
        :rtype: m-type
        """
        if self._scaleSearchDirection or self._initial_H is None:
            H_scale = None
        else:
            H_scale = self._initial_H
        # start the iteration:
        iterCount = 0
        iterCount_last_break_down = self._iterMax*1000
        alpha = 1.

        self._result = m
        args = self.getCostFunction().getArgumentsAndCount(m)
        grad_Fm = self.getCostFunction().getGradientAndCount(m, *args)
        Fm = self.getCostFunction().getValueAndCount(m, *args)
        Fm_0 = Fm
        self.logger.info("Initialization completed.")

        # TODO
        cbargs = {'k': iterCount, 'x': m, 'Fm': Fm, 'grad_Fm': grad_Fm}
        if self._F_tol:
            cbargs.update(norm_dJ=None)
        if self._m_tol:
            cbargs.update(norm_dx=None)
        self._doCallback(**cbargs)

        non_curable_break_down = False
        converged = False
        while not converged and not non_curable_break_down and iterCount < self._iterMax:
            k = 0
            break_down = False
            s_and_y = []
            self.getCostFunction().updateHessianAndCount(m, *args)
            self.logger.debug("Hessian is updated in step %d" % iterCount)
            # initial step length for line search

            while not converged and not break_down and k < self._restart and iterCount < self._iterMax:
                self.logger.info("********** iteration %3d **********" % iterCount)
                self.logger.info("\tF(m) = %g" % Fm)
                self.logger.debug("\tH = %s" % H_scale)

                # determine search direction
                p = -self._twoLoop(H_scale, grad_Fm, s_and_y)
        
                try:
                    alpha, Fm_new, grad_Fm_new, args_new = line_search(self.getCostFunction(), m, p, grad_Fm, Fm, args, alpha,
                                                                        c1=self._c1,
                                                                        c2=self._c2,
                                                                        zoom_interpolationOrder=self._zoom_interpolationOrder,
                                                                        alphaMin=self._alphaMin,
                                                                        alphaMax=self._alphaMax,
                                                                        iterMax=self._linesearch_iterMax,
                                                                        zoom_iterMax=self._zoom_iterMax)
                except LineSearchTerminated as e:
                    self.logger.debug("Line search failed: %s. BFGS is restarted." % str(e))
                    break_down = True
                    break
                # this function returns a scaling alpha for the search
                # direction as well as the cost function evaluation and
                # gradient for the new solution approximation x_new=x+alpha*p
                self.logger.debug("Search direction scaling found as alpha=%g" % alpha)
                # execute the step
                delta_m = alpha * p
                m_new = m + delta_m
                self._result = m_new
                converged = True
                if self._F_tol:
                    dF = abs(Fm_new - Fm)
                    Ftol_abs = self._F_tol * abs(Fm_0)
                    flag = dF <= Ftol_abs
                    if self.logger.isEnabledFor(logging.DEBUG):
                        if flag:
                            self.logger.info("Cost function has converged: dF=%g, F_tol=%g" % (dF, Ftol_abs))
                        else:
                            self.logger.info("Cost function checked: dJ=%g, F_tol=%g" % (dF, Ftol_abs))
                    cbargs.update(norm_dJ=dF)
                    converged = converged and flag

                if self._m_tol:
                    norm_m = self.getCostFunction().getNormAndCount(m_new)
                    norm_dm = self.getCostFunction().getNormAndCount(delta_m)
                    mtol_abs=norm_m * self._m_tol
                    flag = norm_dm <= mtol_abs
                    if flag:
                        self.logger.info("Solution has converged: dm=%g, m*m_tol=%g" % (norm_dm, mtol_abs))
                    else:
                        self.logger.info("Solution checked: dx=%g, x*m_tol=%g" % (norm_dm, mtol_abs))
                    cbargs.update(norm_dx=norm_dm)
                    converged = converged and flag

                if converged:
                    self.logger.info("********** iteration %3d **********" % (iterCount + 1,))
                    self.logger.info("\tF(m) = %g" % Fm_new)
                    break
                # unfortunately there is more work to do!
                if grad_Fm_new is None:
                    self.logger.debug("Calculating missing gradient.")
                    args_new = self.getCostFunction().getArgumentsAndCount(m_new)
                    grad_Fm_new = self.getCostFunction().getGradientAndCount(m_new, *args_new)

                delta_g = grad_Fm_new - grad_Fm
                rho = self.getCostFunction().getDualProductAndCount(delta_m, delta_g)
                if abs(rho) > 0:
                    s_and_y.append((delta_m, delta_g, rho))
                else:
                    self.logger.debug("Break down detected (<dm,dg>=0).")
                    break_down = True
                    break
                # move forward
                m = m_new
                grad_Fm = grad_Fm_new
                Fm = Fm_new
                args = args_new
                k += 1
                iterCount += 1
                cbargs.update(k=iterCount, x=m, Fm=Fm, grad_Fm=grad_Fm)
                self._doCallback(**cbargs)

                # delete oldest vector pair
                if k > self._truncation:
                    s_and_y.pop(0)

                if not break_down:
                    # an estimation of the inverse Hessian if it would be a scalar P=H_scale:
                    # the idea is that
                    #
                    #      H*dm=dg
                    #
                    #  if there is a inverse Hessian approximation provided we use
                    #
                    #     H_scale*|dm|^2 = <dm, dg>=rho ->   H_scale =rho/|dm|^2
                    if H_scale:
                        H_scale = rho/self.getCostFunction().getNormAndCount(delta_m) ** 2
                        self.logger.info("New search direction scale = %g." % H_scale)
            # case handling for inner iteration:
            if break_down:
                if not iterCount > iterCount_last_break_down:
                    non_curable_break_down = True
                    self.logger.warning(">>>>> Incurable break down detected in step %d." % (iterCount + 1) )
                else:
                    iterCount_last_break_down = iterCount
                    self.logger.debug("Break down detected in step %d. Iteration is restarted." % (iterCount + 1) )
            if not k < self._restart:
                self.logger.debug("Iteration is restarted after %d steps." % ( iterCount + 1) )

        # case handling for inner iteration:
        self._result = m
        if iterCount >= self._iterMax:
            self.logger.warning(">>>>>>>>>> Maximum number of iterations reached! <<<<<<<<<<")
            raise MinimizerMaxIterReached("Gave up after %d steps." % iterCount)
        elif non_curable_break_down:
            self.logger.warning(">>>>>>>>>> Incurable breakdown! <<<<<<<<<<")
            raise MinimizerIterationIncurableBreakDown("Gave up after %d steps." % iterCount)
        self.logger.info("Success after %d iterations!" % iterCount)
        return self._result

    def _twoLoop(self, H_scale, grad_Fm, s_and_y):
        """
        Helper for the L-BFGS method.
        See 'Numerical Optimization' by J. Nocedal for an explanation.
        """
        q = grad_Fm
        alpha = []
        for s, y, rho in reversed(s_and_y):
            a = self.getCostFunction().getDualProductAndCount(s, q) / rho
            alpha.append(a)
            q = q - a * y
        p = self.getCostFunction().getInverseHessianApproximationAndCount(q)
        norm_p = self.getCostFunction().getNormAndCount(p)
        if not norm_p > 0:
            raise MinimizerException("approximate Hessian inverse returns zero.")
        # this is if one wants
        if H_scale is not None:
            p_dot_q = self.getCostFunction().getDualProductAndCount(p, q)
            if abs(p_dot_q) > 0:
                #  we  would like to see that  h * norm(p)^2 = < H p, p> = < q, p> with h=norm(H)=H_scale
                # so we rescale p-> a*p ; h * a^2 * norm(p)^2 = a * < q, p> -> a = <q,p>/h/norm(p)**2
                a = p_dot_q * abs(p_dot_q) / H_scale / norm_p ** 2
                p *= a
                self.logger.debug("Search direction scaled : p=%g (scale factor = %g)" % (norm_p * a, a))
        for s, y, rho in s_and_y:
            beta = self.getCostFunction().getDualProductAndCount(p, y) / rho
            a = alpha.pop()
            p += s * (a - beta)
        return p


##############################################################################
class MinimizerNLCG(AbstractMinimizer):
    """
    Minimizer that uses the nonlinear conjugate gradient method
    (Fletcher-Reeves variant).
    """

    def run(self, x):
        """
        The callback function is called with the following arguments:
            k     - iteration number
            x     - current estimate
            Fm    - value of cost function at x
            grad_Fm  - gradient of cost function at x
            gnorm - norm of grad_Fm (stopping criterion)
        """
        # self.logger.setLevel(logging.DEBUG)
        # lslogger.setLevel(logging.DEBUG)
        i = 0
        k = 0
        args = self.getCostFunction().getArgumentsAndCount(x)
        r = -self.getCostFunction().getGradientAndCount(x, *args)
        Fm = self.getCostFunction()(x, *args)
        d = r
        delta = self.getCostFunction().getDualProductAndCount(r, r)
        delta0 = delta
        gnorm0 = Lsup(r)
        gnorm = gnorm0
        alpha = 1
        self._doCallback(k=i, x=x, Fm=Fm, grad_Fm=-r, gnorm=gnorm)

        while i < self._iterMax and gnorm > self._m_tol * gnorm0:
            self.logger.debug("iteration %d" % i)
            self.logger.debug("grad F(m) = %s" % (-r))
            self.logger.debug("d = %s" % d)
            self.logger.debug("x = %s" % x)

            alpha, Fm, Fm_info, args_new = line_search(self.getCostFunction(), x, d, -r, Fm, args, alpha=alpha,
                                                       c1=self._c1,
                                                       c2=self._c2,
                                                       zoom_interpolationOrder=self._zoom_interpolationOrder,
                                                       alphaMin=self._alphaMin,
                                                       alphaMax=self._alphaMax,
                                                       iterMax=self._linesearch_iterMax,
                                                       zoom_iterMax=self._zoom_iterMax)
            self.logger.debug("alpha=%g" % alpha)
            x = x + alpha * d
            args = args_new

            r = -self.getCostFunction().getGradientAndCount(x, *args) if Fm_info is None else -Fm_info
            # Fletcher-Reeves version
            delta_o = delta
            delta = self.getCostFunction().getDualProductAndCount(r, r)
            beta = delta / delta_o
            d = r + beta * d
            k = k + 1
            if self.getCostFunction().getDualProductAndCount(r, d) <= 0:
                d = r
                k = 0
            i += 1
            gnorm = Lsup(r)
            self._doCallback(k=i, x=x, Fm=Fm, grad_Fm=Fm_info, gnorm=gnorm)
            self._result = x

        if i >= self._iterMax:
            self.logger.warning(">>>>>>>>>> Maximum number of iterations %s reached! <<<<<<<<<<" % i)
            raise MinimizerMaxIterReached("Gave up after %d steps." % i)

        self.logger.info("Success after %d iterations! Initial/Final gradient=%e/%e" % (i, gnorm0, gnorm))
        return self._result


if __name__ == "__main__":
    # Example usage with function 'rosen' (minimum=[1,1,...1]):
    import numpy as np
    from scipy.optimize import rosen, rosen_der
    from costfunctions import CostFunction
    import sys

    N = 4
    x0 = np.array([4.] * N)  # initial guess


    class RosenFunc(CostFunction):
        def __init__(self):
            super(RosenFunc, self).__init__()
            self.provides_inverse_Hessian_approximation = False

        def getDualProductAndCount(self, f0, f1):
            return np.dot(f0, f1)

        def getValue(self, x, *args):
            return rosen(x)

        def getGradient(self, x, *args):
            return rosen_der(x)

        def getNorm(self, x):
            return np.norm(x)


    f = RosenFunc()
    m = None
    if len(sys.argv) > 1:
        method = sys.argv[1].lower()
        if method == 'nlcg':
            m = MinimizerNLCG(f)
        elif method == 'bfgs':
            m = MinimizerLBFGS(f)

    if m is None:
        # default
        m = MinimizerLBFGS(f)
        # m.setOptions(historySize=10000)

    logging.basicConfig(format='[%(funcName)s] \033[1;30m%(message)s\033[0m', level=logging.DEBUG)
    m.setTolerance(m_tol=1e-5)
    m.setMaxIterations(3000)
    m.run(x0)
    m.logSummary()
    print("\tLsup(result)=%.8f" % np.amax(abs(m.getResult())))

    # from scipy.optimize import fmin_cg
    # print("scipy ref=%.8f"%np.amax(abs(fmin_cg(rosen, x0, rosen_der, maxiter=10000))))
