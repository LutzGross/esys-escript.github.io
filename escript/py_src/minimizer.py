#############################################################################
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

"""Generic minimization algorithms"""


__copyright__ = """Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__ = "https://github.com/LutzGross/esys-escript.github.io"

__all__ = ['MinimizerException', 'MinimizerIterationIncurableBreakDown',
           'MinimizerMaxIterReached', 'AbstractMinimizer', 'MinimizerLBFGS', 'MinimizerNLCG', 'CostFunction']

import numpy as np
from .costfunctions import *
import logging, sys
EPSILON = sys.float_info.epsilon


class MinimizerException(Exception):
    """
    This is a generic exception thrown by a minimizer.
    """
    pass


class CostFunction1DEvaluationFactory(object):
    """
    This class manages the interface of ``Costfunction`` implementing F to the 1D projection
    Phi(a)=F(m+a*p). The idea is to use a ``EvaluatePhi`` object to hold the values of F, gradF and args
    for a given value.
    """
    def __init__(self, m, p, costfunction=None):
        """
        sets up the factory for offset `m` and direction `p` to use costfunction ``CostFunction``.
        :param m: offset
        :type m: m-type
        :param p: direction
        :type p: m-type
        :param costfunction: a cost function object
        :type costfunction: ``CostFunction``
        """
        self.m = m
        self.p = p
        self.costfunction = costfunction

    def getEvaluation(self, a=1):
        """
        Creates a container for Phi(a)=F(m+a*p), gradF(m+a*p) and the respective arguments.

        :param a: scaling factor of search direction
        :type a: ``float``
        :return: container holding evaluation state
        :rtype: `EvalutedPhi`
        """
        return EvalutedPhi(self, a)

    def __call__(self, a=0.):
        return self.getEvaluation(a)
    
    def getArgs(self, a):
        """

        returns the arguments for m+a*p
        :param a: scaling factor of search direction
        :type a: ``float``
        :returns args: arguments appropriate for ``Costfunction``
        """
        args = self.costfunction.getArgumentsAndCount(self.m + a * self.p)
        return args

    def getValue(self, a, args=None):
        """
        returns value F(m+a*p). If the arguments args==None, args is calculated
        :param a: scaling factor of search direction
        :type a: ``float``
        :param args: arguments appropriate for ``Costfunction``
        :type args: arguments appropriate for ``Costfunction`` or ``None``
        :returns: value of F and corresponding args.
        """
        if args is None:
            args = self.getArgs(a)
        v = self. costfunction.getValueAndCount(self.m + a * self.p, *args)
        return v, args

    def getGrad(self, a, gradF=None, args=None):
        """
        returns value Phi'(a)=<grad F(m+a*p), p>. If the arguments args==None, args is calculated.
        If the gradient gradF==None, the gradient is calculated (using args).
        :param a: scaling factor of search direction
        :type a: ``float``
        :param args: arguments appropriate for ``Costfunction``
        :type args: arguments appropriate for ``Costfunction`` or ``None``
        :param gradF: gradient of ``CostFunction`` for `m+alpha*p`
        :type gradF: appropriate for ``Costfunction`` or None
        :returns: derivate of Phi at value a
        """
        if args is None:
            args = self.getArgs(a)
        if gradF is None:
            gradF = self.costfunction.getGradientAndCount(self.m + a * self.p, *args)
        return self.costfunction.getDualProductAndCount(self.p, gradF), gradF, args


class EvalutedPhi(object):
    """
    this class provides an access to value Phi=F(m+a*p) for a `CostFunction`
    The interface is handeled by a factory class that holds m, p and the costfunction
    where methods of the latter are called via the factory when needed. Main purpose
    of this class is to maintain a simple mechanism to manage the values of Phi and its derivative
    Phi' (=gradPhi) but avoiding a recalculation of the arguments and full gradient once the Phi has been optimized
    via line serach and zoom.
    """
    def __init__(self, factory, alpha=1., args=None, valF=None, gradF=None):
        """
        :param factory: factory class
        :type factory: ``CostFunction1DEvaluationFactory`` or similar
        :param alpha: length factor
        :type alpha: ``float``
        :param args: argument for `m+alpha*p` if known. Is held at self.args
        :type args: appropriate for ``Costfunction``
        :param valF: value of Phi(alpha)=F(m+alpha*p) if known. Is held at self.valF
        :type valF: ``float``
        :param gradF: gradient of ``CostFunction`` for `m+alpha*p` if known. is held at self.gradF
        :type gradF: appropriate for ``Costfunction``
        """
        self.factory = factory
        self.alpha = alpha
        self.args = args
        self.valF = valF
        self.gradF = gradF
        self.gradPhi = None

    def getVal(self):
        """
        Returns the value for Phi(alpha)=F(m+alpha*p).

        In contrast to using self.valF directly, the value is calculated
        if self.valF is None.

        :return: the function value at alpha
        :rtype: ``float``
        """
        if self.valF is None:
            v, args = self.factory.getValue(self.alpha, self.args)
            self.args = args
            self.valF = v
        return self.valF

    def getDiff(self):
        """
        Returns the value for Phi'(alpha)=<p, grad(F(m+alpha*p))>.

        :return: the derivative of Phi at alpha
        :rtype: ``float``
        """
        if self.gradPhi is None:
            gp, g, args = self.factory.getGrad(self.alpha, self.gradF,  self.args)
            self.args = args
            self.gradF = g
            self.gradPhi = gp
        return self.gradPhi


class LineSearchTerminationError(MinimizerException):
    """
    Exception thrown if the line search was unsuccessful
    """
    pass


class LineSearchInterpolationBreakDownError(LineSearchTerminationError):
    """
    Interpolation break down
    """
    pass


def muchLarger(x, y, vareps=0.):
    """
    Tests if x is much larger than y with tolerance.

    Returns True if x > y and |x-y| > vareps * max(|x|,|y|).

    :param x: first value
    :type x: ``float``
    :param y: second value
    :type y: ``float``
    :param vareps: relative tolerance
    :type vareps: ``float``
    :return: True if x is much larger than y
    :rtype: ``bool``
    """
    L = vareps * max(abs(x), abs(y))
    return x > y and abs(x - y) > L


def muchSmaller(x, y, vareps=0.):
    """
    Tests if x is much smaller than y with tolerance.

    Returns True if x < y and |x-y| > vareps * max(|x|,|y|).

    :param x: first value
    :type x: ``float``
    :param y: second value
    :type y: ``float``
    :param vareps: relative tolerance
    :type vareps: ``float``
    :return: True if x is much smaller than y
    :rtype: ``bool``
    """
    L = vareps * max(abs(x), abs(y))
    return x < y and abs(x - y) > L


class LineSearchIterMaxReachedError(LineSearchTerminationError):
    """
    Exception thrown if the line search reaches maximum iteration count
    """
    pass


class LineSearchAlphaMinError(LineSearchTerminationError):
    """
    Exception thrown if the line search if minimum alpha is reached.
    """
    pass


class LineSearchAlphaMaxError(LineSearchTerminationError):
    """
    Exception thrown if the line search if maximum alpha is reached.
    """
    pass


class LineSearchSearchDirectionError(LineSearchTerminationError):
    """
    Exception thrown if The search direction is not a descent direction.
    """
    pass


class LineSearch(object):
    """
    Line search method to minimize F(m+alpha*p)
    The iteration is terminated when the strong Wolfe conditions is fulfilled.

    See Chapter 3 of 'Numerical Optimization' by J. Nocedal for an explanation.

    Note: the line search return a new search direction scale alpha that satisfies the strong Wolfe condition:
          (W1) sufficient decrease condition:  Phi(alpha) <= Phi(0) + c1 * alpha * phi'(0)
          (W2) curveture condition abs(Phi'(alpha)) <= c2 * abs(Phi'(0))

    The first step is to find an alpha_lo and alpha_hi such that the following conditions hold:
          (Z1) [alpha_lo, alpha_hi] contains an alpha fullfilling (W1)+(W2)
          (Z2) alpha_lo is giving the smallest value for Phi amongst alphas generated so far and that
             are satisfying the sufficient decrease condition (W1)
          (Z3) alpha_hi is chosen so that Phi′(alpha_lo)(alpha_hi − alpha_lo ) < 0
                that means alpha_hi  > alpha_lo if Phi′(alpha_lo)<0
                and alpha_hi  < alpha_lo if Phi′(alpha_lo)>0

    Condition (Z1) is fullfilled if one of the following conditions is satisfied:
            (E1) alpha_hi violates the sufficient decrease condition (W1)
            (E2) Phi(αlpha_hi) ≥ Phi(αlpha_lo)
            (E3) Phi′(alpha_lo)(alpha_hi − alpha_lo )  < 0.


    """
    _phiEpsilon = np.sqrt(EPSILON)
    _alphaMin = EPSILON
    _alphaMax = 5000
    _overStepFactor = 2
    _iterMax = 20
    _c1 = 1e-4
    _c2 = 0.9
    _inter_order = 3
    _inter_tol = 1e-6
    _inter_iterMax = 100
    _alphaOffset = 1e-4
    _alphaWidthMin = np.sqrt(EPSILON)
    _zoom_iterMax = 30
    _zoom_reductionMin = 0.66

    def __init__(self, logger=None):
        """
        initialize the line search solver

        :param logger: logger to be used. If not set 'esys.minimizer' is used.
        :type logger: ``logging.Logger``
        """
        if logger is None:
            logger = logging.getLogger('esys.minimizer')
        self.lslogger = logger.getChild('linesearch')
        self.zoomlogger = self.lslogger.getChild('zoom')

    def getOptions(self):
        """
        returns a dictionary of LBFGS options
        rtype: dictionary
        """
        return {'alphaMin': self._alphaMin, 'alphaMax': self._alphaMax, 'overStepFactor': self._overStepFactor,
                'iterMax': self._iterMax, 'c1': self._c1, 'c2': self._c2, 'phiEpsilon': self._phiEpsilon,
                'inter_order': self._inter_order, 'inter_tol': self._inter_tol, 'inter_iterMax': self._inter_iterMax,
                'alphaOffset': self._alphaOffset, 'alphaWidthMin': self._alphaWidthMin,
                'zoom_iterMax': self._zoom_iterMax, 'zoom_reductionMin': self._zoom_reductionMin}

    def setOptions(self, **opts):
        """
        set options for the line search.


        :key alphaMin: minimum search length
        :type alphaMin: float
        :default alphaMin: 1e-20
        :key alphaMax: maximum search length
        :type alphaMax: float
        :default alphaMax: 5000
        :key overStepFactor : factor to increase step size in line search
        :type overStepFactor: float
        :default overStepFactor: 2.
        :key iterMax: maximum number of line search iterations
        :type iterMax: int
        :default iterMax: 25
        :key c1: sufficient decrease condition factor c1
        :type c1: float
        :default c1: 1e-4
        :key c2: curvature condition factor c2
        :type c2: float
        :default c2: 0.9
        :key inter_order: order of the interpolation used for line search
        :type inter_order: 1,2,3
        :default inter_order: 3
        :key inter_iterMax: maximum number of iteration steps to when minimizing interploted cost function
        :type inter_iterMax: int
        :default inter_iterMax: 100
        :key inter_tol: tolerance to when minimizing interploted cost function
        :type inter_tol: float
        :default inter_tol: 1.
        :key zoom_iterMax: maximum number of zoom iterations
        :type zoom_iterMax: int
        :default zoom_iterMax: 20

        :key phiEpsilon : tolerance for `greater than` check of cost function values
        :type phiEpsilon: float
        :default phiEpsilon: ``np.sqrt(EPSILON)``

        :key  alphaOffset : minimal relative distance of new alpha from boundaries
        :type alphaOffset : float
        :default alphaOffset : 1e-4

        :key alphaWidthMin : minimal relative distance of new alpha from boundaries
        :type alphaWidthMin: float
        :default alphaWidthMin: ``np.sqrt(EPSILON)``
        :key zoom_reductionMin: minimal reduction search interval length between zoom steps
        :type zoom_reductionMin : float
        :default zoom_reductionMin :0.66
        """
        self.lslogger.debug("Setting options: %s" % (str(opts)))
        for o in opts:
            if o == 'alphaMin':
                self._alphaMin = max(float(opts[o]), 0.)
            elif o == 'alphaMax':
                self._alphaMax = max(float(opts[o]), 1.)
            elif o == 'iterMax':
                self._iterMax = max(int(opts[o]), 1)
            elif o == 'c1':
                self._c1 = max(float(opts[o]), 0.)
            elif o == 'c2':
                self._c2 = max(float(opts[o]), 0.)
            elif o == 'overStepFactor':
                self._overStepFactor = max(float(opts[o]), 1.)
            elif o == 'inter_tol':
                self._inter_tol = float(opts[o])
            elif o == 'inter_iterMax':
                self._inter_iterMax = max(int(opts[o]), 0)
            elif o == 'inter_order' or o == 'interpolationOrder':
                self._inter_order = max(int(opts[o]), 1)
            elif o == 'phiEpsilon':
                self._phiEpsilon = max(float(opts[o]), 0.)
            elif o == 'alphaWidthMin':
                self._alphaWidthMin = max(float(opts[o]), EPSILON)
            elif o == 'alphaOffset':
                self._alphaOffset = max(float(opts[o]), EPSILON)
            elif o == 'zoom_iterMax':
                self._zoom_iterMax = max(int(opts[o]), 1)
            elif o == 'zoom_reductionMin':
                self._zoom_reductionMin = min(max(float(opts[o]), 0.), 1)
            else:
                raise KeyError("invalid option '%s'" % o)

    def run(self, phi_0, alpha=1.0):
        """
        Line search method to minimize F(m+alpha*p)
        The iteration is terminated when the strong Wolfe conditions is fulfilled.
        See Chapter 3 of 'Numerical Optimization' by J. Nocedal for an explanation.
    
        :param phi_0: ``EvaluatePhi`` for alpha=0
        :param alpha: initial step length. If grad_Fm is properly scaled alpha=1 is a
                      reasonable starting value.
        :return: alpha
    
        """
        # testing parameters:
        assert self._c1 < self._c2, "You need to set c1 < c2."
        assert self._overStepFactor > 1, "overStepFactor (=%g) must be much larger then one." % self._overStepFactor
        assert 0 < self._alphaMin, "alphaMin (=%g) must be positive." % self._alphaMin
        assert self._alphaMin < self._alphaMax, "self._alphaMax (=%g) must be greater alphaMin (=%g)." % (self._alphaMax, self._alphaMin)
        assert alpha < self._alphaMax, "alphaMax (=%g) needs to be greater than initial alpha (=%g)." % (self._alphaMax, alpha)
        assert alpha > self._alphaMin, "alphaMin (=%g) needs to be less than initial alpha (=%g)." % (self._alphaMin, alpha)
        assert self._alphaWidthMin > 0, "alphaWidthMin (=%g) must be positive."  % self._alphaWidthMin
        assert self._alphaOffset > 0 and self._alphaOffset < 1, "alphaOffset (=%g) must be in ]0, 0.5[ " % self._alphaOffset

        Phi = phi_0.factory
        value_phi_0 = phi_0.getVal()
        diff_phi_0 = phi_0.getDiff()


        self.lslogger.info("Initial values: phi(0)=%g, phi'(0)=%g, alpha= %g." % (value_phi_0, diff_phi_0, alpha))
        if not diff_phi_0 < 0.:
            raise LineSearchSearchDirectionError("The search direction is not a descent direction")
    
        success = False
        iterCount = 1

        phi_alpha_old=phi_0
        phi_alpha = Phi(alpha)
        while not success:
            #phi_alpha, phi_alpha_old = phi_alpha_new, phi_alpha
            alpha = phi_alpha.alpha
            value_phi_alpha = phi_alpha.getVal()
            self.lslogger.info("Iteration %d, alpha=%g, phi(alpha)=%g" % (iterCount, alpha, value_phi_alpha))
            # if (value_phi_alpha >  value_phi_0 + self._c1 * alpha *  diff_phi_0) or \
            #        ( (not value_phi_alpha < phi_alpha_old.getVal() ) and (iterCount > 1)):
            if muchLarger(value_phi_alpha, value_phi_0 + self._c1 * alpha * diff_phi_0, vareps=self._phiEpsilon) or \
                    (not muchLarger(phi_alpha_old.getVal(), value_phi_alpha,  vareps=self._phiEpsilon) and iterCount > 1):
                self.lslogger.debug("Sufficient decrease condition or decrease criterium are violated -> start zoom.")
                # [ alpha_old, alpha] meets (E1) or (E2) (alpha_hi=alpha) -> (Z1) is fullfilled.
                # Phi(alpha_old) 
                phi_alpha_new, success = self.zoom(phi_alpha_old, phi_alpha, phi_0)
            else:
                # alpha meets sufficient decrease condition (W1) now!
                diff_phi_alpha = phi_alpha.getDiff()
                self.lslogger.debug("phi'(alpha)=%g" % (diff_phi_alpha,))
                if np.abs(diff_phi_alpha) <= -self._c2 * diff_phi_0:  # check curvature condition (W2)
                    phi_alpha_new, success = phi_alpha, True
                    self.lslogger.debug("Strong Wolfe condition is fulfilled. We are done.")
                elif diff_phi_alpha >= 0:  #
                    # [ alpha, alpha_old ] has alpha with strong wolf condition
                    self.lslogger.debug("Positive phi'(alpha).")
                    phi_alpha_new, success = self.zoom(phi_alpha, phi_alpha_old, phi_0)
            #    else:
            if not success:
                    phi_alpha=phi_alpha_old
                    # the factor is arbitrary as long as there is sufficient increase
                    alpha_new = self._overStepFactor * alpha + (1 - self._overStepFactor) * self._alphaMin
                    self.lslogger.debug("Search interval extended.")
                    if alpha_new > self._alphaMax:
                        raise LineSearchAlphaMaxError(
                            "Line search terminated as maximal alpha (=%d) reached. " % self._alphaMax)
                    phi_alpha_new, success = Phi(alpha_new), False
            iterCount += 1
            phi_alpha, phi_alpha_old = phi_alpha_new, phi_alpha
            if iterCount > self._iterMax:
                raise LineSearchIterMaxReachedError("Maximum number of iteration steps (=%d) reached." % self._iterMax)
    
        self.lslogger.info("Line search completed after %d steps (alpha=%g, phi(alpha)=%g)." % (iterCount-1, phi_alpha.alpha, phi_alpha.getVal()))
        return phi_alpha

    def zoom(self, phi_lo, phi_hi, phi_0):
        """
        Helper function for `line_search` below which tries to tighten the range
        phi_lo.alpha...phi_hi.alpha such that phi_lo, phi_hi fulfill one of these conditions
        (which needs to be fullfilled at entry):
        
        (Z1) [alpha_lo, alpha_hi] contains an alpha fullfilling (W1)+(W2). This condition holds if
            (E1) alpha_hi violates the sufficient decrease condition (W1)
            (E2) Phi(αlpha_hi) ≥ Phi(αlpha_lo)
            (E3) Phi′(alpha_lo)(alpha_hi − alpha_lo )  < 0.
        (Z2) alpha_lo is giving the smallest value for Phi amongst alphas generated so far and that are satisfying
            the sufficient decrease condition (W1)
        (Z3) alpha_hi is chosen so that Phi′(alpha_lo)(alpha_hi − alpha_lo ) < 0
        
        
        See Chapter 3 of 'Numerical Optimization' by
        J. Nocedal for an explanation.
        interpolation options are linear, quadratic, cubic or Newton interpolation.
        """
        def addToHistory(phi):
            if self._inter_order > 1:
                hist_phi.insert(0, phi)
            if self._inter_order <= 2:
                if len(hist_phi) > 2:
                    hist_phi.pop(-2)
            elif self._inter_order <= 3:
                if len(hist_phi) > 3:
                    hist_phi.pop(-2)
            else:
                if len(hist_phi) > self._inter_order + 1:
                    hist_phi.pop(-2)
        assert self._zoom_reductionMin > 0 and self._zoom_reductionMin <= 1., "zoom_reductionMin (=%g) must be in ]0,1]."
        # let the show begin:
        iterCount = 1
        if phi_hi == phi_0:
            hist_phi = [phi_lo, phi_0]
        else:
            hist_phi = [phi_hi, phi_0]
        # note: hist_phi[-1]==phi_0 at all times.
        Phi = phi_0.factory
        # interval of potential alpha's  with (Z1)
        alpha_min = min(phi_lo.alpha, phi_hi.alpha)
        alpha_max = max(phi_lo.alpha, phi_hi.alpha)
        width = alpha_max - alpha_min

        while iterCount <= self._zoom_iterMax:
            self.zoomlogger.debug("Iteration %d: alpha range =[%g, %g] (width= %g)" % (iterCount, alpha_min, alpha_max, width))
            self.zoomlogger.debug("Interpolation nodes: " + str([p.alpha for p in hist_phi]))
            # get a new alpha:
            for k in range(min(len(hist_phi), self._inter_order), 0, -1):
                reduceOrder = False
                try:
                    if k == 1:
                        new_alpha = self._bisection(phi_lo, phi_hi)
                    elif k == 2:
                        new_alpha = self._quadinterpolate(hist_phi)
                    elif k == 3:
                        new_alpha = self._cubicinterpolate(hist_phi)
                    else:
                        raise NotImplementedError("Higher order interpolation is not supported yet.")
                        new_alpha = self._newtoninterpolate(hist_phi)
                except LineSearchInterpolationBreakDownError as e:
                    # if there is any break down order is reduced: (should not happen with bisection)
                    reduceOrder = True
                    self.zoomlogger.debug("Interpolation order %d failed: %s" % (k, repr(e)))
                else:
                    # now we safeguard the new alpha:
                    # none of these should happen with bisection.
                    if new_alpha < alpha_min or new_alpha > alpha_max:
                        reduceOrder = True
                        self.zoomlogger.debug("New alpha %g outside alpha search interval." % new_alpha)
                    elif new_alpha < self._alphaMin:
                        reduceOrder = True
                        self.zoomlogger.debug("alpha = %g is too small." % new_alpha)
                    elif new_alpha  < alpha_min+ self._alphaOffset * width:
                        reduceOrder = True
                        self.zoomlogger.debug("New alpha %g to close to left alpha bound %g." % (new_alpha, alpha_min))
                    elif new_alpha > alpha_max - self._alphaOffset * width :
                        reduceOrder = True
                        self.zoomlogger.debug("New alpha %g to close to right alpha bound %g." % (new_alpha, alpha_max))
                    elif max(new_alpha - alpha_min, alpha_max - new_alpha ) > self._zoom_reductionMin * width:
                        reduceOrder = True
                        self.zoomlogger.debug("Insufficient search interval reduction.")
                if not reduceOrder:
                    break
        
            # get new value of costfunction value:
            new_phi = Phi(new_alpha)
            new_phi_alpha = new_phi.getVal()
            self.zoomlogger.debug(
                "Iteration %d, alpha=%g, phi(alpha)=%g (interpolation order = %d)" % (iterCount, new_alpha, new_phi_alpha, k))
            if muchLarger(new_phi_alpha, hist_phi[-1].getVal() + self._c1 * new_alpha * hist_phi[-1].getDiff(), vareps=self._phiEpsilon) \
                    or muchLarger(new_phi_alpha, phi_lo.getVal(), vareps=self._phiEpsilon):  # new_phi_alpha >= phi_lo.getVal():
                # new alpha violates E1 or E2 with phi_hi=new_phi
                phi_hi = new_phi
                addToHistory(phi_hi)
                self.zoomlogger.debug("Iteration %d: alpha_hi updated, now = %g" % (iterCount, phi_hi.alpha))
            else:
                new_diff_phi_alpha = new_phi.getDiff()
                self.zoomlogger.debug("phi'(alpha)=%g" % (new_diff_phi_alpha,))
                # check curvature condition (W2):
                #if not muchLarger(np.abs(new_diff_phi_alpha), -self._c2 * hist_phi[-1].getDiff(), vareps=self._phiEpsilon):
                if np.abs(new_diff_phi_alpha) <= -self._c2 * hist_phi[-1].getDiff():
                    # now alpha_new has (W1)+(W2)
                    self.zoomlogger.info("Zoom completed after %d steps: alpha=%g, phi(alpha)=%g, phi'(alpha)=%g." %
                                         (iterCount, new_alpha, new_phi_alpha, new_diff_phi_alpha))
                    self.zoomlogger.debug("Strong Wolfe condition is fulfilled. We are done.")
                    return new_phi, True
                # new we have (W1) and we need to set  phi_lo=phi_new
                # To we need to meet (E3) which determines phi_hi:
                if new_diff_phi_alpha * (phi_hi.alpha - phi_lo.alpha) >= 0:
                    phi_lo, phi_hi = new_phi, phi_lo
                else:
                    phi_lo, phi_hi = new_phi, phi_hi
                self.zoomlogger.debug("Iteration %d: alpha_lo updated, now = %g" % (iterCount, phi_lo.alpha))
                addToHistory(new_phi)

            alpha_min = min(phi_lo.alpha, phi_hi.alpha)
            alpha_max = max(phi_lo.alpha, phi_hi.alpha)
            width, width2 = alpha_max - alpha_min, width

            if not width < self._zoom_reductionMin * width2:
                self.zoomlogger.debug(
                    "Zoom is terminated due to insufficient search interval reduction (%g -> %g)." % (width2, width))
                return new_phi, False

            if width < self._alphaWidthMin:
                self.zoomlogger.debug("Zoom is terminated due to small search interval for alpha (=%g)." % width)
                return new_phi, False

            iterCount += 1
        self.zoomlogger.debug("Zoom is terminated (exceeded max iterations).")
        return new_phi, False

    def _bisection(self, phi_lo, phi_hi):
        """
        applies bisection
        """
        self.zoomlogger.debug("Bisection is applied.")
        if phi_lo.alpha < phi_hi.alpha:
            alpha = 0.5 * (max(phi_lo.alpha, self._alphaMin) + phi_hi.alpha)
        else:
            alpha = 0.5 * (max(phi_hi.alpha, self._alphaMin) + phi_lo.alpha)
        return alpha

    def _quadinterpolate(self, hist_phi):
        """
        find next alpha by quadratic  interpolation
        P(0)=phi(0), P'(0)=phi'(0),  P(alpha[0])=phi(alpha[0])
        """
        assert len(hist_phi) > 1
        self.zoomlogger.debug("Quadratic interpolation is applied.")
        denom = 2.0 * (hist_phi[0].getVal() - hist_phi[-1].getVal() - hist_phi[-1].getDiff() * hist_phi[0].alpha)
        if denom == 0:
            self.zoomlogger.debug("Root calculation failed (denom == 0).")
            raise LineSearchInterpolationBreakDownError("Root calculation failed (denom == 0).")
        alpha = - hist_phi[-1].getDiff() * hist_phi[0].alpha * hist_phi[0].alpha / denom
        return alpha

    def _cubicinterpolate(self, hist_phi):
        """
        find next alpha by quadratic  interpolation
        P(0)=phi(0), P'(0)=phi'(0),  P(alpha[0])=phi(alpha[0]), P(alpha[1])=phi(alpha[1])
        """
        assert len(hist_phi) > 2
        self.zoomlogger.debug("Cubic interpolation is applied.")
        a0s, a1s = hist_phi[1].alpha * hist_phi[1].alpha, hist_phi[0].alpha * hist_phi[0].alpha
        denom = a0s * a1s * (hist_phi[0].alpha - hist_phi[1].alpha)
        if denom == 0:
            self.zoomlogger.debug("Root calculation failed (denom == 0).")
            raise LineSearchInterpolationBreakDownError("Root calculation failed (denom == 0).")
        a0c, a1c = a0s * hist_phi[1].alpha, a1s * hist_phi[0].alpha
        tmpA = hist_phi[0].getVal() - hist_phi[-1].getVal() - hist_phi[-1].getDiff() * hist_phi[0].alpha
        tmpB = hist_phi[1].getVal() - hist_phi[-1].getVal() - hist_phi[-1].getDiff() * hist_phi[1].alpha
        a = (a0s * tmpA - a1s * tmpB) / denom
        b = (-a0c * tmpA + a1c * tmpB) / denom
        deter = b * b - 3.0 * a * hist_phi[-1].getDiff()
        if deter < 0:
            self.zoomlogger.debug("Root calculation failed (deter < 0).")
            raise LineSearchInterpolationBreakDownError("Root calculation failed (deter < 0).")
        alpha = (-b + np.sqrt(deter)) / (3.0 * a)
        return alpha

    def _newtoninterpolate(self, hist_phi):
        # Interpolates using a polynomial of increasing order
        # The coefficients of the interpolated polynomial are found using the Newton method
        assert len(hist_phi) > 1
        # calculat Newton DD:
        coefs, alphas = self.__makeNewtonPolynomial(hist_phi)
        alpha = self.__getNewtonPolynomialRoot(coefs, alphas, startingGuess=hist_phi[0].alpha)
        return alpha

    def __makeNewtonPolynomial(self, phis):
        # Returns the coefficients of the newton form polynomial of Phi(alpha)
        # for the points alpha_data and function values phi_data
        m = len(phis)
        x = np.array([p.alpha for p in phis])
        a = np.array([p.getVal() for p in phis])
        for k in range(1, m):
            a[k:m] = (a[k:m] - a[k - 1]) / (x[k:m] - x[k - 1])
        return a, x

    def __getNewtonPolynomialRoot(self, coefs, alpha_data, startingGuess=1.):
        # Solves for the root of a polynomial using the newton method
        dfcoefs = list(range(1, len(coefs)))
        for k in range(0, len(coefs) - 1):
            dfcoefs[k] *= coefs[k + 1]
        point = startingGuess
        counter = 0
        while counter < self._inter_iterMax:
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
                self.zoomlogger.debug("Root calculation failed (denom == 0).")
                raise LineSearchInterpolationBreakDownError("Root calculation failed (denom == 0).")
            newpoint = float(point - numer / denom)
            error = abs(newpoint - point)
            if error < self._inter_tol * max(abs(newpoint), abs(point)):
                self.zoomlogger.debug("Higher order interpolation converged in %d iterations. Got alpha = %g" % (counter, newpoint))
                return newpoint
            else:
                point = newpoint
        self.zoomlogger.debug("Higher order interpolation failed (exceeded max iterations).")
        LineSearchInterpolationBreakDownError("Higher order interpolation failed (exceeded max iterations).")
        return point


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


class AbstractMinimizer(object):
    """
    Base class for function minimization methods.
    """
    __line_search = None
    __F = None
    # see setOption for explanation
    _m_tol = 1e-4
    _grad_tol = 1e-8
    _iterMax = 300
    _truncation = 30
    _restart = 60
    _relAlphaMin = 1.e-8
    _scaleSearchDirection = True
    _initialAlpha = 1.

    # updated during iteration
    __initializeHessian = None

    def __init__(self, F=None, m_tol=1e-4, grad_tol=1e-8, iterMax=300, logger=None):
        """
        Initializes a new minimizer for a given cost function.

        :param F: the cost function to be minimized
        :type F: `CostFunction`
        :param m_tol: relative tolerance for solution `m` for termination of iteration
        :type m_tol: `float`
        :default m_tol: 1e-4
        :param grad_tol: tolerance for gradient relative to initial costfunction value for termination of iteration
        :type grad_tol: `float`
        :default grad_tol: 1e-8
        :parameter iterMax: maximium number of iterations.
        :type iterMax: `int`
        :default iterMax: 300
        """
        if logger is None:
            self.logger = logging.getLogger('esys.minimizer.%s' % self.__class__.__name__)
            self.logger.setLevel(logging.ERROR)
        else:
            self.logger = logger
        self.setCostFunction(F)
        self.__line_search = LineSearch(self.logger)
        self._result = None
        self._callback = None
        self.setOptions(m_tol=m_tol, grad_tol=grad_tol, iterMax=iterMax)

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

    def getLineSearch(self):
        """
        returns the ``LineSearch`` object being used.
        """
        return self.__line_search

    def setTolerance(self, m_tol=1e-4, grad_tol=1e-4):
        """
        Sets the tolerance for the stopping criterion. The minimizer stops
        when an appropriate norm is less than `m_tol`.

        :param m_tol: relative tolerance for solution `m` for termination of iteration
        :type m_tol: `float`
        :default m_tol: 1e-4
        :param grad_tol: tolerance for gradient relative to initial costfunction value for termination of iteration
        :type grad_tol: `float`
        :default grad_tol: 1e-8
        """
        self.setOptions(m_tol=m_tol, grad_tol=grad_tol)

    def setIterMax(self, iterMax=300):
        """
        Sets the maximum number of iterations before the minimizer terminates.

        :parameter iterMax: maximium number of iterations.
        :type iterMax: `int`
        :default iterMax: 300
        """
        self.setOptions(iterMax=iterMax)

    def getOptions(self):
        """
        returns a dictionary of LBFGS options
        rtype: dictionary
        """
        return {'m_tol': self._m_tol, 'grad_tol': self._grad_tol, 'iterMax': self._iterMax, 'truncation': self._truncation,
                'restart': self._restart, 'relAlphaMin': self._relAlphaMin,
                'scaleSearchDirection': self._scaleSearchDirection, 'initialAlpha': self._initialAlpha,
                'line_search': self.getLineSearch().getOptions()}

    def setOptions(self, **opts):
        """
        setOptions for LBFGS.  use       solver.setOptions( key = value)

        :key m_tol: relative tolerance for solution `m` for termination of iteration
        :type m_tol: `float`
        :default m_tol: 1e-4
        :key grad_tol: tolerance for gradient relative to initial cost function value for termination of iteration
        :type grad_tol: `float` or None
        :default grad_tol: 1e-4
        :key truncation: sets the number of previous LBFGS iterations to keep
        :type truncation : `int`
        :default truncation: 30
        :key restart: restart after this many iteration steps.
        :type restart: `int`
        :default restart: 60
        :key iterMax: maximium number of iterations.
        :type iterMax: `int`
        :default iterMax: 300
        :key relAlphaMin: minimal step size relative to serach direction.
                          The value should be chosen such that
                          At any iteration step `F(m + alpha * p)` is just discriminable from
                          `F(m)` for any `alpha > relAlphaMin * |m|/|p|'.
        :type relAlphaMin: ``float``
        :default relAlphaMin: 1e-8
        :key initialAlpha: initial step size alpha in line serach. Typically alpha=1 is a good initial value
                            but a larger or smaller initial value may help to get the iteration started
                            when only an approximation of the Hessian is available.
        :type initialAlpha: ``float``
        :default initialAlpha: 1.
        :key scaleSearchDirection: if set the search direction is rescaled using an estimation of the norm of the Hessian
        :type scaleSearchDirection: ``bool``
        :default scaleSearchDirection: True
        
            Example of usage::
              cf=DerivedCostFunction()
              solver=MinimizerLBFGS(J=cf, m_tol = 1e-5, grad_tol = 1e-5, iterMax=300)
              solver.setOptions(truncation=20)
              solver.getLineSearch().setOptions(zoom_relChangeMin =1e-7)
              solver.run(initial_m)
              result=solver.getResult()
        """
        self.logger.debug("Setting options: %s" % (str(opts)))
        for o in opts:
            if o == "iterMax":
                self._iterMax = max(1, int(opts[o]))
            elif o == "m_tol":
                self._m_tol = max(float(opts[o]), EPSILON)
            elif o == "grad_tol":
                if opts[o] is None:
                    self._grad_tol = None
                else:
                    self._grad_tol = max(float(opts[o]), EPSILON)
            elif o == 'historySize' or o == 'truncation':
                assert opts[o] > 2, "Trancation must be greater than 2."
                self._truncation = max(0, int(opts[o]))
            elif o == 'restart':
                self._restart = max(int(opts[o]), 1)
            elif o == 'relAlphaMin':
                self._relAlphaMin = max(float(opts[o]), EPSILON)
            elif o == 'initialAlpha':
                self._initialAlpha = max(float(opts[o]), EPSILON)
            elif o == 'scaleSearchDirection':
                self._scaleSearchDirection = opts[o]
            else:
                raise KeyError("invalid option '%s'" % o)

    def setCallback(self, callback):
        """
        Sets a callback function to be called after every iteration.

            def callback(iterCount, m, norm_m, dm, Fm, gradFm, norm_grad_Fm, failed, args)

            with iteration count `iterCount`, current approximation `m`, last update `dm`, costfunction value `Fm`,
            gradient `gradFm`, see method ``AbstractMinimizer.doCallback``

        """
        if callback is not None and not callable(callback):
            raise TypeError("Callback function not callable.")
        self._callback = callback

    def doCallback(self, **args):
        """
        The callback function is called with the following arguments:

        :key iterCount: iteration count
        :type iterCount: ``int``
        :key m: current solution
        :type m: m-type (see ``CostFunction``)
        :key dm: last solution incerement
        :type dm: m-type (see ``CostFunction``) or None if iterCount==0
        :norm_m: norm of current solution `m`
        :type norm_m: ``float``
        :key Fm: value of costs function for `m`
        :type Fm: ``float``
        :key gradFm: gradient for `m`
        :type gradFm: g-type (see ``CostFunction``)
        :param args_m: arguments for `m`
        :type args_m: ``tuple``
        :param failed: set if the step was unsuccessful.
        :type failed: ``bool``
        """
        if self._callback is not None:
            self.logger.debug("Callback called.")
            self._callback(**args)

    def getResult(self):
        """
        Returns the result of the minimization.
        :rtype: m-type
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
        :param m: initial guess
        :type m: m-type
        :return: solution
        :rtype: m-type
        """
        assert self._m_tol > 0.
        if not self._grad_tol is None:
            assert self._grad_tol > 0.
        assert self._iterMax > 1
        assert self._relAlphaMin > 0.
        assert self._truncation > 0
        assert self._restart > 0
        # start the iteration:
        iterCount = 0
        alpha = self._initialAlpha
        H_scale = None
        if self._restart < self._iterMax+2:
            self._restart = self._iterMax * 3
            self.logger.info("MinimizerLBFGS: Restart is currently disabled.")
        self._result = m
        args_m = self.getCostFunction().getArgumentsAndCount(m)
        grad_Fm = self.getCostFunction().getGradientAndCount(m, *args_m)
        norm_m = self.getCostFunction().getNormAndCount(m)
        Fm = self.getCostFunction().getValueAndCount(m, *args_m)
        #Fm_old = Fm
        self.logger.info("Initialization completed.")

        self.doCallback(iterCount=0, m=m, dm=None, Fm=Fm, grad_Fm=grad_Fm,
                        norm_m=norm_m, args_m=args_m, failed=False)

        non_curable_break_down = False
        converged = False
        while not converged and not non_curable_break_down and iterCount < self._iterMax:
            k = 0
            break_down = False
            s_and_y = []
            self.__initializeHessian=True
            iterCount_last_break_down = -1

            while not converged and not break_down and k < self._restart and iterCount < self._iterMax:
                self.logger.info("********** iteration %3d **********" % iterCount)
                self.logger.info("\tF(m) = %g" % Fm)
                # determine search direction
                p = -self._twoLoop(H_scale, grad_Fm, s_and_y, m, args_m)
                # Now we call the line search with F(m+alpha*p)
                # at this point we know that grad F(m) is not zero?
                phi=EvalutedPhi(CostFunction1DEvaluationFactory(m, p, costfunction=self.getCostFunction()),
                                alpha=0., args=args_m, valF=Fm, gradF=grad_Fm)
                if norm_m > 0:
                    alphaMin = self._relAlphaMin * norm_m/self.getCostFunction().getNormAndCount(p)
                else:
                    alphaMin = self._relAlphaMin
                self.getLineSearch().setOptions(alphaMin=alphaMin)
                alpha = max(alpha, alphaMin * 1.10)
                alpha_s = self.getCostFunction().getSqueezeFactor(m, p)
                if not alpha_s is None:
                    assert alpha_s >0
                    self.logger.debug("Safe alpha value given as %g." % (alpha_s,))
                    alpha_s*=0.5
                    if alpha_s > alphaMin:
                        alpha = min(alpha, alpha_s)
                self.logger.info("Starting line search with alpha  = %g."% (alpha, ))
                try:
                    phi_new = self.getLineSearch().run(phi, alpha)
                    alpha = phi_new.alpha
                    Fm_new = phi_new.valF
                    grad_Fm_new = phi_new.gradF
                    args_m_new = phi_new.args
                except LineSearchTerminationError as e:
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
                norm_m_new = self.getCostFunction().getNormAndCount(m_new)
                norm_dm = self.getCostFunction().getNormAndCount(delta_m)

                self._result = m_new
                converged = False

                mtol_abs = norm_m_new * self._m_tol
                flag = norm_dm <= mtol_abs
                if flag:
                    self.logger.info("F(m) = %g" % Fm_new)
                    self.logger.info("Solution has converged: |m-m_old|=%g < |m|*m_tol=%g" % (norm_dm, mtol_abs))
                    converged = True
                    break
                else:
                    self.logger.info("Solution checked: |m-m_old|=%g, |m|*m_tol=%g" % (norm_dm, mtol_abs))
                # unfortunately there is more work to do!
                if grad_Fm_new is None:
                    self.logger.debug("Calculating missing gradient.")
                    args_new = self.getCostFunction().getArgumentsAndCount(m_new)
                    grad_Fm_new = self.getCostFunction().getGradientAndCount(m_new, *args_new)

                if not self._grad_tol is None:
                    Ftol_abs = self._grad_tol * abs(max(abs(Fm), abs(Fm_new)))
                    dFm = abs(Fm - Fm_new)
                    flag = dFm <= Ftol_abs
                    if flag:
                        converged = True
                        self.logger.info("F(m) = %g" % Fm_new)
                        self.logger.info("Gradient has converged: |F-Fold|=%g < g_tol*max(|F|,|Fold|)=%g" % (dFm, Ftol_abs))
                        break
                    else:
                        self.logger.info("Gradient checked: |F-Fold|=%g, g_tol*max(|F|,|Fold|)=%g" % (dFm, Ftol_abs))

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
                args_m = args_m_new
                norm_m = norm_m_new
                k += 1
                iterCount += 1
                self.doCallback(iterCount=iterCount, m=m, dm=delta_m, Fm=Fm, grad_Fm=grad_Fm,
                                norm_m=norm_m, args_m=args_m, failed=break_down)

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
                    H_scale = rho/norm_dm ** 2
                    self.logger.info("Scale of Hessian = %g." % H_scale)
            # case handling for inner iteration:
            if break_down:
                if not iterCount > iterCount_last_break_down:
                    non_curable_break_down = True
                    self.logger.warning(">>>>> Incurable break down detected in step %d." % iterCount)
                else:
                    iterCount_last_break_down = iterCount
                    self.logger.debug("Break down detected in step %d. Iteration is restarted." % iterCount)
            if not k < self._restart:
                self.logger.debug("Iteration is restarted after %d steps." % iterCount)

        # case handling for inner iteration:
        if iterCount >= self._iterMax:
            self.logger.warning(">>>>>>>>>> Maximum number of iterations reached! <<<<<<<<<<")
            raise MinimizerMaxIterReached("Gave up after %d steps." % iterCount)
        elif non_curable_break_down:
            self.logger.warning(">>>>>>>>>> Incurable breakdown! <<<<<<<<<<")
            raise MinimizerIterationIncurableBreakDown("Gave up after %d steps." % iterCount)
        self.logger.info("Success after %d iterations!" % iterCount)
        return self._result

    def _twoLoop(self, H_scale, grad_Fm, s_and_y, m, args_m):
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

        if self.__initializeHessian:
            self.logger.debug("Hessian is expected to be updated.")
        p = self.getCostFunction().getInverseHessianApproximationAndCount(q, m,
                                                                            initializeHessian=self.__initializeHessian,
                                                                            *args_m)
        self.__initializeHessian = False
        norm_p = self.getCostFunction().getNormAndCount(p)
        if not norm_p > 0:
            raise MinimizerException("Approximate Hessian inverse returns zero.")
        # this is if one wants
        if H_scale is not None and self._scaleSearchDirection:
            p_dot_q = self.getCostFunction().getDualProductAndCount(p, q)
            if abs(p_dot_q) > 0:
                #  we  would like to see that  h * norm(p)^2 = < H p, p> = < q, p> with h=norm(H)=H_scale
                # so we rescale p-> a*p ; h * a^2 * norm(p)^2 = a * < q, p> -> a = <q,p>/h/norm(p)**2
                a = abs(p_dot_q) / H_scale / norm_p ** 2
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
        # self.lslogger.setLevel(logging.DEBUG)
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
            F = CostFunction1DEvaluationFactory(x, d, costfunction=self.getCostFunction())
            alpha, Fm, Fm_info, args_new = self.getLineSearch(F, x, d, -r, Fm, args, alpha=alpha)
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

