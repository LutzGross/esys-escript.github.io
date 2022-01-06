
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

"""Generic minimization algorithms"""

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['MinimizerException', 'MinimizerIterationIncurableBreakDown',\
           'MinimizerMaxIterReached' , 'AbstractMinimizer', 'MinimizerLBFGS', 'MinimizerNLCG']

import logging
import numpy as np

try:
    from esys.escript import Lsup, sqrt, EPSILON, getMPIRankWorld
except:
    Lsup=lambda x: np.amax(abs(x))
    sqrt=np.sqrt
    EPSILON=1e-18

lslogger=logging.getLogger('inv.minimizer.linesearch')
zoomlogger=logging.getLogger('inv.minimizer.linesearch.zoom')

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


def _zoom(phi, gradphi, phiargs, alpha_lo, alpha_hi, phi_lo, phi_hi, c1, c2,
          phi0, gphi0, interpolationOrder, tol_df, tol_sm, IMAX=25):
    """
    Helper function for `line_search` below which tries to tighten the range
    alpha_lo...alpha_hi. See Chapter 3 of 'Numerical Optimization' by
    J. Nocedal for an explanation.
    interpolation options are linear, quadratic, cubic or Newton interpolation.
    """

    def linearinterpolate(alpha_lo,alpha_hi,old_alpha):
        alpha = 0.5*(alpha_lo+alpha_hi)
        return alpha, 1

    def quadinterpolate(alpha_lo,alpha_hi,old_alpha,old_phi):
        if old_alpha is None:
            return linearinterpolate(alpha_lo,alpha_hi,old_alpha)
        denom=2.0*(old_phi-phi0-gphi0*old_alpha)
        if denom == 0:
            return linearinterpolate(alpha_lo,alpha_hi,old_alpha)
        alpha=-gphi0*old_alpha*old_alpha/denom
        if np.abs(alpha-old_alpha) < tol_df or np.abs(alpha) < tol_sm*np.abs(old_alpha):
            alpha=0.5*old_alpha
        if abs(alpha) < 1e-8:
            return linearinterpolate(alpha_lo,alpha_hi,old_alpha)
        return alpha, 2

    def cubicinterpolate(alpha_lo,alpha_hi,old_alpha,old_phi,very_old_alpha,very_old_phi):
        if very_old_alpha is None:
            return quadinterpolate(alpha_lo,alpha_hi,old_alpha,old_phi)
        a0s,a1s=very_old_alpha*very_old_alpha,old_alpha*old_alpha
        denom=a0s*a1s*(old_alpha-very_old_alpha)
        if denom==0:
            return quadinterpolate(alpha_lo,alpha_hi,old_alpha,old_phi)
        a0c,a1c=a0s*very_old_alpha,a1s*old_alpha
        tmpA=old_phi-phi0-gphi0*old_alpha 
        tmpB=very_old_phi-phi0-gphi0*very_old_alpha
        a=( a0s*tmpA-a1s*tmpB)/denom
        b=(-a0c*tmpA+a1c*tmpB)/denom
        deter=b*b-3.0*a*gphi0
        if deter<0:
            return quadinterpolate(alpha_lo,alpha_hi,old_alpha,old_phi)
        alpha=(-b+sqrt(deter))/(3.0*a)
        if np.abs(alpha-old_alpha) < tol_df or np.abs(alpha) < tol_sm*np.abs(old_alpha):
            alpha=0.5*old_alpha
        if abs(alpha) < 1e-8:
            return quadinterpolate(alpha_lo,alpha_hi,old_alpha,old_phi)
        return alpha, 3

    def newtoninterpolate(alpha_data,phi_data,old_alpha,old_phi):
        # Interpolates using a polynomial of increasing order
        # The coefficients of the interpolated polynomial are found using the Newton method
        if very_old_alpha is None:
            return quadinterpolate(alpha_lo,alpha_hi,old_alpha,old_phi)
        if old_alpha in alpha_data:
            return quadinterpolate(alpha_lo,alpha_hi,old_alpha,old_phi)
        alpha_data.append(old_alpha)
        phi_data.append(old_phi)
        coefs=newton_poly(alpha_data,phi_data)
        alpha, ii=newtonroot(coefs,alpha_data,old_alpha)
        if alpha < alpha_lo or alpha > alpha_hi: # Root is outside the domain
            if getMPIRankWorld() == 0:
                zoomlogger.debug("NOTE: Newton interpolation converged on a root outside of [alpha_lo,alpha_hi]. Falling back on cubic interpolation")
            return cubicinterpolate(alpha_lo,alpha_hi,old_alpha,old_phi,very_old_alpha,very_old_phi)
        if np.abs(alpha-old_alpha) < tol_df or np.abs(alpha) < tol_sm*np.abs(old_alpha):
            alpha=0.5*old_alpha
        if abs(alpha) <= 0:
            if getMPIRankWorld() == 0:
                zoomlogger.debug("NOTE: Newton interpolation returned alpha <= 0. Falling back on cubic interpolation")
            return cubicinterpolate(alpha_lo,alpha_hi,old_alpha,old_phi,very_old_alpha,very_old_phi)
        return alpha, len(alpha_data)

    def newton_poly(alpha_data,phi_data):
        # Returns the coefficients of the newton form polynomial of phi(alpha)
        # for the points alpha_data and function values phi_data
        m=len(alpha_data)
        x=np.copy(alpha_data)
        a=np.copy(phi_data)
        for k in range(1,m):
            a[k:m]=(a[k:m]-a[k-1])/(x[k:m]-x[k-1])
        return a

    def newtonroot(coefs,alpha_data,startingGuess):
        # Solves for the root of a polynomial using the newton method
        dfcoefs=list(range(1,len(coefs)))
        for i in range(0,len(coefs)-1):
            dfcoefs[i]*=coefs[i+1]
        tol=1e-6
        error=100
        maxiterations=100
        point=startingGuess
        counter=0
        for i in range(0,maxiterations):
            counter+=1
            product=1
            numer=coefs[0]
            denom=dfcoefs[0]
            for i in range(1, len(coefs)):
                product*=(point-alpha_data[i-1])
                numer+=product*coefs[i]
                if i < len(coefs)-1:
                    denom+=product*dfcoefs[i]
            if denom == 0:
                if getMPIRankWorld() == 0:
                    zoomlogger.debug("NOTE: The Newton solver failed (denom==0). Falling back on cubic interpolation")                
                return cubicinterpolate(alpha_lo,alpha_hi,old_alpha,old_phi,very_old_alpha,very_old_phi)
            newpoint=float(point-numer/denom)
            error=abs(newpoint-point)
            if error < tol:
                if getMPIRankWorld() == 0:
                    zoomlogger.debug("Newton interpolation converged in %d iterations. Got alpha = %f" % (counter, newpoint))
                return newpoint, counter
            else:
                point=newpoint
        # zoomlogger.debug("NOTE: Newton interpolation failed to converge (exceeded max iterations). Falling back on cubic interpolation" % error)
        return cubicinterpolate(alpha_lo,alpha_hi,old_alpha,old_phi,very_old_alpha,very_old_phi)

    i=0

    # Setup
    old_alpha=None
    old_phi=None
    very_old_alpha=None
    very_old_phi=None
    alpha_data=[0.0]
    phi_data=[phi0]

    while i<=IMAX:
        if interpolationOrder == 1:
            alpha, o=linearinterpolate(alpha_lo,alpha_hi,old_alpha)
        elif interpolationOrder == 2:
            alpha, o=quadinterpolate(alpha_lo,alpha_hi,old_alpha,old_phi)
        elif interpolationOrder == 3:
            alpha, o=cubicinterpolate(alpha_lo,alpha_hi,old_alpha,old_phi,very_old_alpha,very_old_phi)
        elif interpolationOrder == 4:
            alpha, o=newtoninterpolate(alpha_data,phi_data,old_alpha,old_phi)
        else:
            raise TypeError("Invalid interpolation order")

        phiargs(alpha)
        phi_a=phi(alpha)
        if getMPIRankWorld() == 0:
            zoomlogger.debug("iteration %d, alpha=%g, phi(alpha)=%g (interpolation order = %d)"%(i,alpha,phi_a, o))
        if phi_a > phi0+c1*alpha*gphi0 or phi_a >= phi_lo:
            very_old_alpha,very_old_phi=old_alpha,old_phi
            old_alpha,old_phi=alpha_hi,phi_hi
            alpha_hi,phi_hi=alpha,phi_a
        else:
            gphi_a=gradphi(alpha)
            if getMPIRankWorld() == 0:
                zoomlogger.debug("\tphi'(alpha)=%e"%(gphi_a))
            if np.abs(gphi_a) <= -c2*gphi0:
                break
            if gphi_a*(alpha_hi-alpha_lo) >= 0:
                alpha_hi = alpha_lo
            very_old_alpha,very_old_phi=old_alpha,old_phi
            old_alpha,old_phi=alpha_lo,phi_lo
            alpha_lo,phi_lo=alpha,phi_a
        if not alpha_hi > alpha_lo:
            break
        i+=1
    return alpha, phi_a

def line_search(f, x, p, g_Jx, Jx, args_x, alpha=1.0, alpha_truncationax=50.0,
                c1=1e-4, c2=0.9, interpolationOrder=1, tol_df=1e-6, tol_sm=1e-5, IMAX=15, IMAX_ZOOM=25):
    """
    Line search method that satisfies the strong Wolfe conditions.
    See Chapter 3 of 'Numerical Optimization' by J. Nocedal for an explanation.

    :param f: callable objective function f(x)
    :param x: start value for the line search
    :param p: search direction
    :param g_Jx: value for the gradient of f at x
    :param Jx: value of f(x)
    :param args_x: arguments for x
    :param alpha: initial step length. If g_Jx is properly scaled alpha=1 is a
                  reasonable starting value.
    :param alpha_truncationax: algorithm terminates if alpha reaches this value
    :param c1: value for Armijo condition (see reference)
    :param c2: value for curvature condition (see reference)
    :param interpolationOrder: the order of the interpolation used (1, 2 or 3)
    :param tol_df: the interpolated value of alpha, alpha_i, must differ from the previous
    :             value by at least this much i.e. abs(alpha_i-alpha_{i-1}) < tol_df
    :param tol_sm: the interpolated value of alpha must not be less than tol_sm
    :             i.e. abs(alpha_i) < tol_sm*abs(alpha_{i-1})
    :param IMAX: maximum number of iterations to perform
    """
    # this stores the latest gradf(x+a*p) which is returned
    if args_x is None:
        args_x=f.getArguments(x)
        Jx=None
        g_Jx=None
    if Jx is None:
        Jx=f(x, *args_x)
        if getMPIRankWorld() == 0:
            lslogger.debug("initial J(x) calculated= %s"%str(Jx))
    if g_Jx is None:
        g_Jx=f.getGradient(x, *args_x)
        if getMPIRankWorld() == 0:
            lslogger.debug("initial grad J(x) calculated.")
    # this stores the latest gradf(x+a*p) which is returned
    g_Jx_new=[args_x, Jx, g_Jx]
    def phi(a):
        """ phi(a):=f(x+a*p) """
        if g_Jx_new[0] is None:
            phiargs(a)
        g_Jx_new[1]=f(x+a*p, *g_Jx_new[0])
        return g_Jx_new[1]
    def gradphi(a):
        if g_Jx_new[0] is None:
            phiargs(a)
        g_Jx_new[2]=f.getGradient(x+a*p, *g_Jx_new[0])
        if getMPIRankWorld() == 0:    
            lslogger.debug("grad J(x+alpha*p) calculated.")
        return f.getDualProduct(p, g_Jx_new[2])
    
    def phiargs(a):
        try:
            args=f.getArguments(x+a*p)
        except:
            args=()
        g_Jx_new[0]=args
        g_Jx_new[1]=None
        g_Jx_new[2]=None
        return args
    
    old_alpha=0.
    phi0=Jx
    if getMPIRankWorld() == 0:    
        lslogger.debug("phi(0)=%e"%(phi0))
    gphi0=f.getDualProduct(p, g_Jx) #gradphi(0., *args0)
    if getMPIRankWorld() == 0:
        lslogger.debug("phi'(0)=%e"%(gphi0))
    old_phi_a=phi0
    phi_a=phi0
    i=1
    if (interpolationOrder != 1) and (interpolationOrder != 2) and (interpolationOrder != 3) and (interpolationOrder != 4):
        if getMPIRankWorld() == 0:
            lslogger.debug("WARNING invalid interpolation order. Setting interpolationOrder = 1")
        interpolationOrder=1
    #alpha=2.*alpha
    while i< max(IMAX,2) and alpha>0.:
        phiargs(alpha)
        phi_a=phi(alpha)
        if getMPIRankWorld() == 0:
            lslogger.debug("iteration %d, alpha=%e, phi(alpha)=%e"%(i,alpha,phi_a))
        if (phi_a > phi0+c1*alpha*gphi0) or ((phi_a>=old_phi_a) and (i>1)):
            alpha, phi_a = _zoom(phi, gradphi, phiargs, old_alpha, alpha, old_phi_a, phi_a, c1, c2, phi0, gphi0, interpolationOrder, tol_df, tol_sm, IMAX=IMAX_ZOOM)
            break

        gphi_a=gradphi(alpha)
        if getMPIRankWorld() == 0:
            lslogger.debug("phi'(alpha)=%e"%(gphi_a))
        if np.abs(gphi_a) <= -c2*gphi0:
            break
        if gphi_a >= 0:
            alpha, phi_a = _zoom(phi, gradphi, phiargs, alpha, old_alpha, phi_a, old_phi_a, c1, c2, phi0, gphi0, interpolationOrder, tol_df, tol_sm, IMAX=IMAX_ZOOM)
            break

        old_alpha=alpha
        # the factor is arbitrary as long as there is sufficient increase
        alpha=2.*alpha
        old_phi_a=phi_a
        if alpha > alpha_truncationax:
            break
        
        i+=1
    return alpha, phi_a, g_Jx_new[2], g_Jx_new[0] # returns for x+alpha*p (g_Jx_new[2] can be None)


##############################################################################
class AbstractMinimizer(object):
    """
    Base class for function minimization methods.   
    """

    def __init__(self, J=None, m_tol=1e-4, J_tol=None, imax=300):
        """
        Initializes a new minimizer for a given cost function.

        :param J: the cost function to be minimized
        :type J: `CostFunction`
        :param m_tol: terminate interations when relative change of the level set 
                      function is less than or equal m_tol
        :type m_tol: float
 
        """
        self.setCostFunction(J)
        self.setMaxIterations(imax)
        self._result = None
        self._callback = None
        self.logger = logging.getLogger('inv.minimizer.%s'%self.__class__.__name__)
        self.setTolerance(m_tol=m_tol, J_tol=J_tol)
        self.interpolationOrder=1

    def setCostFunction(self, J):
        """
        set the cost function to be minimized

        :param J: the cost function to be minimized
        :type J: `CostFunction`
        """
        self.__J=J

    def getCostFunction(self):
        """
        return the cost function to be minimized

        :rtype: `CostFunction`
        """
        return self.__J

    def setTolerance(self, m_tol=1e-4, J_tol=None):
        """
        Sets the tolerance for the stopping criterion. The minimizer stops
        when an appropriate norm is less than `m_tol`.
        """
        self._m_tol = m_tol
        self._J_tol = J_tol

    def setMaxIterations(self, imax):
        """
        Sets the maximum number of iterations before the minimizer terminates.
        """
        self._imax = imax

    def setCallback(self, callback):
        """
        Sets a callback function to be called after every iteration.
        It is up to the specific implementation what arguments are passed
        to the callback. Subclasses should at least pass the current
        iteration number k, the current estimate x, and possibly f(x),
        grad f(x), and the current error.
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

    def getOptions(self):
        """
        Returns a dictionary of minimizer-specific options.
        """
        return {}

    def setOptions(self, **opts):
        """
        Sets minimizer-specific options. For a list of possible options see
        `getOptions()`.
        """
        raise NotImplementedError

    def run(self, x0):
        """
        Executes the minimization algorithm for *f* starting with the initial
        guess ``x0``.

        :return: the result of the minimization
        """
        raise NotImplementedError

    def logSummary(self):
        """
        Outputs a summary of the completed minimization process to the logger.
        """
        if hasattr(self.getCostFunction(), "Value_calls") and getMPIRankWorld() == 0:
            self.logger.info("Number of cost function evaluations: %d"%self.getCostFunction().Value_calls)
            self.logger.info("Number of gradient evaluations: %d"%self.getCostFunction().Gradient_calls)
            self.logger.info("Number of inner product evaluations: %d"%self.getCostFunction().DualProduct_calls)
            self.logger.info("Number of argument evaluations: %d"%self.getCostFunction().Arguments_calls)
            self.logger.info("Number of norm evaluations: %d"%self.getCostFunction().Norm_calls)

##############################################################################
class MinimizerLBFGS(AbstractMinimizer):
    """
    Minimizer that uses the limited-memory Broyden-Fletcher-Goldfarb-Shanno
    method.
    See Chapter 6 of 'Numerical Optimization' by J. Nocedal for an explanation.
    """ 

    # History size
    _truncation = 30

    # Initial Hessian multiplier
    _initial_H = 1.

    # Restart after this many iteration steps
    _restart = 60

    # maximum number of line search steps
    _max_linesearch_steps = 25

    # maximum number of zoom steps in line search
    _max_zoom_steps = 50

         
    # tolerance for the interpolated value of alpha in line search must not be less than tol_sm
    # i.e. abs(alpha_i) < tol_sm*abs(alpha_{i-1})
    _tol_sm = 1e-5
    
    # the interpolated value of alpha, alpha_i in line search, must differ from the previous 
    # value by at least this much i.e. abs(alpha_i-alpha_{i-1}) < tol_df    
    _tol_df = 1e-6
    
    def getOptions(self):
        """
        returns a dictionary of LBFGS options 
        rtype: dictionary with keys 'truncation', 'initialHessian', 'restart', 'max_linesearch_steps', 'max_zoom_steps'
        'interpolationOrder', 'tol_df', 'tol_sm' 
        """
        return {'truncation':self._truncation,'initialHessian':self._initial_H, 'restart':self._restart,
                 'max_linesearch_steps' : self._max_linesearch_steps, 'max_zoom_steps' : self._max_zoom_steps, 
                'interpolationOrder': self.interpolationOrder,'tol_df': self._tol_df, 'tol.sm': self._tol_sm}



    def setOptions(self, **opts):
        """
        setOptions for LBFGS.  use       solver.setOptions( key = value)
        :key truncation: sets the number of previous LBFGS iterations to keep
        :type truncation : int
        :default truncation: 30
        :key initialHessian: 1 if initial Hessian is given
        :type initialHessian: 1 or 0    
        :default initialHessian: 1       
        :key restart: restart after this many iteration steps
        :type restart: int
        :default restart: 60
        :key max_linesearch_steps: maximum number of line search iterations
        :type max_linesearch_steps: int
        :default max_linesearch_steps: 25  
        :key max_zoom_steps: maximum number of zoom iterations
        :type max_zoom_steps: int
        :default max_zoom_steps: 60
        :key interpolationOrder: order of the interpolation used for line search 
        :type interpolationOrder: 1,2,3 or 4
        :default interpolationOrder: 1      
        :key tol_df: interpolated value of alpha, alpha_i, must differ
                   : from the previous value by at least this much 
                   : abs(alpha_i-alpha_{i-1}) < tol_df  (line search)
        :type tol_df: float
        :default tol_df: 1e-6
        :key tol_sm: interpolated value of alpha must not be less than tol_sm
                   : abs(alpha_i) < tol_sm*abs(alpha_{i-1})
        :type tol_sm: float
        :default tol_sm: 1e-5     
        
            Example of usage::
              cf=DerivedCostFunction()
              solver=MinimizerLBFGS(J=cf, m_tol = 1e-5, J_tol = 1e-5, imax=300)
              solver.setOptions(truncation=20, tol_df =1e-7)
              solver.run(initial_m)
              result=solver.getResult()
        """
        if getMPIRankWorld() == 0:
            self.logger.debug("Setting options: %s"%(str(opts)))
        for o in opts:
            if o=='historySize' or o=='truncation':
                assert opts[o]>2, "Trancation must be greater than 2."
                self._truncation=opts[o]
            elif o=='initialHessian':
                self._initial_H=opts[o]
            elif o=='restart':
                self._restart=opts[o]
            elif o=='max_linesearch_steps':
                self._max_linesearch_steps=opts[o]
            elif o=='max_zoom_steps':
                self._max_zoom_steps=opts[o]
            elif o=='interpolationOrder':
                self.interpolationOrder=opts[o]
            elif o=='tol_df':
                self.tol_df=opts[o]
            elif o=='tol_sm':
                self.tol_sm=opts[o]
            else:
                raise KeyError("Invalid option '%s'"%o)

    def run(self, x):
        """
        The callback function is called with the following arguments:
            k       - iteration number
            x       - current estimate
            Jx      - value of cost function at x
            g_Jx    - gradient of cost function at x
            norm_dJ - ||Jx_k - Jx_{k-1}|| (only if J_tol is set)
            norm_dx - ||x_k - x_{k-1}|| (only if m_tol is set)

        :param x: Level set function representing our initial guess
        :type x: `Data`
        :return: Level set function representing the solution
        :rtype: `Data`
        """
        #if self.getCostFunction().provides_inverse_Hessian_approximation:
        #    self.getCostFunction().updateHessian()
        #    invH_scale = None
        #else:
        #    invH_scale = self._initial_H
        if self._initial_H:
            invH_scale = 1./self._initial_H
        else:
            invH_scale=None
        # start the iteration:
        n_iter = 0
        n_last_break_down=-1
        alpha=1.
        non_curable_break_down = False
        converged = False
        args=self.getCostFunction().getArguments(x)
        g_Jx=self.getCostFunction().getGradient(x, *args)

        if getMPIRankWorld() == 0:
            self.logger.debug("initial grad J(x) calculated.")
        # equivalent to getValue() for Downunder CostFunctions
        Jx=self.getCostFunction()(x, *args)
        Jx_0=Jx

        cbargs = {'k':n_iter, 'x':x, 'Jx':Jx, 'g_Jx':g_Jx}
        if self._J_tol:
            cbargs.update(norm_dJ=None)
        if self._m_tol:
            cbargs.update(norm_dx=None)
        self._doCallback(**cbargs)

        while not converged and not non_curable_break_down and n_iter < self._imax:
          k=0
          break_down = False
          s_and_y=[]
          # initial step length for line search

          while not converged and not break_down and k < self._restart and n_iter < self._imax:
                #self.logger.info("\033[1;31miteration %d\033[1;30m"%n_iter)
                if getMPIRankWorld() == 0:
                    self.logger.info("********** iteration %3d **********"%n_iter)
                    self.logger.info("\tJ(x) = %s"%Jx)
                #self.logger.debug("\tgrad f(x) = %s"%g_Jx)
                if invH_scale and getMPIRankWorld() == 0:
                    self.logger.debug("\tH = %s"%invH_scale)

                # determine search direction
                p = -self._twoLoop(invH_scale, g_Jx, s_and_y, x, *args)

                # determine new step length using the last one as initial value
                # however, avoid using too small steps for too long.
                # FIXME: This is a bit specific to esys.downunder in that the
                # inverse Hessian approximation is not scaled properly (only
                # the regularization term is used at the moment)...
                if invH_scale is None:
                    #if n_iter >1:
                    #if alpha <= 0.5:
                    #alpha=min(2*alpha,1.0)
                    alpha=alpha
                    #else:
                    #alpha=1.
                else:
                    # reset alpha for the case that the cost function does not
                    # provide an approximation of inverse H
                    alpha=1.0*alpha
                alpha, Jx_new, g_Jx_new, args_new = line_search(self.getCostFunction(), x, p, g_Jx, Jx, args, alpha, interpolationOrder=self.interpolationOrder, IMAX=self._max_linesearch_steps, IMAX_ZOOM=self._max_zoom_steps)
                # this function returns a scaling alpha for the search
                # direction as well as the cost function evaluation and
                # gradient for the new solution approximation x_new=x+alpha*p
                if getMPIRankWorld() == 0:
                    self.logger.debug("Search direction scaling alpha=%e"%alpha)

                # execute the step
                delta_x = alpha*p
                x_new = x + delta_x

                converged = True
                if self._J_tol:
                    dJ = abs(Jx_new-Jx)
                    JJtol = self._J_tol * abs(Jx_new-Jx_0)
                    flag = dJ <= JJtol
                    if self.logger.isEnabledFor(logging.DEBUG) and getMPIRankWorld() == 0:
                        if flag:
                            self.logger.debug("Cost function has converged: dJ=%e, J*J_tol=%e"%(dJ,JJtol))
                        else:
                            self.logger.debug("Cost function checked: dJ=%e, J*J_tol=%e"%(dJ,JJtol))
                    cbargs.update(norm_dJ=dJ)
                    converged = converged and flag

                if self._m_tol:
                    norm_x = self.getCostFunction().getNorm(x_new)
                    norm_dx = self.getCostFunction().getNorm(delta_x)
                    flag = norm_dx <= self._m_tol * norm_x
                    if self.logger.isEnabledFor(logging.DEBUG) and getMPIRankWorld() == 0:
                        if flag:
                            self.logger.debug("Solution has converged: dx=%e, x*m_tol=%e"%(norm_dx, norm_x*self._m_tol))
                        else:
                            self.logger.debug("Solution checked: dx=%e, x*m_tol=%e"%(norm_dx, norm_x*self._m_tol))
                    cbargs.update(norm_dx=norm_dx)
                    converged = converged and flag

                
                if converged:
                    if getMPIRankWorld() == 0:
                        self.logger.info("********** iteration %3d **********"%(n_iter+1,))
                        self.logger.info("\tJ(x) = %s"%Jx_new)
                    break
                
                # unfortunately there is more work to do!
                if g_Jx_new is None:
                    if getMPIRankWorld() == 0:
                        self.logger.debug("Calculating missing gradient for x+alpha*p.")
                    args_new=self.getCostFunction().getArguments(x_new)
                    g_Jx_new=self.getCostFunction().getGradient(x_new, *args_new)
                    #self.logger.debug("grad J(x+alpha*p) = %s"%str(g_Jx_new))
                delta_g=g_Jx_new-g_Jx

                rho=self.getCostFunction().getDualProduct(delta_x, delta_g)
                if abs(rho)>0:
                    s_and_y.append((delta_x,delta_g, rho ))
                else:
                    break_down=True

                self.getCostFunction().updateHessian()
                x=x_new
                g_Jx=g_Jx_new
                Jx=Jx_new
                args=args_new

                k+=1
                n_iter+=1
                cbargs.update(k=n_iter, x=x, Jx=Jx, g_Jx=g_Jx)
                self._doCallback(**cbargs)

                # delete oldest vector pair
                if k>self._truncation: s_and_y.pop(0)
                
                if not break_down:
                    # an estimation of the inverse Hessian if it would be a scalar P=invH_scale:
                    # the idea is that 
                    #  
                    #      dx=P*dg
                    # 
                    #  if there is a inverse Hessian approximation provided we use 
                    #
                    #     |dx|^2 = P <dx, dg> ->   invH_scale =|dx|^2 / <dx, dg> 
                    
                    if self.getCostFunction().provides_inverse_Hessian_approximation:
                        # set the new scaling factor (approximation of inverse Hessian)
                        denom=abs(self.getCostFunction().getDualProduct(delta_x, delta_g))
                        if denom > 0:
                            invH_scale=self.getCostFunction().getNorm(delta_x)**2/denom
                        else:
                            if self._initial_H:
                                invH_scale=1./self._initial_H
                                if getMPIRankWorld() == 0:
                                    self.logger.debug("** Break down in H update. Resetting to initial value %s."%(1./self._initial_H))
                    else:
                        # if there no inverse Hessian approximation is provided
                        #
                        #     <dg, dx> = P <dg, dg> ->   invH_scale =<dg, dx> / <dg, dg> 
                        #
                        denom=self.getCostFunction().getDualProduct(delta_g, delta_g)
                        if denom > 0:
                            invH_scale=self.getCostFunction().getDualProduct(delta_x,delta_g)/denom
                        else:
                            invH_scale=1./self._initial_H
                            if getMPIRankWorld() == 0:
                                self.logger.debug("** Break down in H update. Resetting to initial value %s."%(1./self._initial_H))
          # case handling for inner iteration:
          if break_down:
              if n_iter == n_last_break_down+1:
                  non_curable_break_down = True
                  if getMPIRankWorld() == 0:
                    self.logger.debug("** Incurable break down detected in step %d."%n_iter)
              else:
                  n_last_break_down = n_iter
                  if getMPIRankWorld() == 0:
                    self.logger.debug("** Break down detected in step %d. Iteration is restarted."%n_iter)
          if not k < self._restart:
              if getMPIRankWorld() == 0:
                    self.logger.debug("Iteration is restarted after %d steps."%n_iter)

        # case handling for inner iteration:
        self._result=x
        if n_iter >= self._imax:
            if getMPIRankWorld() == 0:
                self.logger.warning(">>>>>>>>>> Maximum number of iterations reached! <<<<<<<<<<")
            raise MinimizerMaxIterReached("Gave up after %d steps."%n_iter)
        elif non_curable_break_down:
            if getMPIRankWorld() == 0:
                self.logger.warning(">>>>>>>>>> Incurable breakdown! <<<<<<<<<<")
            raise MinimizerIterationIncurableBreakDown("Gave up after %d steps."%n_iter)
        if getMPIRankWorld() == 0:
            self.logger.info("Success after %d iterations!"%n_iter)
        return self._result

    def _twoLoop(self, invH_scale, g_Jx, s_and_y, x, *args):
        """
        Helper for the L-BFGS method.
        See 'Numerical Optimization' by J. Nocedal for an explanation.
        """
        q=g_Jx
        alpha=[]
        for s,y, rho in reversed(s_and_y):
            a=self.getCostFunction().getDualProduct(s, q)/rho
            alpha.append(a)
            q=q-a*y
        
        if self.getCostFunction().provides_inverse_Hessian_approximation:
            r = self.getCostFunction().getInverseHessianApproximation(x, q, *args)
            #  we  would like to see that  | r | ^2 = P * < r, q> with P=invH_scale
            
            norm_r=self.getCostFunction().getNorm(r)
            if norm_r > 0:
                if invH_scale is not None:
                    rescale=invH_scale*abs(self.getCostFunction().getDualProduct(r, q))/norm_r**2
                    r*=rescale
                    if getMPIRankWorld() == 0:
                        self.logger.debug("rescale of search direction : r=%g (re-scale = %g, invH = %g)"%(norm_r*rescale, rescale, invH_scale))
            else:
                raise MinimizerException("approximate Hessian inverse returns zero.")
        else:
             # if there no inverse Hessian approximation is provided then set r = P * q with invH_scale=P
             r = invH_scale * q
                        
        for s,y,rho in s_and_y:
            beta = self.getCostFunction().getDualProduct(r, y)/rho
            a = alpha.pop()
            r = r + s * (a-beta)
        return r


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
            Jx    - value of cost function at x
            g_Jx  - gradient of cost function at x
            gnorm - norm of g_Jx (stopping criterion)
        """
        #self.logger.setLevel(logging.DEBUG)
        #lslogger.setLevel(logging.DEBUG)
        i=0
        k=0
        args=self.getCostFunction().getArguments(x)
        r=-self.getCostFunction().getGradient(x, *args)
        Jx=self.getCostFunction()(x, *args)
        d=r
        delta=self.getCostFunction().getDualProduct(r,r)
        delta0=delta
        gnorm0=Lsup(r)
        gnorm=gnorm0
        alpha=1
        self._doCallback(k=i, x=x, Jx=Jx, g_Jx=-r, gnorm=gnorm)

        while i < self._imax and gnorm > self._m_tol * gnorm0:
            self.logger.debug("iteration %d"%i)
            self.logger.debug("grad f(x) = %s"%(-r))
            self.logger.debug("d = %s"%d)
            self.logger.debug("x = %s"%x)
        
            alpha, Jx, g_Jx_new, args_new = line_search(self.getCostFunction(), x, d, -r, Jx, args, alpha=alpha, c2=0.4, interpolationOrder=self.interpolationOrder)
            self.logger.debug("alpha=%e"%(alpha))
            x=x+alpha*d
            args=args_new
            
            r=-self.getCostFunction().getGradient(x, *args) if g_Jx_new is None else -g_Jx_new
            # Fletcher-Reeves version
            delta_o=delta
            delta=self.getCostFunction().getDualProduct(r,r)
            beta=delta/delta_o
            d=r+beta*d
            k=k+1
            if self.getCostFunction().getDualProduct(r,d) <= 0:
                d=r
                k=0
            i+=1
            gnorm=Lsup(r)
            self._doCallback(k=i, x=x, Jx=Jx, g_Jx=g_Jx_new, gnorm=gnorm)
            self._result=x

        if i >= self._imax:
            if getMPIRankWorld() == 0:
                    self.logger.warning(">>>>>>>>>> Maximum number of iterations %s reached! <<<<<<<<<<"%i)
            raise MinimizerMaxIterReached("Gave up after %d steps."%i)

        if getMPIRankWorld() == 0:
           self.logger.info("Success after %d iterations! Initial/Final gradient=%e/%e"%(i,gnorm0, gnorm))
        return self._result


if __name__=="__main__":
    # Example usage with function 'rosen' (minimum=[1,1,...1]):
    from scipy.optimize import rosen, rosen_der
    from esys.downunder import MeteredCostFunction
    import sys
    N=4
    x0=np.array([4.]*N) # initial guess

    class RosenFunc(MeteredCostFunction):
        def __init__(self):
          super(RosenFunc, self).__init__()
          self.provides_inverse_Hessian_approximation=False
        def _getDualProduct(self, f0, f1):
            return np.dot(f0, f1)
        def _getValue(self, x, *args):
            return rosen(x)
        def _getGradient(self, x, *args):
            return rosen_der(x)
        def _getNorm(self,x):
            return Lsup(x)

    f=RosenFunc()
    m=None
    if len(sys.argv)>1:
        method=sys.argv[1].lower()
        if method=='nlcg':
            m=MinimizerNLCG(f)
        elif method=='bfgs':
            m=MinimizerBFGS(f)

    if m is None:
        # default
        m=MinimizerLBFGS(f)
        #m.setOptions(historySize=10000)

    logging.basicConfig(format='[%(funcName)s] \033[1;30m%(message)s\033[0m', level=logging.DEBUG)
    m.setTolerance(m_tol=1e-5)
    m.setMaxIterations(3000)
    m.run(x0)
    m.logSummary()
    print("\tLsup(result)=%.8f"%np.amax(abs(m.getResult())))

    #from scipy.optimize import fmin_cg
    #print("scipy ref=%.8f"%np.amax(abs(fmin_cg(rosen, x0, rosen_der, maxiter=10000))))

