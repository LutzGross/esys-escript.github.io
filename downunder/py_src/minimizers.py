
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

"""Generic minimization algorithms"""

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['AbstractMinimizer', 'MinimizerLBFGS', 'MinimizerBFGS', 'MinimizerNLCG']

import logging
import numpy as np

try:
    from esys.escript import Lsup, sqrt, EPSILON
except:
    Lsup=lambda x: np.amax(abs(x))
    sqrt=np.sqrt
    EPSILON=1e-18

import sys
if sys.version_info[0]>2:
    xrange=range
    
lslogger=logging.getLogger('inv.minimizer.linesearch')
zoomlogger=logging.getLogger('inv.minimizer.linesearch.zoom')

def _zoom(phi, gradphi, phiargs, alpha_lo, alpha_hi, phi_lo, phi_hi, c1, c2, phi0, gphi0, IMAX=25):
    """
    Helper function for `line_search` below which tries to tighten the range
    alpha_lo...alpha_hi. See Chapter 3 of 'Numerical Optimization' by
    J. Nocedal for an explanation.
    """
    i=0
    while True:
        alpha=alpha_lo+.5*(alpha_hi-alpha_lo) # should use interpolation...
        args_a=phiargs(alpha)
        phi_a=phi(alpha, *args_a)
        zoomlogger.debug("iteration %d, alpha=%e, phi(alpha)=%e"%(i,alpha,phi_a))
        if phi_a > phi0+c1*alpha*gphi0 or phi_a >= phi_lo:
            alpha_hi=alpha
        else:
            gphi_a=gradphi(alpha, *args_a)
            zoomlogger.debug("grad(phi(alpha))=%e"%(gphi_a))
            if np.abs(gphi_a) <= -c2*gphi0:
                break
            if gphi_a*(alpha_hi-alpha_lo) >= 0:
                alpha_hi = alpha_lo
            alpha_lo=alpha
            phi_lo=phi_a
        i+=1
        if i>IMAX:
            gphi_a=None
            break
    return alpha, phi_a, gphi_a

def line_search(f, x, p, gf, fx, alpha_max=50.0, c1=1e-4, c2=0.9, IMAX=15):
    """
    Line search method that satisfies the strong Wolfe conditions.
    See Chapter 3 of 'Numerical Optimization' by J. Nocedal for an explanation.

    :param f: callable objective function f(x)
    :param x: start value for the line search
    :param p: search direction
    :param gf: value for the gradient of f at x
    :param fx: value of f(x)
    :param alpha_max: algorithm terminates if alpha reaches this value
    :param c1: value for Armijo condition (see reference)
    :param c2: value for curvature condition (see reference)
    :param IMAX: maximum number of iterations to perform
    """
    # this stores the latest gradf(x+a*p) which is returned
    gf_new=[gf]

    def phi(a, *args):
        """ phi(a):=f(x+a*p) """
        return f(x+a*p, *args)
    def gradphi(a, *args):
        gf_new[0]=f.getGradient(x+a*p, *args)
        return f.getDualProduct(p, gf_new[0])
        #return f.getDirectionalDerivative(x+a*p, p, *args)
    def phiargs(a):
        try:
            args=f.getArguments(x+a*p)
        except:
            args=()
        return args

    old_alpha=0.
    # we assume gf is properly scaled so alpha=1 is a reasonable starting value
    alpha=1.
    if fx is None:
        args0=phiargs(0.)
        phi0=phi(0., *args0)
    else:
        phi0=fx
    lslogger.debug("phi(0)=%e"%(phi0))
    gphi0=f.getDualProduct(p, gf) #gradphi(0., *args0)
    lslogger.debug("grad phi(0)=%e"%(gphi0))
    old_phi_a=phi0
    i=1

    while i<IMAX and alpha>0. and alpha<alpha_max:

        args_a=phiargs(alpha)
        phi_a=phi(alpha, *args_a)
        lslogger.debug("iteration %d, alpha=%e, phi(alpha)=%e"%(i,alpha,phi_a))
        if (phi_a > phi0+c1*alpha*gphi0) or ((phi_a>=old_phi_a) and (i>1)):
            alpha, phi_a, gphi_a = _zoom(phi, gradphi, phiargs, old_alpha, alpha, old_phi_a, phi_a, c1, c2, phi0, gphi0)
            break

        gphi_a=gradphi(alpha, *args_a)
        if np.abs(gphi_a) <= -c2*gphi0:
            break
        if gphi_a >= 0:
            alpha, phi_a, gphi_a = _zoom(phi, gradphi, phiargs, alpha, old_alpha, phi_a, old_phi_a, c1, c2, phi0, gphi0)
            break

        old_alpha=alpha
        # the factor is arbitrary as long as there is sufficient increase
        alpha=2.*alpha
        old_phi_a=phi_a
        i+=1

    return alpha, phi_a, gf_new[0]

##############################################################################
class AbstractMinimizer(object):
    """
    Base class for function minimization methods.
    """

    TOLERANCE_REACHED, MAX_ITERATIONS_REACHED=list(xrange(2))

    def __init__(self, f, tol=1e-5, imax=300):
        """
        Initializes a new minimizer for a given cost function.

        :param f: the cost function to be minimized
        :type f: CostFunction
        """
        self._f=f
        self._tol = tol
        self._imax = imax
        self._result = None
        self._callback = None
        self.logger = logging.getLogger('inv.%s'%self.__class__.__name__)

    def setTolerance(self, tol):
        """
        Sets the tolerance for the stopping criterion. The minimizer stops when
        an appropriate norm is less than `tol`.
        """
        self._tol = tol

    def setMaxIterations(self, imax):
        """
        Sets the maximum number of iterations before the minimizer terminates.
        """
        self._imax = imax

    def setCallback(self, callback):
        """
        Sets a callback function to be called after every iteration.
        The arguments to the function are: (k, x, fx, gfx), where
        k is the current iteration, x is the current estimate, fx=f(x) and
        gfx=grad f(x).
        """
        if callback is not None and not callable(callback):
            raise TypeError("Callback function not callable.")
        self._callback = callback

    def _doCallback(self, *args):
        if self._callback is not None:
            self._callback(*args)

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

        :return: `TOLERANCE_REACHED` or `MAX_ITERATIONS_REACHED`
        """
        raise NotImplementedError

    def logSummary(self):
        """
        Outputs a summary of the completed minimization process to the logger.
        """
        if hasattr(self._f, "Value_calls"):
	  self.logger.warning("Number of function evaluations: %d"%self._f.Value_calls)
	  self.logger.warning("Number of gradient evaluations: %d"%self._f.Gradient_calls)
	  self.logger.warning("Number of inner product evaluations: %d"%self._f.DualProduct_calls)
	  self.logger.warning("Number of argument evaluations: %d"%self._f.Arguments_calls)


##############################################################################
class MinimizerLBFGS(AbstractMinimizer):
    """
    Minimizer that uses the limited-memory Broyden-Fletcher-Goldfarb-Shanno
    method.
    """

    # History size
    _m = 15

    # Initial Hessian multiplier
    _initial_H = 1

    def getOptions(self):
        return {'historySize':self._m,'initialHessian':self._initial_H}

    def setOptions(self, **opts):
        for o in opts:
            if o=='historySize':
                self._m=opts[o]
            elif o=='initialHessian':
                self._initial_H=opts[o]
            else:
                raise KeyError("Invalid option '%s'"%o)

    def run(self, x):
        args=self._f.getArguments(x)
        gf=self._f.getGradient(x, *args)
        fx=self._f(x, *args)
        if self._f.provides_inverse_Hessian_approximation:
            self._f.updateHessian()
            invH_scale = None 
        else:
	    invH_scale = self._initial_H
        k=0
        error=2*self._tol
        s_and_y=[]

        self._doCallback(k, x, fx, gf)

        while error > self._tol and k < self._imax:
            #self.logger.info("\033[1;31miteration %d\033[1;30m, error=%e"%(k,error))
            self.logger.info("iteration %d, error=%e"%(k,error))
            # determine search direction
            p = -self._twoLoop(invH_scale, gf, s_and_y, x, *args)

            if invH_scale: self.logger.debug("H = %s"%invH_scale)
            self.logger.debug("grad f(x) = %s"%gf)
            #self.logger.debug("p = %s"%p)
            #self.logger.debug("x = %s"%x)

            # determine step length
            alpha, fx_new, gf_new = line_search(self._f, x, p, gf, fx)
            self.logger.debug("alpha=%e"%(alpha))
            # execute the step
            x_new = x + alpha*p
            if gf_new is None:
	        args=self._f.getArguments(x_new)
                gf_new=self._f.getGradient(x_new, args)

            delta_x=x_new-x
            delta_g=gf_new-gf
            s_and_y.append((delta_x,delta_g, self._f.getDualProduct(delta_x, delta_g) ))
            
            # this needs to be reviewed
            if fx_new==0.:
                error=fx
            else:
                error=abs(fx_new-fx)/abs(fx_new)
            
            self._f.updateHessian()
            x=x_new
            gf=gf_new
            fx=fx_new
            k+=1
            self._doCallback(k, x, fx, gf)
            if (error<=self._tol): break

            # delete oldest vector pair
            if k>self._m: s_and_y.pop(0)

            if not self._f.provides_inverse_Hessian_approximation:            
		  # set the new scaling factor (approximation of inverse Hessian)
		  denom=self._f.getDualProduct(delta_g, delta_g)
		  if denom > 0:
		      invH_scale=self._f.getDualProduct(delta_x,delta_g)/denom
		  else:
		      invH_scale=self._initial_H
		      self.logger.debug("Break down in H update. Resetting to initial value %s."%self._initial_H) 

        if k >= self._imax:
            reason=self.MAX_ITERATIONS_REACHED
            self.logger.warning("Maximum number of iterations reached!")
        else:
            reason=self.TOLERANCE_REACHED
            self.logger.warning("Success after %d iterations! Final error=%e"%(k,error))
        self._result=x
        return reason

    def _twoLoop(self, invH_scale, gf, s_and_y, x, *args):
        """
        Helper for the L-BFGS method.
        See 'Numerical Optimization' by J. Nocedal for an explanation.
        """
        q=gf
        alpha=[]
        for s,y, rho in reversed(s_and_y):
            a=self._f.getDualProduct(s, q)/rho
            alpha.append(a)
            q=q-a*y

        if self._f.provides_inverse_Hessian_approximation:    
             r= self._f.getInverseHessianApproximation(x, q, *args)
        else:
	     r= invH_scale * q
	     
        for s,y,rho in s_and_y:
            beta=self._f.getDualProduct(r, y)/rho
            a=alpha.pop()
            r=r+s*(a-beta)
        return r

##############################################################################
class MinimizerBFGS(AbstractMinimizer):
    """
    Minimizer that uses the Broyden-Fletcher-Goldfarb-Shanno method.
    """

    # Initial Hessian multiplier
    _initial_H = 1

    def getOptions(self):
        return {'initialHessian':self._initial_H}

    def setOptions(self, **opts):
        for o in opts:
            if o=='initialHessian':
                self._initial_H=opts[o]
            else:
                raise KeyError("Invalid option '%s'"%o)

    def run(self, x):
        args=self._f.getArguments(x)
        gf=self._f.getGradient(x, *args)
        fx=self._f(x, *args)
        k=0
        try:
            n=len(x)
        except:
            n=x.getNumberOfDataPoints()
        I=np.eye(n)
        H=self._initial_H*I
        gnorm=Lsup(gf)
        self._doCallback(k, x, fx, gf)

        while gnorm > self._tol and k < self._imax:
            self.logger.info("iteration %d, gnorm=%e"%(k,gnorm))

            # determine search direction
            d=-self._f.getDualProduct(H, gf)

            self.logger.debug("H = %s"%H)
            self.logger.debug("grad f(x) = %s"%gf)
            self.logger.debug("d = %s"%d)
            self.logger.debug("x = %s"%x)

            # determine step length
            alpha, fx, gf_new = line_search(self._f, x, d, gf, fx)
            self.logger.debug("alpha=%e"%alpha)
            # execute the step
            x_new=x+alpha*d
            delta_x=x_new-x
            x=x_new
            if gf_new is None:
                gf_new=self._f.getGradient(x_new)
            delta_g=gf_new-gf
            gf=gf_new
            k+=1
            self._doCallback(k, x, fx, gf)
            gnorm=Lsup(gf)
            if (gnorm<=self._tol): break

            # update Hessian
            denom=self._f.getDualProduct(delta_x, delta_g)
            if denom < EPSILON * gnorm:
                denom=1e-5
                self.logger.debug("Break down in H update. Resetting.")
            rho=1./denom
            self.logger.debug("rho=%e"%rho)
            A=I-rho*delta_x[:,None]*delta_g[None,:]
            AT=I-rho*delta_g[:,None]*delta_x[None,:]
            H=self._f.getDualProduct(A, self._f.getDualProduct(H,AT)) + rho*delta_x[:,None]*delta_x[None,:]
        if k >= self._imax:
            reason=self.MAX_ITERATIONS_REACHED
            self.logger.warning("Maximum number of iterations reached!")
        else:
            reason=self.TOLERANCE_REACHED
            self.logger.warning("Success after %d iterations! Final gnorm=%e"%(k,gnorm))

        self._result=x
        return reason

##############################################################################
class MinimizerNLCG(AbstractMinimizer):
    """
    Minimizer that uses the nonlinear conjugate gradient method
    (Fletcher-Reeves variant).
    """

    def run(self, x):
        i=0
        k=0
        args=self._f.getArguments(x)
        r=-self._f.getGradient(x, *args)
        fx=self._f(x, *args)
        d=r
        delta=self._f.getDualProduct(r,r)
        delta0=delta
        self._doCallback(i, x, fx, -r)

        while i<self._imax and Lsup(r)>self._tol:
            self.logger.info("iteration %d"%i)
            self.logger.debug("grad f(x) = %s"%(-r))
            self.logger.debug("d = %s"%d)
            self.logger.debug("x = %s"%x)

            alpha, fx, gf_new = line_search(self._f, x, d, -r, fx, c2=0.4)
            self.logger.debug("alpha=%e"%(alpha))
            x=x+alpha*d
            r=-self._f.getGradient(x) if gf_new is None else -gf_new
            delta_o=delta
            delta=self._f.getDualProduct(r,r)
            beta=delta/delta_o
            d=r+beta*d
            k=k+1
            try:
                lenx=len(x)
            except:
                lenx=x.getNumberOfDataPoints()
            if k == lenx or self._f.getDualProduct(r,d) <= 0:
                d=r
                k=0
            i+=1
            self._doCallback(i, x, fx, gf_new)

        if i >= self._imax:
            reason=self.MAX_ITERATIONS_REACHED
            self.logger.warning("Maximum number of iterations reached!")
        else:
            reason=self.TOLERANCE_REACHED
            self.logger.warning("Success after %d iterations! Final delta=%e"%(i,delta))

        self._result=x
        return reason


if __name__=="__main__":
    # Example usage with function 'rosen' (minimum=[1,1,...1]):
    from scipy.optimize import rosen, rosen_der
    from esys.downunder import  MeteredCostFunction
    import sys
    N=10
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
    m.setTolerance(1e-5)
    m.setMaxIterations(600)
    m.run(x0)
    m.logSummary()
    print("\tLsup(result)=%.8f"%np.amax(abs(m.getResult())))

    #from scipy.optimize import fmin_cg
    #print("scipy ref=%.8f"%np.amax(abs(fmin_cg(rosen, x0, rosen_der, maxiter=10000))))

