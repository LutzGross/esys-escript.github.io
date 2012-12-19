
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

__all__ = ['MinimizerException', 'MinimizerIterationIncurableBreakDown', 'MinimizerMaxIterReached' , 'AbstractMinimizer', 'MinimizerLBFGS', 'MinimizerBFGS', 'MinimizerNLCG']

import logging
import numpy as np

try:
    from esys.escript import Lsup, sqrt, EPSILON
except:
    Lsup=lambda x: np.amax(abs(x))
    sqrt=np.sqrt
    EPSILON=1e-18

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
        zoomlogger.debug("Zoom.iteration %d, alpha=%e, phi(alpha)=%e"%(i,alpha,phi_a))
        if phi_a > phi0+c1*alpha*gphi0 or phi_a >= phi_lo:
            alpha_hi=alpha
        else:
            gphi_a=gradphi(alpha, *args_a)
            zoomlogger.debug("Zoom.grad(phi(alpha))=%e"%(gphi_a))
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

def line_search(f, x, p, g_Jx, Jx, alpha_truncationax=50.0, c1=1e-4, c2=0.9, IMAX=15):
    """
    Line search method that satisfies the strong Wolfe conditions.
    See Chapter 3 of 'Numerical Optimization' by J. Nocedal for an explanation.

    :param f: callable objective function f(x)
    :param x: start value for the line search
    :param p: search direction
    :param g_Jx: value for the gradient of f at x
    :param Jx: value of f(x)
    :param alpha_truncationax: algorithm terminates if alpha reaches this value
    :param c1: value for Armijo condition (see reference)
    :param c2: value for curvature condition (see reference)
    :param IMAX: maximum number of iterations to perform
    """
    # this stores the latest gradf(x+a*p) which is returned
    g_Jx_new=[g_Jx]

    def phi(a, *args):
        """ phi(a):=f(x+a*p) """
        return f(x+a*p, *args)
    def gradphi(a, *args):
        g_Jx_new[0]=f.getGradient(x+a*p, *args)
        return f.getDualProduct(p, g_Jx_new[0])

    def phiargs(a):
        try:
            args=f.getArguments(x+a*p)
        except:
            args=()
        return args

    old_alpha=0.
    # we assume g_Jx is properly scaled so alpha=1 is a reasonable starting value
    alpha=1.
    if Jx is None:
        args0=phiargs(0.)
        phi0=phi(0., *args0)
    else:
        phi0=Jx
    lslogger.debug("phi(0)=%e"%(phi0))
    gphi0=f.getDualProduct(p, g_Jx) #gradphi(0., *args0)
    lslogger.debug("grad phi(0)=%e"%(gphi0))
    old_phi_a=phi0
    i=1

    while i<IMAX and alpha>0. and alpha<alpha_truncationax:

        args_a=phiargs(alpha)
        phi_a=phi(alpha, *args_a)
        lslogger.debug("Line Search.iteration %d, alpha=%e, phi(alpha)=%e"%(i,alpha,phi_a))
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

    return alpha, phi_a, g_Jx_new[0]

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
   Exception thrown if the iteration scheme encountered an incurable breakdown.
   """
   pass    


##############################################################################
class AbstractMinimizer(object):
    """
    Base class for function minimization methods.
    """

    TOLERANCE_REACHED, MAX_ITERATIONS_REACHED, INCURABLE_BREAKDOWN = list(range(3))

    def __init__(self, J, x_tol=1e-5, J_tol=None, imax=300):
        """
        Initializes a new minimizer for a given cost function.

        :param f: the cost function to be minimized
        :type f: CostFunction
        """
        self._J=J
        self._x_tol = x_tol
        self._J_tol = J_tol
        self._imax = imax
        self._result = None
        self._callback = None
        self.logger = logging.getLogger('inv.%s'%self.__class__.__name__)
        self.setTolerance()

    def setTolerance(self, x_tol=1e-4, J_tol=None):
        """
        Sets the tolerance for the stopping criterion. The minimizer stops when
        an appropriate norm is less than `tol`.
        """
        if x_tol: self._x_tol = x_tol
        if J_tol: self._J_tol = J_tol

    def setMaxIterations(self, imax):
        """
        Sets the maximum number of iterations before the minimizer terminates.
        """
        self._imax = imax

    def setCallback(self, callback):
        """
        Sets a callback function to be called after every iteration.
        The arguments to the function are: (k, x, Jx, g_Jxx), where
        k is the current iteration, x is the current estimate, Jx=f(x) and
        g_Jxx=grad f(x).
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
        if hasattr(self._J, "Value_calls"):
            self.logger.warning("Number of cost function evaluations: %d"%self._J.Value_calls)
            self.logger.warning("Number of gradient evaluations: %d"%self._J.Gradient_calls)
            self.logger.warning("Number of inner product evaluations: %d"%self._J.DualProduct_calls)
            self.logger.warning("Number of argument evaluations: %d"%self._J.Arguments_calls)
            self.logger.warning("Number of norm evaluations: %d"%self._J.Norm_calls)

##############################################################################
class MinimizerLBFGS(AbstractMinimizer):
    """
    Minimizer that uses the limited-memory Broyden-Fletcher-Goldfarb-Shanno
    method.
    """

    # History size
    _truncation = 30

    # Initial Hessian multiplier
    _initial_H = 1
    
    # restart 
    _restart= 60

    def getOptions(self):
        return {'truncation':self._truncation,'initialHessian':self._initial_H, 'restart':self._restart}

    def setOptions(self, **opts):
        for o in opts:
            if o=='historySize' or 'truncation':
                self._truncation=opts[o]
            elif o=='initialHessian':
                self._initial_H=opts[o]
            elif o=='restart':
                self._restart=opts[o]
            else:
                raise KeyError("Invalid option '%s'"%o)

    def run(self, x):

	if self._J.provides_inverse_Hessian_approximation:
	    self._J.updateHessian()
	    invH_scale = None
	else:
	    invH_scale = self._initial_H
	    
	# start the iteration:
	n_iter = 0
	n_last_break_down=-1
	non_curable_break_down = False
	converged = False
	args=self._J.getArguments(x)
	g_Jx=self._J.getGradient(x, *args)
	Jx=self._J(x, *args)
	Jx_0=Jx	
	
	while not converged and not non_curable_break_down and n_iter < self._imax:
	  	 
	  k=0	
	  break_down = False
	  s_and_y=[]
	  self._doCallback(n_iter, x, Jx, g_Jx)

	  while not converged and not break_down and k < self._restart and n_iter < self._imax:
	    #self.logger.info("\033[1;31miteration %d\033[1;30m, error=%e"%(k,error))
	    self.logger.debug("LBFGS.iteration %d .......... "%n_iter)
	    # determine search direction
	    self.logger.debug("LBFGS.J(x) = %s"%Jx)
	    self.logger.debug("LBFGS.grad f(x) = %s"%g_Jx)
            if invH_scale: self.logger.debug("LBFGS.H = %s"%invH_scale)

	    p = -self._twoLoop(invH_scale, g_Jx, s_and_y, x, *args)
	    # determine step length
	    alpha, Jx_new, g_Jx_new = line_search(self._J, x, p, g_Jx, Jx)
	    # this function returns the a scaling alpha for the serch direction
	    # as well the cost function evaluation and gradient for the new solution 
	    # approximation x_new = x + alpha*p
	    
	    self.logger.debug("LBFGS.search direction scaling alpha=%e"%(alpha))
	    # execute the step
	    delta_x = alpha*p 
	    x_new = x + delta_x
	    self.logger.debug("LBFGS.J(x) = %s"%Jx_new)

	    converged = True
	    if self._J_tol:
		flag= abs(Jx_new-Jx) <= self._J_tol * abs(Jx_new-Jx_0)
		if flag: 
		    self.logger.debug("LBFGS: cost function has converged: dJ, J*J_tol = %e, %e"%(Jx-Jx_new,abs(Jx_new-Jx_0)*self._J_tol))
		else:    
		    self.logger.debug("LBFGS: cost function checked: dJ, J * J_tol = %e, %e"%(Jx-Jx_new,abs(Jx_new)*self._J_tol))

		converged = converged and flag
	    if self._x_tol:	      
	      norm_x=self._J.getNorm(x_new)
	      norm_dx=self._J.getNorm(delta_x)
	      flag= norm_dx <= self._x_tol * norm_x
	      if flag:
		  self.logger.debug("LBFGS: solution has converged: dx, x*x_tol = %e, %e"%(norm_dx,norm_x*self._x_tol))
	      else:
		  self.logger.debug("LBFGS: solution checked: dx, x*x_tol = %e, %e"%(norm_dx,norm_x*self._x_tol))
	      converged = converged and flag
	    
	    x=x_new
	    if converged:
	      break
	    
	    # unfortunatly there is more work to do!
	    if g_Jx_new is None:
		args=self._J.getArguments(x_new)
		g_Jx_new=self._J.getGradient(x_new, args)
	    delta_g=g_Jx_new-g_Jx
	    
	    rho=self._J.getDualProduct(delta_x, delta_g)
	    if abs(rho)>0 :
		s_and_y.append((delta_x,delta_g, rho ))
	    else:
	        break_down=True

	    self._J.updateHessian()	    
	    g_Jx=g_Jx_new
	    Jx=Jx_new
	    
	    k+=1
	    n_iter+=1
	    self._doCallback(k, x, Jx, g_Jx)

	    # delete oldest vector pair
	    if k>self._truncation: s_and_y.pop(0)

	    if not self._J.provides_inverse_Hessian_approximation and not break_down :
		  # set the new scaling factor (approximation of inverse Hessian)
		  denom=self._J.getDualProduct(delta_g, delta_g)
		  if denom > 0:
		      invH_scale=self._J.getDualProduct(delta_x,delta_g)/denom
		  else:
		      invH_scale=self._initial_H
		      self.logger.debug("LBFGS.Break down in H update. Resetting to initial value %s."%self._initial_H)
          # case handeling for inner iteration: 
          if break_down:
	     if n_iter == n_last_break_down +1 :
	        non_curable_break_down = True
	        self.logger.debug("LBFGS. Incurable break down detected in step %d."%(n_iter,))
	     else: 
	        n_last_break_down = n_iter
	        self.logger.debug("LBFGS.Break down detected in step %d. Iteration is restarted."%(n_iter,))
          if not k < self._restart:
	        self.logger.debug("LBFGS.Iteration is restarted after %d steps."%(n_iter,))

	# case handeling for inner iteration:         
	self._result=x    
	if n_iter >= self._imax:
	    self.return_status=self.MAX_ITERATIONS_REACHED
	    raise MinimizerMaxIterReached("Gave up after %d steps."%(n_iter,))
	    self.logger.warning("LBFGS.>>>>>> Maximum number of iterations reached!")
	elif non_curable_break_down:
	    self.return_status=self.INCURABLE_BREAKDOWN	
	    self.logger.warning("LBFGS.>>>>>> Uncurable breakdown!")
	    raise MinimizerIterationIncurableBreakDown("Gave up after %d steps."%(n_iter,))
	else: 
	    self.return_status=self.TOLERANCE_REACHED
	    self.logger.warning("LBFGS.Success after %d iterations!"%k)

	return self.return_status
	
    def _twoLoop(self, invH_scale, g_Jx, s_and_y, x, *args):
        """
        Helper for the L-BFGS method.
        See 'Numerical Optimization' by J. Nocedal for an explanation.
        """
        q=g_Jx
        alpha=[]
        for s,y, rho in reversed(s_and_y):
            a=self._J.getDualProduct(s, q)/rho
            alpha.append(a)
            q=q-a*y

        if self._J.provides_inverse_Hessian_approximation:
             r = self._J.getInverseHessianApproximation(x, q, *args)
        else:
             r = invH_scale * q

        for s,y,rho in s_and_y:
            beta = self._J.getDualProduct(r, y)/rho
            a = alpha.pop()
            r = r + s * (a-beta)
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
        args=self._J.getArguments(x)
        g_Jx=self._J.getGradient(x, *args)
        Jx=self._J(x, *args)
        k=0
        try:
            n=len(x)
        except:
            n=x.getNumberOfDataPoints()
        I=np.eye(n)
        H=self._initial_H*I
        gnorm=Lsup(g_Jx)
        self._doCallback(k, x, Jx, g_Jx)

        while gnorm > self._x_tol and k < self._imax:
            self.logger.debug("iteration %d, gnorm=%e"%(k,gnorm))

            # determine search direction
            d=-self._J.getDualProduct(H, g_Jx)

            self.logger.debug("H = %s"%H)
            self.logger.debug("grad f(x) = %s"%g_Jx)
            self.logger.debug("d = %s"%d)
            self.logger.debug("x = %s"%x)

            # determine step length
            alpha, Jx, g_Jx_new = line_search(self._J, x, d, g_Jx, Jx)
            self.logger.debug("alpha=%e"%alpha)
            # execute the step
            x_new=x+alpha*d
            delta_x=x_new-x
            x=x_new
            if g_Jx_new is None:
                g_Jx_new=self._J.getGradient(x_new)
            delta_g=g_Jx_new-g_Jx
            g_Jx=g_Jx_new
            k+=1
            self._doCallback(k, x, Jx, g_Jx)
            gnorm=Lsup(g_Jx)
            if (gnorm<=self._x_tol): break

            # update Hessian
            denom=self._J.getDualProduct(delta_x, delta_g)
            if denom < EPSILON * gnorm:
                denom=1e-5
                self.logger.debug("Break down in H update. Resetting.")
            rho=1./denom
            self.logger.debug("rho=%e"%rho)
            A=I-rho*delta_x[:,None]*delta_g[None,:]
            AT=I-rho*delta_g[:,None]*delta_x[None,:]
            H=self._J.getDualProduct(A, self._J.getDualProduct(H,AT)) + rho*delta_x[:,None]*delta_x[None,:]
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
        args=self._J.getArguments(x)
        r=-self._J.getGradient(x, *args)
        Jx=self._J(x, *args)
        d=r
        delta=self._J.getDualProduct(r,r)
        delta0=delta
        self._doCallback(i, x, Jx, -r)

        while i<self._imax and Lsup(r)>self._x_tol:
            self.logger.debug("iteration %d"%i)
            self.logger.debug("grad f(x) = %s"%(-r))
            self.logger.debug("d = %s"%d)
            self.logger.debug("x = %s"%x)

            alpha, Jx, g_Jx_new = line_search(self._J, x, d, -r, Jx, c2=0.4)
            self.logger.debug("alpha=%e"%(alpha))
            x=x+alpha*d
            r=-self._J.getGradient(x) if g_Jx_new is None else -g_Jx_new
            delta_o=delta
            delta=self._J.getDualProduct(r,r)
            beta=delta/delta_o
            d=r+beta*d
            k=k+1
            try:
                lenx=len(x)
            except:
                lenx=x.getNumberOfDataPoints()
            if k == lenx or self._J.getDualProduct(r,d) <= 0:
                d=r
                k=0
            i+=1
            self._doCallback(i, x, Jx, g_Jx_new)

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
    m.setTolerance(x_tol=1e-5)
    m.setMaxIterations(600)
    m.run(x0)
    m.logSummary()
    print("\tLsup(result)=%.8f"%np.amax(abs(m.getResult())))

    #from scipy.optimize import fmin_cg
    #print("scipy ref=%.8f"%np.amax(abs(fmin_cg(rosen, x0, rosen_der, maxiter=10000))))

