##############################################################################
#
# Copyright (c) 2014-2015 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from .minimizers import AbstractMinimizer
from esys.escriptcore.splitworld import Job


class SplitMinimizerLBFGS(AbstractMinimizer):
    """
    Minimizer that uses the limited-memory Broyden-Fletcher-Goldfarb-Shanno
    method.
    
    version modified to fit with split world.
    """

    # History size
    _truncation = 30

    # Initial Hessian multiplier
    _initial_H = 1

    # Restart after this many iteration steps
    _restart = 60

    def getOptions(self):
        return {'truncation':self._truncation,'initialHessian':self._initial_H, 'restart':self._restart}

    def setOptions(self, **opts):
        self.logger.debug("Setting options: %s"%(str(opts)))
        for o in opts:
            if o=='historySize' or o=='truncation':
                self._truncation=opts[o]
            elif o=='initialHessian':
                self._initial_H=opts[o]
            elif o=='restart':
                self._restart=opts[o]
            else:
                raise KeyError("Invalid option '%s'"%o)

    def run(self):
        """
        This version relies on the costfunction already having an initial guess loaded.
        It also does not return the result, meaning a job needs to be submitted to
        get the result out.
        """
        if self.getCostFunction().provides_inverse_Hessian_approximation:
            self.getCostFunction().updateHessian()
            invH_scale = None
        else:
            invH_scale = self._initial_H

        # start the iteration:
        n_iter = 0
        n_last_break_down=-1
        non_curable_break_down = False
        converged = False
        

	    
	  
        self.getCostFunction().setPoint()	# Set point to initial guess value (takes the place of a getArgs call)
        #args=self.getCostFunction().getArguments(x)
        
        self.getCostFunction().calculateValue(["Jx","Jx_0"])	#evaluate the function and store the result in the named variables
        self.getCostFunction().calculateGradient("g_Jx")        #compute the gradient and store the result
        
        #g_Jx=self.getCostFunction().getGradient(x, *args)
        #Jx=self.getCostFunction()(x, *args) # equivalent to getValue() for Downunder CostFunctions
        
        
        
        Jx_0=Jx

        while not converged and not non_curable_break_down and n_iter < self._imax:
          k=0
          break_down = False
          s_and_y=[]
          # initial step length for line search
          alpha=1.0
          self._doCallback(n_iter, x, Jx, g_Jx)

          while not converged and not break_down and k < self._restart and n_iter < self._imax:
                #self.logger.info("\033[1;31miteration %d\033[1;30m"%n_iter)
                self.logger.info("********** iteration %3d **********"%n_iter)
                self.logger.info("\tJ(x) = %s"%Jx)
                #self.logger.debug("\tgrad f(x) = %s"%g_Jx)
                if invH_scale:
                    self.logger.debug("\tH = %s"%invH_scale)

                # determine search direction
                p = -self._twoLoop(invH_scale, g_Jx, s_and_y, x, *args)

                # determine new step length using the last one as initial value
                # however, avoid using too small steps for too long.
                # FIXME: This is a bit specific to esys.downunder in that the
                # inverse Hessian approximation is not scaled properly (only
                # the regularization term is used at the moment)...
                if invH_scale is None:
                    if alpha <= 0.5:
                        alpha=2*alpha
                else:
                    # reset alpha for the case that the cost function does not
                    # provide an approximation of inverse H
                    alpha=1.0
                alpha, Jx_new, g_Jx_new = line_search(self.getCostFunction(), x, p, g_Jx, Jx, alpha)
                # this function returns a scaling alpha for the search
                # direction as well as the cost function evaluation and
                # gradient for the new solution approximation x_new=x+alpha*p
                self.logger.debug("\tSearch direction scaling alpha=%e"%alpha)

                # execute the step
                delta_x = alpha*p
                x_new = x + delta_x

                converged = True
                if self._J_tol:
                    flag=abs(Jx_new-Jx) <= self._J_tol * abs(Jx_new-Jx_0)
                    if self.logger.isEnabledFor(logging.DEBUG):
                        if flag:
                            self.logger.debug("Cost function has converged: dJ, J*J_tol = %e, %e"%(Jx-Jx_new,abs(Jx_new-Jx_0)*self._J_tol))
                        else:
                            self.logger.debug("Cost function checked: dJ, J*J_tol = %e, %e"%(Jx-Jx_new,abs(Jx_new)*self._J_tol))

                    converged = converged and flag
                if self._m_tol:
                    norm_x = self.getCostFunction().getNorm(x_new)
                    norm_dx = self.getCostFunction().getNorm(delta_x)
                    flag = norm_dx <= self._m_tol * norm_x
                    if self.logger.isEnabledFor(logging.DEBUG):
                        if flag:
                            self.logger.debug("Solution has converged: dx, x*m_tol = %e, %e"%(norm_dx,norm_x*self._m_tol))
                        else:
                            self.logger.debug("Solution checked: dx, x*m_tol = %e, %e"%(norm_dx,norm_x*self._m_tol))
                    converged = converged and flag

                x=x_new
                if converged:
                    self.logger.info("\tJ(x) = %s"%Jx_new)
                    break

                # unfortunately there is more work to do!
                if g_Jx_new is None:
                    args=self.getCostFunction().getArguments(x_new)
                    g_Jx_new=self.getCostFunction().getGradient(x_new, args)
                delta_g=g_Jx_new-g_Jx

                rho=self.getCostFunction().getDualProduct(delta_x, delta_g)
                if abs(rho)>0:
                    s_and_y.append((delta_x,delta_g, rho ))
                else:
                    break_down=True

                self.getCostFunction().updateHessian()
                g_Jx=g_Jx_new
                Jx=Jx_new

                k+=1
                n_iter+=1
                self._doCallback(n_iter, x, Jx, g_Jx)

                # delete oldest vector pair
                if k>self._truncation: s_and_y.pop(0)

                if not self.getCostFunction().provides_inverse_Hessian_approximation and not break_down:
                    # set the new scaling factor (approximation of inverse Hessian)
                    denom=self.getCostFunction().getDualProduct(delta_g, delta_g)
                    if denom > 0:
                        invH_scale=self.getCostFunction().getDualProduct(delta_x,delta_g)/denom
                    else:
                        invH_scale=self._initial_H
                        self.logger.debug("** Break down in H update. Resetting to initial value %s."%self._initial_H)
          # case handling for inner iteration:
          if break_down:
              if n_iter == n_last_break_down+1:
                  non_curable_break_down = True
                  self.logger.debug("** Incurable break down detected in step %d."%n_iter)
              else:
                  n_last_break_down = n_iter
                  self.logger.debug("** Break down detected in step %d. Iteration is restarted."%n_iter)
          if not k < self._restart:
              self.logger.debug("Iteration is restarted after %d steps."%n_iter)

        # case handling for inner iteration:
        self._result=x
        if n_iter >= self._imax:
            self.logger.warn(">>>>>>>>>> Maximum number of iterations reached! <<<<<<<<<<")
            raise MinimizerMaxIterReached("Gave up after %d steps."%n_iter)
        elif non_curable_break_down:
            self.logger.warn(">>>>>>>>>> Incurable breakdown! <<<<<<<<<<")
            raise MinimizerIterationIncurableBreakDown("Gave up after %d steps."%n_iter)

        self.logger.info("Success after %d iterations!"%n_iter)
        #This version does nor return the result
        #You need to set up a job to extract the result

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
        else:
             r = invH_scale * q

        for s,y,rho in s_and_y:
            beta = self.getCostFunction().getDualProduct(r, y)/rho
            a = alpha.pop()
            r = r + s * (a-beta)
        return r
        
    