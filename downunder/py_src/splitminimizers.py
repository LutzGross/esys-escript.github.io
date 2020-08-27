##############################################################################
#
# Copyright (c) 2014-2018 by The University of Queensland
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

from __future__ import print_function, division

from .minimizers import AbstractMinimizer
from esys.escriptcore.splitworld import Job, FunctionJob
from .splitinversioncostfunctions import SplitInversionCostFunction
from esys.escript.pdetools import ArithmeticTuple
import numpy as np

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

        # This function sets current_point=base_point+alpha*search_direction [m=m+p*a]
    @staticmethod
    def move_point_from_base(self, **kwargs):
            m=self.importValue('base_point')
            p=self.importValue('search_direction')
            a=kwargs['alpha']
            newpoint=m+p*a
            SplitInversionCostFunction.update_point_helper(self, newpoint)

    def run(self):
        """
        This version relies on the costfunction already having an initial guess loaded.
        It also does not return the result, meaning a job needs to be submitted to
        get the result out.
        """
        
        # First we'll define our own versions of the helper functions

        # Updates "g_Jx_new_0" and "g_Jx_new_1" to the gradient at the current point
        # then returns f.dualProduct of search_direction and g_Jx_new
        # Not a splitworld function
        def grad_phi(f, **kwargs):
                f.calculateGradient('g_Jx_0', 'g_Jx_1')
                # need to call dualProduct here
                def dual_p_g_Jx_new(self, **kwargs):
                    p=self.importValue("search_direction")
                    g_Jx_0=self.importValue("g_Jx_0")
                    g_Jx_1=self.importValue("g_Jx_1")
                    reg=self.importValue("regularization")
                        #again, this assumes that only the regularization term is relevant
                    res=reg.getDualProduct(p, (g_Jx_0, g_Jx_1))
                    self.exportValue("dp_result",res)
                # Now we will only run this on one world and rely on getDouble to ship it
                f.splitworld.addJob( FunctionJob, dual_p_g_Jx_new)
                f.splitworld.runJobs()
                res=f.splitworld.getFloatVariable("dp_result")
                return res
        #End of grad_phi

        def _zoom(f, alpha_lo, alpha_hi, phi_lo, phi_hi, c1, c2,
                phi0, gphi0, IMAX=25):
            """
            Helper function for `line_search` below which tries to tighten the range
            alpha_lo...alpha_hi. See Chapter 3 of 'Numerical Optimization' by
            J. Nocedal for an explanation.
            """
            i=0
            while True:
                alpha=alpha_lo+.5*(alpha_hi-alpha_lo) # should use interpolation...
                f.splitworld.addJobPerWorld( FunctionJob, SplitMinimizerLBFGS.move_point_from_base, alpha=alpha, imports=['base_point', 'search_direction'])
                f.splitworld.runJobs()
                f.calculateValue('phi_a')
                phi_a=f.splitworld.getFloatVariable('phi_a')
                #zoomlogger.debug("iteration %d, alpha=%e, phi(alpha)=%e"%(i,alpha,phi_a))
                if phi_a > phi0+c1*alpha*gphi0 or phi_a >= phi_lo:
                    alpha_hi=alpha
                else:
                    gphi_a=grad_phi(f)
                    #zoomlogger.debug("\tgrad(phi(alpha))=%e"%(gphi_a))
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

        def line_search(f, alpha=1.0, alpha_truncationax=50.0,
                        c1=1e-4, c2=0.9, IMAX=15):
            """
            Line search method that satisfies the strong Wolfe conditions.
            See Chapter 3 of 'Numerical Optimization' by J. Nocedal for an explanation.

            This version is converted from the line_search from minimizers.py
            however, it takes fewer parameters because some of the values needed
            by the original version will be available as subworld variables rather
            than as parameters.
            
            
            :param f: callable objective function f(x)
            :param p: search direction
            :param alpha: initial step length. If g_Jx is properly scaled alpha=1 is a
                        reasonable starting value.
            :param alpha_truncationax: algorithm terminates if alpha reaches this value
            :param c1: value for Armijo condition (see reference)
            :param c2: value for curvature condition (see reference)
            :param IMAX: maximum number of iterations to perform
                        
            Removed parameters (now in subworld variables instead):
            x    - The start value for line search: in the variable "current_point".
            p    - search direction: in the variable "search_direction"
            g_Jx - value for the gradient of f at x: in the variables "g_Jx_0" and "g_Jx_1"
            Jx   - value of f(x): in the variable "Jx"
            """
            
            # This will handle subworld side of work
            def line_search_init_worker(self, **kwargs):
                x=self.importValue("current_point")
                p=self.importValue("search_direction")
                g_Jx=ArithmeticTuple(self.importValue("g_Jx_0"), self.importValue("g_Jx_1"))
                Jx=self.importValue("Jx")
                regular=self.importValue("regularization")
                phi0=Jx
               
                # In the original, this part calls getDualProduct on f
                # However, since that only ends up referring to the 
                # regularisation term, I've called that directly
                # If your dual product operation requires access to
                # the other models, then  this step needs
                # a rethink since not all models are local
                gphi0=regular.getDualProduct(p, g_Jx)
            
                #Still need to decide what this worker will do
                old_phi_a=phi0
                phi_a=phi0
                self.exportValue("old_phi_a", old_phi_a)
                self.exportValue("phi_a", phi_a)
                self.exportValue("gphi0", gphi0)
                self.exportValue("phi0", phi0)
                self.exportValue("base_point", x)       # To ensure we can revert if needed
            #End of line_search_init_worker
            
            old_alpha=0
            i=1
            f.splitworld.addJobPerWorld( FunctionJob, line_search_init_worker, imports=['search_direction', 'g_Jx_0', 'g_Jx_1', 'Jx', 'regularization'])       
            f.splitworld.runJobs()

        
            # Updates "g_Jx_new_0" and "g_Jx_new_1" to the gradient at the current point
            # then returns f.dualProduct of search_direction and g_Jx_new
            # Not a splitworld function
            def grad_phi(f, **kwargs):
                f.calculateGradient('g_Jx_0', 'g_Jx_1')
                # need to call dualProduct here
                def dual_p_g_Jx_new(self, **kwargs):
                    p=self.importValue("search_direction")
                    g_Jx_0=self.importValue("g_Jx_0")
                    g_Jx_1=self.importValue("g_Jx_1")
                    reg=self.importValue("regularization")
                        #again, this assumes that only the regularization term is relevant
                    res=reg.getDualProduct(p, (g_Jx_0, g_Jx_1))
                    self.exportValue("dp_result",res)
                # Now we will only run this on one world and rely on getDouble to ship it
                f.splitworld.addJob( FunctionJob, dual_p_g_Jx_new)
                f.splitworld.runJobs()
                res=f.splitworld.getFloatVariable("dp_result")
                return res
            #End of grad_phi

            while i<IMAX and alpha>0. and alpha<alpha_truncationax:
                alpha_at_loop_start=alpha
                f.splitworld.addJobPerWorld( FunctionJob, SplitMinimizerLBFGS.move_point_from_base, alpha=alpha, imports=['current_point', 'search_direction'])
                f.splitworld.runJobs()
                f.calculateValue('phi_a')
                #lslogger.debug("iteration %d, alpha=%e, phi(alpha)=%e"%(i,alpha,phi_a))
                phi_a=f.splitworld.getFloatVariable('phi_a')
                phi0=f.splitworld.getFloatVariable('phi0')
                gphi0=f.splitworld.getFloatVariable('gphi0')
                old_phi_a=f.splitworld.getFloatVariable('old_phi_a')
                if (phi_a > phi0+c1*alpha*gphi0) or ((phi_a>=old_phi_a) and (i>1)):
                    alpha, phi_a, gphi_a = _zoom(f, old_alpha, alpha, old_phi_a, phi_a, c1, c2, phi0, gphi0)
                    break

                   # Need to check if alpha has changed. If it has, we need to move the point again
                if alpha_at_loop_start!=alpha:
                   f.splitworld.addJobPerWorld( FunctionJob, SplitMinimizerLBFGS.move_point_from_base, alpha=alpha, imports=['current_point', 'search_direction'])
                   f.splitworld.runJobs()

                gphi_a=grad_phi(f)
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
            #return alpha, phi_a, g_Jx_new[0]
            return alpha, phi_a
        #End of line_search
        
        
        
        
        
        
        
        splitworld=self.getCostFunction().splitworld 
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
        

            
          
        self.getCostFunction().setPoint()       # Set point to initial guess value (takes the place of a getArgs call)
        #args=self.getCostFunction().getArguments(x)
        
        self.getCostFunction().calculateValue("Jx")    #evaluate the function and store the result in the named variable
                      # note that call sets Jx=Jx_original
        splitworld.copyVariable("Jx", "Jx_original")
                      
        self.getCostFunction().calculateGradient("g_Jx_0","g_Jx_1")        #compute the gradient and store the result
        
        while not converged and not non_curable_break_down and n_iter < self._imax:
          k=0
          break_down = False
          reset_s_and_y = True
          # initial step length for line search
          alpha=1.0
          #self._doCallback(n_iter, x, Jx, g_Jx)

          while not converged and not break_down and k < self._restart and n_iter < self._imax:
                #self.logger.info("\033[1;31miteration %d\033[1;30m"%n_iter)
                self.logger.info("********** iteration %3d **********"%n_iter)
                #self.logger.info("\tJ(x) = %s"%Jx)
                #self.logger.debug("\tgrad f(x) = %s"%g_Jx)
                if invH_scale:
                    self.logger.debug("\tH = %s"%invH_scale)

                splitworld.copyVariable("g_Jx_0", "old_g_Jx_0")
                splitworld.copyVariable("g_Jx_1", "old_g_Jx_1")
                splitworld.copyVariable("Jx", "Jx_old")

                # determine search direction
                self._twoLoop(self.getCostFunction().splitworld, reset_s_and_y)
                reset_s_and_y = False

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
                alpha, phi_a = line_search(self.getCostFunction(), alpha)
                # this function returns a scaling alpha for the search
                # direction as well as the cost function evaluation and
                # gradient for the new solution approximation x_new=x+alpha*p
                self.logger.debug("\tSearch direction scaling alpha=%e"%alpha)

                # execute the step
                # This update operation has already been done in the line_search
                #delta_x = alpha*p
                #x_new = x + delta_x

                Jx=splitworld.getFloatVariable("Jx")
                converged = True
                if self._J_tol:
                    Jx_old=splitworld.getFloatVariable("Jx_old")
                    Jx_original=splitworld.getFloatVariable("Jx_original")
                    flag=abs(Jx-Jx_old) <= self._J_tol * abs(Jx-Jx_original)
                    #flag=abs(Jx_new-Jx) <= self._J_tol * abs(Jx_new-Jx_0)
                    if self.logger.isEnabledFor(logging.DEBUG):
                        if flag:
                            self.logger.debug("Cost function has converged: dJ, J*J_tol = %e, %e"%(Jx_old-Jx,abs(Jx-Jx_original)*self._J_tol))
                        else:
                            self.logger.debug("Cost function checked: dJ, J*J_tol = %e, %e"%(Jx_old-Jx,abs(Jx)*self._J_tol))

                    converged = converged and flag
                if self._m_tol:
                    def converged_check(self, **kwargs):
                        alpha=kwargs["alpha"]
                        m_tol=kwargs["m_tol"]
                        reg=self.importValue("regularization")
                        p=self.importValue("search_direction")
                        delta_x=alpha*p
                        x=self.importValue("current_point")
                        norm_x = reg.getNorm(x)
                        norm_dx = reg.getNorm(delta_x)
                        flag = norm_dx <= m_tol * norm_x
                        #if self.logger.isEnabledFor(logging.DEBUG):
                            #if flag:
                                #self.logger.debug("Solution has converged: dx, x*m_tol = %e, %e"%(norm_dx,norm_x*self._m_tol))
                            #else:
                                #self.logger.debug("Solution checked: dx, x*m_tol = %e, %e"%(norm_dx,norm_x*self._m_tol))
                        self.exportValue('conv_flag', flag)
                    # End of converged_check 
                    self.getCostFunction().splitworld.addJobPerWorld(FunctionJob, converged_check, alpha=alpha, m_tol=self._m_tol, imports=["regularization", "search_direction", "current_point"])
                    self.getCostFunction().splitworld.runJobs()
                    flag=self.getCostFunction().splitworld.getFloatVariable("conv_flag")>0.001
                    converged = converged and flag

                #Already done in the line_search call
                #x=x_new
                if converged:
                    self.logger.info("\tJ(x) = %s"%Jx)
                    break

                # unfortunately there is more work to do!
                def run_worker(self, **kwargs):
                    break_down=False
                    need_trunc=kwargs["need_trunc"]
                    # Need to do some imports here
                    alpha=kwargs["alpha"]
                    g_Jx=ArithmeticTuple(self.importValue("g_Jx_0"), self.importValue("g_Jx_1"))
                    old_g_Jx=ArithmeticTuple(self.importValue("old_g_Jx_0"), self.importValue("old_g_Jx_1"))
                    p=self.importValue("search_direction")
                    reg=self.importValue("regularization")
                    s_and_y=self.importValue("s_and_y")

                    ##The original code had this check
                    #if g_Jx_new is None:
                        #args=self.getCostFunction().getArguments(x_new)
                        #g_Jx_new=self.getCostFunction().getGradient(x_new, args)
                    #delta_g=g_Jx_new-g_Jx
                    delta_g=g_Jx-old_g_Jx
                    delta_x=alpha*p;
                    rho=reg.getDualProduct(delta_x, delta_g)
                    if abs(rho)>0:
                        s_and_y.append((delta_x,delta_g, rho ))
                    else:
                        break_down=True
                    if need_trunc:
                        s_and_y.pop(0)
                    self.exportValue("break_down", break_down)
                    self.exportValue("s_and_y",s_and_y)
                # End run_worker
                # Only one world has s_and_y - so its important that only single jobs be run that manipulate
                # s_and_y
                self.getCostFunction().splitworld.addJob( FunctionJob, run_worker, alpha=alpha, need_trunc=(k>=self._truncation))
                self.getCostFunction().splitworld.runJobs()

                break_down=(splitworld.getLocalObjectVariable("break_down")>0.001)
                self.getCostFunction().updateHessian()
                #g_Jx=g_Jx_new
                #Jx=Jx_new

                k+=1
                n_iter+=1
                #self._doCallback(n_iter, x, Jx, g_Jx)

                # delete oldest vector pair
                # Already done in the run_worker
                #if k>self._truncation: s_and_y.pop(0)

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
        #self._result=x
        self._result=None
        if n_iter >= self._imax:
            self.logger.warn(">>>>>>>>>> Maximum number of iterations reached! <<<<<<<<<<")
            raise MinimizerMaxIterReached("Gave up after %d steps."%n_iter)
        elif non_curable_break_down:
            self.logger.warn(">>>>>>>>>> Incurable breakdown! <<<<<<<<<<")
            raise MinimizerIterationIncurableBreakDown("Gave up after %d steps."%n_iter)

        self.logger.info("Success after %d iterations!"%n_iter)
        #This version does nor return the result
        #You need to set up a job to extract the result

    def _twoLoop(self, splitworld, reset):
        """
        Helper for the L-BFGS method.
        See 'Numerical Optimization' by J. Nocedal for an explanation.

        This has been converted to use splitworld.
          As such it doesn't return a result, instead it stores
          The "p" result in the "search_direction" variable

        Depends on the following splitworld variables:
          g_Jx_0, g_Jx_1
          current_point
          s_and_y 

        Other conversion notes:
          The current cost function's inverseHessianApproximation and
          dualProduct only
          depend on the regularization term so this function can be 
          done entirely by one world.   If the Hessian/dp ever needed
          other worlds, this will need to be reworked.
          Implicit in the above is that the overall cost function
          has cf.provides_inverse_Hessian_approximation==True
        """
        # Make a fn to push the work into subworld
        def two_loop_worker(self, **kwargs):
            reset=kwargs['reset']
            x=self.importValue("current_point")
            reg=self.importValue("regularization")
            g_Jx=ArithmeticTuple(self.importValue("g_Jx_0"), self.importValue("g_Jx_1"))
            if reset:
                s_and_y = []
                self.exportValue('s_and_y', list())
            else:
                s_and_y=self.importValue("s_and_y")
            q=g_Jx
            alpha=[]
            for s,y, rho in reversed(s_and_y):
                a=reg.getDualProduct(s, q)/rho
                alpha.append(a)
                q=q-a*y

            r=reg.getInverseHessianApproximationAtPoint(q)

            for s,y,rho in s_and_y:
                beta = reg.getDualProduct(r, y)/rho
                a = alpha.pop()
                r = r + s * (a-beta)
            # In the original version, the caller negated the result
            self.exportValue("search_direction", -r)
        # Note: this is only running on world world (of possibly many)
        # Any jobs which also need the variables exported here
        # need to be shared or the jobs need to run on the same world
        splitworld.addJob(FunctionJob, two_loop_worker, reset=reset, imports=["current_point", "regularization", "g_Jx_0", "g_Jx_1"])
        splitworld.runJobs() 


