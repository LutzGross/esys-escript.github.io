
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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

#from esys.escript import *
from esys.escript.modelframe import Model,IterationDivergenceError
from esys.escript.linearPDEs import LinearPDE

class Mechanics(Model):
      """
      base class for mechanics models in updated lagrangean framework

      :note: Instance variable domain - domain (in)
      :note: Instance variable internal_force - =Data()
      :note: Instance variable external_force - =Data()
      :note: Instance variable prescribed_velocity - =Data()
      :note: Instance variable location_prescribed_velocity - =Data()
      :note: Instance variable temperature -  = None
      :note: Instance variable expansion_coefficient -  = 0.
      :note: Instance variable bulk_modulus - =1.
      :note: Instance variable shear_modulus - =1.
      :note: Instance variable rel_tol - =1.e-3
      :note: Instance variable abs_tol - =1.e-15
      :note: Instance variable max_iter - =10
      :note: Instance variable displacement - =None
      :note: Instance variable stress - =None
      """
      SAFTY_FACTOR_ITERATION=1./100.
      def __init__(self,**kwargs):
         """
         set up the model
         
         :keyword debug: debug flag
         :type debug: ``bool``
         """
         super(Mechanics, self).__init__(self,**kwargs)
         self.declareParameter(domain=None, \
                               displacement=None, \
                               stress=None, \
                               velocity=None, \
                               internal_force=None, \
                               external_force=None, \
                               prescribed_velocity=None, \
                               location_prescribed_velocity=None, \
                               temperature = None, \
                               expansion_coefficient = 0., \
                               bulk_modulus=2., \
                               shear_modulus=1., \
                               rel_tol=1.e-3,abs_tol=1.e-15,max_iter=10)
         self.__iter=0

      def doInitialization(self):
           """
           initialize model
           """
           if not self.displacement: self.displacement=Vector(0.,ContinuousFunction(self.domain))
           if not self.velocity: self.velocity=Vector(0.,ContinuousFunction(self.domain))
           if not self.stress: self.stress=Tensor(0.,ContinuousFunction(self.domain))
           if not self.internal_force: self.internal_force = Data()
           if not self.external_force: self.external_force = Data()
           if not self.prescribed_velocity: self.prescribed_velocity = Data()
           if not self.location_prescribed_velocity: self.location_prescribed_velocity =Data()
           # save the old values:
           self.__stress_safe=self.stress
           self.__temperature_safe=self.temperature
           self.__displacement_safe=self.displacement
           self.__velocity_safe=self.velocity
           self.__velocity_old=None
           self.__old_dt=None
           self.__very_old_dt=None
           # get node cooridnates and apply initial displacement
           self.__x=self.domain.getX()
           self.domain.setX(self.__x+self.displacement)
           # open PDE:
           self.__pde=LinearPDE(self.domain)
           self.__pde.setSolverMethod(self.__pde.DIRECT)
           self.__solver_options=self.__pde.getSolverOptions()
           self.__solver_options.setSolverMethod(self.__solver_options.DIRECT)
           self.__solver_options.setVerbosity(self.debug)

           # self.__pde.setSymmetryOn()

      def doStepPreprocessing(self,dt):
            """
            step up pressure iteration

            if run within a time dependend problem extrapolation of pressure from previous time steps is used to
            get an initial guess (that needs some work!!!!!!!)
            """
            # reset iteration counters:
            self.__iter=0
            self.__diff=self.UNDEF_DT
            # set initial guesses for the iteration:
            self.displacement=self.__displacement_safe
            self.stress=self.__stress_safe
            self.velocity=self.__velocity_safe
            # update geometry
            self.domain.setX(self.__x+self.displacement)

      def doStep(self,dt):
          """
          """
          self.__iter+=1
          k3=kronecker(self.domain)
          # set new thermal stress increment
          if self.temperature == None:
             self.deps_th=0.
          else:
             self.deps_th=self.expansion_coefficient*(self.temperature-self.__temperature_safe)
          # set PDE coefficients:
          self.__pde.setValue(A=self.S)
          self.__pde.setValue(X=-self.stress-self.bulk_modulus*self.deps_th*k3)
          if self.internal_force: self.__pde.setValue(Y=self.internal_force)
          if self.external_force: self.__pde.setValue(y=self.external_force)
          self.__pde.setValue(q=self.location_prescribed_velocity, \
                              r=Data())
          if not self.prescribed_velocity.isEmpty() and self.__iter==1:
               self.__pde.setValue(r=dt*self.prescribed_velocity)
          # solve the PDE:
          self.__solver_options.setTolerance(self.rel_tol**2)
          self.du=self.__pde.getSolution()
          # update geometry
          self.displacement=self.displacement+self.du
          self.domain.setX(self.__x+self.displacement)
          self.velocity=(self.displacement-self.__displacement_safe)/dt

          if self.debug:
             for i in range(self.domain.getDim()):
                self.trace("du %d range %e:%e"%(i,inf(self.du[i]),sup(self.du[i])))
             for i in range(self.domain.getDim()):
                self.trace("displacement %d range %e:%e"%(i,inf(self.displacement[i]),sup(self.displacement[i])))
             for i in range(self.domain.getDim()):
                self.trace("velocity %d range %e:%e"%(i,inf(self.velocity[i]),sup(self.velocity[i])))
          self.__stress_last=self.stress

      def terminateIteration(self):
          """iteration is terminateIterationd if relative pressure change is less than rel_tol"""
          if self.__iter>self.max_iter:
              raise IterationDivergenceError("Maximum number of iterations steps reached")
          if self.__iter==0:
             self.__diff=self.UNDEF_DT
          else:
             self.__diff,diff_safe=Lsup(self.stress-self.__stress_last),self.__diff
             s_sup=Lsup(self.stress)
             self.trace("stress max and increment :%e, %e"%(s_sup,self.__diff))
             if self.__iter>2 and diff_safe<self.__diff:
                 raise IterationDivergenceError("no improvement in stress iteration")
             return self.__diff<=self.rel_tol*self.SAFTY_FACTOR_ITERATION*s_sup+self.abs_tol

      def doStepPostprocessing(self,dt):
           """
           accept all the values:
           """
           self.__displacement_safe=self.displacement
           self.__temperature_safe=self.temperature
           self.__stress_safe=self.stress
           self.__velocity_safe=self.velocity

      def getSafeTimeStepSize(self,dt):
           """
           returns new step size
           """
           a=sup(length(self.velocity)/self.domain.getSize())
           if a>0:
              return 1./a
           else:
              return self.UNDEF_DT



class DruckerPrager(Mechanics):
      """

      """

      def __init__(self,**kwargs):
           """
           set up model
           """
           super(DruckerPrager, self).__init__(**kwargs)
           self.declareParameter(plastic_stress=0.,
                                 hardening=0.,
                                 friction_parameter=0.,
                                 dilatancy_parameter=0.,
                                 shear_length=1.e15)
      def doInitialization(self):
          """
          """
          super(DruckerPrager, self).doInitialization()
          self.__plastic_stress_safe=self.plastic_stress
          self.__shear_length_safe=self.shear_length
          self.__hardening_safe=self.hardening
          self.__chi_safe=0
          self.__tau_safe=0

      def doStepPreprocessing(self,dt):
          """
          """
          super(DruckerPrager, self).doStepPreprocessing(dt)
          # set initial guess for iteration:
          self.shear_length=self.__shear_length_safe
          self.plastic_stress=self.__plastic_stress_safe
          self.hardening=self.__hardening_safe
          self.__chi=self.__chi_safe
          self.__tau=self.__tau_safe

      def doStep(self,dt):
          # set new tangential operator:
          self.setTangentialTensor()
          # do the update step:
          super(DruckerPrager, self).doStep(dt)
          # update stresses:
          self.setStress()

      def doStepPostprocessing(self,dt):
          super(DruckerPrager, self).doStepPostprocessing(dt)
          self.__plastic_stress_safe=self.plastic_stress
          self.__shear_length_safe=self.shear_length
          self.__hardening_safe=self.hardening
          self.__chi_safe=self.__chi
          self.__tau_safe=self.__tau

      def setStress(self):
           d=self.domain.getDim()
           G=self.shear_modulus
           K=self.bulk_modulus
           alpha=self.friction_parameter
           beta=self.dilatancy_parameter
           h=self.hardening
           k3=kronecker(self.domain)
           # elastic trial stress:
           g=grad(self.du)
           D=symmetric(g)
           W=nonsymmetric(g)
           s_e=self.stress+K*self.deps_th*k3+ \
                      2*G*D+(K-2./3*G)*trace(D)*k3 \
                      +2*symmetric(matrix_mult(W,self.stress))
           p_e=-1./d*trace(s_e)
           s_e_dev=s_e+p_e*k3
           tau_e=sqrt(1./2*inner(s_e_dev,s_e_dev))
           # yield conditon for elastic trial stress:
           F=tau_e-alpha*p_e-self.shear_length
           self.__chi=whereNonNegative(F+(self.rel_tol*(self.SAFTY_FACTOR_ITERATION)**2)*self.shear_length)
           # plastic stress increment:
           l=self.__chi*F/(h+G+beta*K)
           self.__tau=tau_e-G*l
           self.stress=self.__tau/(tau_e+self.abs_tol*whereZero(tau_e,self.abs_tol))*s_e_dev-(p_e+l*beta*K)*k3
           self.plastic_stress=self.plastic_stress+l
           # update hardening
           self.hardening=(self.shear_length-self.__shear_length_safe)/(l+self.abs_tol*whereZero(l))

      def setTangentialTensor(self):
           d=self.domain.getDim()
           G=self.shear_modulus
           K=self.bulk_modulus
           alpha=self.friction_parameter
           beta=self.dilatancy_parameter
           tau_Y=self.shear_length
           chi=self.__chi
           tau=self.__tau
           h=self.hardening
           k3=kronecker(Function(self.domain))

           sXk3=outer(self.stress,k3)
           k3Xk3=outer(k3,k3)
           s_dev=self.stress-trace(self.stress)*(k3/d)
           tmp=G*s_dev/(tau+self.abs_tol*whereZero(tau,self.abs_tol))

           self.S=G*(swap_axes(k3Xk3,0,3)+swap_axes(k3Xk3,1,3)) \
                 + (K-2./3*G)*k3Xk3 \
                 + (sXk3-swap_axes(swap_axes(sXk3,1,2),2,3))   \
                 + 1./2*(swap_axes(swap_axes(sXk3,0,2),2,3)   \
                        -swap_axes(swap_axes(sXk3,0,3),2,3)   \
                        -swap_axes(sXk3,1,2)                  \
                        +swap_axes(sXk3,1,3)                ) \
                 - outer(chi/(h+G+alpha*beta*K)*(tmp+beta*K*k3),tmp+alpha*K*k3)
