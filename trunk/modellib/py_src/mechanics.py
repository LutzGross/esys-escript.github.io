# $Id:$

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""

from esys.escript import *
from esys.escript.modelframe import Model,IterationDivergenceError
from esys.escript.linearPDEs import LinearPDE

class Mechanics(Model):
      """
      base class for mechanics models in updated lagrangean framework

      @ivar domain: domain (in)
      @ivar internal_force: =Data()
      @ivar external_force: =Data()
      @ivar prescribed_velocity: =Data()
      @ivar location_fixed_velocity: =Data()
      @ivar temperature:  = None
      @ivar expansion_coefficient:  = 0.
      @ivar bulk_modulus: =1.
      @ivar shear_modulus: =1.
      @ivar rel_tol: =1.e-3
      @ivar abs_tol: =1.e-15
      @ivar max_iter: =10
      @ivar displacement: =None
      @ivar stress: =None
      @ivar velocity: =None
      """
      def __init__(self,debug=False):
         """
         set up the model
         
         @param debug: debug flag
         @type debug: C{bool}
         """
         super(Mechanics, self).__init__(self,debug=debug)
         self.declareParameter(domain=None, \
                               displacement=None, \
                               stress=None, \
                               velocity=None, \
                               internal_force=Data(), \
                               external_force=Data(), \
                               location_fixed_velocity=Data(), \
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
           # open PDE:
           self.stress_old=self.stress
           self.__velocity_old=self.velocity
           self.__temperature_old=self.temperature
           self.__displacement_old=self.displacement
           # get node cooridnates:
           self.__x=self.domain.getX()
           self.domain.setX(self.__x+self.displacement)
           # open PDE:
           self.__pde=LinearPDE(self.domain)
           self.__pde.setSymmetryOn()
      def doStepPreprocessing(self,dt):
            """
            step up pressure iteration

            if run within a time dependend problem extrapolation of pressure from previous time steps is used to
            get an initial guess (that needs some work!!!!!!!)
            """
            self.__iter=0
            self.__diff=self.UNDEF_DT
            # set initial values for the iteration:
            self.velocity=self.__velocity_old
            self.temperature=self.__temperature_old
            self.displacement=self.__displacement_old+self.velocity*dt

            # set temperature strain:
            if self.temperature:
               self.deps_therm=self.self.expansion_coefficient*(self.temperature-self.__temperature_old)
            else:
               self.deps_therm=0.

            # update geometry
            self.domain.setX(self.__x+self.displacement)
            self.stress=self.stress_old

      def doStep(self,dt):
          """

          performs an iteration step of the penalty method.
          IterationDivergenceError is raised if pressure error cannot be reduced or max_iter is reached.
     
          requires self.S to be set 
          updates the thermal stress increment

          """
          k3=kronecker(self.domain)
          # set new thermal stress increment
          if self.temperature:
             self.deps_therm=self.self.expansion_coefficient*(self.temperature-self.__temperature_old)
          else:
             self.deps_therm=0.
          # set PDE coefficients:
          self.__pde.setValue(A=self.S)
          self.__pde.checkSymmetry()
          self.__pde.setValue(X=-self.stress+self.stress_old-self.bulk_modulus*self.deps_therm*k3)
          if self.internal_force: self.__pde.setValue(Y=self.internal_force)
          if self.external_force: self.__pde.setValue(y=self.external_force)
          self.__pde.setValue(q=self.location_fixed_velocity)
          # solve the PDE:
          self.__pde.setTolerance(self.rel_tol**2)
          du=self.__pde.getSolution(verbose=True)
          self.velocity+=du/dt
          # update geometry
          self.displacement=self.__displacement_old+dt*self.velocity
          self.domain.setX(self.__x+self.displacement)

          self.__iter+=1
          if self.debug:
             for i in range(self.domain.getDim()):
                self.trace("velocity %d range %e:%e"%(i,inf(self.velocity[i]),sup(self.velocity[i])))
             for i in range(self.domain.getDim()):
                self.trace("displacement %d range %e:%e"%(i,inf(self.displacement[i]),sup(self.displacement[i])))

      def terminateIteration(self):
          """iteration is terminateIterationd if relative pressure change is less then rel_tol"""
          if self.__iter>2 and self.diff_last<self.diff:
              raise IterationDivergenceError,"no improvement in stress iteration"
          if self.__iter>self.max_iter:
              raise IterationDivergenceError,"Maximum number of iterations steps reached"
          return self.diff<=self.rel_tol*1.e-2*Lsup(self.stress)+self.abs_tol

      def doStepPostprocessing(self,dt):
           """
           accept all the values:
           """
           self.__displacement_old=self.displacement
           self.__velocity_old,self.__velocity_old=self.velocity,self.__velocity_old
           self.__temperature_old=self.temperature
           self.stress_old=self.stress

      def getSafeTimeStepSize(self,dt):
           """
           returns new step size
           """
           d=Lsup(self.velocity-self.__velocity_old)*Lsup(self.displacement)
           if d>0:
              return dt/d*self.rel_tol
           else:
              return self.UNDEF_DT



class DruckerPrager(Mechanics):
      """

      """

      def __init__(self,debug=False):
           """
           set up model
           """
           super(DruckerPrager, self).__init__(debug=debug)
           self.declareParameter(plastic_stress=0.,
                                 friction_parameter=0.,
                                 dilatancy_parameter=0.,
                                 shear_length=1.e15)
      def doInitialization(self):
          """
          """
          super(DruckerPrager, self).doInitialization()

      def doStepPreprocessing(self):
          """
          """
          super(DruckerPrager, self).doStepPreprocessing(dt)
          # set initial guess for iteration:
          self.shear_length=self.__shear_length_old
          self.plastic_stress=self.__plastic_stress_old
          # set initial stress:
          self.updateStress()

 

#================================
           
          # set new tangential operator:
          self.S=self.getTangentialTensor()
          # do the update step:

          super(DruckerPrager, self).doStep(dt)
          # update stresses:
          self.updateStress()

      def doStepPostprocessing(self,dt):
          super(DruckerPrager, self).doStepPostprocessing(dt)
          self.__plastic_stress_old=self.plastic_stress
          self.__shear_length_old=self.shear_length

      def setStress(self):
           G=self.shear_modulus
           K=self.bulk_modulus
           alpha=self.friction_parameter
           beta=self.dilatancy_parameter
           h=self.hardening
           k3=kronecker(self.domain)
           # elastic trial stress:
           g=grad(self.velocity)
           D=symmetric(g)
           W=nonsymmetric(g)
           s_e=self.stress_old+K*self.deps_therm*k3+ \
                      dt*(2*G*D+(K-2./3*G)*trace(D)*k3 \
                          +2*symmetric(matrix_mult(W,self.stress_old)))
           p_e=-1./3*trace(s_e)
           s_e_dev=s_e+p_e*k3
           tau_e=sqrt(1./2*inner(s_e_dev,s_e_dev))
           # yield conditon for elastic trial stress:
           F=tau_e-alpha*p_e-self.shear_length
           self.__chi=whereNonNegative(F)
           # plastic stress increment:
           l=self.__chi*F/(h+G+beta*K)
           self.__tau=tau_e-G*l
           self.stress=self.__tau/(tau_e+self.abs_tol*whereZero(tau_e,self.abs_tol))*s_e_dev+(p_e+l*beta*K)*k3
           self.plastic_stress=self.plastic_stress_old+l
           # update hardening
           self.hardening=(self.shear_length-self.__shear_length_old)/(l+self.abs_tol*whereZero(l))

      def setTangentialTensor(self):
           G=self.shear_modulus
           K=self.bulk_modulus
           alpha=self.friction_parameter
           beta=self.dilatancy_parameter
           tau_Y=self.shear_length
           chi=self.__chi
           h=self.hardening
           tmp=G*s_dev/(tau+self.abs_tol*whereZero(tau,self.abs_tol))

           k3=kronecker(Function(self.domain))
           s_dev=self.stress+trace(self.stress)*(k3/d)
           sXk3=outer(self.stress_old,k3)
           k3Xk3=outer(k3,k3)
           S=G*(swap_axes(k3Xk3,1,2)+swap_axes(k3Xk3,1,3)) \
               + (K-2./3*G)*k3Xk3 \
               + (sXk3-swap_axes(sXk3,1,3)) \
               + 1./2*(swap_axes(sXk3,0,3)-swap_axes(sXk3,1,2) \
                      -swap_axes(sXk3,1,3)+swap_axes(sXk3,0,2)) \
               - outer(chi/(h+G+alpha*beta*K)*(tmp+beta*K*k3),tmp+alpha*K*k3)
           return S
