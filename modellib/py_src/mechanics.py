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
      @ivar displacement: current displacements

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
                               prescribed_velocity=Data(), \
                               location_prescribed_velocity=Data(), \
                               temperature = None, \
                               expansion_coefficient = 0., \
                               bulk_modulus=1., \
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
           self.__pde=LinearPDE(self.domain)
           self.__displacement_old=self.displacement
           self.stress_old=self.stress
           self.__velocity_old=self.velocity
           self.__temperature_old=self.temperature
           
      def doStepPreprocessing(self,dt):
            """
            step up pressure iteration

            if run within a time dependend problem extrapolation of pressure from previous time steps is used to
            get an initial guess (that needs some work!!!!!!!)
            """
            self.__iter=0
            self.__diff=self.UNDEF_DT
            # set new values:
            self.displacement=self.__displacement_old
            self.stress=self.stress_old
            self.velocity=self.__velocity_old
            self.temperature=self.__temperature_old
            self.__velocity_last=self.velocity

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
             self.dthermal_stress=self.bulk_modulus*self.self.expansion_coefficient*(self.temperature-self.__temperature_old)
          else:
             self.dthermal_stress=0.
          # set PDE coefficients:
          self.__pde.setValue(A=self.S)
          self.__pde.setValue(X=self.stress_old-self.dthermal_stress*k3)
          if self.internal_force: self.__pde.setValue(Y=self.internal_force)
          if self.external_force: self.__pde.setValue(y=self.external_force)
          self.__pde.setValue(r=self.prescribed_velocity, \
                              q=self.location_prescribed_velocity)
          # solve the PDE:
          self.__pde.setTolerance(self.rel_tol/100.)
          self.velocity=self.__pde.getSolution()
          # calculate convergence indicators:
          self.__diff,diff_old=Lsup(self.velocity-self.__velocity_last),self.__diff
          self.__velocity_last=self.velocity
          self.displacement=self.__displacement_old+dt*self.velocity
          self.__iter+=1
          self.trace("velocity range %e:%e"%(inf(self.velocity),sup(self.velocity)))
          if self.__iter>2 and diff_old<self.__diff:
              raise IterationDivergenceError,"no improvement in stress iteration"
          if self.__iter>self.max_iter:
              raise IterationDivergenceError,"Maximum number of iterations steps reached"

      def terminateIteration(self):
          """iteration is terminateIterationd if relative pressure change is less then rel_tol"""
          return self.__diff<=self.rel_tol*Lsup(self.velocity)+self.abs_tol

      def doStepPostprocessing(self,dt):
           """
           accept all the values:
           """
           self.displacement=self.__displacement_old
           self.stress=self.stress_old
           self.velocity=self.__velocity_old
           self.temperature=self.__temperature_old

      def getSafeTimeStepSize(self,dt):
           """
           returns new step size
           """
           d=Lsup(self.velocity-self.__velocity_old)/dt
           if d>0:
               return Lsup(self.displacement)/d
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
                                 shear_length=0.)

      def doInitialization(self):
           """
           initialize model
           """
           super(DruckerPrager, self).doInitialization()
           self.__plastic_stress_old=self.plastic_stress
           self.__tau_y_old=self.shear_length

      def doStepPreprocessing(self,dt):
            """
            step up pressure iteration

            if run within a time dependend problem extrapolation of pressure from previous time steps is used to
            get an initial guess (that needs some work!!!!!!!)
            """
            super(DruckerPrager, self).doStepPreprocessing(dt)
            self.plastic_stress=self.__plastic_stress_old

      def doStep(self,dt):
           G=self.shear_modulus
           K=self.bulk_modulus
           alpha=self.friction_parameter
           beta=self.dilatancy_parameter
           tau_Y=self.shear_length
           if self.__plastic_stress_old:
              dps=self.plastic_stress-self.__plastic_stress_old
              h=(tau_Y-self.__tau_y_old)/(dps+self.abs_tol*whereZero(dps))
           else:
              h=0
           # set new tangential operator:
           self.S=self.getTangentialTensor(self.stress,
                                           tau_Y,G,K,alpha,beta,h)
           # do the update step:
           super(DruckerPrager, self).doStep(dt)

           # update stresses:
           self.stress,self.plastic_stress=self.getNewStress(self.stress_old,self.__plastic_stress_old,
                                                        self.velocity*dt,
                                                        self.dthermal_stress,tau_Y,G,K,alpha,beta,h)

      def doStepPostprocessing(self,dt):
          super(DruckerPrager, self).doStepPostprocessing(dt)
          self.plastic_stress=self.__plastic_stress_old=self.plastic_stress

      def getNewStress(self,s,gamma_p,du,ds_therm,tau_Y,G,K,alpha,beta,h):
            k3=kronecker(self.domain)
            dt=1.
            g=grad(du)
            D=symmetric(g)
            W=nonsymmetric(g)
            s_e=s+ds_therm+dt*(2*G*D+(K-2./3*G)*trace(D)*k3 \
                               +2*nonsymmetric(matrix_mult(W,s)))
            p_e=-1./3*trace(s_e)
            s_e_dev=s_e+p_e*k3
            tau_e=sqrt(1./2*inner(s_e_dev,s_e_dev))
            F=tau_e-alpha*p_e-tau_Y
            chi=whereNonNegative(F)
            l=chi*F/(h+G+beta*K)
            s=(1.-l*G/tau_e)*s_e_dev+(p_e+l*beta*K)*k3
            gamma_p=gamma_p+l
            return s, gamma_p
            

      def getTangentialTensor(self,s,tau_Y,G,K,alpha,beta,h):
           d=self.domain.getDim()
           k3=kronecker(Function(self.domain))
           p=-1./d*trace(s)
           s_dev=s+p*k3 
           tau=sqrt(1./2*inner(s_dev,s_dev))
           chi=whereNonNegative(tau-alpha*p-tau_Y)
           sXk3=outer(s,k3)
           k3Xk3=outer(k3,k3)
           tmp=G*s_dev/tau
           S=G*(swap_axes(k3Xk3,1,2)+swap_axes(k3Xk3,1,3)) \
               + (K-2./3*G)*k3Xk3 \
               + sXk3-swap_axes(sXk3,1,3) \
               + 1./2*(swap_axes(sXk3,0,3)+swap_axes(sXk3,1,2) \
                      -swap_axes(sXk3,1,3)-swap_axes(sXk3,0,2))
               # - chi/(h+G+alpha*beta*K)*outer(tmp+beta*K*k3,tmp+alpha*K*k3)\
           return S
