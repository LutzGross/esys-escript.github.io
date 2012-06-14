
########################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import Data, kronecker, Lsup, div, inf, sup
from esys.escript.modelframe import Model,IterationDivergenceError
from esys.escript.linearPDEs import LameEquation

class SteadyIncompressibleFlow(Model):
       """

       *-\left(\eta\left(v_{i,j}+v_{j,i}\right)\right)_{,j}+p_{,i}=F_i*
       
       *\sigma_{ij}=2\eta D_{ij}-p\,\delta_{ij}*

       *D_{ij}=\frac{1}{2}\left( v_{j,i} + v_{i,j }\right)*
       
       *v_{k,k} = 0*

       """

       def __init__(self,**kwargs):
           """
           set up model
           """
           Model.__init__(self,**kwargs)
           self.declareParameter(domain=None, \
                                 velocity=0., \
                                 pressure=0., \
                                 viscosity=1., \
                                 internal_force=0., \
                                 location_prescribed_velocity=None, \
                                 prescribed_velocity=None, \
                                 rel_tol=1.e-3,abs_tol=0.,max_iter=10,relaxation=0.0001)
           self.__iter=0

       def doInitialization(self):
           """
           initialize model
           """
           self.__p_old=None
           self.__p_very_old=None
           self.__dt_old=None
           self.__pde=LameEquation(self.domain)
           self.__pde.getSolverOptions().setSolverMethod(self.__pde.getSolverOptions().DIRECT)
           if self.location_prescribed_velocity == None: self.location_prescribed_velocit=Data()
           if self.prescribed_velocity == None: self.prescribed_velocity=Data()

       def stress(self):
           """
           returns current stress"""
           return 2*self.viscosity*self.stretching-self.pressure*kronecker(self.domain)

       def stretching(self):
           """
           returns stertching tensor
           """
           g=grad(self.velocity)
           return (g+transpose(g))/2

       def doStepPreprocessing(self,dt):
            """
            step up pressure iteration

            if run within a time dependend problem extrapolation of pressure from previous time steps is used to
            get an initial guess (that needs some work!!!!!!!)
            """
            self.__iter=0
            self.__diff=1.e40
            if not self.__p_old==None:
               if self.__p_very_old==None:
                  self.pressure=self.__p_old
               else:
                  self.pressure=(1.+dt/self.__dt_old)*self.__p_old-dt/self.__dt_old*self.__p_very_old

       def doStep(self,dt):
          """

          performs an iteration step of the penalty method.
          IterationDivergenceError is raised if pressure error cannot be reduced or max_iter is reached.

          """
          penalty=self.viscosity/self.relaxation
          self.__pde.setValue(lame_mu=self.viscosity, \
                              lame_lambda=penalty, \
                              F=self.internal_force, \
                              sigma=self.pressure*kronecker(self.__pde.getDomain()), \
                              r=self.prescribed_velocity, \
                              q=self.location_prescribed_velocity)
          self.__pde.getSolverOptions().setTolerance(self.rel_tol/10.)
          self.velocity=self.__pde.getSolution()
          update=penalty*div(self.velocity)
          self.pressure=self.pressure-update
          self.__diff,diff_old=Lsup(update),self.__diff
          self.__iter+=1
          self.trace("velocity range %e:%e"%(inf(self.velocity),sup(self.velocity)))
          self.trace("pressure range %e:%e"%(inf(self.pressure),sup(self.pressure)))
          self.trace("pressure correction: %e"%self.__diff)
          if self.__iter>2 and diff_old<self.__diff:
              self.trace("Pressure iteration failed!")
              raise IterationDivergenceError("no improvement in pressure iteration")
          if self.__iter>self.max_iter:
              raise IterationDivergenceError("Maximum number of iterations steps reached")

       def terminateIteration(self):
          """iteration is terminateIterationd if relative pressure change is less then rel_tol"""
          return self.__diff<=self.rel_tol*Lsup(self.pressure)+self.abs_tol

       def doStepPostprocessing(self,dt):
          self.__dt_old=dt
          self.__p_old,self.__p_very_old=self.pressure,self.__p_old
