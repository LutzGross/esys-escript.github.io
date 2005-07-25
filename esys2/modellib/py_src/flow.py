# $Id$


from esys.escript import *
from esys.modelframe import Model,IterationDivergenceError
from esys.linearPDEs import LameEquation

class SteadyIncompressibleFlow(Model):
       """ 

                   \f[
                   -\left(\eta\left(v_{i,j}+v_{j,i}\right)\right)_{,j}+p_{,i}=F_i
                   \f]

                   \f[
                   \sigma_{ij}=2\eta D_{ij}-p\,\delta_{ij}
                   \f[
                   D_{ij}=\frac{1}{2}\left( v_{j,i} + v_{i,j }\right) \; .
                   \f]

                   \f[
                   -v_{k,k} = 0 \; .
                   \f]
                   
       """

       def __init__(self,debug=False):
           """set up model"""
           Model.__init__(self,debug=debug)
           self.declareParameter(domain=None, \
                                 velocity=0., \
                                 pressure=0., \
                                 viscocity=1., \
                                 internal_force=0., \
                                 location_prescribed_velocities=Data(), \
                                 prescribed_velocities=Data(), \
                                 rel_tol=1.e-8,abs_tol=1.e-15,max_iter=10,relaxation=0.01) 
           self.__iter=0

       def doInitialization(self,t):
           """initialize model"""
           self.__pressure_history=None
           self.__dt_history=None
           self.__pde=LameEquation(self.domain)

       def stress(self):
           """returns current stress"""
           return 2*self.viscocity*self.stretching-self.pressure

       def stretching(self):
           """returns stertching tensor"""
           g=grad(self.velocity)
           return (g+transpose(g))/2
           
       def doIterationInitialization(self,dt):
            """
               step up pressure iteration

               if run within a time dependend problem extrapolation of pressure from previous time steps is used to
               get an initial guess 
            """
            if self.__dt_history==None:
               self.__pressure_history=self.pressure
               self.__iter=0
            else:
               self.__pressure_history,self.pressure=self.pressure,(1.+dt/self.__dt_history)*self.pressure+dt/self.__dt_history*self.__pressure_history
               self.__iter=1
            self.diff=1.e400

       def doIterationStep(self,dt):
          """

             performs an iteration step of the penalty method.
             IterationDivergenceError is raised if pressure error cannot be reduced or max_iter is reached.

          """
          penalty=self.viscosity/self.relaxation
          self.__pde.setValue(lame_lambda=self.viscosity, \
                              lame_my=penalty, \
                              F=F, \
                              sigma0=self.pressure*kronecker(self.__pde.getDomain()), \
                              r=self.location_prescribed_velocities, \
                              q=self.location_prescribed_velocities) 
          self.__pde.setTolerance(self.tol*1.e-2)
          self.velocity=self.pde.getSolution()
          update=penalty*div(self.velocity)
          self.pressure=self.pressure-update
          self.diff,diff=Lsup(update),self.diff
          print "Pressure iteration: step %d: correction %e"%(self.__iter,self.diff/Lsup(self.pressure))
          if diff<=self.diff:
              raise IterationDivergenceError,"no improvement in pressure iteration"

       def terminate(self):
          """iteration is terminated if relative pressure change is less then rel_tol"""
          if self.iter<2:
              return False
          else:
             return self.diff<self.rel_tol*Lsup(self.pressure)+self.abs_tol

       def doIterationFinalization(self,dt):
          self.__dt_history=dt
if __name__=="__main__":



   from esys.modelframe import Link,Simulation,ExplicitSimulation
   from esys.geometry import RectangularDomain,VectorConstrainer
   from esys.probe import Probe
   from esys.input import InterpolatedTimeProfile,GausseanProfile

   dom=RectangularDomain()
   constraints=VectorConstrainer()
   constraints.domain=Link(dom)
   constraints.left=[1,1,1]
   constraints.right=[1,1,1]
   constraints.top=[1,1,1]
   constraints.bottom=[1,1,1]
   constraints.front=[1,1,1]
   constraints.back=[1,1,1]

   flow=SteadyIncompressibleFlow()
   flow.domain=Link(dom)
   flow.velocity=0.
   flow.pressure=0. 
   flow.viscocity=1.
   flow.internal_force=[1.,1.,0.]
   flow.location_prescribed_velocities=Link(constraints,"location_of_constraint")
   flow.prescribed_velocities=Data(), \

   ptest=Probe()
   ptest.reference("x[0]+x[1]")
   ptest.value=Link(flow,"pressure")

   s=ExplicitSimulation([dom,constraints,flow,ptest])
   s.writeXML()
   s.run()
