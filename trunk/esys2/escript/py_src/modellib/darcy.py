# $Id$


from esys.escript import *
from esys.modelframe import Model,IterationDivergenceError
from esys.linearPDEs import AdvectivePDE,LinearPDE

class DarcyFlow(Model):
       """ """

       def __init__(self,debug=False):
           Model.__init__(self,debug=debug)
           self.declareParameter(domain=None, \
                                 tend=0., \
                                 dt=0.1, \
                                 temperature=1., \
                                 density=1., \
                                 c_p=1., \
                                 thermal_permabilty=1., \
                                 reference_temperature=0., \
                                 radiation_coefficient=0., \
                                 thermal_source=0., \
                                 location_fixed_temperature=Data(), \
                                 iterate=True, \
                                 tol=1.e-8, \
                                 implicit=True)
           self.iter=0

       def doInitialization(self,t):
           self.tn=t
           self.pde=LinearPDE(self.domain)
           self.pde.setSymmetryOn()

       def doIterationInitialization(self,dt):
            self.iter=0
            self.T_last=self.temperature
            self.diff=1.e400

       def doIterationStep(self,dt):
          T=self.temperature
          diff=self.diff
          dim=self.pde.getDim()
          self.iter+=1
          rhocp=self.density*self.c_p
          self.pde.setValue(A=self.thermal_permabilty*kronecker(dim), \
                              D=rhocp/dt, \
                              Y=self.thermal_source+rhocp/dt*self.T_last, \
                              d=self.radiation_coefficient, \
                              y=self.radiation_coefficient*self.reference_temperature, \
                              q=self.location_fixed_temperature, \
                              r=self.T_last)
          if isinstance(self,TemperatureAdvection): self.pde.setValue(C=self.velocity[:dim]*rhocp)
          self.pde.setTolerance(self.tol*1.e-2)
          self.temperature=self.pde.getSolution()
          self.diff=Lsup(T-self.temperature)
          if diff<=self.diff:
              raise IterationDivergenceError,"no improvement in the temperature iteration"

       def terminate(self):
          if not self.implicit:
              return True
          elif self.iter<1:
              return False
          else:
             return self.diff<self.tol*Lsup(self.temperature)

       def doIterationFinalization(self,dt):
          self.tn+=dt

       def getSafeTimeStepSize(self,dt):
           return self.dt

       def finalize(self):
            return self.tn>=self.tend

if __name__=="__main__":
   from esys.modelframe import Link,Simulation,ExplicitSimulation
   from esys.visualization import WriteVTK
   from esys.materials import SimpleMaterialTable
   from esys.geometry import RectangularDomain,ScalarConstrainer
   from esys.input import InterpolatedTimeProfile,GausseanProfile

   dom=RectangularDomain()
   constraints=ScalarConstrainer()
   constraints.domain=Link(dom)
   constraints.top=1
   constraints.bottom=1
   
   mt=MaterialTable()
   
   pf=InterpolatedTimeProfile()
   pf.t=[0.,0.25,0.5,0.75]
   pf.values=[0.,1.,1.,0.]
   
   q=GausseanProfile()
   q.domain=Link(dom)
   q.width=0.05
   q.x_c=numarray.array([0.5,0.5,0.5])
   q.r=0.01
   q.A=Link(pf,"out")
   
   tt=DarcyFlow()
   tt.domain=Link(dom)
   tt.tend=1.
   tt.dt=0.1
   tt.temperature=0.
   tt.density=Link(mt)
   tt.c_p=Link(mt)
   tt.thermal_permabilty=Link(mt)
   tt.reference_temperature=0.
   tt.radiation_coefficient=Link(mt)
   tt.thermal_source=Link(q,"out")
   tt.location_fixed_temperature=Link(constraints,"location_of_constraint")
   tt.implicit=True
   
   vis=WriteVTK()
   vis.scalar=Link(tt,"temperature")
   
   s=ExplicitSimulation([dom,constraints,pf,q,Simulation([mt,tt],debug=True),vis],debug=True)
   # s=Simulation([dom,constraints,pf,q,Simulation([mt,tt]),vis])
   s.writeXML()
   s.run()
