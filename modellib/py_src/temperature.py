
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

from __future__ import division, print_function

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import Data, inf, sup, length, grad, inner
from esys.escript.modelframe import Model,IterationDivergenceError
import numpy
import esys.escript.linearPDEs as lpde


class TemperatureAdvection(Model):
       """

       The conservation of internal heat energy is given by

       *rho c_p ( dT/dt+v[j] * grad(T)[j])-grad(\kappa grad(T)_{,i}=Q*

       *n_i \kappa T_{,i}=0*

       it is assummed that *\rho c_p* is constant in time.

       solved by Taylor Galerkin method

       """
       def __init__(self,**kwargs):
           super(TemperatureAdvection, self).__init__(**kwargs)
           self.declareParameter(domain=None, \
                                 temperature=1., \
                                 velocity=numpy.zeros([3]),
                                 density=1., \
                                 heat_capacity=1., \
                                 thermal_permabilty=1., \
                                 # reference_temperature=0., \
                                 # radiation_coefficient=0., \
                                 thermal_source=0., \
                                 fixed_temperature=0.,
                                 location_fixed_temperature=Data(),
                                 safety_factor=0.1)

       def doInitialization(self):
           self.__pde=lpde.LinearPDE(self.domain)
           self.__pde.setSymmetryOn()
           self.__pde.setReducedOrderOn()
           self.__pde.getSolverOptions().setSolverMethod(lpde.SolverOptions.LUMPING)
           self.__pde.setValue(D=self.heat_capacity*self.density)

       def getSafeTimeStepSize(self,dt):
           """
           returns new step size
           """
           h=self.domain.getSize()
           return self.safety_factor*inf(h**2/(h*abs(self.heat_capacity*self.density)*length(self.velocity)+self.thermal_permabilty))

       def G(self,T,alpha):
           """
           tangential operator for taylor galerikin
           """
           g=grad(T)
           self.__pde.setValue(X=-self.thermal_permabilty*g, \
                               Y=self.thermal_source-self.__rhocp*inner(self.velocity,g), \
                               r=(self.__fixed_T-self.temperature)*alpha,\
                               q=self.location_fixed_temperature)
           return self.__pde.getSolution()
           

       def doStepPostprocessing(self,dt):
           """
           perform taylor galerkin step
           """
           T=self.temperature
           self.__rhocp=self.heat_capacity*self.density
           self.__fixed_T=self.fixed_temperature
           self.temperature=dt*self.G(dt/2*self.G(T,1./dt)+T,1./dt)+T
           self.trace("Temperature range is %e %e"%(inf(self.temperature),sup(self.temperature)))
