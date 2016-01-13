
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript.modelframe import Model,ParameterSet
from esys.escript.util import exp
import numarray

class GravityForce(ParameterSet):
       """
       Sets a gravity force of given direction in given domain:

       @ivar domain: domain of interest (in).
       @type domain: L{esys.escript.Domain}
       @ivar density: density, default 1 (in).
       @ivar gravity: the gravity constant, default 9.81 (in).
       @ivar direction: the direction of gravity, default [1.,0.,0.] (in).
       """
       def __init__(self,**kwargs):
           """
           initializes the set
           """
           super(GravityForce, self).__init__(**kwargs)
           self.declareParameter(domain=None,
                                 gravity=9.81, \
                                 density=1., \
                                 direction=[1.,0.,0.])

       def gravity_force(self):
           """
           return the gravity force as L{density}*L{gravity}*L{direction}
           """
           if isinstance(self.direction,list): 
               dir=numarray.array(self.direction[:self.domain.getDim()])
           else:
               dir=self.direction[:self.domain.getDim()]
           return self.gravity*self.density*dir
     

    
class MaterialTable(ParameterSet):
       """
       A simple matrial table which allows setting physical ivar of a model

       @ivar density: density(in/out)
       @ivar heat_capacity:  heat_capacity(in/out)
       @ivar thermal_permabilty: permabilty(in/out)
       @ivar viscosity: viscosity(in/out)
       @ivar radiation_coefficient: (in/out)
       """
       def __init__(self,**kwargs):
           super(MaterialTable, self).__init__(**kwargs)
           self.declareParameter(gravity=9.81, \
                                 density=1., \
                                 heat_capacity=1., \
                                 thermal_permabilty=1., \
                                 radiation_coefficient=0.)

class SimpleEarthModel(ParameterSet):
       """
       A simple matrial table run convection models::

           density=density0*(1-expansion_coefficient*(temperature-reference_temperature))
           viscocity=viscocity0*(exp(alpha*(1/reference_temperature - 1/temperature))

       @ivar gravity: gravity constants (9.81) (in)
       @ivar reference_temperature: reference temperature (in)
       @ivar density0: density at reference temperature  (in)
       @ivar viscosity0: viscosity0 at reference temperature  (in)
       @ivar alpha: viscosity contrast  (in)
       @ivar expansion_coefficient: Raleigh number  (in)
       @ivar heat_capacity: heat capacity  (in)
       @ivar thermal_permabilty: permabilty  (in)
       @ivar temperature: temperature  (in)
       @ivar viscosity: viscosity  (out)
       @ivar density: density  (out)
       """
       def __init__(self,**kwargs):
           super(SimpleEarthModel, self).__init__(**kwargs)
           self.declareParameter(reference_temperature=1.,
                                 gravity=9.81, \
                                 density0=1., \
                                 expansion_coefficient=0., \
                                 viscosity0=1., \
                                 alpha=0., \
                                 temperature=0.,\
                                 heat_capacity=1., \
                                 thermal_permabilty=1.)
      
       def density(self):
           return self.density0*(1-self.expansion_coefficient*(self.temperature-self.reference_temperature))

       def viscosity(self):
           return self.viscosity0*exp(self.alpha*(1/self.reference_temperature - 1/(self.temperature+1.e-15)))

class SimpleSolidMaterial(MaterialTable):
       """
       A simple matrial table which allows setting physical parameters of 
       a model.

       @ivar density: density  (in/out)
       @ivar heat_capacity: heat_capacity  (in/out)
       @ivar thermal_permabilty: permabilty  (in/out)
       @ivar viscosity: viscosity  (in/out)
       @ivar radiation_coefficient:  (in/out)
       """
       def __init__(self,**kwargs):
           super(MaterialTable, self).__init__(**kwargs)
           self.declareParameter(lame_lambda=1.,\
                                 lame_my=1.)

class SimpleFluidMaterial(MaterialTable):
       """
       A simple matrial table which allows setting physical ivar of a model.

       @ivar density: density(in/out)
       @ivar heat_capacity: heat_capacity(in/out)
       @ivar thermal_permabilty: permabilty (in/out)
       @ivar viscosity: viscosity (in/out)
       @ivar radiation_coefficient:(in/out)
       """
       def __init__(self,**kwargs):
           super(MaterialTable, self).__init__(**kwargs)
           self.declareParameter(viscosity=1., \
                                 hydraulic_conductivity=1.e-4)

# vim: expandtab shiftwidth=4:
