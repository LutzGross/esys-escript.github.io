#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""

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

       @ivar density (in/out): density
       @ivar heat_capacity (in/out):  heat_capacity
       @ivar thermal_permabilty (in/out): permabilty
       @ivar viscosity (in/out): viscosity
       @ivar radiation_coefficient (in/out):
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
       B simple matrial table run convection models::

           density=density0*(1-expansion_coefficient*(temperature-reference_temperature))
           viscocity=viscocity0*(exp(alpha*(1/reference_temperature - 1/temperature))

       @ivar gravity (in): gravity constants (9.81)
       @ivar reference_temperature (in): reference temperature
       @ivar density0 (in): density at reference temperature
       @ivar viscosity0 (in): viscosity0 at reference temperature
       @ivar alpha (in): viscosity contrast
       @ivar expansion_coefficient (in): Raleigh number
       @ivar heat_capacity (in): heat capacity
       @ivar thermal_permabilty (in): permabilty
       @ivar temperature (in): temperature
       @ivar viscosity (out): viscosity
       @ivar density (out): density
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

       @ivar density (in/out): density
       @ivar heat_capacity (in/out): heat_capacity
       @ivar thermal_permabilty (in/out): permabilty
       @ivar viscosity (in/out): viscosity
       @ivar radiation_coefficient (in/out):
       """
       def __init__(self,**kwargs):
           super(MaterialTable, self).__init__(**kwargs)
           self.declareParameter(lame_lambda=1.,\
                                 lame_my=1.)

class SimpleFluidMaterial(MaterialTable):
       """
       A simple matrial table which allows setting physical ivar of a model.

       @ivar density (in/out): density
       @ivar heat_capacity (in/out): heat_capacity
       @ivar thermal_permabilty (in/out): permabilty
       @ivar viscosity (in/out): viscosity
       @ivar radiation_coefficient (in/out):
       """
       def __init__(self,**kwargs):
           super(MaterialTable, self).__init__(**kwargs)
           self.declareParameter(viscosity=1., \
                                 hydraulic_conductivity=1.e-4)

# vim: expandtab shiftwidth=4:
