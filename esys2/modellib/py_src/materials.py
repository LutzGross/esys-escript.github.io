# $Id$

from escript.modelframe import Model,ParameterSet
from escript.util import exp
import numarray

class GravityForce(ParameterSet):
       """@brief sets a gravity force of given direction in given domain:

              @param domain (in) - domain of interest
              @param density (in) - density
              @param direction (in) - density
              @param gravity_force(out)  - gravity force

       """
       def __init__(self,debug=False):
           ParameterSet.__init__(self,debug=debug)
           self.declareParameter(domain=None,
                                 gravity=9.81, \
                                 density=1., \
                                 direction=[1.,0.,0.])

       def gravity_force(self):
            if isinstance(self.direction,list): 
               dir=numarray.array(self.direction[:self.domain.getDim()])
            else:
               dir=self.direction[:self.domain.getDim()]
            return self.gravity*self.density*dir
     

    
class MaterialTable(ParameterSet):
       """@brief a simple matrial table which allows setting physical parameters of a model

              @param density (in/out) - density
              @param heat_capacity (in/out)  - heat_capacity
              @param thermal_permabilty (in/out) - permabilty
              @param viscosity (in/out)  - viscosity
              @param radiation_coefficient (in/out) -

       """
       def __init__(self,debug=False):
           ParameterSet.__init__(self,debug=debug)
           self.declareParameter(gravity=9.81, \
                                 density=1., \
                                 heat_capacity=1., \
                                 thermal_permabilty=1., \
                                 radiation_coefficient=0.)

class SimpleEarthModel(ParameterSet):
       """@brief a simple matrial table run convection models:

              density=density0*(1-rayleigh_number*(temperature-reference_temperature))
              viscocity=viscocity0*(exp(alpha*(1/reference_temperature - 1/temperature))

              @param gravity  (in) - gravity constants (9.81)
              @param reference_temperature  (in) - reference temperature
              @param density0  (in) - density at reference temperature
              @param viscosity0  (in) - viscosity0 at reference temperature
              @param alpha  (in) - viscosity contrast
              @param rayleigh_number (in)  - Raleigh number
              @param heat_capacity (in)  - heat capacity
              @param thermal_permabilty (in) - permabilty
              @param temperature (in) - temperature
              @param viscosity (out)  - viscosity
              @param density (out)  - density

       """
       def __init__(self,debug=False):
           ParameterSet.__init__(self,debug=debug)
           self.declareParameter(reference_temperature=1.,
                                 gravity=9.81, \
                                 density0=1., \
                                 rayleigh_number=0., \
                                 viscosity0=1., \
                                 alpha=0., \
                                 heat_capacity=1., \
                                 thermal_permabilty=1.)
      
       def density(self):
           return self.density0*(1-self.rayleigh_number*(self.temperature-self.reference_temperature))

       def viscosity(self):
           return self.viscosity0*exp(self.alpha*(1/self.reference_temperature - 1/(self.temperature+1.e-15)))

class SimpleSolidMaterial(MaterialTable):
       """@brief a simple matrial table which allows setting physical parameters of a model

              @param density (in/out) - density
              @param heat_capacity (in/out)  - heat_capacity
              @param thermal_permabilty (in/out) - permabilty
              @param viscosity (in/out)  - viscosity
              @param radiation_coefficient (in/out) -

       """
       def __init__(self,debug=False):
           MaterialTable.__init__(self,debug=debug)
           self.declareParameter(lame_lambda=1.,\
                                 lame_my=1.)

class SimpleFluidMaterial(MaterialTable):
       """@brief a simple matrial table which allows setting physical parameters of a model

              @param density (in/out) - density
              @param heat_capacity (in/out)  - heat_capacity
              @param thermal_permabilty (in/out) - permabilty
              @param viscosity (in/out)  - viscosity
              @param radiation_coefficient (in/out) -

       """
       def __init__(self,debug=False):
           MaterialTable.__init__(self,debug=debug)
           self.declareParameter(viscosity=1., \
                                 hydraulic_conductivity=1.e-4)

