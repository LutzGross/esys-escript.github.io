# $Id$


from esys.modelframe import Model

class MaterialTable(Model):
       """@brief a simple matrial table which allows setting physical parameters of a model

              @param density (in/out) - density
              @param c_p (in/out)  - c_p
              @param thermal_permabilty (in/out) - permabilty
              @param viscosity (in/out)  - viscosity
              @param radiation_coefficient (in/out) -

       """
       def __init__(self,debug=False):
           Model.__init__(self,debug=debug)
           self.declareParameter(g=9.81, \
                                 density=1., \
                                 c_p=1., \
                                 thermal_permabilty=1., \
                                 radiation_coefficient=0.)

class SimpleSolidMaterial(MaterialTable):
       """@brief a simple matrial table which allows setting physical parameters of a model

              @param density (in/out) - density
              @param c_p (in/out)  - c_p
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
              @param c_p (in/out)  - c_p
              @param thermal_permabilty (in/out) - permabilty
              @param viscosity (in/out)  - viscosity
              @param radiation_coefficient (in/out) -

       """
       def __init__(self,debug=False):
           MaterialTable.__init__(self,debug=debug)
           self.declareParameter(viscosity=1., \
                                 hydraulic_conductivity=1.e-4)

