__copyright__ = "Copyright (c) 2020 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Lutz Gross, Andrea Codd"


from esys.escript import *
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
import esys.escript.unitsSI as U
import numpy as np

class GravityModel(object):
    """
    This class is a simple wrapper for a 2D or 3D gravity PDE model.  
    It solves PDE
        - div (grad u) = -4 pi G rho
    where 
       G  is the gravitational constant and 
       rho is density
       u  is computed anomaly potential

    Possible boundary conditions included in the class are Dirichlet on 
       - top (default) or  
       - top and base.
       
    It has a function to set ground property 
       - setDensity
    get ground property
       - getDensity
    and functions to output solutions
       - getgravityPotential : u
       - getzGravity : -grad(u)[2] (3D) or -grad(u)[1] (2D)
       - getGravityVector : grad (u) .
    """
    
    def __init__(self, domain, fixBase = False):
        """
        Initialise the class with domain and boundary conditions.
        Setup PDE and density.
        : param domain: the domain 
        : type domain: `Domain`
        : param fixBase: if true the gravitational potential at the bottom is set to zero.
        : type fixBase: `bool`
        : param fixVert: if true the magnetic field on all vertical sudes is set to zero.
        : type fixBase: `bool`
        : if fixBase is True then gravity field is set to zero at base and top surfaces.
        """
        self.domain = domain
        self.fixBase=fixBase
        assert self.domain.getDim() == 2 or self.domain.getDim() == 3
        self.pde=self.__createPDE()
        self.setDensity()

    def __createPDE(self):
        """
        Create the PDE and set boundary conditions.
        """
        pde = LinearSinglePDE(self.domain, isComplex=False)
        domdim = self.domain.getDim() 
        zdim = domdim - 1       
        optionsG = pde.getSolverOptions()
        optionsG.setSolverMethod(SolverOptions.PCG)
        pde.setSymmetryOn()
        pde.setValue(A = kronecker(domdim))
        x = self.domain.getX()
        q = whereZero(x[zdim]-sup(x[zdim]))
        if self.fixBase:
            q += whereZero(x[zdim]-inf(x[zdim])) 
        pde.setValue(q = q)
        if hasFeature('trilinos'):
            optionsG.setPackage(SolverOptions.TRILINOS)
            optionsG.setPreconditioner(SolverOptions.AMG)
        return pde
             
    def setDensity(self, rho=0):
        """
        set density
        : param rho: density
        : type rho: `Data` or `float` 
        """
        self.rho = rho
        self.reset = True

    def getDensity(self):
        """
        returns density
        : returns: rho
        """
        return self.rho
        
    def getGravityPotential(self):
        """
        get the potential of the density anomaly
        """
        if self.reset:
            self.pde.setValue(Y= -4.0*np.pi*U.Gravitational_Constant*self.rho)
        return self.pde.getSolution()

    def getzGravity(self):
        """
        get Bouger gravity in -z direction (vertical)
        """
        zdim= self.domain.getDim() -1
        return -self.getGravityVector()[zdim] 

    def getGravityVector(self):
        """
        get the Bouger gravity vector
        """
        return grad(self.getGravityPotential(), ReducedFunction(self.pde.getDomain()))

