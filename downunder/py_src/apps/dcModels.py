__copyright__ = "Copyright (c) 2020 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Lutz Gross, Andrea Codd"

from esys.escript import *
import esys.escript.unitsSI as U
import numpy as np
import argparse, sys, os
from esys.finley import ReadMesh
from esys.weipa import saveVTK, saveSilo
from esys.escript.pdetools import Locator
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions

class DCResistivityModel(object):
    """
    This class is a simple wrapper for 3D Dirac function direct current sources PDE model.  
    It solves PDE
        - div (sigma grad u) = sum (I_s d_dirac(x_s))
    where 
       I_s  is the applied cirrent at x_s and 
       we sum over all sources s. 
    Possible boundary conditions included in the class are Dirichlet on 
       - top (default) or  
       - top and base
       
    It solves just the one load condition, either with 
       - sources as list of points with corresponding list of currents (analytic primary solution)   
       - or sources as one function = sum (I_s d_dirac(x_s) (FE primary solution)

    It has a function to set ground property 
       - setConductivity
    get ground property
       - getConductivity
    
    It uses primary and secondary potential in solution
       - setPrimaryPotentialForHalfSpace (analytic solution - sources is list of points) 
       - setPrimaryPotential    (finite element computed solution - sources is a Dirac Delta)
       - getPrimaryPotential
       - getPrimaryField

    Output solution
       - getSecondaryPotential
       - getSecondaryField
       - getPotential
       - getField
    """

    def __init__(self, domain, sigma0=0.001, fixAllFaces=True, useFastSolver=False):
        """
        Initialise the class with domain and boundary conditions.
        Setup PDE and conductivity.
        :param domain: the domain 
        :type domain: `Domain`
        :param sigma0: background electric conductivity for the primary potential
        :type sigma0: typically `float`
        :param fixAllFaces: if true the potential at all faces except the top are set to zero. 
            Otherwise only the base is set to zero.
        :type fixAllFaces: `bool`
        :param useFastSolver: use multigrid solver. This may fail.  (Changing the mesh could stop this fail.) 
        :type useFastSolver: `bool`
        """
        self.domain = domain
        self.sigma0 = sigma0
        self.fixAllFaces = fixAllFaces
        self.useFastSolver = useFastSolver
        self.primary_potential = None
        self.primary_E = None
        self.secondary_potential = None
        self.sigma = None
        self.pde=self.__createPDE()

    def __createPDE(self):
        """
        Create the PDE and set boundary conditions.
        """
        pde = LinearSinglePDE(self.domain, isComplex=False)
        optionsG = pde.getSolverOptions()
        optionsG.setSolverMethod(SolverOptions.PCG)
        pde.setSymmetryOn()
        if self.useFastSolver and hasFeature('trilinos'):
            optionsG.setPackage(SolverOptions.TRILINOS)
            optionsG.setPreconditioner(SolverOptions.AMG)
        if self.fixAllFaces:
            x=pde.getDomain().getX()[0]
            y=pde.getDomain().getX()[1]
            z=pde.getDomain().getX()[2]
            pde.setValue(q=whereZero(x-inf(x))+whereZero(x-sup(x))+ whereZero(y-inf(y))+whereZero(y-sup(y))+whereZero(z-inf(z)))
        else:
            z=pde.getDomain().getX()[2]
            pde.setValue(q=whereZero(z-inf(z)))
        return pde
        
    def setPrimaryPotentialForHalfSpace(self, sources= [], charges=[]):
        """
        sets the primary potential for the infinite half space with constant conductivity `sigma0`
        uses analytic solution.  Method of images used for sources at depth)  
        :param sources: list of source locations (x,y,z). Need to be defined on the surface of the domain.
        :type sources: list of tuples.
        :param charges: list of charges (in [A])
        :type sources: list of floats        
        """
        assert len(sources) == len(charges)
        surf = sup(self.domain.getX()[2])
        x=Function(self.domain).getX()
        E=Vector(0., x.getFunctionSpace())
        u=Scalar(0., x.getFunctionSpace())
        # u=I/(2*pi*sigma*r)
        for s,I in zip(sources, charges):
            if s[2] == surf:  
                r=length(x-s)
                u+= I/(2*np.pi*self.sigma0*r)
                E+= I/(2*np.pi*self.sigma0*r**3)*(x-s) 
            else:
                s2 = (s[0],s[1],2*surf-s[2])
                r1 = length(x-s)
                r2 = length(x-s2)
                u+= I/(4*np.pi*self.sigma0*r1)+I/(4*np.pi*self.sigma0*r2)
                E+= I/(4*np.pi*self.sigma0*r1**3)*(x-s)+I/(4*np.pi*self.sigma0*r2**3)*(x-s2) 
        self.primary_potential  = u
        self.primary_E  = E
        return self


    
    def setPrimaryPotential(self, source = None):
        """
        set the primary potential using the source term `source`.
        :param source: source of charges as `DiracDeltaFunctions` object
        :type source: `Data` object
        """
        assert source.getFunctionSpace() == DiracDeltaFunctions(self.domain)
        
        self.pde.setValue(A=self.sigma0*kronecker(self.domain), y_dirac=source, X=Data())
        u=self.pde.getSolution()
        self.primary_potential=interpolate(u, Function(self.domain))
        self.primary_E  = -grad(u, Function(self.domain))
        return self
            
    def getPrimaryPotential(self):
        """
        returns the primary potential
        """
        if self.primary_potential is None:
            raise ValueError("no primary potential set.")

        return self.primary_potential

    def getPrimaryField(self):
        """
        returns the primary electric field. This is -grad(getPrimaryPotential())
        """
        if self.primary_E is None:
            raise ValueError("no primary field set.")

        return self.primary_E
    
    def setConductivity(self, sigma):
        """
        sets the conductivity. This solves a `LinearSinglePDE`.
        :param sigma: conductivity distribution. The value at source locations should be sigma0.
        :type sigma: `Data`
        """
        self.sigma=interpolate(sigma, Function(self.domain))
        self.pde.setValue(A=self.sigma*kronecker(self.domain.getDim()),
                          X=(self.sigma-self.sigma0)*self.getPrimaryField(), y_dirac=Data())
        self.secondary_potential=self.pde.getSolution()
        return self

    def getConductivity(self):
        """
        returns the conductivity
        """
        return self.sigma
    
    def getSecondaryPotential(self):
        """
        returns the secondary potential
        """
        if self.secondary_potential is None:
            raise ValueError("no secondary potential set.")

        return self.secondary_potential

    def getSecondaryField(self):
        """
        returns secondary electric field. This is -grad(getSecondaryPotential())
        """
        if self.secondary_potential is None:
            raise ValueError("no secondary potential set.")

        return -grad(self.secondary_potential)

    
    def getPotential(self):
        """
        returns total potential
        """
        return self.getPrimaryPotential()+self.getSecondaryPotential()

    def getField(self):
        """
        returns the total potential. This is -grad(getPotential())
        """
        return self.getPrimaryField()+self.getSecondaryField()  
    
class DCResistivityModelNoPrimary(object):
    """
    This class is a simple wrapper for 3D Dirac function direct current sources PDE model.  
    It solves PDE
        - div (sigma grad u) = sum (I_s d_dirac(x_s))
    where 
       I_s  is the applied cirrent at x_s and 
       we sum over all sources s. 
    Possible boundary conditions included in the class are Dirichlet on 
       - top (default) or  
       - top and base
       
    It solves just the one load condition, either with 
       - sources as one function = sum (I_s d_dirac(x_s) (FE primary solution)

    It has a function to set ground property 
       - setConductivity
    get ground property
       - getConductivity
    
    It just solves PDE without splitting into primary and secondary.
    
    Output solution
       - getSecondaryPotential
       - getSecondaryField
       - getPotential
       - getField
    """

    def __init__(self, domain, source, sigma=0.001, fixAllFaces=True, useFastSolver=False):
        """
        Initialise the class with domain and boundary conditions.
        Setup PDE and conductivity.
        :param domain: the domain 
        :type domain: `Domain`
        :param sigma0: background electric conductivity for the primary potential
        :type sigma0: typically `float`
        :param fixAllFaces: if true the potential at all faces except the top are set to zero. 
            Otherwise only the base is set to zero.
        :type fixAllFaces: `bool`
        :param useFastSolver: use multigrid solver. This may fail.  (Changing the mesh could stop this fail.) 
        :type useFastSolver: `bool`
        """
        self.domain = domain
        self.fixAllFaces = fixAllFaces
        self.useFastSolver = useFastSolver
        self.sigma = sigma
        self.potential = None
        self.source=source
        self.pde=self.__createPDE()


    def __createPDE(self):
        """
        Create the PDE and set boundary conditions.
        """
        pde = LinearSinglePDE(self.domain, isComplex=False)
        optionsG = pde.getSolverOptions()
        optionsG.setSolverMethod(SolverOptions.PCG)
        pde.setSymmetryOn()
        assert self.source.getFunctionSpace() == DiracDeltaFunctions(self.domain)
        sigma=interpolate(self.sigma, Function(self.domain))
        pde.setValue(A=sigma*kronecker(self.domain), y_dirac=self.source)
        if self.useFastSolver and hasFeature('trilinos'):
            optionsG.setPackage(SolverOptions.TRILINOS)
            optionsG.setPreconditioner(SolverOptions.AMG)
        if self.fixAllFaces:
            x=pde.getDomain().getX()[0]
            y=pde.getDomain().getX()[1]
            z=pde.getDomain().getX()[2]
            pde.setValue(q=whereZero(x-inf(x))+whereZero(x-sup(x))+ whereZero(y-inf(y))+whereZero(y-sup(y))+whereZero(z-inf(z)))
        else:
            z=pde.getDomain().getX()[2]
            pde.setValue(q=whereZero(z-inf(z)))
        return pde
        
    def getPotential(self, source = None):
        """
        get the potential using the source term `source`.
        :param source: source of charges as `DiracDeltaFunctions` object
        :type source: `Data` object
        """
        u = self.pde.getSolution()
        self.potential= u  #interpolate(u, Function(self.domain))
        return self.potential
            
    def getPrimaryField(self):
        """
        returns the primary electric field. This is -grad(getPrimaryPotential())
        """
        if self.potential is None:
            raise ValueError("no potential yet.")
        self.field  = -grad(self.potential, Function(self.domain))
        return self.field
    
    def setConductivity(self,sig):
        """
        returns the conductivity
        """
        self.sigma=sig        
        sigma=interpolate(self.sigma, Function(self.domain))
        self.pde.setValue(A=sigma*kronecker(self.domain.getDim()))
        return self
    
    def getConductivity(self):
        """
        returns the conductivity
        """
        return self.sigma
    
   
 
