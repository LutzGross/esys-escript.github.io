__copyright__ = "Copyright (c) 2020 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Lutz Gross, Andrea Codd"

from esys.escript import *
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions

class MagneticModel2D(object):
    """
    This class is a simple wrapper for a 2D magnetic PDE model.  
    It solves PDE
        - div (grad u) = -div(k Bh)
    where 
       k  is magnetic susceptibility and 
       Bh is background magnetic field.
       u  is computed anomaly potential

    Possible boundary conditions are Dirichlet on 
       - one corner (bottom left),
       - vertical sides, or 
       - base.

    Requires domain and possibly boundary condition choice.  
    
    Default boundary conditions 
       - fix the left bottom corner.  
    Otherwise add 
       - fixVert = True (for all vertical surfaces fixed) or 
       - fixBase = True (for base fixed).  

    It has functions 
       - setSusceptibility
       - getSusceptibility
       - setBackgroundMagneticField
       - getAnomalyPotential
       - getMagneticFieldAnomaly.
    """

    def __init__(self, domain, fixVert=False, fixBase=False):
        """
        Initialise the class with domain and boundary conditions.
        Setup PDE, susceptibility and background magnetic field.
        : param domain: the domain 
        : type domain: `Domain`
        : param fixBase: if true the magnetic field at the bottom is set to zero. 
        : type fixBase: `bool`
        : param fixVert: if true the magnetic field on all vertical sudes is set to zero.
        : type fixVert: `bool`
        : if fixBase and fixVert are False then magnetic field is set to zero at bottom, front, left corner.
        """
        self.domain  = domain
        self.fixVert = fixVert
        self.fixBase = fixBase
        assert self.domain.getDim() == 2
        self.pde=self.__createPDE()
        self.setSusceptibility()
        self.setBackgroundMagneticField()

    def __createPDE(self):
        """
        Create the PDE and set boundary conditions.
        """
        pde=LinearSinglePDE(self.domain, isComplex=False)
        optionsG=pde.getSolverOptions()
        optionsG.setSolverMethod(SolverOptions.PCG)
        pde.setSymmetryOn()
        pde.setValue(A=kronecker(self.domain.getDim()))
        x=self.domain.getX()
        pde.setValue(A=kronecker(self.domain))
        if self.fixVert:
            pde.setValue(q = whereZero(x[0]-inf(x[0])) + whereZero(x[0]-sup(x[0]))) 
        elif self.fixBase:
            pde.setValue(q=whereZero(x[1]-inf(x[1]))) 
        else:
            pde.setValue(q = whereZero(x[0]-inf(x[0]))*whereZero(x[1]-inf(x[1]))) 
        if hasFeature('trilinos'):
            optionsG.setPackage(SolverOptions.TRILINOS)
            optionsG.setPreconditioner(SolverOptions.AMG)
        return pde
             
    def setSusceptibility(self, k=0):
        """
        set susceptibility
        : param k: susceptibility
        : type k: `Data` or `float`
        """
        self.k=k
        self.reset=True

    def getSusceptibility(self):
        """
        returns susceptibility
        : returns: k
        """
        return self.k
        
    def setBackgroundMagneticField(self, Bh=[ 0., 45000.0] ):
        """
        sets background magnetic field in nT

        """
        self.Bh=Bh
        self.reset=True
    
    def getAnomalyPotential(self):
        """
        get the potential of the the magnetic anomaly
        """
        if self.reset:
            self.pde.setValue(X = self.k*self.Bh)
        return self.pde.getSolution()

    def getMagneticFieldAnomaly(self):
        """
        get the total Magnetic field
        """
        return -grad(self.getAnomalyPotential(), ReducedFunction(self.pde.getDomain()))

class MagneticModel3D(object):
    """
    This class is a simple wrapper for a 3D magnetic forward model.  
    It solves PDE
        - div (grad u) = -div(k Bh)
    where 
       k  is magnetic susceptibility and 
       Bh is background magnetic field.
    Possible boundary conditions are Dirichlet on 
       - one corner (bottom left), 
       - vertical sides, or 
       - base.

    Input is domain and boundary condition choice.  
    Default boundary conditions fix the left front bottom corner.  
    Otherwise add 
       - fixVert = True or 
       - fixBase = True.  

    It has functions 
       - setSusceptibility
       - getSusceptibility
       - setBackgroundMagneticField
       - getAnomalyPotential
       - getMagneticFieldAnomaly.
    """

    def __init__(self, domain, fixVert=False, fixBase=False):
        """
        Initialise the class with domain and boundary conditions.
        Setup PDE, susceptibility and background magnetic field.
        :param domain: the domain 
        :type domain: `Domain`
        :param fixBase: if true the magnetic field at the bottom is set to zero. .
        :type fixBottom: `bool`
        :param fixVert: if true the magnetic field on all vertical sudes is set to zero.
        :type fixBottom: `bool`
        :if fixBottom and fixBase are False then magnetic field is set to zero at bottom, front, left corner.
        """
        self.domain  = domain
        self.fixVert = fixVert
        self.fixBase = fixBase
        assert self.domain.getDim() == 3
        self.pde=self.__createPDE()
        self.setSusceptibility()
        self.setBackgroundMagneticField()

    def __createPDE(self):
        """
        Create the PDE and set boundary conditions.
        """
        pde=LinearSinglePDE(self.domain, isComplex=False)
        optionsG=pde.getSolverOptions()
        optionsG.setSolverMethod(SolverOptions.PCG)
        pde.setSymmetryOn()
        pde.setValue(A=kronecker(self.domain.getDim()))
        x=self.domain.getX()
        pde.setValue(A=kronecker(self.domain))
        if self.fixVert:
            pde.setValue(q = whereZero(x[0]-inf(x[0])) + whereZero(x[0]-sup(x[0])) + whereZero(x[1]-inf(x[1])) + whereZero(x[1]-sup(x[1])))
        elif self.fixBase:
            pde.setValue(q=whereZero(x[2]-inf(x[2]))) 
        else:
            pde.setValue(q = whereZero(x[0]-inf(x[0]))*whereZero(x[1]-inf(x[1])) * whereZero(x[2]-inf(x[2]))) 
        if hasFeature('trilinos'):
            optionsG.setPackage(SolverOptions.TRILINOS)
            optionsG.setPreconditioner(SolverOptions.AMG)
        return pde
             
    def setSusceptibility(self, k=0):
        """
        sets susceptibility
        : param k: susceptibility
        : type k: `Data` or `float`
        """
        self.k=k
        self.reset=True

    def getSusceptibility(self):
        """
        returns the susceptibility
        : returns: k
        """
        return self.k
        
    def setBackgroundMagneticField(self, Bh=[ 0., 45000.0, 0.] ):
        """
        sets background magnetic field in nT

        """
        self.Bh=Bh
        self.reset=True
    
    def getAnomalyPotential(self):
        """
        get the potential of the the magnetic anomaly
        """
        if self.reset:
            self.pde.setValue(X = self.k*self.Bh)
        return self.pde.getSolution()

    def getMagneticFieldAnomaly(self):
        """
        get the total Magnetic field
        """
        return -grad(self.getAnomalyPotential(), ReducedFunction(self.pde.getDomain()))


