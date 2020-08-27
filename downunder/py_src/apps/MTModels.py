__copyright__ = "Copyright (c) 2020 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Lutz Gross, Andrea Codd"

from esys.escript import *
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.escript.pdetools import Locator
import cmath
import numpy as np

class MT2DTEModel(object):
    """
    This class is a simple wrapper for 2D MT PDE model in the TE mode.
    MT  
        curl ((1/sigma) curl H) + i omega H = 0
        curl ((1/mu) curl E) + i omega sigma E = 0
    
    2D reduces to 
        -div (1/mu grad u) + i omega sigma u = 0
    where 
        u = Ex is transverse component of electric field 
        mu is magnetic permeability
        sigma is electrical conductivity
        omega is angular frequency
        i = sqrt(-1)   
 
    Domain typically includes air and ground layers.
    Conductivity sigma = 0 in the air layer.
        
    Boundary conditions included in the class are  
       - Ex is set to one at the top of the domain, typically at the top of an air layer.
       - At the bottom of the domain Ex=0 (set `fixBottom`=True) 
         or radiation condition dEx/dn+k*Ex=0 with k^2=2*pi*f*mu*sigma is set

    It has a function to set ground property
       - setConductivity
    and functions to output solutions
       - getImpedance
       - getApparentResitivity
       - getPhase. 

    """

    def __init__(self, domain, fixBottom=False, useFastSolver=False,  mu=4*np.pi*1e-7):
        """
        :param domain: the domain 
        :type domain: `Domain`
        :param fixBottom: if true the electric field at the bottom is set to zero. 
            Otherwise radiation condition is set.
        :type fixBottom: `bool`
        :param useFastSolver: use multigrid solver. (not supported yet)
        :type useFastSolver: `bool`
        """
        self.domain=domain
        self.mu=mu
        self.sigma=None
        self.sigma_boundary=None
        self.fixBottom=fixBottom
        self.useFastSolver=False
        self.pde=self.__createPDE()

    def __createPDE(self):
        """
        Create the PDE and set boundary conditions.
        """
        pde=LinearSinglePDE(self.domain, isComplex=True,)
        optionsG=pde.getSolverOptions()
        optionsG.setSolverMethod(SolverOptions.DIRECT)
        pde.setSymmetryOn()
        if self.useFastSolver and hasFeature('trilinos'): # ignored for now!
            optionsG.setPackage(SolverOptions.TRILINOS)
            optionsG.setPreconditioner(SolverOptions.AMG)
        pde.setValue(A=kronecker(self.domain.getDim()))
        z=self.domain.getX()[self.domain.getDim()-1]
        t=sup(z)
        if self.fixBottom:
            b=inf(z)
            pde.setValue(q=whereZero(z-t)+whereZero(z-b), r=(z-b)/(t-b))
        else:
            pde.setValue(q=whereZero(z-t), r=1.)
        return pde
    
    def setConductivity(self, sigma, sigma_boundary=None):
        """
        sets the conductivity `sigma`. 
        
        :param sigma: conductivity distribution. 
        :type sigma: `Data` or `float`
        :param sigma_boundary: conductivity distribution on bottom boundary. 
           Only required if fixBottom is not set and `sigma` cannot be 
           interpolated to boundary.
        :type sigma: `Data` or `float`
        """
        self.sigma=interpolate(sigma, Function(self.domain))
        if not self.fixBottom:
            if sigma_boundary:
                self.sigma_boundary=interpolate(sigma_boundary, FunctionOnBoundary(self.domain))
            else:
                self.sigma_boundary=interpolate(sigma, FunctionOnBoundary(self.domain))
        return self
        
    def getImpedance(self, f=1.):
        """
        return the impedance Zxy for frequency `f` in [Hz]. The electric field
        and magnetic field cane be accessed as attributes `Ex` and `Hy` after 
        completion.
        
        :param f: frequency in [Hz]
        :type f: `float`
        :returns: Zxy
        """
        o=2*np.pi*f
        self.pde.setValue(D=1j*o*self.mu*self.sigma)
        if not self.fixBottom:
            z=FunctionOnBoundary(self.domain).getX()[self.domain.getDim()-1]
            b=inf(z)
            k=(1+1j)*sqrt(o*self.mu*self.sigma_boundary/2)
            self.pde.setValue(d=k*whereZero(z-b))
        Ex=self.pde.getSolution()
        g=grad(Ex, ReducedFunction(self.domain))
        self.Hy=-1./(1j*o*self.mu)*g[self.domain.getDim()-1]
        self.Ex=interpolate(Ex, self.Hy.getFunctionSpace())
        Zxy=self.Ex/self.Hy
        return Zxy
    
    def getApparentResitivity(self, f, Zxy):
        """
        return the apparent resistivity from a given frequency `f` and impedance `Zxy`

        :param f: frequency in [Hz]
        :type f: `float`
        :param Zxy: impedance
        :type Zxy: `Data` or `np.array`
        """
        o=2*np.pi*f
        return abs(Zxy)**2/(self.mu*o)

    def getPhase(self, f, Zxy):
        """
        return the phase in [deg] from a given frequency `f` and impedance `Zxy`
        
        :param f: frequency in [Hz]
        :type f: `float`
        :param Zxy: impedance
        :type Zxy: `Data` or `np.array`
        """
        return atan2(Zxy.imag(),Zxy.real())/np.pi*180

class MT2DTMModel(object):
    """
    This a class for solving the 2D MT model in the TM mode.
    MT  
        curl ((1/sigma)curl H) + i omega H = 0
        curl ((1/mu)curl E) + i omega sigma E = 0
    
    2D 
        -div (1/sigma grad u) + i omega mu u = 0
    where 
        u = Hx is transverse component of magnetic field
        mu is magnetic permeability
        sigma is electrical conductivity
        omega is angular frequency
        i = sqrt(-1)    
    
    Hx is set to one in the air layer and the interface of the air layer and the subsurface.
    At the bottom of the domain Ex=0 (set `fixBottom`=True) 
       or radiation condition dEx/dn+k*Ex=0 with k^2=2*pi*f*mu*sigma is set.
    
    """
    def __init__(self, domain, fixBottom=False, airLayer=None, useFastSolver=False,  mu=4*np.pi*1e-7):
        """
        :param domain: the domain
        :type domain: `Domain`
        :param fixBottom: if true the potential at all faces except the top is set to zero. 
            Otherwise on the the bottom is set to zero.
        :type fixBottom: `bool`
        :param airLayer: defines the air layer including the interface between air layer and subsurface. 
             If set to `None` then just the (plane) top surface is used. 
             If set to a `float` then the region above `airlayer` (including the interface) 
                is defined as air layer. 
             Otherwise `airlayer` needs to be defined as `Data` with value `1` marking 
                 the air layer and its interface.  
        :type airLayer: `None`, `float`, `Data`
        :param useFastSolver: use multigrid solver. This may fail.
        :type useFastSolver: `bool`
        
        :note: the attribute `airLayer` gives the mask for the air layer 
               including the interface between the air layer and the subsurface.
        """
        self.domain=domain
        self.mu=mu
        self.rho=None
        self.rho_boundary=None

        self.fixBottom=fixBottom
        self.useFastSolver=False
        self.pde=self.__createPDE(airLayer)

    def __createPDE(self, airLayer=None):
        pde=LinearSinglePDE(self.domain, isComplex=True,)
        optionsG=pde.getSolverOptions()
        optionsG.setSolverMethod(SolverOptions.DIRECT)
        pde.setSymmetryOn()
        if self.useFastSolver and hasFeature('trilinos'): # ignored for now!
            optionsG.setPackage(SolverOptions.TRILINOS)
            optionsG.setPreconditioner(SolverOptions.AMG)
            optionsG.setTrilinosParameter("problem: symmetric", True)
        z=self.domain.getX()[self.domain.getDim()-1]
        b=inf(z)        
        if airLayer is None:
            self.airLayer=whereNonNegative(z-sup(z))
        elif isinstance(airLayer, float) or  isinstance(airLayer, int):
            self.airLayer=whereNonNegative(z-airLayer)
        else:
            self.airLayer=wherePositive(interpolate(airLayer, Solution(self.domain)))
        if self.fixBottom:
            pde.setValue(q=self.airLayer+whereZero(z-b), r=self.airLayer)
        else:
            pde.setValue(q=self.airLayer, r=self.airLayer)
        return pde
    
    def setResistivity(self, rho, rho_boundary=None):
        """
        sets the resistivity. 
        :param rho: resistivity distribution. 
        :type rho: `Data`
        :param rho_boundary: rho distribution on bottom boundary. Only required if fixBottom is set and
                               `rho` cannot be interpolated to boundary.
        :type sigma: `Data`
        """
        self.rho=interpolate(rho, Function(self.domain))
        if not self.fixBottom:
            if rho_boundary:
                self.rho_boundary=interpolate(rho_boundary, FunctionOnBoundary(self.domain))
            else:
                self.rho_boundary=interpolate(rho, FunctionOnBoundary(self.domain))

        self.pde.setValue(A=self.rho*kronecker(self.domain.getDim()))        
        return self
        
    def getImpedance(self, f=1.):
        """
        return the impedance Zyx and the electric field Ex for frequency `f` 
        
        :param f: frequency in [Hz]
        :type f: `float`
        :returns: Zyx
        """
        o=2*np.pi*f
        self.pde.setValue(D=1j*o*self.mu)
        if not self.fixBottom:
            z=FunctionOnBoundary(self.domain).getX()[self.domain.getDim()-1]
            b=inf(z)
            k=(1+1j)*sqrt(o*self.mu*self.rho_boundary/2)
            self.pde.setValue(d=k*whereZero(z-b))
        Hx=self.pde.getSolution()
        g=grad(Hx, ReducedFunction(self.domain))
        
        self.Ey=self.rho*g[self.domain.getDim()-1]
        self.Hxi=interpolate(Hx, self.Ey.getFunctionSpace())
        Zyx=self.Ey/self.Hxi
        return Zyx
    
    def getApparentResitivity(self, f, Zyx):
        """
        return the apparent resistivity from a given frequency `f` and impedance `Zyx`

        :param f: frequency in Hz
        :type f: `float`
        :param Zyx: impedance
        :type Zyx: `Data` or `np.array`
        """
        o=2*np.pi*f
        return abs(Zyx)**2/(self.mu*o)

    def getPhase(self, f, Zyx):
        """
        return the phase in [deg] from a given frequency f [Hz] and impedance Zyx

        :param f: frequency in Hz
        :type f: `float`
        :param Zyx: impedance
        :type Zyx: `Data` or `np.array`
        """
        return atan2(Zyx.imag(),Zyx.real())/np.pi*180

