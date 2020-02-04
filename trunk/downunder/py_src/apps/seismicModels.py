from esys.escript import *
from esys.escript.linearPDEs import LinearPDE, SolverOptions, LinearSinglePDE
import numpy as np




class SonicWaveInFrequencyDomain(object):
    """  
    This class is a simple wrapper for the solution of the 2D or 3D sonic wave equtaion in the frequency domain.  
    It solves complex PDE
        div (grad p) - k^2 p = Q
    where 
        k = omega/c
        Q = source term : Dirac function
 
    Boundary conditions implemented are   NEED TO SET PML CONDITION FIRST
       - perfectly matched layer : see class PMLCondition
       - fixed BC - all but the top surface

    It has functions
       - getDomain
       - setFrequency
       - setVp
       
    Output solution
       - getWave
    """

    def __init__(self, domain, pml_condition=None, frequency=None, vp=None, fix_boundary=True):
        """
        Initialise the class with domain, boundary conditions and frequency.
        Setup PDE 
        :param domain: the domain : setup with Dirac points
        :type domain: `Domain`
        :param pml_condition: 
        :type pml_condition:  
        :param frequency: 
        :type frequency: 
        :param vp: velocity field 
        :type vp: `Data`
        :param fix_boundary: if true fix all the boundaries except the top
        :type fix_boundary: `bool`
        
        """
        Dim=domain.getDim()
        self.pde=LinearPDE(domain, numEquations=1, numSolutions=1, isComplex=True)
        self.pde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
        self.pml=pml_condition
        if not self.pml:
            self.pde.setValue(A=kronecker(domain))
            self.J=1.
        if fix_boundary:
            x=domain.getX()
            q=whereZero(x[Dim-1]-inf(x[Dim-1]))
            for i in range(Dim-1):
                q+=whereZero(x[i]-inf(x[i]))+whereZero(x[i]-sup(x[i]))
            self.pde.setValue(q=q)
        self.setFrequency(frequency)        
        self.setVp(vp)
        
    def getDomain(self):
        """
        returns the domain
        """
        return self.pde.getDomain()
    
    def setFrequency(self, frequency):
        """
        sets the frequency 
        """
        self.frequency=frequency
        self.newFreq=True
        return self
        
    def setVp(self, vp):
        """
        sets the velocity for the rock
        """
        self.newVp=True
        if vp is None:
            self.sigma2=None
        else:
            self.sigma2=1./vp**2
        return self
    
    def getWave(self, source):
        """
        solve the PDE 
        """
        omega=self.frequency*2*np.pi
        if self.pml and self.newFreq:
            Dim=self.getDomain().getDim()
            A=kronecker(Function(self.getDomain()))
            alpha, J = self.pml.getPMLWeights(self.getDomain(), omega=omega)
            for d in range(Dim):
                A[d,d]=J/alpha[d]**2
            self.J=J
            self.pde.setValue(A=A)

        if self.newFreq or self.newVp:
            if self.frequency is None or self.sigma2 is None:
                raise ValueError("freqency/propagation speed vp is not set.")
            omega=self.frequency*2*np.pi
            self.pde.setValue(D=-self.J*omega**2*self.sigma2)
            self.newFreq=False
            self.newVp=False
        self.pde.setValue(y_dirac=source)
        return self.pde.getSolution()

class PMLCondition(object):
    """
    this defines the PML weights over a domain
    :Lleft: thicknesses of PML to the left, bottom, front (None -> no PML)
    :Lright: thicknesses of PML to the right, top, back (None -> no PML)
    :return: return a mask where PML is applied. 
    """
    def __init__(self, sigma0=1., Lleft=[None, None, None], Lright=[None, None, None], m=3):
        """
        initializes the PML pml_condition
        
        :param sigma0: maximum PML damping 
        :type sigma0: `Data`
        :param  Lleft: thicknesses of PML to the left, bottom, front (None -> no PML)
        :type Lleft: list of `floats`, length is domain dimension
        :param Lright: thicknesses of PML to the right, top, back (None -> no PML)
        :type Lright: list of `floats`, length is domain dimension
        :param m: exponent of increase over PML layer (default 3 )
        :type m: `int` or `float`
        """
        self.sigma0=sigma0
        self.Lleft=Lleft
        self.Lright=Lright
        self.m=m
        

    def getPMLWeights(self,domain, omega):
        """
        this defines the PML weights over a domain
        :return: weighting alphas, J=product of alphas and J/alpha 
        """
        alpha=[1. ] * domain.getDim()
        J=1.
        X=Function(domain).getX()
        for d in range(domain.getDim()):
            x=domain.getX()[d]
            Q=0
            if self.Lleft[d] is not None:
                xmin=inf(x)
                Q=wherePositive(xmin+self.Lleft[d]-X[d])*((xmin+self.Lleft[d]-X[d])/self.Lleft[d])**self.m
            if self.Lright[d] is not None:
                xmax=sup(x)
                Q=Q+wherePositive(X[d]-xmax+self.Lright[d])*((X[d]-xmax+self.Lright[d])/self.Lright[d])**self.m
            alpha[d]=1-1j*self.sigma0/omega*Q
            J=J*alpha[d]

        return alpha, J


    def getPMLMask(self,domain):
        """
        returns the mask for PML layer 
        """
        x=domain.getX()
        mask=Scalar(0., x.getFunctionSpace())
        for d in range(domain.getDim()):
            xx=x[d]
            if self.Lleft[d] is not None:
                xmin=inf(xx)
                mask+=wherePositive(xmin+self.Lleft[d]-xx)
            if self.Lright[d] is not None:
                xmax=sup(xx)
                mask+=wherePositive(xx-xmax+self.Lright[d])

        return wherePositive(mask) 

