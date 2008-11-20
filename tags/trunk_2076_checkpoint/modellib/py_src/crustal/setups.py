
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
__url__="http://www.uq.edu.au/esscc/escript-finley"

"""
set up of mining activities in modelframe

@var __author__: name of author
@var __licence__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

from esys.escript.modelframe import Model, ParameterSet
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from mines import parse
import numarray

class MiningHistory(Model):
    """
    manages the history of mines. It mandels the time steping according to the available 
    data and a dictionary of mass changes per year for all the mines. 

    @ivar history: mine site history file.
    @type history: C{DataSource}
    @ivar mass_changes: current mass change per year.
    @type mass_changes: C{DataSource}

    """
    def __init__(self,**kwargs):
        """
        """
        super(MiningHistory,self).__init__(**kwargs)
        self.declareParameter(t=0.,
                              history=None,
                              mass_changes=None)

    def doInitialization(self):
        """
        initialize time integration
        """
        self.__minesite=parse(open(self.history.getLocalFileName(),'r'))
        self.mass_changes=self.__minesite.getMassChanges(self.t)

    def doStepPreprocessing(self, dt):
        self.mass_changes=self.__minesite.getMassChanges(self.t)

    def getSafeTimeStepSize(self, dt):
        return self.__minesite.getNextTimeMarker(self.t)-self.t


class DensityChange(Model):
    """
    translates local mass extraction into density changes.
    "local" refers a taged region.

    @ivar domain: mining region
    @type domain: L{Domian}
    @ivar mass_rate: rate of change mass in a tagged region
    @type mass_rate: C{dict}
    @ivar density_rate: density in each region
    @type density_rate:  C{Data}
    """
    def __init__(self,**kwargs):
        """
        """
        super(DensityChange,self).__init__(**kwargs)
        self.declareParameter(domain=None,
                              mass_rate={},
                              density_rate=None)

    def doInitialization(self):
        """
        initialize time integration
        """
        self.__volumes={}
        for i in getTagNames(self.domain):
            c=integrate(insertTaggedValues(Function(self.domain),**{i : 1.}))
            if c>0: self.__volumes[i] = c

    def doStepPreprocessing(self, dt):
         d={}
         for i in self.__volumes:
            if self.mass_rate.has_key(i): d[i]=-self.mass_rate[i]/self.__volumes[i]
         self.density_rate=inserTaggedValues(Scalar(0.,Function(self.domain)),**d)

class LinearElasticStressChange(Model):
    """
    calculates the stress according to an initial garvity field and a consecutive 
    chenge of density  based on linar elastic material model. It is assumed that 
    the lame coefficients don't change over time, the direction of gravity is
    pointing into the x_{dim} direction. 
  
    @note: add stress changes due to tectonic loading and slip

    @ivar domain: mining region
    @type domain: L{Domian}
    @ivar displacement: displacement field
    @type displacement: L{Vector} or C{None}
    @ivar stress: displacement field 
    @type stress: L{Vector}  or C{None}
    @ivar density: initial density distribution
    @type density: C{float} or L{Scalar}
    @ivar density_rate: density rate by tag (may be changed of time)
    @type density_rate:  C{dict}
    @ivar lame_lambda: elasticity coefficient lambda (assumed to be constant over time)
    @type lame_lambda: C{float} or L{Scalar}
    @ivar lame_mu: elasticity coefficient mu (assumed to be constant over time)
    @type lame_mu: C{float} or L{Scalar}
    @ivar location_of_fixed_displacement: mask of locations and component with zero displacements
    @type location_of_fixed_displacement: L{Vector} or C{None}
    """
    def __init__(self,**kwargs):
        """
        """
        super(LinearElasticStressChange,self).__init__(**kwargs)
        self.declareParameter(domain=None,
                              displacement=None, \
                              stress=None, \
                              density=1., \
                              lame_lambda=2., \
                              lame_mu=1., \
                              location_of_fixed_displacement=None, \
                              density_rate=None)

    def doInitialization(self):
        """
        initialize time integration
        """
        if self.stress == None: self.stress=Tensor(0,Function(self.domain))
        if self.displacement==None: self.displacement=Scalar(0,Solution(self.domain))
        self.__pde=LinearPDE(self.domain)
        self.__pde.setSymmetryOn()
        A =Tensor4(0.,Function(self.domain))
        for i in range(self.domain.getDim()):
           for j in range(self.domain.getDim()):
              A[i,j,j,i] += self.lame_mu
              A[i,j,i,j] += self.lame_mu
              A[i,i,j,j] += self.lame_lambda
        self.__pde.setValue(A=A,q=self.location_of_fixed_displacement)

        self.__stress=self.stress
        self.__displacement=self.displacement
        self.__first=True
        
    def doInitialStep(self):
        """
        """
        d=self.domain.getDim()
        self.__pde.setValue(Y=-kronecker(Function(self.domain))[d-1]*9.81*self.density)
        ddisp=self.__pde.getSolution()
        deps=symmetric(grad(ddisp))
        self.stress=self.stress+self.lame_mu*deps+self.lame_lambda*trace(deps)*kronecker(Function(self.domain))
        self.displacement=self.displacement+ddisp
        self.__first=False

    def terminateInitialIteration(self):
        return not self.__first

    def terminateIteration(self):
        return not self.__first

    def doInitialPostprocessing(self):
        self.__stress=self.stress
        self.__displacement=self.displacement
        self.trace("initial displacement range : %s: %s"%(inf(self.displacement[2]),sup(self.displacement[2])))
        self.trace("initial stress range : %s: %s"%(inf(self.stress),sup(self.stress)))

    def doStepPreprocessing(self,dt):
        self.stress=self.__stress
        self.displacement=self.__displacement
        self.__first=True
         
    def doStep(self,dt):
        """
        """
        d=self.domain.getDim()
        if not self.density_rate == None:
           self.__pde.setValue(Y=-kronecker(Function(self.domain))[d-1]*9.81*self.density_rate)
        ddisp=self.__pde.getSolution()
        deps=symmetric(grad(ddisp))
        self.stress=self.stress+dt*(self.lame_mu*deps+self.lame_lambda*trace(deps)*kronecker(Function(self.domain)))
        self.displacement=self.displacement+dt*ddisp
        self.__first=False

    def doStepPostprocessing(self, dt):
        self.__stress=self.stress
        self.__displacement=self.displacement
        self.trace("displacement range : %s: %s"%(inf(self.displacement[2]),sup(self.displacement[2])))
        self.trace("stress range : %s: %s"%(inf(self.stress),sup(self.stress)))

class CoulombFailureStress(ParameterSet):
    """
    calculates the Coulomb failure stress an planes of a given orientation.

    @ivar domain: mining region
    @type domain: L{Domian}
    @ivar displacement: displacement field
    @type displacement: L{Vector} or C{None}
    @ivar stress: displacement field 
    @type stress: L{Vector}  or C{None}
    @ivar density: initial density distribution
    @type density: C{float} or L{Scalar}
    @ivar density_rate: density rate by tag (may be changed of time)
    @type density_rate:  C{dict}
    @ivar lame_lambda: elasticity coefficient lambda (assumed to be constant over time)
    @type lame_lambda: C{float} or L{Scalar}
    @ivar lame_mu: elasticity coefficient mu (assumed to be constant over time)
    @type lame_mu: C{float} or L{Scalar}
    @ivar location_of_fixed_displacement: mask of locations and component with zero displacements
    @type location_of_fixed_displacement: L{Vector} or C{None}
    """
    def __init__(self,**kwargs):
        """
        """
        super(CoulombFailureStress,self).__init__(**kwargs)
        self.declareParameter(stress=numarray.zeros((3,3)),
                              friction_coefficient=0.,
                              normal=numarray.array([1.,0.,0.]))
    def cfs(self):
        """
        returns  Coulomb failure stress
        """
        sn=matrixmult(self.stress,self.normal)
        nsn=inner(self.normal,sn)
        nssn=inner(sn,sn)
        return (sqrt(nssn-nsn**2)-self.friction_coefficient*nsn)/length(self.normal)

