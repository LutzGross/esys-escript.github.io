
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
Some models for heat advection-diffusion

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

from . import escriptcpp as escore
from . import linearPDEs as lpe
from . import util

class TemperatureCartesian(lpe.TransportPDE):
    """
    Represents and solves the temperature advection-diffusion problem

    *rhocp(T_{,t} + v_i T_{,i} - ( k T_{,i})_i = Q*

    *k T_{,i}*n_i=surface_flux* and *T_{,t} = 0* where ``given_T_mask``>0.

    If surface_flux is not given 0 is assumed.

    Typical usage::

        sp = TemperatureCartesian(domain)
        sp.setTolerance(1.e-4)
        t = 0
        T = ...
        sp.setValues(rhocp=...,  v=..., k=..., given_T_mask=...)
        sp.setInitialTemperature(T)
        while t < t_end:
            sp.setValue(Q=...)
            T = sp.getTemperature(dt)
            t += dt
    """
    def __init__(self,domain,**kwargs):
        """
        Initializes the temperature advection-diffusion problem.

        :param domain: domain of the problem
        :note: the approximation order is switched to reduced if the approximation order is nnot linear (equal to 1).
        """
        lpe.TransportPDE.__init__(self,domain,numEquations=1, **kwargs)
        order=escore.Solution(domain).getApproximationOrder()
        if order>1:
            if escore.ReducedSolution(domain).getApproximationOrder()>1: raise ValueError("Reduced order needs to be equal to 1.")
            self.setReducedOrderOn()
        else:
            self.setReducedOrderOff()
        self.__rhocp=None
        self.__v=None

    def setInitialTemperature(self,T):
        """
        Same as `setInitialSolution`.
        """
        self.setInitialSolution(T)

    def setValue(self,rhocp=None,v=None,k=None,Q=None,surface_flux=None,given_T_mask=None):
        if rhocp is not None:
            self.__rhocp=rhocp
        if v is not None:
            self.__v=v
        if rhocp is not None:
            super(TemperatureCartesian,self).setValue(M=self.__rhocp)
        if (rhocp is not None or v is not None) and self.__rhocp is not None and self.__v is not None:
            super(TemperatureCartesian,self).setValue(C=-self.__rhocp*self.__v)
        if k is not None:
            super(TemperatureCartesian,self).setValue(A=-k*util.kronecker(self.getDomain()))
        if Q is not None:
            super(TemperatureCartesian,self).setValue(Y=Q)
        if surface_flux is not None:
            super(TemperatureCartesian,self).setValue(y=surface_flux)
        if given_T_mask is not None:
            super(TemperatureCartesian,self).setValue(q=given_T_mask)

    def getTemperature(self,dt,**kwargs):
        """
        Same as `getSolution`.
        """
        return self.getSolution(dt,**kwargs)


class Tracer(lpe.TransportPDE):
    """
    Represents and solves the tracer problem

    *C_{,t} + v_i C_{,i} - ( k T_{,i})_i) = 0*

    *C_{,t} = 0* where ``given_C_mask``>0.
    *C_{,i}*n_i=0* 

    Typical usage::

        sp = Tracer(domain)
        sp.setTolerance(1.e-4)
        t = 0
        T = ...
        sp.setValues(given_C_mask=...)
        sp.setInitialTracer(C)
        while t < t_end:
            sp.setValue(v=...)
            dt.getSaveTimeStepSize()
            C = sp.getTracer(dt)
            t += dt
    """
    def __init__(self,domain,useBackwardEuler=False,**kwargs):
        """
        Initializes the Tracer advection problem

        :param domain: domain of the problem
        :param useBackwardEuler: if set the backward Euler scheme is used. Otherwise the Crank-Nicholson scheme is applied. Not that backward Euler scheme will return a safe time step size which is practically infinity as the scheme is unconditional unstable. So other measures need to be applied to control the time step size. The Crank-Nicholson scheme provides a higher accuracy but requires to limit the time step size to be stable.
        :type useBackwardEuler: ``bool``
        :note: the approximation order is switched to reduced if the approximation order is nnot linear (equal to 1).
        """
        lpe.TransportPDE.__init__(self,domain,numEquations=1,useBackwardEuler=useBackwardEuler,**kwargs)
        order=escore.Solution(domain).getApproximationOrder()
        if order>1:
            if escore.ReducedSolution(domain).getApproximationOrder()>1: raise ValueError("Reduced order needs to be equal to 1.")
            self.setReducedOrderOn()
        else:
            self.setReducedOrderOff()
        super(Tracer,self).setValue(M=1.)

    def setInitialTracer(self,C):
        """
        Same as `setInitialSolution`.
        """
        self.setInitialSolution(C)

    def setValue(self, v=None, given_C_mask=None, k=None):
        if v is not None:
            super(Tracer,self).setValue(C=-v)
        if k is not None:
            super(Tracer,self).setValue(A=-k*util.kronecker(self.getDomain()))
        if given_C_mask is not None:
            super(Tracer,self).setValue(q=given_C_mask)

    def getTracer(self,dt,**kwargs):
        """
        Same as `getSolution`.
        """
        return self.getSolution(dt,**kwargs)

