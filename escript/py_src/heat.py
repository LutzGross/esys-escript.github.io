
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
Some models for heat advection-diffusion

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

# from escript import *
import util
from linearPDEs import TransportPDE

class TemperatureCartesian(TransportPDE):
      """
      solves

          rhocp(T_{,t} + v_i T_{,i} - ( k T_{,i})_i = Q

                 k T_{,i}*n_i=surface_flux

                 T_{,t} = 0  where given_T_mask>0

      if surface_flux is not given 0 is assumed.

      typical usage:

            sp=TemperatureCartesian(domain)
            sp.setTolerance(1.e-4)
            t=0
            T=...
            sp.setValues(rhocp = ..,  v=.., k=.., given_T_mask=..)
            sp.setInitialTemperature(T)
            while t < t_end:
                sp.setValue(Q= ...)
                T=sp.getTemperature(dt)
                t+=dt
      """
      def __init__(self,domain,theta=0.5,**kwargs):
        """
        initializes temperature advection-diffuision problem

        @param domain: domain of the problem
        @param theta: method control: theta=1. backward Euler
                                      theta=0.5: Crank-Nicholson scheme
                                      theta=0.: forward Euler (not recommendable)
        """
        TransportPDE.__init__(self,domain,numEquations=1,theta=theta,**kwargs)
        self.setReducedOrderOn()
        self.__rhocp=None
        self.__v=None

      def setInitialTemperature(self,T):
        self.setInitialSolution(T)

      def setValue(self,rhocp=None,v=None,k=None,Q=None,surface_flux=None,given_T_mask=None):
           if rhocp!=None:
               self.__rhocp=rhocp
           if v!=None:
               self.__v=v
           if rhocp!=None:
               super(TemperatureCartesian,self).setValue(M=self.__rhocp)
           if (rhocp!=None or v!=None) and self.__rhocp!=None and self.__v!=None:
               super(TemperatureCartesian,self).setValue(C=-self.__rhocp*self.__v)
           if k!=None:
               super(TemperatureCartesian,self).setValue(A=-k*util.kronecker(self.getDomain()))
           if Q!=None:
               super(TemperatureCartesian,self).setValue(Y=Q)
           if surface_flux!=None:
               super(TemperatureCartesian,self).setValue(y=surface_flux)
           if given_T_mask!=None:
               super(TemperatureCartesian,self).setValue(q=given_T_mask)

      def getTemperature(self,dt,**kwargs):
          return self.getSolution(dt,**kwargs)

