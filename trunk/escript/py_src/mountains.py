
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript.linearPDEs import LinearPDE, SolverOptions
from esys import finley
import sys

class Mountains:
  """
  The Mountains class is defined by the following equations:
  
  (1) w_i,aa+(1/eps)*w_i,33=0 where 0<eps<<1 and a=1,2 and w_i is the extension of the surface velocity where w_i(x_3=1)=v_i.
  
  (2) Integration of topography PDE using Taylor-Galerkin upwinding to stabilize the advection terms
      H^(t+dt)=H^t+dt*w_3+w_hat*dt*[(div(w_hat*H^t)+w_3)+(dt/2)+H^t],
      where w_hat=w*[1,1,0], dt<0.5*d/max(w_i), d is a characteristic element size; H(x_3=1)=lambda (?) 
      
  """
  def __init__(self,domain,eps=0.01):
    """
    Sets up the level set method.

    @param domain: the domain where the mountains is used
    @param eps: the smoothing parameter for (1)
    """
    if eps<=0:
        raise ValueError("Smmoting parameter eps must be positive.")
    self.__domain = domain
    self.__eps = eps
    self.__DIM=domain.getDim()
    z=domain.getX()[self.__DIM-1]

    self.__PDE_W = LinearPDE(domain)
    self.__PDE_W.setSymmetryOn()
    A=kronecker(domain)
    A[self.__DIM-1,self.__DIM-1]=1/self.__eps
    self.__PDE_W.setValue(A=A, q=whereZero(sup(z)-z)+whereZero(inf(z)-z))

    self.__PDE_H = LinearPDE(domain)
    self.__PDE_H.setSymmetryOn()
    self.__PDE_H.setValue(D=1.0)
    self.__PDE_H.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
    self.__PDE_H.setValue(q=whereZero(inf(z)-z))

    self.setVelocity()
    self.setTopography()
  def getSolverOptionsForSmooting(self):
     """
     returns the solver options for the smoothing/extrapolation
     """
     return self.__PDE_W.getSolverOptions()

  def getSolverOptionsForUpdate(self):
     """
     returns the solver options for the topograthy update
     """
     return self.__PDE_H.getSolverOptions()
  def getDomain(self):
      """
      Returns the domain.
      """
      return self.__domain

  def setVelocity(self,v=None):
      """
      set a new velocity. v is define on the entire domain but only the surface values are used.

      @param v: velocity field. If None zero is used.
      @type v: vector
      """
      self.__dt=None
      self.__v=Vector(0.,Solution(self.getDomain()))
      if not v == None:
        z=self.getDomain().getX()[self.__DIM-1]
        z_min=inf(z)
        z_max=sup(z)
        f=(z-z_min)/(z_max-z_min)
        for d in range(self.__DIM):
           self.__PDE_W.setValue(r=v[d]*f)
           self.__v[d]=self.__PDE_W.getSolution()
  def getVelocity(self):
      """
      returns the smoothed/extrapolated velocity
      @rtype: vector L{Data} 
      """
      return self.__v

  def setTopography(self,H=None):
    """
    set the topography to H where H defines the vertical displacement. H is defined for the entire domain.

    @param H: the topography.  If None zero is used.
    @type H: scalar
    """

    if H==None: 
       self.__H=Scalar(0.0, Solution(self.getDomain()))
    else:
       self.__H=interpolate(H, Solution(self.getDomain()))
       
  def getTopography(self):
     """
     returns the current topography.
     @rtype: scalar L{Data}
     """
     return self.__H

  def getSafeTimeStepSize(self):
      """
      Returns the time step value.

      @rtype: C{float}
      """
      if self.__dt == None:
           self.__dt=0.5*inf(self.getDomain().getSize()/length(self.getVelocity()))
      return self.__dt
  def update(self,dt=None):
      """
      Sets a new W and updates the H function.

      @param dt: time step forward. If None the save time step size is used.
      @type dt: positve C{float} which is less or equal than the safe time step size.
      
      """
      if dt == None: 
            dt = self.getSafeTimeStepSize()
      if dt<=0:
           raise ValueError("Time step size must be positive.")
      if dt>self.getSafeTimeStepSize():
           raise ValueError("Time step must be less than safe time step size = %e."%self.getSafeTimeStepSize())
      
      H=self.getTopography()
      w=self.getVelocity()
      w_tilda=1.*w
      w_tilda[self.__DIM-1]=0
      w_z=w[self.__DIM-1]

      self.__PDE_H.setValue(X=((div(w_tilda*H)+w_z)*dt/2+H)*w_tilda*dt, Y=w_z*dt+H, r=H)
      self.setTopography(self.__PDE_H.getSolution())

      return self.getTopography()

