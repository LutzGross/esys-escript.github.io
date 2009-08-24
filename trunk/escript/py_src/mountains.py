
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
import math
import sys

class SubSteppingException(Exception):
   """
   Thrown if the L{Mountains} class uses substepping.
   """
   pass


class Mountains:
  """
  The Mountains class is defined by the following equations:
  
  (1) eps*w_i,aa+w_i,33=0 where 0<=eps<<1 and a=1,2 and w_i is the extension of the surface velocity where w_i(x_3=1)=v_i.
  
  (2) Integration of topography PDE using Taylor-Galerkin upwinding to stabilize the advection terms
      H^(t+dt)=H^t+dt*w_3+w_hat*dt*[(div(w_hat*H^t)+w_3)+(dt/2)+H^t],
      where w_hat=w*[1,1,0], dt<0.5*d/max(w_i), d is a characteristic element size; H(x_3=1)=lambda (?) 
      
  """
  def __init__(self,domain,eps=0.01, reduced=False):
    """
    Sets up the level set method.

    :param domain: the domain where the mountains is used
    :param eps: the smoothing parameter for (1)
    :param reduced: if True topography reduced order is used for reconstruction of internal velocity and topography
    """
    if eps<0:
        raise ValueError("Smooting parameter eps must be non-negative.")
    self.__domain = domain
    self.__reduced=reduced
    self.__DIM=domain.getDim()
    z=domain.getX()[self.__DIM-1]

    self.__PDE_W = LinearPDE(domain)
    self.__PDE_W.setSymmetryOn()
    A=kronecker(domain)*eps*0
    A[self.__DIM-1,self.__DIM-1]=(0.3*(sup(z)-inf(z))/log(2.))**2
    # A[self.__DIM-1,self.__DIM-1]=(sup(FunctionOnBoundary(self.__domain).getSize())/log(2.))**2
    self.__PDE_W.setValue(D=1, A=A, q=whereZero(sup(z)-z)) 

    self.__PDE_H = LinearPDE(domain)
    self.__PDE_H.setSymmetryOn()
    if reduced: self.__PDE_H.setReducedOrderOn()
    # A=kronecker(domain)*0
    # A[self.__DIM-1,self.__DIM-1]=0.1
    self.__PDE_H.setValue(D=1.0)
    # self.__PDE_H.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)

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

      :param v: velocity field. If None zero is used.
      :type v: vector
      """
      self.__dt=None
      self.__v=Vector(0.,Solution(self.getDomain()))
      if not v == None:
        for d in range(self.__DIM):
           self.__PDE_W.setValue(r=v[d])
           self.__v[d]=self.__PDE_W.getSolution()
  def getVelocity(self):
      """
      returns the smoothed/extrapolated velocity
      :rtype: vector `Data` 
      """
      return self.__v

  def setTopography(self,H=None):
    """
    set the topography to H where H defines the vertical displacement. H is defined for the entire domain.

    :param H: the topography.  If None zero is used.
    :type H: scalar
    """
    if self.__reduced:
         fs=ReducedSolution(self.getDomain())
    else:
         fs=Solution(self.getDomain())

    if H==None: 
       self.__H=Scalar(0.0, fs)
    else:
       self.__H=interpolate(H, fs)
       
  def getTopography(self):
     """
     returns the current topography.
     :rtype: scalar `Data`
     """
     return self.__H

  def getSafeTimeStepSize(self):
      """
      Returns the time step value.

      :rtype: ``float``
      """
      if self.__dt == None:
           h=self.getDomain().getSize()
           self.__dt=0.5*inf(h/length(interpolate(self.getVelocity(),h.getFunctionSpace())))
      return self.__dt
  def update(self,dt=None, allow_substeps=True):
      """
      Sets a new W and updates the H function.

      :param dt: time step forward. If None the save time step size is used.
      :type dt: positve ``float`` which is less or equal than the safe time step size.
      
      """
      if dt == None: 
            dt = self.getSafeTimeStepSize()
      if dt<=0:
           raise ValueError("Time step size must be positive.")
      dt_safe=self.getSafeTimeStepSize()
      n=max(int(math.ceil(dt/dt_safe)+0.5),1)
      if n>1 and not allow_substeps:
         raise SubSteppingException,"Substepping required."
      dt/=n
 
      H=self.getTopography()
      w=self.getVelocity()
      w_tilda=1.*w
      w_tilda[self.__DIM-1]=0
      w_z=w[self.__DIM-1]

      t=0
      for i in range(n):
         L=integrate(w_z*dt+H)/vol(self.__PDE_H.getDomain())
         self.__PDE_H.setValue(X=(inner(w_tilda,grad(H))*dt/2+H)*w_tilda*dt, Y=w_z*dt+H-L)
         H=self.__PDE_H.getSolution()
         print "DDD : ava = ",integrate(H)
         t+=dt
      self.setTopography(H)

      return self.getTopography()
