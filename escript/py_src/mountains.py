
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
from esys.escript.linearPDEs import LinearPDE
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
  def __init__(self,domain,v,eps=0.01,z=1):
    """
    Sets up the level set method.

    @param domain: the domain where the mountains is used
    @param eps: the parameter for (1)
    @param z: the height of the box
    """
    self.__domain = domain
    self.__eps = eps
    self.__PDE_W = LinearPDE(domain)
    self.__PDE_W.setSymmetryOn()
    self.__PDE_H = LinearPDE(domain)
    self.__PDE_H.setSymmetryOn()
    self.__PDE_H.setValue(D=1.0)
    self.__PDE_H.setSolverMethod(solver=LinearPDE.LUMPING)
    self.__h = inf(domain.getSize())
    self.__v=v
    self.__z=z
    self.__x=domain.getX()
    self.__DIM=domain.getDim()
    self.__fixed_w_mask=whereZero(self.__z-self.__x[self.__DIM-1])+whereZero(self.__x[self.__DIM-1])
    self.__fixed_H_mask=whereZero(self.__x[self.__DIM-1])
    self.__A=kronecker(domain)
    self.__A[self.__DIM-1,self.__DIM-1]=1/self.__eps
    self.__PDE_W.setValue(A=self.__A, q=self.__fixed_w_mask)
    self.__PDE_H.setValue(q=self.__fixed_H_mask)
    self.__H_t=Scalar(0.0, Solution(domain))
    self.__dt=0.

  def update(self,u=None,H_t=None,dt=None,verbose=False):
      """
      Sets a new W and updates the H function.

      @param dt: time step forward
      """
      if H_t==None:
        H_t=self.__H_t
        
      if u==None:
         u=self.__v
         
      for d in range(self.__DIM):
        self.__PDE_W.setValue(r=u[d]*self.__x[self.__DIM-1]/self.__z)
        u[d]=self.__PDE_W.getSolution(verbose=verbose)

      if dt==None:
        dt=0.5*self.__h/sup(u)
        self.__dt=dt
      else:
        self.__dt=dt
 
      w=u
      w_tilda=1.*w
      w_tilda[self.__DIM-1]=0

      self.__PDE_H.setValue(X=((div(w_tilda*H_t)+w[self.__DIM-1])*dt/2+H_t)*w_tilda*dt, Y=w[self.__DIM-1]*dt+H_t)
      self.__H_t=self.__PDE_H.getSolution()
      
      self.__v=w

      return w,self.__H_t

  def getH(self):
      """
      Returns the mesh size.
      """
      return self.__h

  def getDt(self):
      """
      Returns the time step value.
      """
      return self.__dt
    

  def getDomain(self):
      """
      Returns the domain.
      """
      return self.__domain

  def setTolerance(self,tolerance=1e-3):
    self.__PDE_W.setTolerance(tolerance)
    self.__PDE_H.setTolerance(tolerance)


