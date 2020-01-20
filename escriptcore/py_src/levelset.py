
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

from . import escriptcpp as esc
import esys.escriptcore.linearPDEs as lpe
import esys.escriptcore.util as es
import math

class LevelSet(object):
  """
  The level set method tracking an interface defined by the zero contour of the
  level set function phi which defines the signed distance of a point x from the
  interface. The contour phi(x)=0 defines the interface.

  It is assumed that phi(x)<0 defines the volume of interest,
  phi(x)>0 the outside world. 
  """
  def __init__(self,phi,reinit_max=10,reinitialize_after=1,smooth=2., useReducedOrder=False):
    """
    Sets up the level set method.

    :param phi: the initial level set function
    :param reinit_max: maximum number of reinitialization steps
    :param reinitialize_after: ``phi`` is reinitialized every ``reinit_after`` step
    :param smooth: smoothing width
    """
    self.__domain = phi.getDomain()
    self.__phi = phi


    self.__transport=lpe.SingleTransportPDE(self.__domain)
    if useReducedOrder: self.__transport.setReducedOrderOn()
    self.__transport.setValue(M=1.0)
    self.__transport.setInitialSolution(phi)


    self.__reinitPDE = lpe.LinearPDE(self.__domain, numEquations=1)
    self.__reinitPDE.getSolverOptions().setSolverMethod(lpe.SolverOptions.LUMPING)
    if useReducedOrder: self.__reinitPDE.setReducedOrderOn()
    self.__reinitPDE.setValue(D=1.0)

    # revise:
    self.__reinit_max = reinit_max
    self.__reinit_after = reinitialize_after
    self.__h = es.inf(self.__domain.getSize())
    self.__smooth = smooth
    self.__n_step=0

  def getH(self):
      """
      Returns the mesh size.
      """
      return self.__h

  def getDomain(self):
      """
      Returns the domain.
      """
      return self.__domain

  def getAdvectionSolverOptions(self):
      """
      Returns the solver options for the interface advective.
      """
      return self.__transport.getSolverOptions()

  def getReinitializationSolverOptions(self):
      """
      Returns the options of the solver for the reinitialization
      """
      return self.__reinitPDE.getSolverOption()

  def getLevelSetFunction(self):
      """
      Returns the level set function.
      """
      return self.__phi

  def getTimeStepSize(self,flux):
       """
       Returns a new ``dt`` for a given ``flux`` using the Courant condition.

       :param flux: flux field
       """
       self.__transport.setValue(C=-flux)
       dt=self.__transport.getSafeTimeStepSize()
       return dt

  def update(self,dt):
      """
      Updates the level set function.

      :param dt: time step forward
      """
      self.__phi = self.__transport.getSolution(dt)
      self.__n_step+=1
      if self.__n_step%self.__reinit_after ==0: 
          self.__phi = self.__reinitialise(self.__phi)
          self.__transport.setInitialSolution(self.__phi)
      return self.__phi


  #==============================================================================================
  def __reinitialise(self, phi):
    """
    Reinitializes the level set.

    It solves the PDE...

    :return: reinitialized level set
    """
    fs=esc.ReducedFunction(self.__domain)
    dtau = 0.2*self.__h
    n =0
    #g=grad(phi,fs)
    #error = Lsup((1-length(g))*whereNegative(abs(phi.interpolate(fs))-5.0*self.__h))
    #print(("LevelSet:reinitialization: iteration :", n, " error:", error))
    #mask = whereNegative(abs(phi)-1.*self.__h)
    #self.__reinitPDE.setValue(q=mask, r=phi)
    s= es.sign(phi.interpolate(fs)) 
    while (n<=self.__reinit_max):
      g_phi=es.grad(phi, fs)
      self.__reinitPDE.setValue(Y =  dtau* s * (1 - es.length(g_phi) ), X = - dtau**2/2 * g_phi) 
      phi = phi + self.__reinitPDE.getSolution()
      

      n +=1
    #g=grad(phi,fs)
    #error = Lsup((1-length(g))*whereNegative(abs(phi.interpolate(fs))-5.0*self.__h))
    #print(("LevelSet:reinitialization: iteration :", n, " error:", error))
    return phi


  def getVolume(self):
      """
      Returns the volume of the *phi(x)<0* region.
      """
      return integrate(es.whereNegative(self.__phi.interpolate(esc.Function(self.__domain))))


  def getJumpingParameter(self, param_neg=-1, param_pos=1, phi=None):
      """
      Creates a function with ``param_neg`` where ``phi<0`` and ``param_pos``
      where ``phi>0`` (no smoothing).

      :param param_neg: value of parameter on the negative side (phi<0)
      :param param_pos: value of parameter on the positive side (phi>0)
      :param phi: level set function to be used. If not present the current
                  level set is used.
      """
      mask_neg = es.whereNegative(self.__phi)
      mask_pos = es.whereNonNegative(self.__phi)
      param = param_pos*mask_pos + param_neg*mask_neg
      return param

  def getSmoothedParameter(self, param_neg=-1, param_pos=1, phi=None, smoothing_width=None):
      """
      Creates a smoothed function with ``param_neg`` where ``phi<0`` and
      ``param_pos`` where ``phi>0`` which is smoothed over a length
      ``smoothing_width`` across the interface.

      :param smoothing_width: width of the smoothing zone relative to mesh size.
                              If not present the initial value of ``smooth`` is
                              used.
      """
      if smoothing_width is None: smoothing_width = self.__smooth
      if phi is None: phi=self.__phi
      s=self.getSmoothedJump(phi,smoothing_width)
      return ((param_pos-param_neg)*s+param_pos+param_neg)/2

  def getSmoothedJump(self,phi=None,smoothing_width=None):
      """
      Creates a smooth interface from -1 to 1 over the length
      *2*h*smoothing_width* where -1 is used where the level set is negative
      and 1 where the level set is 1.
      """
      if smoothing_width is None: smoothing_width = self.__smooth
      if phi is None: phi = self.__phi
      s=smoothing_width*self.__h
      phi_on_h=es.interpolate(phi,esc.Function(self.__domain))
      mask_neg = es.whereNonNegative(-s-phi_on_h)
      mask_pos = es.whereNonNegative(phi_on_h-s)
      mask_interface = 1.-mask_neg-mask_pos
      interface=phi_on_h/s
      return - mask_neg + mask_pos + mask_interface * interface

  def getInterface(self,phi=None,smoothing_width=None):
      """
      creates a characteristic function which is 1 over the over the length
      *2*h*smoothing_width*  around the interface and zero elsewhere
      """
      if smoothing_width is None: smoothing_width = self.__smooth
      if phi is None: phi = self.__phi
      s=smoothing_width*self.__h 
      phi_on_h=es.interpolate(phi,esc.Function(self.__domain))
      return es.whereNegative(abs(phi_on_h)-s)
      
  def makeCharacteristicFunction(self, contour=0, phi=None, positiveSide=True, smoothing_width=None):
      """
      Makes a smooth characteristic function of the region ``phi(x)>contour`` if
      ``positiveSide`` and ``phi(x)<contour`` otherwise.

      :param phi: level set function to be used. If not present the current
                  level set is used.
      :param smoothing_width: width of the smoothing zone relative to mesh size.
                              If not present the initial value of ``smooth`` is
                              used.
      """
      if phi is None: phi=self.__phi
      if smoothing_width is None: smoothing_width=self.__smooth
      s=self.getSmoothedJump(phi=phi-contour,smoothing_width=smoothing_width)
      if positiveSide:
          return (1+s)/2
      else:
          return (1-s)/2
