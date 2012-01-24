
########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

#from esys.escript import *
from esys.escript.linearPDEs import LinearPDE, SingleTransportPDE
from esys.escript.pdetools import Projector

import math

USE_OLD_VERSION=False
USE_OLD_VERSION_REINIT=False


class LevelSet:
  """
  The level set method tracking an interface defined by the zero contour of the
  level set function phi which defines the signed distance of a point x from the
  interface. The contour phi(x)=0 defines the interface.

  It is assumed that phi(x)<0 defines the volume of interest,
  phi(x)>0 the outside world. 
  """
  def __init__(self,phi,reinit_max=10,reinitialize_after=20,smooth=2., useReducedOrder=False):
    """
    Sets up the level set method.

    :param phi: the initial level set function
    :param reinit_max: maximum number of reinitialization steps
    :param reinitialize_after: ``phi`` is reinitialized every ``reinit_after`` step
    :param smooth: smoothing width
    """
    self.__domain = phi.getDomain()
    self.__phi = phi

    if USE_OLD_VERSION:
       self.__PDE = LinearPDE(self.__domain)
       if useReducedOrder: self.__PDE.setReducedOrderOn()
       self.__PDE.setSolverMethod(solver=LinearPDE.PCG)
       self.__PDE.setValue(D=1.0)
    else:
       self.__PDE=SingleTransportPDE(self.__domain)
       if useReducedOrder: self.__PDE.setReducedOrderOn()
       self.__PDE.setValue(M=1.0)

    if USE_OLD_VERSION_REINIT:
       self.__reinitPDE = LinearPDE(self.__domain, numEquations=1)
       self.__reinitPDE.getSolverOption().setSolverMethod(solver=LinearPDE.LUMPING)
       if useReducedOrder: self.__reinitPDE.setReducedOrderOn()
       self.__reinitPDE.setValue(D=1.0)
    else:
       self.__reinitPDE=SingleTransportPDE(self.__domain)
       if useReducedOrder: self.__reinitPDE.setReducedOrderOn()
       self.__reinitPDE.setValue(M=1.0)

    # revise:
    self.__reinit_max = reinit_max
    self.__reinit_after = reinitialize_after
    self.__h = inf(self.__domain.getSize())
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
      return self.__PDE.getSolverOptions()

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

  def getTimeStepSize(self,velocity):
       """
       Returns a new ``dt`` for a given ``velocity`` using the Courant condition.

       :param velocity: velocity field
       """
       if USE_OLD_VERSION:
          self.__velocity=velocity
          dt=0.5*self.__h/sup(length(velocity))
       else:
          self.__PDE.setValue(C=-velocity)
          dt=self.__PDE.getSafeTimeStepSize()

       return dt

  def update(self,dt):
      """
      Sets a new velocity and updates the level set function.

      :param dt: time step forward
      """
      self.__phi=self.__advect(dt)
      self.__n_step+=1
      if self.__n_step%self.__reinit_after ==0: self.__phi = self.__reinitialise()
      return self.__phi

  def update_phi(self, velocity, dt):
      """
      Updates ``phi`` under the presence of a velocity field.

      If dt is small this call is equivalent to::

          dt=LevelSet.getTimeStepSize(velocity)
          phi=LevelSet.update(dt)

      otherwise substepping is used.

      :param velocity: velocity field
      :param dt: time step forward
      """
      dt2=self.getTimeStepSize(velocity)
      n=int(math.ceil(dt/dt2))
      dt_new=dt/n
      for i in range(n):
           phi=self.update(dt_new)
      return phi

  def __advect(self,  dt):
    """
    Advects the level set function in the presence of a velocity field.

    :param dt: time increment
    :return: the advected level set function
    """
    if USE_OLD_VERSION:
       velocity=self.__velocity
       Y = self.__phi-dt/2.*inner(velocity,grad(self.__phi))
       self.__PDE.setValue(Y=Y)
       phi_half = self.__PDE.getSolution()
       Y = self.__phi-dt*inner(velocity,grad(phi_half))
       self.__PDE.setValue(Y=Y)
       phi = self.__PDE.getSolution()
    else:
       phi = self.__PDE.getSolution(dt, self.__phi)
    print("LevelSet: Advection step done")
    return phi

  #==============================================================================================
  def __reinitialise(self):
    """
    Reinitializes the level set.

    It solves the PDE...

    :return: reinitialized level set
    """
    phi=self.__phi
    s = sign(phi.interpolate(Function(self.__domain)))
    g=grad(phi)
    w = s*g/(length(g)+1e-10)
    dtau = 0.5*self.__h
    iter =0
    mask = whereNegative(abs(phi)-1.*self.__h)
    self.__reinitPDE.setValue(q=mask, r=phi)
    g=grad(phi)
    while (iter<=self.__reinit_max):
      Y = phi+(dtau/2.)*(s-inner(w,g))
      self.__reinitPDE.setValue(Y = Y)
      phi_half = self.__reinitPDE.getSolution()
      Y = phi+dtau*(s-inner(w,grad(phi_half)))
      self.__reinitPDE.setValue(Y = Y)
      phi = self.__reinitPDE.getSolution()
      g=grad(phi)
      error = Lsup(length(g)*whereNegative(abs(phi.interpolate(g.getFunctionSpace()))-3.0*self.__h))
      print(("LevelSet:reinitialization: iteration :", iter, " error:", error))
      iter +=1
    return phi


  def getVolume(self):
      """
      Returns the volume of the *phi(x)<0* region.
      """
      return integrate(whereNegative(self.__phi.interpolate(Function(self.__domain))))


  def update_parameter_sharp(self, param_neg=-1, param_pos=1, phi=None):
      """
      Creates a function with ``param_neg`` where ``phi<0`` and ``param_pos``
      where ``phi>0`` (no smoothing).

      :param param_neg: value of parameter on the negative side (phi<0)
      :param param_pos: value of parameter on the positive side (phi>0)
      :param phi: level set function to be used. If not present the current
                  level set is used.
      """
      mask_neg = whereNegative(self.__phi)
      mask_pos = whereNonNegative(self.__phi)
      param = param_pos*mask_pos + param_neg*mask_neg
      return param

  def update_parameter(self, param_neg=-1, param_pos=1, phi=None, smoothing_width=None):
      """
      Creates a smoothed function with ``param_neg`` where ``phi<0`` and
      ``param_pos`` where ``phi>0`` which is smoothed over a length
      ``smoothing_width`` across the interface.

      :param smoothing_width: width of the smoothing zone relative to mesh size.
                              If not present the initial value of ``smooth`` is
                              used.
      """
      if smoothing_width==None: smoothing_width = self.__smooth
      if phi==None: phi=self.__phi
      s=self.__makeInterface(phi,smoothing_width)
      return ((param_pos-param_neg)*s+param_pos+param_neg)/2

  def __makeInterface(self,phi,smoothing_width):
      """
      Creates a smooth interface from -1 to 1 over the length
      *2*h*smoothing_width* where -1 is used where the level set is negative
      and 1 where the level set is 1.
      """
      s=smoothing_width*self.__h
      phi_on_h=interpolate(phi,Function(self.__domain))
      mask_neg = whereNonNegative(-s-phi_on_h)
      mask_pos = whereNonNegative(phi_on_h-s)
      mask_interface = 1.-mask_neg-mask_pos
      interface=phi_on_h/s
      return - mask_neg + mask_pos + mask_interface * interface

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
      if phi==None: phi=self.__phi
      if smoothing_width == None: smoothing_width=self.__smooth
      s=self.__makeInterface(phi=phi-contour,smoothing_width=smoothing_width)
      if positiveSide:
          return (1+s)/2
      else:
          return (1-s)/2


class LevelSet2(object):
     def __init__(self,phi,reinit_max=10,reinit_after=2,smooth=2.):
         """
         Initializes the model.
         """
         self.__domain = phi.getDomain()
         x=self.__domain.getX()
         diam=0
         for i in range(self.__domain.getDim()):
             xi=x[i]
             diam+=(inf(xi)-sup(xi))**2
         self.__diam=sqrt(diam)
         self.__h = sup(Function(self.__domain).getSize())
         self.__reinit_after=reinit_after
         self.__reinit_max = reinit_max
         self.__smooth = smooth
         self.__phi = phi
         self.__update_count=0
         self.velocity = None

         #==================================================
         self.__FC=True
         if self.__FC:
            self.__fc=TransportPDE(self.__domain,num_equations=1,useBackwardEuler=False)
            # self.__fc.setReducedOrderOn()
            self.__fc.setValue(M=Scalar(1.,Function(self.__domain)))
            self.__fc.setInitialSolution(phi+self.__diam)
         else:
            self.__fcpde = LinearPDE(self.__domain,numEquations=1, numSolutions=1)
            self.__fcpde.setSolverMethod(solver=LinearPDE.LUMPING)
            self.__fcpde.setReducedOrderOn()
            self.__fcpde.setSymmetryOn()
            self.__fcpde.setValue(D=1.)

         #=======================================
         self.__reinitFC=False
         if self.__reinitFC:
            self.__reinitfc=TransportPDE(self.__domain,num_equations=1,useBackwardEuler=True)
            self.__reinitfc.setValue(M=Scalar(1.,Function(self.__domain)))
            self.__reinitfc.setTolerance(1.e-5)
         else:
            self.__reinitPde = LinearPDE(self.__domain,numEquations=1, numSolutions=1)
            # self.__reinitPde.setSolverMethod(solver=LinearPDE.LUMPING)
            self.__reinitPde.setReducedOrderOn()
            self.__reinitPde.setSymmetryOn()
            self.__reinitPde.setValue(D=1.)
            self.__reinitPde.setTolerance(1.e-2)
            # self.__reinitPde.setSolverMethod(preconditioner=LinearPDE.ILU0)
         #=======================================
         self.__updateInterface()
         print(("phi range:",inf(phi), sup(phi)))

     def setFixedRegion(self,mask,contour=0):
         q=whereNonPositive(abs(self.__phi-contour)-self.__h)*mask
         # self.__fc.setValue(q=q)
         self.__reinitPde.setValue(q=q)


     def __updateInterface(self):
         self.__smoothed_char=self.__makeInterface(self.__phi)


     def update(self,dt):
         """
         Sets a new velocity and updates the level set function.

         :param dt: time step forward
         """
         if self.__FC:
            self.__phi=self.__fc.solve(dt, verbose=False)-self.__diam
         else:
            self.__fcpde.setValue(Y = self.__phi-(dt/2.)*inner(self.velocity,grad(self.__phi)))
            phi_half = self.__fcpde.getSolution()
            self.__fcpde.setValue(Y = self.__phi-dt*inner(self.velocity,grad(phi_half)))
            self.__phi= self.__fcpde.getSolution()
         self.__update_count += 1
         if self.__update_count%self.__reinit_after == 0:
            self.__phi=self.__reinitialise(self.__phi)
            if self.__FC:
                self.__fc.setInitialSolution(self.__phi+self.__diam)
         self.__updateInterface()

     def setTolerance(self,tolerance=1e-3):
         if self.__FC:
            self.__fc.setTolerance(tolerance)

     def __reinitialise(self,phi):
         print("reintialization started:")
         s=self.__makeInterface(phi,1.)
         print(("phi range:",inf(phi), sup(phi)))
         g=grad(phi)
         w=s*g/(length(g)+EPSILON)
         #=======================================================
         # # positive part:
         # phi_p=wherePositive(phi)*phi
         # self.__reinitfc.setInitialSolution(phi_p)
         # self.__reinitfc.setValue(C=-wherePositive(s)*w,Y=wherePositive(s)*s,q=whereNonPositive(phi))
         # dtau=self.__h
         # print "step size: dt (pos)= ",dtau
         # print "phi_p range:",inf(phi_p), sup(phi_p)
         # iter=0
         # while (iter<=self.__reinit_max):
         # phi_p=self.__reinitfc.solve(dtau)
         # print "phi_p range:",inf(phi_p), sup(phi_p)
         # iter+=1
         # # negative part:
         # phi_n=-whereNegative(phi)*phi
         # self.__reinitfc.setInitialSolution(phi_n)
         # self.__reinitfc.setValue(C=-whereNegative(s)*w,Y=-whereNegative(s)*s,q=whereNonNegative(phi))
         # # dtau=self.__reinitfc.getSafeTimeStepSize()
         # dtau=self.__h
         # print "step size: dt (neg)= ",dtau
         # print "phi_n range:",inf(phi_n), sup(phi_n)
         # iter=0
         # while (iter<=self.__reinit_max):
         # phi_n=self.__reinitfc.solve(dtau)
         # print "phi_n range:",inf(phi_n), sup(phi_n)
         # iter+=1
         # phi=phi_p-phi_n
         # print "phi range:",inf(phi), sup(phi)
         # print "reintialization completed."
         # return phi

         #=======================================================
         if self.__reinitFC:
             self.__reinitfc.setValue(C=-w,Y=s)
             self.__reinitfc.setInitialSolution(phi+self.__diam)
             dtau=self.__reinitfc.getSafeTimeStepSize()
             # dtau = 0.5*inf(Function(self.__domain).getSize())
         else:
             dtau = 0.5*inf(Function(self.__domain).getSize())
         print(("step size: dt = ",dtau,inf(abs(phi.interpolate(Function(self.__domain)))/abs(s-inner(w,grad(phi))))))
         iter =0
         # self.__reinitPde.setValue(q=whereNegative(abs(phi)-2*self.__h), r=phi)
         # self.__reinitPde.setValue(r=phi)
         while (iter<=self.__reinit_max):
                 phi_old=phi
                 if self.__reinitFC:
                   phi = self.__reinitfc.solve(dtau)-self.__diam
                 else:
                   self.__reinitPde.setValue(Y = phi+(dtau/2.)*(s-inner(w,grad(phi))))
                   phi_half = self.__reinitPde.getSolution()
                   self.__reinitPde.setValue(Y = phi+dtau*(s-inner(w,grad(phi_half))))
                   # g=grad(phi)
                   # S=inner(w,grad(phi))
                   # self.__reinitPde.setValue(Y = phi+dtau*(s-S),X=dtau**2/2*w*(s-S))
                   phi = self.__reinitPde.getSolution()
                 change = Lsup(phi-phi_old)/self.__diam
                 print(("phi range:",inf(phi), sup(phi)))
                 print(("iteration :", iter, " change:", change))
                 iter +=1
         print("reintialization completed.")
         return phi

     def createParameter(self,value_negative=-1.,value_positive=1):
         out = (value_negative+value_positive)/2. + self.__smoothed_char * ((value_positive-value_negative)/2.)
         return out


     #=========== things from here onwards are not used nor tested: ===========


     def getCharacteristicFunction(self):
         return self.__smoothed_char


     def RK2(self,L):
           k0=L(phi)
           phi_1=phi+dt*k0
           k1=L(phi_1)
           phi=phi+dt/2*(k0+k1)

     def RK3(self,L):
           k0=L(phi)
           phi_1=phi+dt*L(phi)
           k1=L(phi_1)
           phi_2=phi+dt/4*(k0+k1)
           k2=L(phi_2)
           phi=phi+dt/6*(k0+4*k2+K1)

     def TG(self,L):
           k0=L(phi)
           phi_1=phi+dt/2*k0
           k1=L(phi_1)
           phi=phi+dt*k1


     def __reinitialise_old(self):
        #=============================================
        f=0.1
        f2=0.3
        h=self.__h
        fs=h.getFunctionSpace()
        vol=self.getVolumeOfNegativeDomain()
        self.__reinitPde.setValue(D=1.0)
        #============
        grad_phi=grad(self.__phi,fs)
        len_grad_phi=length(grad_phi)
        w=grad_phi/len_grad_phi

        # s=wherePositive(self. __makeInterface(2.))-whereNegative(self. __makeInterface(2.))
        s=self. __makeInterface(2.)
        # s2=whereNegative(len_grad_phi-f2)*(1-(1-len_grad_phi/f2)**2)+(1-whereNegative(len_grad_phi-f2))
        s2=1.
        s=s*s2
        dtau=f*inf(h/abs(s))

        #========================
        diff=1.
        d=1.
        c =0
        #==============================================
        TVD=integrate(length(grad_phi))
        print(("initial range ",inf(self.__phi),sup(self.__phi),"error:",Lsup(1.-len_grad_phi),"volume =",vol,TVD))
        # saveVTK("test.%s.vtu"%c,l=length(grad(self.__phi,fs))-1,s=s,phi=self.__phi)

        dtau=f*inf(h/abs(s))
        while c < self.__reinit_max: # and abs(diff) >= 0.01:
          #
          grad_phi=grad(self.__phi,fs)
          len_grad_phi=length(grad_phi)
          #
          # s=self.__makeInterface(1.)
          # s2=whereNegative(len_grad_phi-f2)*(1-(1-len_grad_phi/f2)**2)+(1-whereNegative(len_grad_phi-f2))
          # s=s*s2
          # phi_on_h=interpolate(self.__phi,fs)
          # self.__reinitPde.setValue(Y =self.__phi+dtau/2*s*(1.-len_grad_phi))
          # phi_half=self.__reinitPde.getSolution()
          # self.__reinitPde.setValue(Y =self.__phi+dtau*s*(1.-length(grad(phi_half,fs))))
          # self.__reinitPde.setValue(Y =self.__phi, Y_reduced=dtau*s*(1.-inner(w,grad(phi_half,ReducedFunction(self.__domain)))))
          # self.__reinitPde.setValue(Y =self.__phi+dtau*(s-inner(w,grad_phi)),X=dtau/2*h*(s-inner(w,grad_phi))*w)
          # self.__reinitPde.setValue(Y =self.__phi+dtau*s*(1.-len_grad_phi),X=-dtau/2*h*s**2*grad_phi)
          self.__reinitPde.setValue(Y=self.__phi+dtau*s*(1.-len_grad_phi),X=f*dtau/2*h*abs(s)*(1.-len_grad_phi)*grad_phi/len_grad_phi)

          # self.__reinitPde.setValue(Y=self.__phi+dtau*s*(1.-len_grad_phi),X=f*dtau/2*h*abs(s)*(1.-len_grad_phi)*grad_phi/len_grad_phi)
          self.__phi, previous = self.__reinitPde.getSolution(), self.__phi
          # self.__updateInterface()

          vol,vol_old=self.getVolumeOfNegativeDomain(),vol
          diff=(vol-vol_old)
          r=Lsup(length(grad(self.__phi))-1.)
          TVD=integrate(length(grad(self.__phi,fs)))
          print(("iteration :", c, "range ",inf(self.__phi),sup(self.__phi),"error :",r,"volume change:",diff,TVD))
          # saveVTK("test.%s.vtu"%(c+1),l=length(grad(self.__phi,fs)),s=s,phi=self.__phi,v=grad(self.__phi,fs))
          c += 1
        return
        #==============================================
        f2=0.7
        h=self.__h
        fs=h.getFunctionSpace()
        vol=self.getVolumeOfNegativeDomain()
        s=abs(self. __makeInterface(2.))
        grad_phi=grad(self.__phi,fs)
        len_grad_phi=length(grad_phi)

        #----
        # aphi_on_h=abs(interpolate(self.__phi,self.__h.getFunctionSpace()))
        # q=2*self.__h
        # mask_neg = whereNegative(aphi_on_h-q)
        # mask_pos = wherePositive(aphi_on_h-q-self.__h)
        # mask_interface = 1.-mask_neg-mask_pos
        # s2=mask_neg + mask_interface * (q+self.__h-aphi_on_h)/self.__h
        #----
        m=whereNonPositive(len_grad_phi-f2)
        s2=m*len_grad_phi/f2+(1.-m)
        #----
        # s2=1.
        #----
        self.__reinitPde.setValue(D=(1-s)/h**2,Y=(1-s)/h**2*self.__phi)
        self.__reinitPde.setValue(A=(s2*len_grad_phi+(1.-s2))*kronecker(3),X=grad_phi)
        # self.__reinitPde.setValue(A=(len_grad_phi-1)/len_grad_phi*kronecker(3))
        # self.__reinitPde.setValue(A=kronecker(3), X=grad_phi/len_grad_phi)
        self.__phi = self.__reinitPde.getSolution()
        r=Lsup(length(grad(self.__phi))-1.)
        vol,vol_old=self.getVolumeOfNegativeDomain(),vol
        diff=(vol-vol_old)/vol
        print(("iteration :", inf(self.__phi),sup(self.__phi),r,diff))
        # saveVTK("test.%s.vtu"%0,l=length(grad(self.__phi,fs)),s=s,phi=self.__phi,v=grad(self.__phi,fs),s2=s2)
        return
        #=============================================

        #============


     def getVolumeOfNegativeDomain(self):
         """
         Returns the current volume of domain with phi<0.
         """
         return integrate((1.-self.__makeInterface(1.))/2.)


     def getBoundingBoxOfNegativeDomain(self):
         """
         Returns the height of the region with phi<0.
         """
         fs=self.__h.getFunctionSpace()
         mask_phi1=wherePositive(interpolate(self.__phi,fs))
         mask_phi2=wherePositive(self.__phi)
         x1=fs.getX()
         x2=self.__domain.getX()
         out=[]
         for i in range(fs.getDim()):
             x1_i=x1[i]
             x2_i=x2[i]
             d=2*(sup(x2_i)-inf(x2_i))
             offset1=d*mask_phi1
             offset2=d*mask_phi2
             out.append((min(inf(x1_i+offset1),inf(x2_i+offset2)),max(sup(x1_i-offset1),sup(x2_i-offset2))))
         return tuple(out)

