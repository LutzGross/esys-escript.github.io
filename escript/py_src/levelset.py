
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

from esys.escript import *
from esys.escript.linearPDEs import LinearPDE, TransportPDE
from esys.escript.pdetools import Projector
from esys import finley
import math



from esys.escript import *
from esys.finley import finley
from esys.escript.linearPDEs import LinearPDE
from esys.escript.pdetools import Projector
import sys


class LevelSet:
  """
  The level set method tracking an interface defined by the zero contour of the level set function phi. 
  it is assumed the phi(x)<0 defines the volume of interest.

  """
  def __init__(self,domain,phi,reinit_max=10,reinit_each=2,smooth=2.):
    """
    set up the level set method

    @param domain: the domain where the level set is used
    @param phi: the initial level set function
    @param reinit_max: maximum number of reinitalization steps
    @param reinit_each: phi is reinitialized every reinit_each step
    @param smooth: smoothing width
    """
    self.__domain = domain
    self.__phi = phi
    self.__reinit_max = reinit_max
    self.__reinit_each = reinit_each
    self.__PDE = LinearPDE(domain)
    self.__PDE.setReducedOrderOn()
    self.__PDE.setValue(D=1.0)
    self.__PDE.setSolverMethod(solver=LinearPDE.PCG)
    self.__reinitPDE = LinearPDE(domain, numEquations=1)
    self.__reinitPDE.setReducedOrderOn()
    self.__reinitPDE.setValue(D=1.0)
    self.__reinitPDE.setSolverMethod(solver=LinearPDE.LUMPING)
    self.__h = inf(domain.getSize())
    self.__smooth = smooth
    self.__n_step=0
  
  def __advect(self, velocity, dt):
    """
    advects the level set function in the presense of a velocity field.

    This implementation uses the 2-step Taylor-Galerkin method 
    @param velocity: velocity field
    @param dt: dime increment
    @return: the advected level set function
    """
    Y = self.__phi-dt/2.*inner(velocity,grad(self.__phi))
    self.__PDE.setValue(Y=Y)    
    phi_half = self.__PDE.getSolution()
    Y = self.__phi-dt*inner(velocity,grad(phi_half))
    self.__PDE.setValue(Y=Y)    
    phi = self.__PDE.getSolution()
    print "LevelSet: Advection step done"
    return phi

  def __reinitialise(self):
    """
    reinializes the level set 

    It solves the 

    @return: reinitalized level set
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
      print "LevelSet:reinitialization: iteration :", iter, " error:", error
      iter +=1
    return phi

  def getTimeStepSize(self,velocity):
       """
       returns a new dt for a given velocity using the courant coundition

       @param velocity: velocity field
       """
       self.__velocity=velocity
       dt=0.5*self.__h/sup(length(velocity))
       return dt
   
  def update(self,dt):
      """
      sets a new velocity and updates the level set fuction

      @param dt: time step forward
      """
      self.__phi=self.__advect(self.__velocity, dt)
      self.__n_step+=1
      if self.__n_step%self.__reinit_each ==0: self.__phi = self.__reinitialise()
      return self.__phi

  def update_phi(self, velocity, dt):
      """
      updates phi under the presense of a velocity field 

      If dt is small this call is equivalent to call 
  
      dt=LevelSet.getTimeStepSize(velocity)
      phi=LevelSet.update(dt)
                                  
      otherwise substepping is used.
      @param velocity: velocity field
      @param dt: time step forward
      """
      dt2=self.getTimeStepSize(velocity)
      n=math.ceil(dt/dt2)
      dt_new=dt/n
      for i in range(n):
           phi=self.update(dt_new)
           t+=dt_new
      return phi


  def getVolume(self):
    """
    return the volume of the phi(x)<0 region
    """
    return integrate(whereNegative(self.__phi.interpolate(Function(self.__domain))))

  def getSurface(self,rel_width_factor=0.5):
    """
    return a mask for phi(x)=1 region 
  
    @param rel_width_factor: relative wideth of region around zero contour.
    """
    return whereNegative(abs(self.__phi)-rel_width_factor*self.__h)
    
  def getH(self):
     """
     returns mesh size
     """
     return self.__h

  def getDomain(self):
     """
     returns domain
     """
     return self.__domain

  def getLevelSetFunction(self):
      """
      returns the level set function
      """
      return self.__phi

  def update_parameter_sharp(self, param_neg=-1, param_pos=1, phi=None):
    """
    creates a function whith param_neg where phi<0 and param_pos where phi>0 (no smoothing)

    @param param_neg: value of parameter on the negative side (phi<0)
    @param param_pos: value of parameter on the positve side (phi>0)
    @param phi: level set funtion to be used. if not present the current level set is used.
    """
    mask_neg = whereNegative(self.__phi)
    mask_pos = whereNonNegative(self.__phi)
    param = param_pos*mask_pos + param_neg*mask_neg
    return param

  def update_parameter(self, param_neg=-1, param_pos=1, phi=None, smoothing_width=None):
    """
    creates a smoothed function whith param_neg where phi<0 and param_pos where phi>0 which is smoothed over a length
    smoothing_width accross the interface

    @param smoothing_width: width of the smoothing zone relative to mesh size. If not present the initial value of C{smooth} is used.
    """
    if smoothing_width==None: smoothing_width = self.__smooth
    if phi==None: phi=self.__phi
    s=self.__makeInterface(phi,smoothing_width)
    return ((param_pos-param_neg)*s+param_pos+param_neg)/2

  def __makeInterface(self,phi,smoothing_width):
      """
      creates a smooth interface from -1 to 1 over the length 2*h*smoothing_width where -1 is used where the level set is negative
      and 1 where the level set is 1
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
      makes a smooth charateristic function of the region phi(x)>contour if positiveSide and phi(x)<contour otherwise.

      @param phi: level set funtion to be used. if not present the current level set is used.
      @param smoothing_width: width of the smoothing zone relative to mesh size. If not present the initial value of C{smooth} is used.
      """
      if phi==None: phi=self.__phi
      if smoothing_width == None: smoothing_width=self.__smooth
      s=self.__makeInterface(phi=phi-contour,smoothing_width=smoothing_width)
      if positiveSide:
          return (1+s)/2
      else:
          return (1-s)/2

  def setTolerance(self,tolerance=1e-3):
    self.__PDE.setTolerance(tolerance)
    self.__reinitPDE.setTolerance(tolerance)


class LevelSet2(object):
     def __init__(self,phi,reinit_max=10,reinit_each=2,smooth=2.):
         """
         initialize model
         """
         self.__domain = phi.getDomain() 
         x=self.__domain.getX()
         diam=0
         for i in range(self.__domain.getDim()):
             xi=x[i]
             diam+=(inf(xi)-sup(xi))**2
         self.__diam=sqrt(diam)
         self.__h = sup(Function(self.__domain).getSize())
         self.__reinit_each=reinit_each
         self.__reinit_max = reinit_max
         self.__smooth = smooth
         self.__phi = phi
         self.__update_count=0
         self.velocity = None

         #==================================================
         self.__FC=True
         if self.__FC:
            self.__fc=TransportPDE(self.__domain,num_equations=1,theta=0.5)
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
            self.__reinitfc=TransportPDE(self.__domain,num_equations=1,theta=1.0)
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
         print "phi range:",inf(phi), sup(phi)

     def setFixedRegion(self,mask,contour=0):
         q=whereNonPositive(abs(self.__phi-contour)-self.__h)*mask
         # self.__fc.setValue(q=q)
         self.__reinitPde.setValue(q=q)


     def __updateInterface(self):
         self.__smoothed_char=self.__makeInterface(self.__phi)




     def update(self,dt):
         """
         sets a new velocity and updates the level set fuction

         @param dt: time step forward
         """
         if self.__FC:
            self.__phi=self.__fc.solve(dt, verbose=False)-self.__diam
         else:
            self.__fcpde.setValue(Y = self.__phi-(dt/2.)*inner(self.velocity,grad(self.__phi)))
            phi_half = self.__fcpde.getSolution()
            self.__fcpde.setValue(Y = self.__phi-dt*inner(self.velocity,grad(phi_half)))
            self.__phi= self.__fcpde.getSolution()
         self.__update_count += 1
         if self.__update_count%self.__reinit_each == 0:
            self.__phi=self.__reinitialise(self.__phi)
            if self.__FC:
                self.__fc.setInitialSolution(self.__phi+self.__diam)
         self.__updateInterface()

     def setTolerance(self,tolerance=1e-3):
         if self.__FC:
            self.__fc.setTolerance(tolerance)

     def __reinitialise(self,phi):
         print "reintialization started:"
         s=self.__makeInterface(phi,1.)
         print "phi range:",inf(phi), sup(phi)
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
         print "step size: dt = ",dtau,inf(abs(phi.interpolate(Function(self.__domain)))/abs(s-inner(w,grad(phi))))
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
                 print "phi range:",inf(phi), sup(phi)
                 print "iteration :", iter, " change:", change
                 iter +=1
         print "reintialization completed."
         return phi

     def createParameter(self,value_negative=-1.,value_positive=1):
         out = (value_negative+value_positive)/2. + self.__smoothed_char * ((value_positive-value_negative)/2.)
         return out


     #================ things from here onwards are not used nor tested: ==========================================
     

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
        print "initial range ",inf(self.__phi),sup(self.__phi),"error:",Lsup(1.-len_grad_phi),"volume =",vol,TVD
        # saveVTK("test.%s.xml"%c,l=length(grad(self.__phi,fs))-1,s=s,phi=self.__phi)

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
          print "iteration :", c, "range ",inf(self.__phi),sup(self.__phi),"error :",r,"volume change:",diff,TVD
          # saveVTK("test.%s.xml"%(c+1),l=length(grad(self.__phi,fs)),s=s,phi=self.__phi,v=grad(self.__phi,fs))
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
        print "iteration :", inf(self.__phi),sup(self.__phi),r,diff
        # saveVTK("test.%s.xml"%0,l=length(grad(self.__phi,fs)),s=s,phi=self.__phi,v=grad(self.__phi,fs),s2=s2)
        return
        #=============================================

        #============


     def getVolumeOfNegativeDomain(self):
         """
         return the current volume of domain with phi<0.
         """
         return integrate((1.-self.__makeInterface(1.))/2.)

              
     def getBoundingBoxOfNegativeDomain(self):
         """
         get the height of the region with  phi<0
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

