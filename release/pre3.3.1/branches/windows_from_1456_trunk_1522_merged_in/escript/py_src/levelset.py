from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.escript.pdetools import Projector
from esys import finley
import math



class LevelSet(object):
     def __init__(self,phi,reinitialization_steps_max=10,reinitialization_each=2,relative_smoothing_width=2., location_fixed_phi=Data()):
         """
         initialize model
         """
         self.domain = phi.getDomain() 
         self.__h = Function(self.domain).getSize()
         self.__location_fixed_phi=location_fixed_phi

         self.__relative_smoothing_width = relative_smoothing_width
         self.__pde = LinearPDE(self.domain)         
         self.__pde.setReducedOrderOn()
         self.__pde.setValue(D=1.0,q=self.__location_fixed_phi)
         # self.__pde.setTolerance(1.e-10)
         #self.__pde.setSolverMethod(solver=LinearPDE.DIRECT)
                     
         self.__reinitPde = LinearPDE(self.domain,numEquations=1)
         # self.__reinitPde.setSolverMethod(solver=LinearPDE.LUMPING)       
         self.__reinitPde.setReducedOrderOn()
         self.__pde.setValue(q=self.__location_fixed_phi)
         self.__reinitPde.setSymmetryOn()
         self.__reinitPde.setTolerance(1.e-8)
         # self.__reinitPde.setSolverMethod(preconditioner=LinearPDE.ILU0)
        
         self.__reinitialization_each=reinitialization_each
         self.__reinitialization_steps_max = reinitialization_steps_max
         print self.__reinitialization_steps_max
         self.__update_count=0


         self.__phi = phi
         self.__updateInterface()
         self.velocity = None

     def getDomain(self):
         return self.domain
     
    
     def __updateInterface(self):
         self.__smoothed_char=self.__makeInterface(self.__relative_smoothing_width)

     def __makeInterface(self,smoothing_width=1.):
         # s=smoothing_width*self.__h*length(grad(self.__phi,self.__h.getFunctionSpace()))
         s=smoothing_width*self.__h
         phi_on_h=interpolate(self.__phi,self.__h.getFunctionSpace())         
         mask_neg = whereNonNegative(-s-phi_on_h)
         mask_pos = whereNonNegative(phi_on_h-s)
         mask_interface = 1.-mask_neg-mask_pos
         interface=1.-(phi_on_h-s)**2/(2.*s**3)*(phi_on_h+2*s)  # function f with f(s)=1, f'(s)=0, f(-s)=-1, f'(-s)=0, f(0)=0, f'(0)=
         # interface=phi_on_h/s
         return - mask_neg + mask_pos + mask_interface * interface

     def getCharacteristicFunction(self):
         return self.__smoothed_char

     def update(self,dt):
         """
         sets a new velocity and updates the level set fuction

         @param dt: time step forward
         """
         self.__advectPhi(dt)
         self.__updateInterface()
         if self.__update_count%self.__reinitialization_each == 0:
            self.__reinitialise()
            self.__updateInterface()
         self.__update_count += 1

     def __advectPhi(self,dt):
        """
        advects the level set function 
        """
        fdf=0.01
        grad_phi = grad(self.__phi)   
        len_grad_phi=length(grad_phi)
        f=whereNonNegative(len_grad_phi-fdf)
        s=(1.-f)*len_grad_phi/fdf+f
        
        # self.__pde.setValue(Y=self.__phi-dt/2.*s*inner(self.velocity,grad_phi),r=self.__phi)   
        # phi_half = self.__pde.getSolution()
        # self.__pde.setValue(Y=self.__phi-dt*s*inner(self.velocity,grad(phi_half)),r=self.__phi)   
        # self.__phi = self.__pde.getSolution()

        self.__pde.setValue(Y=s*inner(self.velocity,grad_phi))
        phi_half = self.__phi-dt/2*self.__pde.getSolution()
        self.__pde.setValue(Y=s*inner(self.velocity,grad(phi_half)))
        self.__phi =self.__phi-dt*self.__pde.getSolution()
        print("level set advection done")

     def getTimeStepSize(self,velocity,safety_factor=0.6):
         """
         returns a new dt for a given velocity using the courant coundition
         """
         self.velocity=velocity
         return safety_factor*inf(self.__h/length(interpolate(velocity,self.__h.getFunctionSpace())))

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
           
           
     def __reinitialise(self):
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
        while c < self.__reinitialization_steps_max: # and abs(diff) >= 0.01:
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
          # self.__reinitPde.setValue(Y =self.__phi, Y_reduced=dtau*s*(1.-inner(w,grad(phi_half,ReducedFunction(self.domain)))))
          # self.__reinitPde.setValue(Y =self.__phi+dtau*(s-inner(w,grad_phi)),X=dtau/2*h*(s-inner(w,grad_phi))*w)
          # self.__reinitPde.setValue(Y =self.__phi+dtau*s*(1.-len_grad_phi),X=-dtau/2*h*s**2*grad_phi)
          self.__reinitPde.setValue(Y=self.__phi+dtau*s*(1.-len_grad_phi),X=f*dtau/2*h*abs(s)*(1.-len_grad_phi)*grad_phi/len_grad_phi)

          # self.__reinitPde.setValue(Y=self.__phi+dtau*s*(1.-len_grad_phi),X=f*dtau/2*h*abs(s)*(1.-len_grad_phi)*grad_phi/len_grad_phi)
          self.__phi, previous = self.__reinitPde.getSolution(), self.__phi
          self.__updateInterface()

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

     def getLevelSetFunction(self):
         """
         returns the level set function
         """
         return self.__phi
              
     def getBoundingBoxOfNegativeDomain(self):
         """
         get the height of the region with  phi<0
         """
         fs=self.__h.getFunctionSpace()
         mask_phi1=wherePositive(interpolate(self.__phi,fs))
         mask_phi2=wherePositive(self.__phi)
         x1=fs.getX()
         x2=self.domain.getX()
         out=[]
         for i in range(fs.getDim()):
             x1_i=x1[i]
             x2_i=x2[i]
             d=2*(sup(x2_i)-inf(x2_i))
             offset1=d*mask_phi1
             offset2=d*mask_phi2
             out.append((min(inf(x1_i+offset1),inf(x2_i+offset2)),max(sup(x1_i-offset1),sup(x2_i-offset2))))
         return tuple(out)

     def createParameter(self,value_negative=-1.,value_positive=1):
         out = (value_negative+value_positive)/2. + self.__smoothed_char * ((value_positive-value_negative)/2.)
         return out
