"""
seismic wave propagation
                                                                                                                                                                                                     
@var __author__: name of author
@var __licence__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""
                                                                                                                                                                                                     
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Lutz Gross, l.gross@uq.edu.au"
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision$"
__date__="$Date$"

from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Brick
import time

output=False
n_end=100

resolution=4000.  # number of elements per m
l=80000.          # width and length m
h=30000.          # height in m

rho_bedrock=8e3
mu_bedrock=1.7e11
lambda_bedrock=1.7e11

l_x_water=10000    # length of water in x
l_y_water=10000    # length of water in y
h_water=max(2000,resolution)       # depth of water region

water_tag=2
rho_water=1e3
mu_water=0.
lambda_water=1.e9

d0_sand=10000      # thickness of sand region on surface
d_sand=max(2000,resolution)        # thickness of sand layer under the water

sand_tag=3
rho_sand=5e3
mu_sand=1.5e10
lambda_sand=1.5e10


# location and geometrical size of event:
xc=[30000.,40000.,10000.]
src_radius = 0.1*h
# direction of event:
event=numarray.array([1.,0.,0.])*1.e6
# time and length of the event
tc_length=2.
tc=sqrt(5.*tc_length)
if output:
   print "event location = ",xc
   print "radius of event = ",src_radius
   print "time of event = ",tc
   print "length of event = ",tc_length
   print "direction = ",event

t_end=20.

def getDomain():
    """
    this defines a dom as a brick of length and width l and hight h

      
    """
    global netotal
    ne_l=int(l/resolution+0.5)
    ne_h=int(h/resolution+0.5)
    netotal=ne_l*ne_l*ne_h
    if output: print "grid : %s x %s x %s"%(ne_l,ne_l,ne_h)
    dom=Brick(ne_l,ne_l,ne_h,l0=l,l1=l,l2=h,order=1)

    fs=Function(dom)
    x=Function(dom).getX()
    fs.setTags(sand_tag,wherePositive(x[0]-(l-l_x_water-d0_sand)) \
                       *wherePositive(x[1]-(l-l_y_water-d0_sand)) \
                       *wherePositive(x[2]-(h-h_water-d_sand)))
    fs.setTags(water_tag,wherePositive(x[0]-(l-l_x_water)) \
                        *wherePositive(x[1]-(l-l_y_water)) \
                        *wherePositive(x[2]-(h-h_water)))
    return dom

def getMaterialProperties(dom):
   rho        =Scalar(rho_bedrock,Function(dom))
   rho.setTaggedValue(sand_tag,rho_sand)
   rho.setTaggedValue(water_tag,rho_water)

   lame_mu    =Scalar(mu_bedrock,Function(dom))
   lame_mu.setTaggedValue(sand_tag,mu_sand)
   lame_mu.setTaggedValue(water_tag,mu_water)

   lame_lambda=Scalar(lambda_bedrock,Function(dom))
   lame_lambda.setTaggedValue(sand_tag,lambda_sand)
   lame_lambda.setTaggedValue(water_tag,lambda_water)

   return rho,lame_mu,lame_lambda


def wavePropagation(dom,rho,lame_mu,lame_lambda):
   x=Function(dom).getX()
   # ... open new PDE ...
   mypde=LinearPDE(dom)
   mypde.setSolverMethod(LinearPDE.LUMPING)
   k=kronecker(Function(dom))
   mypde.setValue(D=k*rho)

   v_p=sqrt((2*lame_mu+lame_lambda)/rho)
   if output: print "v_p=",v_p
   v_s=sqrt(lame_mu/rho)
   if output: print "v_s=",v_s
   dt=(1./5.)*inf(dom.getSize()/v_p)
   if output: print "time step size = ",dt
   # ... set initial values ....
   n=0
   t=0
   # initial value of displacement at point source is constant (U0=0.01)
   # for first two time steps
   u     =Vector(0.,Solution(dom))
   u_last=Vector(0.,Solution(dom))

   starttime = time.clock()
   while t<t_end and n<n_end:
     if output: print n+1,"-th time step t ",t+dt," max u and F: ",Lsup(u),
     # ... get current stress ....
     eps=symmetric(grad(u))
     stress=lame_lambda*trace(eps)*k+2*lame_mu*eps
     # ... force due to event:
     F=exp(-((t-tc)/tc_length)**2)*exp(-(length(x-xc)/src_radius)**2)*event
     if output: print Lsup(F)
     # ... get new acceleration ....
     mypde.setValue(X=-stress,Y=F)
     a=mypde.getSolution()
     # ... get new displacement ...
     u_new=2*u-u_last+dt**2*a
     # ... shift displacements ....
     u_last,u=u,u_new
     # ... save current acceleration in units of gravity and displacements 
     if output:
          if n%10==0: saveVTK("disp.%i.vtu"%(n/10),displacement=u, amplitude=length(u))

     t+=dt
     n+=1

   endtime = time.clock()
   totaltime = endtime-starttime
   global netotal
   print ">>number of elements: %s, total time: %s, per time step: %s <<"%(netotal,totaltime,totaltime/n)
if __name__ =="__main__":
   dom=getDomain()
   rho,lame_mu,lame_lambda=getMaterialProperties(dom)
   wavePropagation(dom,rho,lame_mu,lame_lambda)
