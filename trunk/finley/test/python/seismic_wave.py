
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
seismic wave propagation

:var __author__: name of author
:var __licence__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Brick
from esys.weipa import saveVTK
import time

WORKDIR="./waves/"
output=True
n_end=10000

resolution=1000.  # number of elements per m in the finest region
resolution=400.  # number of elements per m in the finest region
o=1              # element order

l=100000.           # width and length m (without obsorber)
h=30000.            # height in m        (without obsorber)
d_absorber=l*0.10   # thickness of absorbing layer

l_sand=20000.          # thickness of sand region on surface
h_sand=5000.           # thickness of sand layer under the water

l_x_water=10000.       # length of water in x
l_y_water=10000.       # length of water in y
h_water=2000.          # depth of water region

x_sand=l/2-l_x_water/2-l_sand # x coordinate of location of sand region (without obsorber)
y_sand=l/2-l_y_water/2-l_sand # y coordinate of location of sand region (without obsorber)


# origin
origin={"x": -d_absorber, "y" : -d_absorber , "z" : -h-d_absorber }
# location and geometrical size of event reltive to origin:
xc=[l*0.2,l*0.3,-h*0.7]
src_radius  = 2*resolution
# direction of event:
event=numpy.array([0.,0.,1.])*1.e6
# time and length of the event
tc=2.
tc_length=0.5

# material properties:
bedrock=0
absorber=1
water=2
sand=3

rho_tab={}
rho_tab[bedrock]=8e3
rho_tab[absorber]=rho_tab[bedrock]
rho_tab[water]=1e3
rho_tab[sand]=5e3

mu_tab={}
mu_tab[bedrock]=1.7e11
mu_tab[absorber]=mu_tab[bedrock]
mu_tab[water]=0.
mu_tab[sand]=1.5e10

lmbd_tab={}
lmbd_tab[bedrock]=1.7e11
lmbd_tab[absorber]=lmbd_tab[bedrock]
lmbd_tab[water]=1.e9
lmbd_tab[sand]=1.5e10

eta_tab={}
eta_tab[absorber]=-log(0.05)*sqrt(rho_tab[absorber]*(lmbd_tab[absorber]+2*mu_tab[absorber]))/d_absorber
eta_tab[sand]=eta_tab[absorber]/40.
eta_tab[water]=eta_tab[absorber]/40.
eta_tab[bedrock]=eta_tab[absorber]/40.


# material properties:
bedrock=0
absorber=1
water=2
sand=3

rho={}
rho[bedrock]=8e3
rho[absorber]=rho[bedrock]
rho[water]=1e3
rho[sand]=5e3

mu={}
mu[bedrock]=1.7e11
mu[absorber]=mu[bedrock]
mu[water]=0.
mu[sand]=1.5e10

lmbd={}
lmbd[bedrock]=1.7e11
lmbd_absorber=lmbd[bedrock]
lmbd[water]=1.e9
lmbd[sand]=1.5e10

eta={}
eta[absorber]=-log(0.05)*sqrt(rho[absorber]*(lmbd_absorber+2*mu[absorber]))/d_absorber
eta[sand]=eta[absorber]/40.
eta[water]=eta[absorber]/40.
eta[bedrock]=eta[absorber]/40.

if output:
   print("event location = ",xc)
   print("radius of event = ",src_radius)
   print("time of event = ",tc)
   print("length of event = ",tc_length)
   print("direction = ",event)

t_end=30.
dt_write=0.1


def getDomain():
    """
    this defines a dom as a brick of length and width l and hight h

      
    """
    global netotal
    
    v_p={}
    for tag in sorted(rho_tab.keys()):
       v_p[tag]=sqrt((2*mu_tab[tag]+lmbd_tab[tag])/rho_tab[tag])
    v_p_ref=min(v_p.values())
    print("velocities: bedrock = %s, sand = %s, water =%s, absorber =%s, reference =%s"%(v_p[bedrock],v_p[sand],v_p[water],v_p[absorber],v_p_ref))

    sections={}
    sections["x"]=[d_absorber, x_sand, l_sand, l_x_water, l_sand, l-x_sand-2*l_sand-l_x_water, d_absorber]
    sections["y"]=[d_absorber, y_sand, l_sand, l_y_water, l_sand, l-y_sand-2*l_sand-l_y_water, d_absorber]
    sections["z"]=[d_absorber,h-h_water-h_sand,h_sand,h_water]
    if output:
      print("sections x = ",sections["x"])
      print("sections y = ",sections["y"])
      print("sections z = ",sections["z"])

    mats= [ 
            [ [absorber, absorber, absorber, absorber, absorber, absorber, absorber],
              [absorber, absorber, absorber, absorber, absorber, absorber, absorber],
              [absorber, absorber, absorber, absorber, absorber, absorber, absorber],
              [absorber, absorber, absorber, absorber, absorber, absorber, absorber],
              [absorber, absorber, absorber, absorber, absorber, absorber, absorber],
              [absorber, absorber, absorber, absorber, absorber, absorber, absorber],
              [absorber, absorber, absorber, absorber, absorber, absorber, absorber] ],

            [ [absorber, absorber, absorber, absorber, absorber, absorber, absorber],
              [absorber, bedrock , bedrock , bedrock , bedrock , bedrock , absorber],
              [absorber, bedrock , bedrock , bedrock , bedrock , bedrock , absorber],
              [absorber, bedrock , bedrock , bedrock , bedrock , bedrock , absorber],
              [absorber, bedrock , bedrock , bedrock , bedrock , bedrock , absorber],
              [absorber, bedrock , bedrock , bedrock , bedrock , bedrock , absorber],
              [absorber, absorber, absorber, absorber, absorber, absorber, absorber] ],

            [ [absorber, absorber, absorber, absorber, absorber, absorber, absorber],
              [absorber, bedrock , bedrock , bedrock , bedrock , bedrock , absorber],
              [absorber, bedrock , sand    , sand    , sand    , bedrock , absorber],
              [absorber, bedrock , sand    , sand    , sand    , bedrock , absorber],
              [absorber, bedrock , sand    , sand    , sand    , bedrock , absorber],
              [absorber, bedrock , bedrock , bedrock , bedrock , bedrock , absorber],
              [absorber, absorber, absorber, absorber, absorber, absorber, absorber] ],

            [ [absorber, absorber, absorber, absorber, absorber, absorber, absorber],
              [absorber, bedrock , bedrock , bedrock , bedrock , bedrock , absorber],
              [absorber, bedrock , sand    , sand    , sand    , bedrock , absorber],
              [absorber, bedrock , sand    , water   , sand    , bedrock , absorber],
              [absorber, bedrock , sand    , sand    , sand    , bedrock , absorber],
              [absorber, bedrock , bedrock , bedrock , bedrock , bedrock , absorber],
              [absorber, absorber, absorber, absorber, absorber, absorber, absorber] ] ]
    
    num_elem={}
    for d in sections:
       num_elem[d]=[]
       for i in range(len(sections[d])):
           if d=="x":
              v_p_min=v_p[mats[0][0][i]]
              for q in range(len(sections["y"])):
                 for r in range(len(sections["z"])):
                    v_p_min=min(v_p[mats[r][q][i]],v_p_min)
           elif d=="y":
              v_p_min=v_p[mats[0][i][0]]
              for q in range(len(sections["x"])):
                 for r in range(len(sections["z"])):
                    v_p_min=min(v_p[mats[r][i][q]],v_p_min)
           elif d=="z":
              v_p_min=v_p[mats[i][0][0]]
              for q in range(len(sections["x"])):
                 for r in range(len(sections["y"])):
                    v_p_min=min(v_p[mats[i][r][q]],v_p_min)
           num_elem[d].append(max(1,int(sections[d][i] * v_p_ref/v_p_min /resolution+0.5)))
       
    ne_x=sum(num_elem["x"])
    ne_y=sum(num_elem["y"])
    ne_z=sum(num_elem["z"])
    netotal=ne_x*ne_y*ne_z
    if output: print("grid : %s x %s x %s (%s elements)"%(ne_x,ne_y,ne_z,netotal))
    dom=Brick(ne_x,ne_y,ne_z,l0=o*ne_x,l1=o*ne_y,l2=o*ne_z,order=o)
    x_old=dom.getX()
    x_new=0

    for d in sections:
        if d=="x": 
            i=0
            f=[1,0,0]
        if d=="y": 
            i=1
            f=[0,1,0]
        if d=="z": 
            i=2
            f=[0,0,1]
        x=x_old[i]

        p=origin[d]
        ne=0
        s=0.
 
        for i in range(len(sections[d])-1):
            msk=whereNonPositive(x-o*ne+0.5)
            s=s*msk + (sections[d][i]/(o*num_elem[d][i])*(x-o*ne)+p)*(1.-msk)
            ne+=num_elem[d][i]
            p+=sections[d][i]
        x_new=x_new + s * f
    dom.setX(x_new)

    fs=Function(dom)
    x=Function(dom).getX()
    x0=x[0]
    x1=x[1]
    x2=x[2]
    p_z=origin["z"]
    for i in range(len(mats)):
       f_z=wherePositive(x2-p_z)*wherePositive(x2-p_z+sections["z"][i])
       p_y=origin["y"]
       for j in range(len(mats[i])):
         f_y=wherePositive(x1-p_y)*wherePositive(x1-p_z+sections["y"][j])
         p_x=origin["x"]
         for k in range(len(mats[i][j])):
             f_x=wherePositive(x0-p_x)*wherePositive(x0-p_x+sections["x"][k]) 
             fs.setTags(mats[i][j][k],f_x*f_y*f_z)
             p_x+=sections["x"][k]
         p_y+=sections["y"][j]
       p_z+=sections["z"][i]
    return dom

def getMaterialProperties(dom):
   rho =Scalar(rho_tab[bedrock],Function(dom))
   eta =Scalar(eta_tab[bedrock],Function(dom))
   mu  =Scalar(mu_tab[bedrock],Function(dom))
   lmbd=Scalar(lmbd_tab[bedrock],Function(dom))
   tags=Scalar(bedrock,Function(dom))
   
   for tag in sorted(rho_tab.keys()):
      rho.setTaggedValue(tag,rho_tab[tag])
      eta.setTaggedValue(tag,eta_tab[tag])
      mu.setTaggedValue(tag,mu_tab[tag])
      lmbd.setTaggedValue(tag,lmbd_tab[tag])
      tags.setTaggedValue(tag,tag)
   return rho,mu,lmbd,eta

def wavePropagation(dom,rho,mu,lmbd,eta):
   x=Function(dom).getX()
   # ... open new PDE ...
   mypde=LinearPDE(dom)
   mypde.setSolverMethod(LinearPDE.LUMPING)
   k=kronecker(Function(dom))
   mypde.setValue(D=k*rho)

   dt=(1./5.)*inf(dom.getSize()/sqrt((2*mu+lmbd)/rho))
   if output: print("time step size = ",dt)
   # ... set initial values ....
   n=0
   t=0
   t_write=0.
   n_write=0
   # initial value of displacement at point source is constant (U0=0.01)
   # for first two time steps
   u=Vector(0.,Solution(dom))
   v=Vector(0.,Solution(dom))
   a=Vector(0.,Solution(dom))
   a2=Vector(0.,Solution(dom))
   v=Vector(0.,Solution(dom))

   if not os.path.isdir(WORKDIR): os.mkdir(WORKDIR)

   starttime = time.clock()
   while t<t_end and n<n_end:
     if output: print(n+1,"-th time step t ",t+dt," max u and F: ",Lsup(u), end=' ')
     # prediction:
     u_pr=u+dt*v+(dt**2/2)*a+(dt**3/6)*a2
     v_pr=v+dt*a+(dt**2/2)*a2
     a_pr=a+dt*a2
     # ... get current stress ....
     eps=symmetric(grad(u_pr))
     stress=lmbd*trace(eps)*k+2*mu*eps
     # ... force due to event:
     if abs(t-tc)<5*tc_length:
        F=exp(-((t-tc)/tc_length)**2)*exp(-(length(x-xc)/src_radius)**2)*event
        if output: print(Lsup(F))
     else:
        if output: print(0.)
     # ... get new acceleration ....
     mypde.setValue(X=-stress,Y=F-eta*v_pr)
     a=mypde.getSolution()
     # ... get new displacement ...
     da=a-a_pr
     u=u_pr+(dt**2/12.)*da
     v=v_pr+(5*dt/12.)*da
     a2+=da/dt
     # ... save current acceleration in units of gravity and displacements 
     if output:
          if t>=t_write: 
             saveVTK(os.path.join(WORKDIR,"disp.%i.vtu"%n_write),displacement=u, amplitude=length(u))
             t_write+=dt_write
             n_write+=1
     t+=dt
     n+=1

   endtime = time.clock()
   totaltime = endtime-starttime
   global netotal
   print(">>number of elements: %s, total time: %s, per time step: %s <<"%(netotal,totaltime,totaltime/n))
if __name__ =="__main__":
   dom=getDomain()
   rho,mu,lmbd,eta=getMaterialProperties(dom)
   wavePropagation(dom,rho,mu,lmbd,eta)
