
########################################################
#
# Copyright (c) 2009-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2009-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Author: Antony Hallam antony.hallam@uqconnect.edu.au
"""


############################################################FILE HEADER
# example09m.py
# Create a simple 3D model for use in example09.
#
#######################################################EXTERNAL MODULES
from esys.pycad import * #domain constructor
from esys.pycad.gmsh import Design #Finite Element meshing package
from esys.finley import MakeDomain #Converter for escript
from esys.escript import mkDir, getMPISizeWorld
from esys.escript.unitsSI import *
import os
from math import *
import pylab as pl
import numpy as np
########################################################MPI WORLD CHECK
if getMPISizeWorld() > 1:
	import sys
	print "This example will not run in an MPI world."
	sys.exit(0)

# make sure path exists 
save_path= os.path.join("data","example09") 
mkDir(save_path)

################################################ESTABLISHING PARAMETERS
#Model Parameters
xwidth=2000.0*m   #x width of model
ywidth=2000.0*m   #y width of model
depth=500.0*m   #depth of model

####################################################DOMAIN CONSTRUCTION
# Domain Corners
p0=Point(0.0,    0.0,      0.0)
p1=Point(xwidth, 0.0,      0.0)
p2=Point(xwidth, ywidth,   0.0)
p3=Point(0.0,    ywidth,   0.0)
p4=Point(0.0,    ywidth, depth)
p5=Point(0.0,    0.0,    depth)
p6=Point(xwidth, 0.0,    depth)
p7=Point(xwidth, ywidth, depth)
# Join corners in anti-clockwise manner.
l01=Line(p0, p1)
l12=Line(p1, p2)
l23=Line(p2, p3)
l30=Line(p3, p0)
l56=Line(p5, p6)
l67=Line(p6, p7)
l74=Line(p7, p4)
l45=Line(p4, p5)
l05=Line(p0, p5)
l16=Line(p1, p6)
l27=Line(p2, p7)
l34=Line(p3, p4)

# Join line segments to create domain boundaries and then surfaces
ctop=CurveLoop(l01, l12, l23, l30);     stop=PlaneSurface(ctop)
cbot=CurveLoop(-l67, -l56, -l45, -l74); sbot=PlaneSurface(cbot)
cx0 =CurveLoop(l05, l56, -l16, -l01);   sx0=PlaneSurface(cx0)
cy0 =CurveLoop(-l30, l34, l45, -l05);   sy0=PlaneSurface(cy0)
cxy =CurveLoop(-l34, -l23, l27, l74);   sxy=PlaneSurface(cxy)
cyx =CurveLoop(-l12, l16, l67, -l27);   syx=PlaneSurface(cyx)

# Create the volume.
sloop=SurfaceLoop(stop,sbot,sx0,sy0,sxy,syx)
model=Volume(sloop)


#create interface
sspl=51 #number of discrete points in spline
xdsp=xwidth/(sspl-1) #dx of spline steps for xwidth
ydsp=ywidth/(sspl-1) #dy of spline steps for ywidth
dep_sp=250.0*m #avg depth of spline
h_sp=25.0*m #heigh of spline
orit=-1.0 #orientation of spline 1.0=>up -1.0=>down

# Generate Material Boundary
#x=[ Point(i*xdsp\
 #   ,-dep_sp+orit*h_sp*cos(pi*i*xdsp/dep_sp+pi))\
  #   for i in range(0,sspl)\
  #]
#function along x
x=[1+orit*cos(2*sspl*pi*i*xdsp/xwidth)\
     for i in range(0,sspl)]
#function along y
y=[1+orit*cos(2*sspl*pi*i*ydsp/ywidth)\
     for i in range(0,sspl)]
#containers for surface data
xs=np.zeros(sspl,'float')
ys=np.zeros(sspl,'float')
zs=np.zeros([sspl,sspl],'float')

for i in range(0,sspl):
    xs[i]=i*xdsp
    for j in range(0,sspl):
        if (i == 0):
            ys[j]=j*ydsp
        zs[i,j]=x[i]*y[j]*h_sp+dep_sp

pl.plot(x); pl.plot(y,'ro')
pl.savefig(os.path.join(save_path,"interface.png"))
pl.clf()
pl.contourf(xs,ys,zs)
pl.savefig(os.path.join(save_path,"interface_cf.png"))

spl_ar=[] #interface spline array
linex0_ar=[] #interface line array along x
linexy_ar=[] #interface line array along x,ymax
sintf_ar=[]  #interface surface array
nsintf_ar=[]
#loop through x
for i in range(0,sspl):
    #create all the points required
    point_ar=[Point(xs[i],ys[j],zs[i,j]) for j in range(0,sspl)]
    #create the spline and store it
    spl_ar.append(Spline(*tuple(point_ar)))
    #create the lines along the x axis and x axis at ymax
    if (i < sspl-1): 
        #print i,xs[i],ys[0],zs[i,0]
        tp0=Point(xs[i],ys[0],zs[i,0])
        tp1=Point(xs[i+1],ys[0],zs[i+1,0])
        tp2=Point(xs[i],ys[-1],zs[i,-1])
        tp3=Point(xs[i+1],ys[-1],zs[i+1,-1])
        tl0=Line(tp0,tp1); tl1=Line(tp2,tp3)
        linex0_ar.append(tl0); linexy_ar.append(tl1)
        
for i in range(0,sspl):
    #create the mini surfaces via curveloops and then ruledsurfaces
    if (i < sspl-1):
        #print 'i',i
        tc0=CurveLoop(linex0_ar[i],spl_ar[i+1],-linexy_ar[i],-spl_ar[i])
        ntc0=CurveLoop(spl_ar[i],linexy_ar[i],-spl_ar[i+1],-linex0_ar[i])
        ts0=RuledSurface(tc0); sintf_ar.append(ts0)
        nts0=RuledSurface(ntc0); nsintf_ar.append(nts0)
#create the interface using the surface loop constructor        
sintf=SurfaceLoop(*tuple(sintf_ar))

# for each side
ip0=Point(xs[0],ys[0],zs[0,0])
ip1=Point(xs[-1],ys[0],zs[-1,0])
ip2=Point(xs[-1],ys[-1],zs[-1,-1])
ip3=Point(xs[0],ys[-1],zs[0,-1])

linte_ar=[]; #lines for vertical edges
linhe_ar=[]; #lines for horizontal edges
linte_ar.append(Line(p0,ip0))
linte_ar.append(Line(ip0,p5))
linte_ar.append(Line(p1,ip1))
linte_ar.append(Line(ip1,p6))
linte_ar.append(Line(p2,ip2))
linte_ar.append(Line(ip2,p7))
linte_ar.append(Line(p3,ip3))
linte_ar.append(Line(ip3,p4))

linhe_ar.append(Line(ip0,ip1))
linhe_ar.append(Line(ip1,ip2))
linhe_ar.append(Line(ip2,ip3))
linhe_ar.append(Line(ip3,ip0))

cintfa_ar=[]; cintfb_ar=[] #curveloops for above and below interface on sides
cintfa_ar.append(CurveLoop(linte_ar[0],linhe_ar[0],-linte_ar[2],-l01))
cintfa_ar.append(CurveLoop(linte_ar[2],linhe_ar[1],-linte_ar[4],-l12))
cintfa_ar.append(CurveLoop(linte_ar[4],linhe_ar[2],-linte_ar[6],-l23))
cintfa_ar.append(CurveLoop(linte_ar[6],linhe_ar[3],-linte_ar[0],-l30))

cintfb_ar.append(CurveLoop(linte_ar[1],l56,-linte_ar[3],-linhe_ar[0]))
cintfb_ar.append(CurveLoop(linte_ar[3],l67,-linte_ar[5],-linhe_ar[1]))
cintfb_ar.append(CurveLoop(linte_ar[5],l74,-linte_ar[7],-linhe_ar[2]))
cintfb_ar.append(CurveLoop(linte_ar[7],l45,-linte_ar[1],-linhe_ar[3]))

sintfa_ar=[PlaneSurface(cintfa_ar[i]) for i in range(0,4)]
sintfb_ar=[PlaneSurface(cintfb_ar[i]) for i in range(0,4)]

vintfa=Volume(SurfaceLoop(stop,*tuple(sintfa_ar+nsintf_ar)))
vintfb=Volume(SurfaceLoop(sbot,*tuple(sintfb_ar+sintf_ar)))


#############################################EXPORTING MESH FOR ESCRIPT
# Create a Design which can make the mesh
d=Design(dim=3, element_size=200.0*m)
# Add the subdomains and flux boundaries.
d.addItems(model,stop,sbot,sx0,sy0,sxy,syx)
#d.addItems(PropertySet('intf',sintf))

d.addItems(PropertySet('vintfa',vintfa),PropertySet('vintfb',vintfb))


d.setScriptFileName(os.path.join(save_path,"example09m.geo"))

#d.setMeshFileName(os.path.join(save_path,"example09m.msh"))
#
#  make the finley domain:
#
domain=MakeDomain(d)
# Create a file that can be read back in to python with
# mesh=ReadMesh(fileName)
#domain.write(os.path.join(save_path,"example09m.fly"))


