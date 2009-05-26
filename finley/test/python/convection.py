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
"""
this is a convection simulation over a domain [0,L] X [0,L] x [0,H]

"""
__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript.models import TemperatureCartesian, StokesProblemCartesian
from esys.finley import Rectangle, Brick, LoadMesh
from optparse import OptionParser
from math import pi, ceil
import sys
import time

# ======================= Default Values ==================================================
DIM=2                           # spatial dimension
H=1.                            # height
L=4*H                           # length
NE=10                           # number of elements in H-direction. 
PERTURBATION=0.1                # initial temperature perturbation
DT=1.e-4                        # initial time step size
USE_BACKWARD_EULER=False        # use backward Euler in time integartion of temperature otherwise Crank-Nicholson is used
CREATE_TOPGRAPHY=False          # create topgraphy
TOL=1.e-4                       # tolerance 
DT_MIN=1.e-5                    # minumum time step size
T_END=0.1                       # end time

DT_OUT=T_END/500
Dn_OUT=2
VERBOSE=True
CREATE_RESTARTFILES_EVERY_STEP=10
if True:
   # this is a simple linear Stokes model:
   RA=1.e6 # Rayleigh number
   A=22 # Arenious number 
   DI = 0.  # dissipation number 
   MU=None
   ETA0=1.
   TAU0=None
   N=None
   NPL=None
   ETAP0=ETA0
   TAUY=None
   useJAUMANNSTRESS=False
   # this is a simple linear Stokes model:
   RA=1.e5 # Rayleigh number
   A=11 # Arenious number 
   DI = 0.  # dissipation number 
   MU=None
   ETA0=1.
   TAU0=250.
   N=None
   NPL=None
   ETAP0=ETA0
   TAUY=TAU0
   useJAUMANNSTRESS=False
else:
   RA=1.e4 # Rayleigh number
   A=22 # Arenious number 
   DI = 0.  # dissipation number 
   MU=1.e4
   ETA0=1.
   TAU0=0.866*10**2.5 
   N=3
   NPL=14
   TAUY=TAU0
   ETAP0=ETA0
   useJAUMANNSTRESS=True


# =========================================================================================

def removeRestartDirectory(dir_name):
   if dom.onMasterProcessor() and os.path.isdir(dir_name):
       for root, dirs, files in os.walk(dir_name, topdown=False):
           for name in files: os.remove(os.path.join(root,name))
           for name in dirs: os.remove(os.path.join(root,name))
       os.rmdir(dir_name)
       print "Restart files %s have been removed."%dir_name
   dom.MPIBarrier()
# =========================================================================================
#
# read options:
#
parser = OptionParser(usage="%prog [Options]")
parser.add_option("-r", "--restart", dest="restart", help="restart from latest directory. It will be deleted after a new directory has been created.", default=False, action="store_true")
parser.add_option("-d", "--dir", dest="restart_dir", help="restart from directory DIR. The directory will not be deleted but new restart directories are created.",metavar="DIR", default=None)
parser.add_option("-p", "--param", dest="param", help="name of file to be imported ",metavar="PARAM", default=None)
(options, args) = parser.parse_args()
restart=options.restart or (options.restart_dir !=None)
#
#  overwrite the default options:
#
print "Execution started ",time.asctime()
if options.param !=None: 
     exec(open(options.param,'r'))
     print "Parameters are imported from file ",options.param
print "Parameters:"
print "\tDimension\tDIM=\t",DIM
print "\tHeight:\tH=\t",H
print "\tLength:\tL=\t",L
print "\telements in H:\tNE=\t",NE
print "\ttotal #element\t\t=\t",NE**DIM*int(L/H)**(DIM-1)
print "\ttolerance\tTOL=\t",TOL
print "\ttemperature perturbation\tPERTURBATION=\t",PERTURBATION
print "\tinitial time step size\tDT=\t",DT
print "\tminimum time step size\tDT_MIN=\t",tDT_MIN
print "\tend time\tT_END=\t",T_END
print "\tbackward Euler?\tUSE_BACKWARD_EULER=\t",USE_BACKWARD_EULER
print "\ttopgraphy?\tCREATE_TOPOGRAPHY=\t",CREATE_TOPOGRAPHY
#  some control variables (will be overwritten in  case of a restart:
#
t=0         # time stamp
n=0         # time step counter
dt=DT       # current time step size

t_out=0     # 
n_out=0
out_count=0
a=None
dt_a=None
#
#   set up domain or read restart file
#
if restart:
   if options.restart_dir ==None:
      restart_files=[]
      for f in os.listdir("."):
          if f.startswith("restart"): restart_files.append(f)
      if len(restart_files)==0:
          raise IOError,"no restart files"
      restart_files.sort()
      f=restart_files[-1]
   else:
       f=options.restart_dir
   print ">>>Restart from directory ",f
   try:
      dom=LoadMesh("mesh.nc")
   except:
      pass
   FF=open(os.path.join(f,"stamp.%d"%dom.getMPIRank()),"r").read().split(";")
   t=float(FF[0])          # time stamp
   t_out=float(FF[1])      # 
   n_out=int(FF[2])
   n=int(FF[3])
   out_count=int(FF[4])
   dt=float(FF[5])
   stress=load(os.path.join(f,"stress.nc"),dom)
   v=load(os.path.join(f,"v.nc"),dom)
   p=load(os.path.join(f,"p.nc"),dom)
   T=load(os.path.join(f,"T.nc"),dom)
   if n>1:
      dt_a=float(FF[6])
      a=load(os.path.join(f,"a.nc"),dom)
   else:
      dt_a=None
      a=None
   if dom.onMasterProcessor(): nusselt_file=open("nusselt.csv","a")
else:
  if DIM==2:
    dom=Rectangle(int(ceil(L*NE/H)),NE,l0=L,l1=H,order=2, useFullElementOrder=True,optimize=True)
  else:
    dom=Brick(int(ceil(L*NE/H)),int(ceil(L*NE/H)),NE,l0=L,l1=L,l2=H,order=2, useFullElementOrder=True,optimize=True)
  try:
     dom.dump("mesh.nc")
  except:
     pass
  x=dom.getX()
  T=Scalar(1,ReducedSolution(dom))
  for d in range(DIM):
      if d == DIM-1: 
         T*=sin(x[d]/H*pi)
      else:
         T*=cos(x[d]/L*pi)

  T=1.-x[DIM-1]+PERTURBATION*T
  v=Vector(0,Solution(dom))
  stress=Tensor(0,Function(dom))
  x2=ReducedSolution(dom).getX()
  p=-RA*(x2[DIM-1]-0.5*x2[DIM-1]**2)
  if dom.onMasterProcessor(): nusselt_file=open("nusselt.csv","w")

vol=integrate(1.,Function(dom))
p-=integrate(p)/vol
x=dom.getX()
#
#   set up heat problem:
#
heat=TemperatureCartesian(dom,useBackwardEuler=USE_BACKWARD_EULER)
heat.setTolerance(TOL)
fixed_T_at=whereZero(x[DIM-1])+whereZero(H-x[DIM-1])
print "initial Temperature range ",inf(T),sup(T)
heat.setInitialTemperature(clip(T,0))
heat.setValue(rhocp=Scalar(1.,Function(dom)),k=Scalar(1.,Function(dom)),given_T_mask=fixed_T_at)
#
#   velocity constraints:
#
x2=ReducedSolution(dom).getX()
fixed_v_mask=Vector(0,Solution(dom))
for d in range(DIM):
    if d == DIM-1: 
       ll = H
    else:
       ll = L
    if CREATE_TOPOGRAPHY and d==DIM-1:
       fixed_v_mask+=whereZero(x[d])*unitVector(d,DIM)
    else:
       fixed_v_mask+=(whereZero(x[d])+whereZero(x[d]-ll))*unitVector(d,DIM)
1/0
#
#   set up velovity problem
# ????????????????????
sp=StokesProblemCartesian(dom)
# ,stress,v,p,t,useJaumannStress=useJAUMANNSTRESS)
# sp=PlateMantelModel(dom,stress,v,p,t,useJaumannStress=useJAUMANNSTRESS)
# sp.initialize(mu=MU, tau_0=TAU0, n=N, eta_Y=ETAP0, tau_Y=TAU0, n_Y=NPL, q=fixed_v_mask)
sp.setTolerance(TOL*extratol)
# sp.setToleranceReductionFactor(TOL)
#
#  let the show begin:
#
t1 = time.time()
while t<T_END:
    v_last=v*1
    print "============== solve for v ========================"
    FF=exp(A*(1./(1+T.interpolate(Function(dom)))-1./2.))
    print "viscosity range :", inf(FF)*ETA0, sup(FF)*ETA0
    # sp.initialize(eta_N=ETA0*FF, eta_0=ETA0*FF, F=T*(RA*unitVector(DIM-1,DIM)))
    sp.initialize(eta=ETA0*FF, f=T*(RA*unitVector(DIM-1,DIM)),fixed_u_mask=fixed_v_mask)
    v,p=sp.solve(v,p,max_iter=20, verbose=VERBOSE, show_details=False)

    for d in range(DIM):
         print "range %d-velocity"%d,inf(v[d]),sup(v[d])
    if t>=t_out or n>n_out:
      saveVTK("state.%d.vtu"%out_count,T=T,v=v)
      print "visualization file %d for time step %e generated."%(out_count,t)
      out_count+=1
      t_out+=DT_OUT
      n_out+=Dn_OUT
    # calculation of nusselt number:
    se=ETA0*FF*length(symmetric(grad(v)))**2 # this propably wroing!
    print "Xse:",inf(se),sup(se)
    Nu=1.+integrate(se)/(RA*vol)
    if dom.onMasterProcessor(): nusselt_file.write("%e %e\n"%(t,Nu))
    heat.setValue(v=interpolate(v,ReducedSolution(dom)),Q=DI/RA*se)
    print "Xnusselt number = ",Nu, "dt =",dt
    if n>0:
        a,a_alt = (v_last-v)/dt, a
        dt_a,dt_a_alt = dt, dt_a
    if n>1:
       z=(a-a_alt)/((dt_a+dt_a_alt)/2)
       f=Lsup(z)/Lsup(v)
       print "estimated error ",f*dt**2
       # dt_new=min(2*dt,max(dt/2,sqrt(0.05/f)))
       dt_new=sqrt(0.05/f)
       dt=min(dt_new,heat.getSafeTimeStepSize())
    else:
       dt=heat.getSafeTimeStepSize()
    dt=max(DT_MIN,dt)
    print n,". time step t=",t," step size ",dt
    print "============== solve for T ========================"
    T=heat.getSolution(dt, verbose=VERBOSE)	
    print "Temperature range ",inf(T),sup(T)
    n+=1
    t+=dt
    # =========================
    #
    #    create restart files:
    #
    if create_restartfiles_every_step>0:
       if (n-1)%create_restartfiles_every_step == 0:
         c=(n-1)/create_restartfiles_every_step
         old_restart_dir="restart_%s_"%(c-1)
         new_restart_dir="restart_%s_"%c

         print "Write restart files to ",new_restart_dir
         if dom.onMasterProcessor() and not os.path.isdir(new_restart_dir): os.mkdir(new_restart_dir)
	 dom.MPIBarrier()
         # sp.getStress().dump(os.path.join(new_restart_dir,"stress.nc"))
         v.dump(os.path.join(new_restart_dir,"v.nc"))
         p.dump(os.path.join(new_restart_dir,"p.nc"))
         T.dump(os.path.join(new_restart_dir,"T.nc"))
         if n>1:
             file(os.path.join(new_restart_dir,"stamp.%d"%dom.getMPIRank()),"w").write("%e; %e; %s; %s; %s; %e; %e;\n"%(t, t_out, n_out, n, out_count, dt, dt_a))
             a.dump(os.path.join(new_restart_dir,"a.nc"))
         else:
             file(os.path.join(new_restart_dir,"stamp.%d"%dom.getMPIRank()),"w").write("%e; %e; %s; %s; %s; %e;\n"%(t, t_out, n_out, n, out_count, dt))
         removeRestartDirectory(old_restart_dir)
elapsed = time.time() - t1
print "plot","\t",NE,"\t",elapsed

