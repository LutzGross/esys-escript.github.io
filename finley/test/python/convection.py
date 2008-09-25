
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
from esys.escript.models import TemperatureCartesian, PlateMantelModel
from esys.finley import Rectangle, Brick, LoadMesh
from optparse import OptionParser
from math import pi, ceil


def removeRestartDirectory(dir_name):
   if os.path.isdir(dir_name):
       for root, dirs, files in os.walk(dir_name, topdown=False):
           for name in files: os.remove(os.path.join(root,name))
           for name in dirs: os.remove(os.path.join(root,name))
       os.rmdir(dir_name)
       print "Restart files %s have been removed."%dir_name


import sys
import time
t1 = time.time()

extratol=1

# read options:
parser = OptionParser(usage="%prog [-r [DIR]] [-e [NE]] [-s [solver]]")
parser.add_option("-s", "--solver", dest="solver", help="solver to be used for saddle point problem. The possible options are PCG, GMRES, NewtonGMRES, MINRES and TFQMR.", metavar="solver", default="PCG")
parser.add_option("-e", "--elements", dest="NE", help="number of elements in one direction.",metavar="NE", default=16)

parser.add_option("-r", "--restart", dest="restart", help="restart from latest directory. It will be deleted after a new directory has been created.", default=False, action="store_true")
parser.add_option("-d", "--dir", dest="restart_dir", help="restart from directory DIR. The directory will not be deleted but new restart directories are created.",metavar="DIR", default=None)
(options, args) = parser.parse_args()
restart=options.restart or (options.restart_dir !=None)

solver=options.solver
NE=int(options.NE)

DIM=2
H=1.
L=4*H
THETA=0.5 # time stepping THETA=0.5 cranck nicolson
TOL=1.e-4
PERTURBATION=0.1
DT=1.e-4
DT_MIN=1.e-5
T_END=0.1
DT_OUT=T_END/500
Dn_OUT=2
VERBOSE=True
SUPG=False
create_restartfiles_every_step=10
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
   A=0 # Arenious number 
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

print "total number of elements = ",NE**DIM*int(L/H)**(DIM-1)

#
#   set up domain:
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
   print "restart from directory ",f
   try:
      dom=LoadMesh("mesh.nc")
   except:
      pass
   FF=open(os.path.join(f,"stamp.%d"%dom.getMPIRank()),"r").read().split(";")
   t=float(FF[0])
   t_out=float(FF[1])
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
   if dom.getMPIRank()==0: nusselt_file=open("nusselt.csv","a")
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
  if dom.getMPIRank() ==0: nusselt_file=open("nusselt.csv","w")
  t=0
  t_out=0
  n_out=0
  n=0
  out_count=0
  dt=DT
  a=None
  dt_a=None

vol=integrate(1.,Function(dom))
p-=integrate(p)/vol
x=dom.getX()
#
#   set up heat problem:
#
heat=TemperatureCartesian(dom,theta=THETA,useSUPG=SUPG)
heat.setTolerance(TOL*extratol)

fixed_T_at=whereZero(x[DIM-1])+whereZero(H-x[DIM-1])
heat.setInitialTemperature(T)
print "initial Temperature range ",inf(T),sup(T)
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
    fixed_v_mask+=(whereZero(x[d])+whereZero(x[d]-ll))*unitVector(d,DIM)
#
#   set up velovity problem
#
sp=PlateMantelModel(dom,stress,v,p,t,useJaumannStress=useJAUMANNSTRESS)
sp.initialize(mu=MU, tau_0=TAU0, n=N, eta_Y=ETAP0, tau_Y=TAU0, n_Y=NPL, q=fixed_v_mask)
sp.setTolerance(TOL*extratol)
sp.setToleranceReductionFactor(TOL)
#
#  let the show begin:
#
while t<T_END:
    v_last=v*1
    print "============== solve for v ========================"
    FF=exp(A*(1./(1+T.interpolate(Function(dom)))-1./2.))
    print "viscosity range :", inf(FF)*ETA0, sup(FF)*ETA0
    sp.initialize(eta_N=ETA0*FF, eta_0=ETA0*FF, F=T*(RA*unitVector(DIM-1,DIM)))
    sp.update(dt,max_inner_iter=20, verbose=VERBOSE, show_details=False, tol=10., solver=solver)
    v=sp.getVelocity()

    for d in range(DIM):
         print "range %d-velocity"%d,inf(v[d]),sup(v[d])
    if t>=t_out or n>n_out:
      saveVTK("state.%d.vtu"%out_count,T=T,v=v)
      print "visualization file %d for time step %e generated."%(out_count,t)
      out_count+=1
      t_out+=DT_OUT
      n_out+=Dn_OUT
    # calculation of nusselt number:
    se=sp.getMechanicalPower()
    print "Xse:",inf(se),sup(se)
    Nu=1.+integrate(se)/(RA*vol)
    if dom.getMPIRank() ==0: nusselt_file.write("%e %e\n"%(t,Nu))
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
    T=heat.solve(dt, verbose=VERBOSE)	
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
         if not os.path.isdir(new_restart_dir): os.mkdir(new_restart_dir)
         sp.getStress().dump(os.path.join(new_restart_dir,"stress.nc"))
         sp.getVelocity().dump(os.path.join(new_restart_dir,"v.nc"))
         sp.getPressure().dump(os.path.join(new_restart_dir,"p.nc"))
         T.dump(os.path.join(new_restart_dir,"T.nc"))
         if n>1:
             file(os.path.join(new_restart_dir,"stamp.%d"%dom.getMPIRank()),"w").write("%e; %e; %s; %s; %s; %e; %e;\n"%(t, t_out, n_out, n, out_count, dt, dt_a))
             a.dump(os.path.join(new_restart_dir,"a.nc"))
         else:
             file(os.path.join(new_restart_dir,"stamp.%d"%dom.getMPIRank()),"w").write("%e; %e; %s; %s; %s; %e;\n"%(t, t_out, n_out, n, out_count, dt))
         removeRestartDirectory(old_restart_dir)
elapsed = time.time() - t1
print "plot","\t",NE,"\t",elapsed

