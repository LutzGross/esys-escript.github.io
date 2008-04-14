from esys.escript import *
from esys.escript.models import TemperatureCartesian, StokesProblemCartesian
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

if (len(sys.argv)>=3):
 NE=int(sys.argv[2])
else:
 NE=20

if (len(sys.argv)>=2):
 solver=sys.argv[1]
else:
 solver='PCG'

if solver!='PCG':
 extratol=0.001
else:
 extratol=1

DIM=2
H=1.
L=4*H
THETA=0.5
TOL=1.e-3
PERTURBATION=0.1
T_END=0.3
DT_OUT=T_END/500
VERBOSE=False
RA=1.e5 # Rayleigh number
A=22.  # Arenious number 
DI = 0.  # dissipation number
SUPG=False
create_restartfiles_every_step=10

print "total number of elements = ",NE**DIM*int(L/H)**(DIM-1)

# read options:
parser = OptionParser(usage="%prog [-r [DIR]]")
parser.add_option("-r", "--restart", dest="restart", help="restart from latest directory. It will be deleted after a new directory has been created.", default=False, action="store_true")
parser.add_option("-d", "--dir", dest="restart_dir", help="restart from directory DIR. The directory will not be deleted but new restart directories are created.",metavar="DIR", default=None)
(options, args) = parser.parse_args()
restart=options.restart or (options.restart_dir !=None)
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
   ff=open(os.path.join(f,"stamp.%d"%dom.getMPIRank()),"r").read().split(";")
   t=float(ff[0])
   t_out=float(ff[1])
   n=int(ff[2])
   n_out=int(ff[3])
   dt=float(ff[4])
   v=load(os.path.join(f,"v.nc"),dom)
   p=load(os.path.join(f,"p.nc"),dom)
   T=load(os.path.join(f,"T.nc"),dom)
   if n>1:
      dt_a=float(ff[5])
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
         T*=cos(x[d]*((d+1)/L*pi))

  T=1.-x[DIM-1]+PERTURBATION*T
  v=Vector(0,Solution(dom))
  if dom.getMPIRank() ==0: nusselt_file=open("nusselt.csv","w")
  t=0
  t_out=0
  n=0
  n_out=0
  dt=None
  a=None
  dt_a=None

vol=integrate(1.,Function(dom))
x=dom.getX()
#
#   set up heat problem:
#
heat=TemperatureCartesian(dom,theta=THETA,useSUPG=SUPG)
heat.setTolerance(TOL)

fixed_T_at=whereZero(x[DIM-1])+whereZero(H-x[DIM-1])
heat.setInitialTemperature(T)
print "initial Temperature range ",inf(T),sup(T)
heat.setValue(rhocp=Scalar(1.,Function(dom)),k=Scalar(1.,Function(dom)),given_T_mask=fixed_T_at)
#
#   set up velovity problem
#
sp=StokesProblemCartesian(dom)
sp.setTolerance(TOL*extratol)
sp.setToleranceReductionFactor(TOL)
x2=ReducedSolution(dom).getX()
p=-RA*(x2[DIM-1]-0.5*x2[DIM-1]**2)

fixed_v_mask=Vector(0,Solution(dom))
for d in range(DIM):
    if d == DIM-1: 
       ll = H
    else:
       ll = L
    fixed_v_mask+=(whereZero(x[d])+whereZero(x[d]-ll))*unitVector(d,DIM)

#
#  let the show begin:
#
while t<T_END:
    v_last=v*1
    print "============== solve for v ========================"
    viscosity=exp(A*(1./(1+T.interpolate(Function(dom)))-2./3.))
    print "viscosity range :", inf(viscosity), sup(viscosity)
    sp.initialize(f=T*(RA*unitVector(DIM-1,DIM)),eta=viscosity,fixed_u_mask=fixed_v_mask)
    #v,p=sp.solve(v,p,show_details=VERBOSE, verbose=True,max_iter=500,solver='PCG')
    #v,p=sp.solve(v,p,show_details=VERBOSE, verbose=True,max_iter=500,solver='GMRES')
    #v,p=sp.solve(v,p,show_details=VERBOSE, verbose=True,max_iter=500,solver='MINRES')
    v,p=sp.solve(v,p,show_details=VERBOSE, verbose=True,max_iter=500,solver=solver)

    for d in range(DIM):
         print "range %d-velocity"%d,inf(v[d]),sup(v[d])

    if t>=t_out:
      saveVTK("state.%d.vtu"%n_out,T=T,v=v,p=p)
      print "visualization file %d for time step %e generated."%(n_out,t)
      n_out+=1
      t_out+=DT_OUT
    Nu=1.+integrate(viscosity*length(grad(v))**2)/(RA*vol)
    if dom.getMPIRank() ==0: nusselt_file.write("%e %e\n"%(t,Nu))
    heat.setValue(v=v,Q=DI/RA*viscosity*length(symmetric(grad(v)))**2)
    print "nusselt number = ",Nu,n
    if n>0:
        a,a_alt = (v_last-v)/dt, a
        dt_a,dt_a_alt = dt, dt_a
    if n>1:
       z=(a-a_alt)/((dt_a+dt_a_alt)/2)
       f=Lsup(z)/Lsup(v)
       print "estimated error ",f*dt**2
       dt_new=min(2*dt,max(dt/2,sqrt(0.05/f)))
       dt=min(dt_new,heat.getSafeTimeStepSize())
    else:
       dt=heat.getSafeTimeStepSize()
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
         v.dump(os.path.join(new_restart_dir,"v.nc"))
         p.dump(os.path.join(new_restart_dir,"p.nc"))
         T.dump(os.path.join(new_restart_dir,"T.nc"))
         if n>1:
             file(os.path.join(new_restart_dir,"stamp.%d"%dom.getMPIRank()),"w").write("%e; %e; %s; %s; %e; %e;\n"%(t, t_out, n, n_out, dt, dt_a))
             a.dump(os.path.join(new_restart_dir,"a.nc"))
         else:
             file(os.path.join(new_restart_dir,"stamp.%d"%dom.getMPIRank()),"w").write("%e; %e; %s; %s; %e;\n"%(t, t_out, n, n_out, dt))
         removeRestartDirectory(old_restart_dir)
elapsed = time.time() - t1
print "plot","\t",NE,"\t",elapsed

