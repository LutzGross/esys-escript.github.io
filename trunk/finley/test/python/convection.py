from esys.escript import *
from esys.escript.models import TemperatureCartesian, StokesProblemCartesian
from esys.finley import Rectangle, Brick
from optparse import OptionParser
from math import pi, ceil

NE=20
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
print "total number of elements = ",NE**DIM*int(L/H)**(DIM-1)

# read options:
parser = OptionParser(usage="%prog [-r [DIR]]")
parser.add_option("-r", "--restart", dest="restart", help="restart from latest directory", default=False, action="store_true")
parser.add_option("-d", "--dir", dest="restart_dir", help="restart from directory DIR",metavar="DIR", default=None)
(options, args) = parser.parse_args()
restart=options.restart or (options.restart_dir !=None)
print restart
#
#   set up domain:
#
if restart:
   if options.restart_dir ==None:
      restart_files=[]
      for f in os.listdir("."):
          if f.startswith("restart"): restart_files.append(f)
      if len(restart_files)==0:
          raise RunTimeError,"no restart files"
      restart_files.sort()
      f=restart_files[-1]
   else:
       f=options.restart_dir
   print "restart from ",f
   dom=finley.load("mesh.nc")
   velocity=load(os.path.join(f,"v.nc"),dom)
   pressure=load(os.path.join(f,"p.nc"),dom)
   ff=open(os.path.join(f,"s.%d"%dom.getMPIRank()),"r").read().split(";")
   t=float(ff[0])
   t_out=float(ff[1])
   n=int(ff[2])
   n_out=int(ff[3])
   dt=float(ff[4])
   if dom.getMPIRank()==0: nusselt_file=open("nusselt.csv","a")
else:
  if DIM==2:
    dom=Rectangle(int(ceil(L*NE/H)),NE,l0=L,l1=H,order=2, useFullElementOrder=True,optimize=True)
  else:
    dom=Brick(int(ceil(L*NE/H)),int(ceil(L*NE/H)),NE,l0=L,l1=L,l2=H,order=2, useFullElementOrder=True,optimize=True)
  dom.dump("mesh.nc")
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
sp.setTolerance(TOL)
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
dt_new=1.
a=None
a_alt=None
if dom.getMPIRank() ==0: 
     nusselt_file=open("nusselt.csv","w")
while t<T_END:
    v_last=v*1
    print "============== solve for v ========================"
    viscosity=exp(A*(1./(1+T.interpolate(Function(dom)))-2./3.))
    print "viscosity range :", inf(viscosity), sup(viscosity)
    sp.initialize(f=T*(RA*unitVector(DIM-1,DIM)),eta=viscosity,fixed_u_mask=fixed_v_mask)
    v,p=sp.solve(v,p,show_details=VERBOSE, verbose=True,max_iter=500,solver='PCG')
    # v,p=sp.solve(v,p,show_details=VERBOSE, verbose=True,max_iter=500,solver='GMRES')

    for d in range(DIM):
         print "range %d-velocity"%d,inf(v[d]),sup(v[d])

    if t>=t_out:
      saveVTK("state.%d.vtu"%n_out,T=T,v=v,p=p)
      print "visualization file %d for time step %e generated."%(n_out,t)
      n_out+=1
      t_out+=DT_OUT
    Nu=1.+integrate(viscosity*length(grad(v))**2)/(RA*vol)
    if dom.getMPIRank() ==0: nusselt_file.write("%e %e\n"%(t,Nu))
    print "nusselt number = ",Nu
    if a==None:
       a=0
    else:
       a,a_alt = (v_last-v)/dt, a
       dt_alt=dt
       if a_alt!=None:
          z=(a-a_alt)/((dt+dt_alt)/2)
          f=Lsup(z)/Lsup(v)
          print "estimated error ",f*dt**2
          dt_new, dt_alt =min(2*dt,max(dt/2,sqrt(0.05/f))), dt
          # dt_new, dt_alt =sqrt(0.05/f), dt
    heat.setValue(v=v,Q=DI/RA*viscosity*length(symmetric(grad(v)))**2)
    dt=min(dt_new,heat.getSafeTimeStepSize())
    print n,". time step t=",t," step size ",dt
    print "============== solve for T ========================"
    T=heat.solve(dt, verbose=VERBOSE)	
    print "Temperature range ",inf(T),sup(T)
    n+=1
    t+=dt
