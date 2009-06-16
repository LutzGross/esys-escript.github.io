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

It is solved in dimensionless form

"""
__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript.models import TemperatureCartesian, IncompressibleIsotropicFlowCartesian
from esys.finley import Rectangle, Brick, LoadMesh
from optparse import OptionParser
from math import pi, ceil
import sys
import time

# ======================= Default Values ==================================================
DIM=2                           # spatial dimension
H=1.                            # height
L=H                           # length
NE=50                           # number of elements in H-direction. 
PERT=0.05               # initial temperature perturbation
DT=1.e-7                        # initial time step size
CREATE_TOPOGRAPHY=False          # create topgraphy
TOPOGRAPHY_LOAD=False        # considers the loading of the topography in the flow
DT_MIN=1.e-10                    # minumum time step size
T_END=0.1                       # end time

RHO_0=1.                        # surface density
G=1.                            # gravitational constant
ALPHA_0=4.e3                      # thermal expansion coefficient
T_0=0.                          # surface temperature
T_1=1.                          # bottom temperature
C_P=1                           # heat capacity 
K=1.                            # thermal conductivity
CHI=0.                          # Taylor Quinny coefficient
MUE=None                        # elastic shear modulus
TAU_Y=5*10**(2.5)               # Drucker-Prager cohesion factor
BETA=0                          # Drucker-Prager friction factor
TAU_0=2*10**(2.5)               # transition stress
N=3                             # power for power law

E=23                            # activation energy
V=18*0                        # activation volume 
T_OFFSET=1                      # temperature offset on surface (dimensionless formulation T_OFFSET=1 otherwise =0)
R=1                             # gas constant
ETA_N0=1.                       # viscosity at surface 

T_TOL=1.e-4                     # tolerance temperature transport
TOL2=0.1                        # tolerance for flow update at a timestep. (large value will do only one correction step)
FLOW_TOL=1.e-4                  # tolerance for inconcompressible flow solver
FLOW_SUB_TOL=1.e-8             # sub-tolerance for inconcompressible flow solver
NUSSELT_FN="nusselt.csv"
VERBOSE=True
DT_VIS=T_END/500                # time distane between two visulaization files
DN_VIS=5                        # maximum counter increment between two visulaization files
VIS_DIR="results"               # name of the director for vis files
DN_RESTART=1000000              # create a restart file every DN_RESTART step
PREFIX_RESTART="restart"        # name prefix for restart directories
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
print "<%s> Execution started."%time.asctime()
if options.param !=None: 
     exec(open(options.param,'r'))
     print "Parameters are imported from file ",options.param

CREATE_TOPOGRAPHY=CREATE_TOPOGRAPHY or TOPOGRAPHY_LOAD

print "Input Parameters:"
print "\tDimension                     DIM\t\t=\t",DIM
print "\tHeight                        H\t\t\t=\t",H
print "\tLength                        L\t\t\t=\t",L
print "\telements in H                 NE\t\t=\t",NE
print "\ttemperature perturbation      PERT\t\t=\t",PERT
print "\tinitial time step size        DT\t\t=\t",DT
print "\tminimum time step size        DT_MIN\t\t=\t",DT_MIN
print "\tend time                      T_END\t\t=\t",T_END
print "\tcreate topgraphy              CREATE_TOPOGRAPHY\t=\t",CREATE_TOPOGRAPHY
print "\tloading by topography         TOPOGRAPHY_LOAD\t=\t",TOPOGRAPHY_LOAD
print "\tsurface density               RHO_0\t\t=\t",RHO_0
print "\tgravitational constant        G\t\t\t=\t",G
print "\tthermal expansion coefficient ALPHA_0\t\t=\t",ALPHA_0
print "\tsurface temperature           T_0\t\t=\t",T_0
print "\tbottom temperature            T_1\t\t=\t",T_1
print "\theat capacity                 C_P\t\t=\t",C_P
print "\tthermal conductivity          K\t\t\t=\t",K
print "\tTaylor-Qinny coefficient      CHI\t\t=\t",CHI
print "\telastic shear modulus         MUE\t\t=\t",MUE
print "\tcohesion factor               TAU_Y\t\t=\t",TAU_Y
print "\tfriction factor               BETA\t\t=\t",BETA
print "\ttransition stress             TAU_0\t\t=\t",TAU_0
print "\tpower for power law           N\t\t\t=\t",N
print "\tviscosity at surface          ETA_N0\t\t=\t",ETA_N0
print "\tactivation energy             E\t\t\t=\t",E
print "\tactivation volume             V\t\t\t=\t",V
print "\ttemperature offset            T_OFFSET\t\t=\t",T_OFFSET
print "\tgas constant                  R\t\t\t=\t",R

print "\ttransport tolerance           T_TOL\t\t=\t",T_TOL
print "\ttolerance for flow updates    TOL2\t\t=\t",TOL2
print "\tflow tolerance                FLOW_TOL\t\t=\t",FLOW_TOL
print "\tflow sub-tolerance            FLOW_SUB_TOL\t=\t",FLOW_SUB_TOL
print "\tfile for Nusselt number       NUSSELT_FN\t=\t",NUSSELT_FN
print "\tmin. time incr. for vis file  DT_VIS\t\t=\t",DT_VIS
print "\tmin. count incr. for vis file DN_VIS\t\t=\t",DN_VIS
print "\tvisualization dir             VIS_DIR\t\t=\t",VIS_DIR
print "\trestart counter incerement    DN_RESTART\t\t=\t",DN_RESTART
print "\tprefix for restart dirs       PREFIX_RESTART\t\t=\t",PREFIX_RESTART
print "\tverbosity                     VERBOSE\t\t=\t",VERBOSE

print "Control Parameters:"
t_REF=RHO_0*C_P*H**2/K
P_REF=ETA_N0/t_REF
Ar=E/R/(T_1-T_0)
if E>0 and V>0:
   V_REF=P_REF*V/E
else:
  V_REF=0
T_OFFSET_REF=T_OFFSET/(T_1-T_0)
Ra=RHO_0*G*H*(T_1-T_0)*ALPHA_0/P_REF
Di=ALPHA_0*G*H/C_P
CHI_REF=CHI*K*ETA_N0/(RHO_0**2*C_P**2*(T_1-T_0)*H**2)
if MUE == None:
  De=None
else:
  De=ETA_N0/MUE/t_REF
ETA_BOT=exp(Ar*((1.+V_REF)/(T_OFFSET_REF+1)-1./T_OFFSET_REF))*ETA_N0
print "\ttotal #element                \t\t\t=\t",NE**DIM*int(L/H)**(DIM-1)
print "\treference time                t_REF\t\t=\t%e"%t_REF
print "\treference pressure            P_REF\t\t=\t%e"%P_REF
print "\treference Taylor-Qinny        CHI_REF\t\t=\t%e"%CHI_REF
print "\tDissipation number            DI\t\t=\t%e"%Di
print "\tRayleigh number surface       Ra\t\t=\t%e"%Ra
if MUE == None:
   print "\tDebora number surface         De\t\t=\t",None
else:
   print "\tDebora number surface         De\t\t=\t%e"%De
     
print "\tbottom viscosity              \t\t\t=\t%e"%ETA_BOT
print "\tRayleigh number bottom        \t\t\t=\t%e"%(RHO_0*G*H*(T_1-T_0)*ALPHA_0*t_REF/ETA_BOT)
if MUE == None:
   print "\tDebora number bottom          \t\t\t=\t",None
else:
   print "\tDebora number bottom          \t\t\t=\t%d"%(ETA_BOT/MUE/t_REF)
print "\tArrhenius                     Ar\t\t=\t%e"%Ar
print "\tscaled activation volume      V_REF\t\t=\t%e"%V_REF
#  some control variables (will be overwritten in  case of a restart:
#
t=0         # time stamp
n=0         # time step counter
dt=DT       # current time step size
t_vis=0     
n_vis=0
counter_vis=0
if getMPIRankWorld()==0:
    if not os.path.isdir(VIS_DIR): os.mkdir(VIS_DIR)
#=========================
#
#   set up domain or read restart file
#
if restart:
   if options.restart_dir ==None:
      restart_files=[]
      for f in os.listdir("."):
          if f.startswith(PREFIX_RESTART): restart_files.append(f)
      if len(restart_files)==0:
          raise IOError,"no restart files"
      restart_files.sort()
      f=restart_files[-1]
   else:
       f=options.restart_dir
   try:
      dom=LoadMesh("mesh.nc")
   except:
      pass
   FF=open(os.path.join(f,"stamp.%d"%dom.getMPIRank()),"r").read().split(";")
   t=float(FF[0])          # time stamp
   n=int(FF[3])            # time step counter
   t_vis=float(FF[1])      # 
   n_vis=int(FF[2])
   counter_vis=int(FF[4])
   dt=float(FF[5])
   stress=load(os.path.join(f,"stress.nc"),dom)
   v=load(os.path.join(f,"v.nc"),dom)
   p=load(os.path.join(f,"p.nc"),dom)
   T=load(os.path.join(f,"T.nc"),dom)
   
   nusselt_file=FileWriter(NUSSELT_FN,append=True)
   print "<%s> Restart from file %s at time step %s (t=%s) completed."%(time.asctime(),f,t,n)
else:
  if DIM==2:
    dom=Rectangle(int(ceil(L*NE/H)),NE,l0=L/H,l1=1,order=2, useFullElementOrder=True,optimize=True)
  else:
    dom=Brick(int(ceil(L*NE/H)),int(ceil(L*NE/H)),NE,l0=L/H,l1=L/H,l2=1,order=2, useFullElementOrder=True,optimize=True)
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

  T=(1.-x[DIM-1])+PERT*T
  v=Vector(0,Solution(dom))
  stress=Tensor(0,Function(dom))
  x2=ReducedSolution(dom).getX()
  p=Ra*(x2[DIM-1]-0.5*x2[DIM-1]**2-0.5)
  nusselt_file=FileWriter(NUSSELT_FN,append=False)

p_last=p
x=dom.getX()
#
#   set up heat problem:
#
heat=TemperatureCartesian(dom,useBackwardEuler=False)
print "<%s> Temperature transport has been set up."%time.asctime()
heat.getSolverOptions().setTolerance(T_TOL)
heat.getSolverOptions().setVerbosity(VERBOSE)
fixed_T_at=whereZero(x[DIM-1])+whereZero(H-x[DIM-1])
heat.setInitialTemperature(clip(T,T_0))
heat.setValue(rhocp=1,k=1,given_T_mask=fixed_T_at)
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
    if TOPOGRAPHY_LOAD and d==DIM-1:
       fixed_v_mask+=whereZero(x[d])*unitVector(d,DIM)
    else:
       fixed_v_mask+=(whereZero(x[d])+whereZero(x[d]-ll))*unitVector(d,DIM)
#
#   set up velovity problem
#
flow=IncompressibleIsotropicFlowCartesian(dom, stress=stress, v=v, p=p, t=t, numMaterials=2, verbose=VERBOSE)
flow.setDruckerPragerLaw(tau_Y=TAU_Y/P_REF+BETA*(1.-Function(dom).getX()[DIM-1]))

flow.setElasticShearModulus(MUE)
flow.setTolerance(TOL2)
flow.setFlowTolerance(FLOW_TOL)
flow.setEtaTolerance(FLOW_TOL)
flow.setExternals(fixed_v_mask=fixed_v_mask)
print "<%s> Flow solver has been set up."%time.asctime()
#
#  let the show begin:
#
t1 = time.time()
print "<%s> Start time step %s (t=%s)."%(time.asctime(),n,t)
while t<T_END:
    #======= solve for velovity ====================================================================
    eta_N=exp(Ar*((1.+V_REF*(1-Function(dom).getX()[DIM-1]))/(T_OFFSET_REF+interpolate(T,Function(dom)))-1./T_OFFSET_REF))
    print "viscosity range :", inf(eta_N), sup(eta_N)
    flow.setPowerLaws([eta_N, eta_N ], [ 1., TAU_0],  [1,N])
    flow.setExternals(F=Ra*T*unitVector(DIM-1,DIM)) #f= add the mountanins here
    flow.update(dt, iter_max=100, inner_iter_max=200, verbose=False)
    v=flow.getVelocity()
    for d in range(DIM):
         print "range %d-velocity"%d,inf(v[d]),sup(v[d])
    print "<%s> flow solver completed."%time.asctime()
    if t>=t_vis or n>n_vis:
      saveVTK(os.path.join(VIS_DIR,"state.%d.vtu"%counter_vis),T=T,v=v,eta=flow.getCurrentEtaEff())
      print "<%s> Visualization file %d for time step %e generated."%(time.asctime(),counter_vis,t)
      counter_vis+=1
      t_vis+=DT_VIS
      n_vis+=DN_VIS
    n+=1
    t+=dt
    print "<%s> Time step %s (t=%s) completed."%(time.asctime(),n,t)
    #======= set Temperature problem ====================================================================
    #
    heat.setValue(v=v,Q=CHI_REF*flow.getTau()**2/flow.getCurrentEtaEff())
    dt=heat.getSafeTimeStepSize()
    print "<%s> New time step size is %e"%(time.asctime(),dt)
    print "<%s> Start time step %s (t=%s)."%(time.asctime(),n+1,t+dt)
    #
    #  solve temperature:
    #
    T=heat.getSolution(dt)
    print "Temperature range ",inf(T),sup(T)
    print "<%s> temperature update completed."%time.asctime()
    #======= anaysis ==================================================================================
    #
    #  .... Nusselt number
    #
    dTdz=grad(T)[DIM-1]
    Nu=1.-integrate(v[DIM-1]*T)/integrate(dTdz)
    print "nusselt number = ",Nu
    nusselt_file.write("%e %e\n"%(t,Nu))
    eta_bar=integrate(flow.getTau())/integrate(flow.getTau()/flow.getCurrentEtaEff())
    Ra_eff= (t_REF*RHO_0*G*H*(T_1-T_0)*ALPHA_0)/eta_bar
    print "av. eta = %e"%eta_bar
    print "effective Rayleigh number = %e"%Ra_eff
    # =========================
    #
    #    create restart files:
    #
    if DN_RESTART>0:
       if (n-1)%DN_RESTART == 0:
         c=(n-1)/DN_RESTART
         old_restart_dir="%s_%s_"%(PREFIX_RESTART,c-1)
         new_restart_dir="%s_%s_"%(PREFIX_RESTART,c)

         if dom.onMasterProcessor() and not os.path.isdir(new_restart_dir): os.mkdir(new_restart_dir)
	 dom.MPIBarrier()
         file(os.path.join(new_restart_dir,"stamp.%d"%dom.getMPIRank()),"w").write("%e; %e; %s; %s; %s; %e;\n"%(t, t_vis, n_vis, n, counter_vis, dt))
         v.dump(os.path.join(new_restart_dir,"v.nc"))
         p.dump(os.path.join(new_restart_dir,"p.nc"))
         T.dump(os.path.join(new_restart_dir,"T.nc"))
         removeRestartDirectory(old_restart_dir)
         print "<%s> Restart files written to %s."%(time.asctime(),new_restart_dir)
print "<%s> Calculation finalized after %s sec."%(time.asctime(),time.time() - t1)

