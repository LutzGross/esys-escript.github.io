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
"""
this is a localization a simulation over domain [0,L] X [0,L] x [0,H]
with a plastic layer above a viscous layer of thickness H_VISC.
The yield condition is perturbed along a line at the boundary between
viscous and plastic layer to trigger localization.
"""

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript.unitsSI import DEG
from esys.escript.models import StokesProblemCartesian
from esys.finley import Rectangle, Brick, LoadMesh
from esys.weipa import saveVTK
from math import pi, ceil
import sys
import time

# ======================= Default Values ==================================================
DIM=2                           # spatial dimension
H=1.                            # height
L=4*H                           # length
NE=10                           # number of elements in H-direction. 
H_VISC=0.2*H                    # height of viscous zone (must aline with element boundary)

ETA_TOP=10.
ETA_BOTTOM=ETA_TOP/10.
TRANSITION_WIDTH = H/NE          # half width of the transition zone between top and bottom (typically = H/NE)
FRICTION_ANGLE=70*DEG
C=None                       # =None interesting value is calculated.
V0=None                        # =None interesting value is calculated.
COMPRESSION=False              
RHO=1.
G=1.
W_WEAK=0.03*L                   # width of weak zone (>H/NE)
ALPHA_WEAK=65*DEG               # angle of week zone against x-axis
OFFSET_X_WEAK=L*0.3              # offset of weak zone 
OFFSET_Z_WEAK=H_VISC              # offset of weak zone 
WEAK_FRACTION=1.0              # reduction factor for weak zone
STOP_FRAC=0.2                   # iteratiun stops when height ~ (1-STOP_FRAC) * H or height ~ (1.+STOP_FRAC)*H 
RESTART = False
VIS_DIR="./data"
TOL=1.e-4
MAX_ITER=100
VERBOSE=True

#
#   derived values:
#
if C==None:
     C=RHO*G*H*(1-sin(FRICTION_ANGLE))/cos(FRICTION_ANGLE)
if V0==None:
     V0=RHO*G*H**2/ETA_TOP * 0.1
if COMPRESSION:
    DIRECTION=-1.
    T_END=(STOP_FRAC+2)/(STOP_FRAC+1)*L/V0
else:
    DIRECTION=1.
    T_END=STOP_FRAC/(1-STOP_FRAC)*L/V0
FRICTION_ANGLE_WEAK=FRICTION_ANGLE*WEAK_FRACTION
C_WEAK=C*WEAK_FRACTION
NE_L=int(ceil(L/H*NE))
H_VISC=int((H_VISC*NE)/H+0.5)*(H/NE)
OFFSET_X_WEAK=int((OFFSET_X_WEAK*NE_L)/L+0.5)*(L/NE_L)
OFFSET_Z_WEAK=int((OFFSET_Z_WEAK*NE)/H+0.5)*(H/NE)
#
#   print input
#
print("spatial dimenstion DIM = ",DIM)
print("height H =",H)
print("length L = ",L)
print("vertical number of element NE = ", NE)
print("horizontal number of element NE_L = ", NE_L)
print("total number of element = ", NE_L**(DIM-1)*NE)
print("height of viscosity layer H_VISC = ",H_VISC)
print("viscosity in top layer ETA_TOP = ",ETA_TOP)
print("viscosity in bottom/viscous layer ETA_BOTTOM = ",ETA_BOTTOM)
print("friction angle FRICTION_ANGLE =",FRICTION_ANGLE/DEG)
print("C =",C)
print("width of weak zone W_WEAK = ",W_WEAK)
print("elements in weak zone = ",W_WEAK/H*NE)
print("strike of weak zone ALPHA_WEAK = ",ALPHA_WEAK/DEG) 
print("x-offset of weak zone OFFSET_X_WEAK = ",OFFSET_X_WEAK)
print("z-offset of weak zone OFFSET_Z_WEAK = ",OFFSET_Z_WEAK)
print("friction angle in weak zone FRICTION_ANGLE_WEAK =",FRICTION_ANGLE_WEAK/DEG)
print("C in weak zone C_WEAK =",C_WEAK)
print("density RHO = ",RHO)
print("Gravity G = ",G)
if COMPRESSION:
     print("Direction: compression")
else:
     print("Direction: expansion")
print("surface velocity = ",V0)
print("end time T_END = ",T_END)
mkDir(VIS_DIR)
#
#   set up geomtry:
#
#
if RESTART:
   dom=LoadMesh("mesh.nc")
   if not dom.getDim() == DIM:
         raise ValueError("Expected DIM and dimension of mesh in restart file don't match.")

   FF=open("stamp.csv","r").read().split(",")
   t=float(FF[0])          # time stamp
   t_vis=float(FF[1])      # 
   n_vis=int(FF[2])
   n=int(FF[3])            # time step counter
   counter_vis=int(FF[4])
else:
    if DIM==2:
        dom=Rectangle(NE_L,NE,l0=L,l1=H,order=-1,optimize=True)
    else:
        dom=Brick(NE_L,NE_L,NE,l0=L,l1=L,l2=H,order=-1,optimize=True)
    t=0         
    n=0        
    t_vis=0     
    n_vis=0
    counter_vis=0
#
#   material properties:
#
x=Function(dom).getX()
mask_visc=whereNegative(x[DIM-1]-H_VISC)
mask_visc=clip((H_VISC+TRANSITION_WIDTH-x[DIM-1])/(2.*TRANSITION_WIDTH),maxval=1.,minval=0)

if DIM == 3:
   offset=[OFFSET_X_WEAK, 0, OFFSET_Z_WEAK]
   strike=numpy.array([cos(ALPHA_WEAK),sin(ALPHA_WEAK),0.])
   d2=length(x-offset)**2-inner(x-offset,strike)**2
else:
   offset=[OFFSET_X_WEAK, OFFSET_Z_WEAK]
   d2=length(x-offset)**2
factor_weak=exp(-2*d2/W_WEAK**2)
c=(C_WEAK-C)*factor_weak+C
friction_angle=(FRICTION_ANGLE_WEAK-FRICTION_ANGLE)*factor_weak+FRICTION_ANGLE
#
#   velocity constraints:
#
x=dom.getX()
x0=x[0]
v=DIRECTION*V0/2*(1.-2*(sup(x0)-x0)/(sup(x0)-inf(x0)))*unitVector(0,DIM)
fixed_v_mask=(whereZero(x0-inf(x0))+whereZero(x0-sup(x0)))*unitVector(0,DIM)  + \
             whereZero(x[DIM-1])*unitVector(DIM-1,DIM)
if DIM==3: fixed_v_mask+=(whereZero(x[1])+whereZero(x[1]-L))*unitVector(1,DIM)

p=Scalar(0.,ReducedSolution(dom))
gamma_dot=sqrt(2.)*length(deviatoric(symmetric(grad(v))))
tau=0.
tau_lsup=0

flow=StokesProblemCartesian(dom)
flow.setTolerance(1.e-5)

while n<1:

  print("========= Time step %d ======="%( n+1,))
  m=0
  dtau_lsup =1.
  while dtau_lsup > TOL * tau_lsup:
     print("--- iteration step %d ---- "%m)
     if n==0 and m==0:
        eta_top=ETA_TOP
     else:
        tau_Y=clip(c*cos(friction_angle)+p*sin(friction_angle), minval=0.)
        eta_top=clip(safeDiv(tau_Y,gamma_dot),maxval=ETA_TOP)
        print("eta_top=",eta_top)
     eta_eff=eta_top*(1-mask_visc) + ETA_BOTTOM*mask_visc
     print("eta_eff=",eta_eff)
     flow.initialize(fixed_u_mask=fixed_v_mask,eta=eta_eff,f=-RHO*G*unitVector(DIM-1,DIM))
     v,p=flow.solve(v,-3.*p,max_iter=MAX_ITER,verbose=VERBOSE,usePCG=True)
     print("p=",p)
     p*=-(1./3.)
     gamma_dot=sqrt(2.)*length(deviatoric(symmetric(grad(v))))
     tau, tau_old = eta_eff*gamma_dot, tau
     dtau_lsup=Lsup(tau-tau_old)
     tau_lsup=Lsup(tau)
     print("increment tau = ",dtau_lsup,tau_lsup, dtau_lsup > TOL * tau_lsup, TOL)
     print("flux balance = ",integrate(inner(v,dom.getNormal())))
     m+=1
     if m>MAX_ITER:
        raise ValueError("no convergence.")
  print("iteration complted after %d steps"%(m,))
  print(p)
  saveVTK("test.vtu",v=v, p=p, eta=eta_eff, tau=tau);

  print("Max velocity =", Lsup(v))
  #Courant condition
  #dt=0.4*h/(Lsup(velocity))
  n+=1

1/0


     
#
#   set up heat problem:
#
heat=TemperatureCartesian(dom,)
print("<%s> Temperature transport has been set up."%time.asctime())
heat.getSolverOptions().setTolerance(T_TOL)
heat.getSolverOptions().setVerbosity(VERBOSE)
fixed_T_at=whereZero(x[DIM-1])+whereZero(H-x[DIM-1])
heat.setInitialTemperature(clip(T,T_0))
heat.setValue(rhocp=1,k=1,given_T_mask=fixed_T_at)
#
#   velocity constraints:
#
fixed_v_mask=Vector(0,Solution(dom))
faces=Scalar(0.,Solution(dom))
for d in range(DIM):
    if d == DIM-1: 
       ll = H
    else:
       ll = L
    if CREATE_TOPOGRAPHY and d==DIM-1:
       fixed_v_mask+=whereZero(x[d])*unitVector(d,DIM)
    else:
       s=whereZero(x[d])+whereZero(x[d]-ll)
       faces+=s
       fixed_v_mask+=s*unitVector(d,DIM)
#
#   set up velovity problem
#
flow=IncompressibleIsotropicFlowCartesian(dom, stress=stress, v=v, p=p, t=t, numMaterials=2, verbose=VERBOSE)
flow.setDruckerPragerLaw(tau_Y=TAU_Y/P_REF+BETA*(1.-Function(dom).getX()[DIM-1]))

flow.setElasticShearModulus(MUE)
flow.setTolerance(FLOW_TOL)
flow.setEtaTolerance(FLOW_TOL)
flow.setExternals(fixed_v_mask=fixed_v_mask)
print("<%s> Flow solver has been set up."%time.asctime())
#
# topography set-up
#
boundary=FunctionOnBoundary(dom).getX()[DIM-1]
top_boundary_mask=whereZero(boundary-sup(boundary))
surface_area=integrate(top_boundary_mask)
if CREATE_TOPOGRAPHY:
    mts=Mountains(dom,eps=TOPO_SMOOTH)
    mts.setTopography(topography)
    print("<%s> topography has been set up."%time.asctime())
#
#  let the show begin:
#
t1 = time.time()
print("<%s> Start time step %s (t=%s)."%(time.asctime(),n,t))
while t<T_END:
    if CREATE_TOPOGRAPHY: topography_old=topography
    v_old, p_old, stress_old=v, p, stress
    T_old=T
    #======= solve for velovity ====================================================================
    eta_N=exp(Ar*((1.+V_REF*(1-Function(dom).getX()[DIM-1]))/(T_OFFSET_REF+interpolate(T,Function(dom)))-1./T_OFFSET_REF))
    print("viscosity range :", inf(eta_N), sup(eta_N))
    flow.setPowerLaws([eta_N, eta_N ], [ 1., TAU_0],  [1,N])
    flow.setExternals(F=Ra*T*unitVector(DIM-1,DIM))
    # if dt<=0 or not CREATE_TOPOGRAPHY:
    if not CREATE_TOPOGRAPHY:
            flow.update(dt, iter_max=100, verbose=False)
    else:
        topography_last=topography
        Topo_norm, error_Topo=1,1
        i=0
        print("DDDDD : ====",dt)
        while error_Topo > TOPO_TOL * Topo_norm:
            flow.setStatus(t, v_old, p_old, stress_old)
            flow.setExternals(f=-SURFACE_LOAD*(topography-dt*v)*unitVector(DIM-1,DIM)*top_boundary_mask, restoration_factor=SURFACE_LOAD*dt*top_boundary_mask) 
            flow.update(dt, iter_max=100, verbose=False)
            v=flow.getVelocity()
            mts.setTopography(topography_old)
            mts.setVelocity(v)
            topography_last, topography=topography, mts.update(dt, allow_substeps=True)
            error_Topo=sqrt(integrate(((topography-topography_last)*top_boundary_mask)**2))
            Topo_norm=sqrt(integrate((topography*top_boundary_mask)**2))
            print("DDDDD :", "input=",integrate(v*top_boundary_mask)[DIM-1],"output=", integrate(topography*top_boundary_mask)/Lsup(topography), error_Topo, Topo_norm)
            print("topography update step %s error = %e, norm = %e."%(i, error_Topo, Topo_norm), Lsup(v))
            i+=1
            if i > TOPO_ITER_MAX: 
               raise RuntimeError("topography did not converge after %s steps."%TOPO_ITER_MAX)
    v=flow.getVelocity()
    for d in range(DIM):
         print("range %d-velocity"%d,inf(v[d]),sup(v[d]))
    print("Courant = ",inf(dom.getSize()/(length(v)+1e-19)), inf(dom.getSize()**2))
    print("<%s> flow solver completed."%time.asctime())
    n+=1
    t+=dt
    # print "influx= ",integrate(inner(v,dom.getNormal())), sqrt(integrate(length(v)**2,FunctionOnBoundary(dom))), integrate(1., FunctionOnBoundary(dom))
    print("<%s> Time step %s (t=%s) completed."%(time.asctime(),n,t))
    #======= setup Temperature problem ====================================================================
    #
    heat.setValue(v=v,Q=CHI_REF*flow.getTau()**2/flow.getCurrentEtaEff())
    dt=heat.getSafeTimeStepSize()
    print("<%s> New time step size is %e"%(time.asctime(),dt))
    if n == 10: 1/0
    #======= set-up topography ==================================================================================
    if CREATE_TOPOGRAPHY:
        dt=min(mts.getSafeTimeStepSize()*0.5,dt)
        print("<%s> New time step size is %e"%(time.asctime(),dt))
    print("<%s> Start time step %s (t=%s)."%(time.asctime(),n+1,t+dt))
    #
    #  solve temperature:
    #
    T=heat.getSolution(dt)
    print("Temperature range ",inf(T),sup(T))
    print("<%s> temperature update completed."%time.asctime())
    #======= anaysis ==================================================================================
    #
    #  .... Nusselt number
    #
    dTdz=grad(T)[DIM-1]
    Nu=1.-integrate(v[DIM-1]*T)/integrate(dTdz)
    eta_bar=integrate(flow.getTau())/integrate(flow.getTau()/flow.getCurrentEtaEff())
    Ra_eff= (t_REF*RHO_0*G*H*(T_1-T_0)*ALPHA_0)/eta_bar
    print("nusselt number = %e"%Nu)
    print("av. eta = %e"%eta_bar)
    print("effective Rayleigh number = %e"%Ra_eff)
    if CREATE_TOPOGRAPHY:
       topo_level=integrate(topography*top_boundary_mask)/surface_area
       valleys_deep=inf(topography)
       mountains_heigh=sup(topography)
       print("topography level = ",topo_level)
       print("valleys deep = ", valleys_deep)
       print("mountains_heigh = ", mountains_heigh)
       diagnostics_file.write("%e %e %e %e %e %e %e\n"%(t,Nu, topo_level, valleys_deep, mountains_heigh, eta_bar, Ra_eff))
    else:
       diagnostics_file.write("%e %e %e %e\n"%(t,Nu, eta_bar, Ra_eff))
    #
    #  .... visualization
    #
    if t>=t_vis or n>n_vis:
      if CREATE_TOPOGRAPHY:
         saveVTK(os.path.join(VIS_DIR,"state.%d.vtu"%counter_vis),T=T,v=v,eta=flow.getCurrentEtaEff(), topography=topography_old, vex=mts.getVelocity())
      else:
         saveVTK(os.path.join(VIS_DIR,"state.%d.vtu"%counter_vis),T=T,v=v,eta=flow.getCurrentEtaEff())
      print("<%s> Visualization file %d for time step %e generated."%(time.asctime(),counter_vis,t))
      counter_vis+=1
      t_vis+=DT_VIS
      n_vis+=DN_VIS
    # =========================
    #
    #    create restart files:
    #
    if DN_RESTART>0:
       if (n-1)%DN_RESTART == 0:
         c=(n-1)/DN_RESTART
         old_restart_dir="%s_%s_"%(PREFIX_RESTART,c-1)
         new_restart_dir="%s_%s_"%(PREFIX_RESTART,c)
         mkDir(new_restart_dir)
         dom.MPIBarrier()
         file(os.path.join(new_restart_dir,"stamp.%d"%dom.getMPIRank()),"w").write("%e; %e; %s; %s; %s; %e;\n"%(t, t_vis, n_vis, n, counter_vis, dt))
         v.dump(os.path.join(new_restart_dir,"v.nc"))
         p.dump(os.path.join(new_restart_dir,"p.nc"))
         T.dump(os.path.join(new_restart_dir,"T.nc"))
         if CREATE_TOPOGRAPHY: topography.dump(os.path.join(new_restart_dir,"topo.nc"))
         removeRestartDirectory(old_restart_dir)
         print("<%s> Restart files written to %s."%(time.asctime(),new_restart_dir))
print("<%s> Calculation finalized after %s sec."%(time.asctime(),time.time() - t1))

