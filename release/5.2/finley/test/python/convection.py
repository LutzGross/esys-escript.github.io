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
this is a convection simulation over a domain [0,L] X [0,L] x [0,H]

It is solved in dimensionless form

"""

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript import DataManager as DM
from esys.escript.models import TemperatureCartesian, IncompressibleIsotropicFlowCartesian, Mountains, SubSteppingException
from esys.finley import Rectangle, Brick, LoadMesh
from optparse import OptionParser
from math import pi, ceil
import sys
import time

# ============================= Default Values ===============================

DIM=2                     # spatial dimension
H=1.                      # height
L=2*H                     # length
NE=30                     # number of elements in H-direction
PERT=0.15                 # initial temperature perturbation
DT=1.e-4                  # initial time step size
CREATE_TOPO=False         # create topography
DT_MIN=1.e-10             # minimum time step size
T_END=10.                 # end time

RHO_0=100.                # surface density (lf ~ RHO_0**2) 
G=1.                      # gravitational constant
ALPHA_0=0.1               # thermal expansion coefficient (Ra ~ RHO_0**2 * ALPHA_0 = lf * ALPHA_0)
T_0=0.                    # surface temperature
T_1=1.                    # bottom temperature
C_P=1                     # heat capacity
K=1.                      # thermal conductivity
CHI=0.                    # Taylor-Quinny coefficient
MUE=None                  # elastic shear modulus
TAU_Y=5*10**(2.5)         # Drucker-Prager cohesion factor
BETA=0                    # Drucker-Prager friction factor
TAU_0=2*10**(2.5)         # transition stress
N=3                       # power for power law

E=23*0                    # activation energy
V=18*0                    # activation volume 
T_OFFSET=1                # temperature offset on surface (dimensionless formulation T_OFFSET=1 otherwise =0)
R=1                       # gas constant
ETA_N0=1.                 # viscosity at surface 

TOPO_SMOOTH=1e-5          # smoothing factor for extrapolation of surface velocity to interior
T_TOL=1.e-4               # tolerance temperature transport
FLOW_TOL=1.e-3            # tolerance for inconcompressible flow solver
TOPO_TOL=0.1              # tolerance for update of topography
DIAGNOSTICS_FN="diag.csv" # filename for diagnostics
VERBOSE=False             # output more messages from solvers
DT_VIS=T_END/500          # time difference between visualization files
DN_VIS=5                  # max. number of steps between visualization files
DN_RESTART=1000           # create restart files every DN_RESTART steps
PREFIX_RESTART="restart"  # name prefix for restart directories
TOPO_ITER_MAX=20          # max. number of iteration steps to update topography

# ============================================================================

#
# read options:
#
parser = OptionParser(usage="%prog [Options]")
parser.add_option("-r", "--restart", dest="restart",
        help="restart from latest checkpoint directory. It will be deleted after new data is exported.", default=False, action="store_true")
parser.add_option("-d", "--dir", dest="restart_dir",
        help="locate/create restart directories under DIR.", metavar="DIR", default='.')
parser.add_option("-p", "--param", dest="param",
        help="name of file to be imported ", metavar="PARAM", default=None)
(options, args) = parser.parse_args()
restart=options.restart

#
#  overwrite the default options:
#
print(("<%s> Execution started."%time.asctime()))
if options.param !=None: 
    exec(open(options.param,'r'))
    print(("Parameters imported from file ",options.param))

print("Input Parameters:")
print(("\tDimension                     DIM\t\t= %d"%DIM))
print(("\tHeight                        H\t\t\t= %s"%H))
print(("\tLength                        L\t\t\t= %s"%L))
print(("\tElements in H                 NE\t\t= %d"%NE))
print(("\tTemperature perturbation      PERT\t\t= %s"%PERT))
print(("\tInitial time step size        DT\t\t= %s"%DT))
print(("\tMinimum time step size        DT_MIN\t\t= %s"%DT_MIN))
print(("\tEnd time                      T_END\t\t= %s"%T_END))
print(("\tCreate topography             CREATE_TOPO\t= %s"%CREATE_TOPO))
print(("\tSurface density               RHO_0\t\t= %s"%RHO_0))
print(("\tGravitational constant        G\t\t\t= %s"%G))
print(("\tThermal expansion coefficient ALPHA_0\t\t= %s"%ALPHA_0))
print(("\tSurface temperature           T_0\t\t= %s"%T_0))
print(("\tBottom temperature            T_1\t\t= %s"%T_1))
print(("\tHeat capacity                 C_P\t\t= %s"%C_P))
print(("\tThermal conductivity          K\t\t\t= %s"%K))
print(("\tTaylor-Quinny coefficient     CHI\t\t= %s"%CHI))
print(("\tElastic shear modulus         MUE\t\t= %s"%MUE))
print(("\tCohesion factor               TAU_Y\t\t= %s"%TAU_Y))
print(("\tFriction factor               BETA\t\t= %s"%BETA))
print(("\tTransition stress             TAU_0\t\t= %s"%TAU_0))
print(("\tPower for power law           N\t\t\t= %s"%N))
print(("\tViscosity at surface          ETA_N0\t\t= %s"%ETA_N0))
print(("\tActivation energy             E\t\t\t= %s"%E))
print(("\tActivation volume             V\t\t\t= %s"%V))
print(("\tTemperature offset            T_OFFSET\t\t= %s"%T_OFFSET))
print(("\tGas constant                  R\t\t\t= %s"%R))
print(("\tTopography smoothing          TOPO_SMOOTH\t= %s"%TOPO_SMOOTH))
print(("\tTolerance for topography      TOPO_TOL\t\t= %s"%TOPO_TOL))
print(("\tTransport tolerance           T_TOL\t\t= %s"%T_TOL))
print(("\tFlow tolerance                FLOW_TOL\t\t= %s"%FLOW_TOL))
#print("\tFile for diagnostics          DIAGNOSTICS_FN\t= %s"%DIAGNOSTICS_FN)
print(("\tRestart counter increment     DN_RESTART\t= %d"%DN_RESTART))
print(("\tPrefix for restart dirs       PREFIX_RESTART\t= %s"%PREFIX_RESTART))
print(("\tVerbosity                     VERBOSE\t\t= %s"%VERBOSE))

print("Control Parameters:")
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
if CREATE_TOPO:
    SURFACE_LOAD=RHO_0*G*H/P_REF
else:
    SURFACE_LOAD=0.
if MUE == None:
    De=None
else:
    De=ETA_N0/MUE/t_REF
ETA_BOT=exp(Ar*((1.+V_REF)/(T_OFFSET_REF+1)-1./T_OFFSET_REF))*ETA_N0
print(("\tTotal #elements               \t\t\t= %d"%(NE**DIM*int(L/H)**(DIM-1))))
print(("\tReference time                t_REF\t\t= %s"%t_REF))
print(("\tReference pressure            P_REF\t\t= %s"%P_REF))
print(("\tReference Taylor-Quinny       CHI_REF\t\t= %s"%CHI_REF))
print(("\tDissipation number            DI\t\t= %s"%Di))
print(("\tRayleigh number surface       Ra\t\t= %s"%Ra))
print(("\tDebora number surface         De\t\t= %s"%De))
print(("\tBottom viscosity              \t\t\t= %s"%ETA_BOT))
print(("\tRayleigh number bottom        \t\t\t= %s"%(RHO_0*G*H*(T_1-T_0)*ALPHA_0*t_REF/ETA_BOT)))
if MUE == None:
   print("\tDebora number bottom          \t\t\t= None")
else:
   print(("\tDebora number bottom          \t\t\t= %s"%(ETA_BOT/MUE/t_REF)))
print(("\tArrhenius                     Ar\t\t= %s"%Ar))
print(("\tSurface load factor           SURFACE_LOAD\t= %s"%SURFACE_LOAD))
print(("\tScaled activation volume      V_REF\t\t= %s"%V_REF))
print()

# some control variables (will be overwritten in case of a restart:
t=0          # time stamp
n=0          # time step counter
t_vis=DT_VIS # time of next visualization file export
dt=DT        # current time step size
#=========================
#
#   set up domain and data from scratch or from restart files
#
dataMgr=DM(formats=[DM.RESTART,DM.VTK], work_dir=options.restart_dir, restart_prefix=PREFIX_RESTART, do_restart=restart)
dataMgr.setCheckpointFrequency(DN_RESTART/DN_VIS)
if dataMgr.hasData():
    dom=dataMgr.getDomain()
    t=dataMgr.getValue('t')
    t_vis=dataMgr.getValue('t_vis')
    n=dataMgr.getValue('n')
    dt=dataMgr.getValue('dt')
    stress=dataMgr.getValue('stress')
    v=dataMgr.getValue('v')
    p=dataMgr.getValue('p')
    T=dataMgr.getValue('T')
    if CREATE_TOPO:
        topography=dataMgr.getValue('topography')
   
    #diagnostics_file=FileWriter(DIAGNOSTICS_FN,append=True)
    print(("<%s> Restart at time step %d (t=%e) completed."%(time.asctime(),n,t)))
else:
    if DIM==2:
        dom=Rectangle(int(ceil(L*NE/H)),NE,l0=L/H,l1=1,order=-1,optimize=True)
    else:
        dom=Brick(int(ceil(L*NE/H)),int(ceil(L*NE/H)),NE,l0=L/H,l1=L/H,l2=1,order=-1,optimize=True)
    x=dom.getX()
    T=Scalar(1,Solution(dom))
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

    if CREATE_TOPO:
        topography=Scalar(0.,Solution(dom))
    #diagnostics_file=FileWriter(DIAGNOSTICS_FN,append=False)
    #diagnostics_file.write("Ra = %e Lambda= %e\n"%(Ra, SURFACE_LOAD))

p_last=p
x=dom.getX()
#
#   set up heat problem:
#
heat=TemperatureCartesian(dom)
print(("<%s> Temperature transport has been set up."%time.asctime()))
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
    if CREATE_TOPO and d==DIM-1:
        fixed_v_mask+=whereZero(x[d])*unitVector(d,DIM)
    else:
        s=whereZero(x[d])+whereZero(x[d]-ll)
        faces+=s
        fixed_v_mask+=s*unitVector(d,DIM)
#
#   set up velocity problem
#
flow=IncompressibleIsotropicFlowCartesian(dom, stress=stress, v=v, p=p, t=t, numMaterials=2, verbose=VERBOSE)
flow.setDruckerPragerLaw(tau_Y=TAU_Y/P_REF+BETA*(1.-Function(dom).getX()[DIM-1]))

flow.setElasticShearModulus(MUE)
flow.setTolerance(FLOW_TOL)
flow.setEtaTolerance(FLOW_TOL)
flow.setExternals(fixed_v_mask=fixed_v_mask)
print(("<%s> Flow solver has been set up."%time.asctime()))
#
#   topography setup
#
boundary=FunctionOnBoundary(dom).getX()[DIM-1]
top_boundary_mask=whereZero(boundary-sup(boundary))
surface_area=integrate(top_boundary_mask)
if CREATE_TOPO:
    mts=Mountains(dom,eps=TOPO_SMOOTH)
    mts.setTopography(topography)
    print(("<%s> topography has been set up."%time.asctime()))

#
#   let the show begin:
#
t1 = time.time()
print(("<%s> Start time step %d (t=%s)."%(time.asctime(),n,t)))
while t<T_END:
    if CREATE_TOPO: topography_old=topography
    v_old, p_old, stress_old=v, p, stress
    T_old=T
    #========================= solve for velocity ============================
    eta_N=exp(Ar*((1.+V_REF*(1-Function(dom).getX()[DIM-1]))/(T_OFFSET_REF+interpolate(T,Function(dom)))-1./T_OFFSET_REF))
    #print("viscosity range :", inf(eta_N), sup(eta_N))
    flow.setPowerLaws([eta_N, eta_N ], [ 1., TAU_0],  [1,N])
    flow.setExternals(F=Ra*T*unitVector(DIM-1,DIM))
    # if dt<=0 or not CREATE_TOPO:
    if not CREATE_TOPO:
        flow.update(dt, iter_max=100, verbose=False)
    else:
        topography_last=topography
        Topo_norm, error_Topo=1,1
        i=0
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
            #print("topography update step %d error = %e, norm = %e."%(i, error_Topo, Topo_norm), Lsup(v))
            i+=1
            if i > TOPO_ITER_MAX: 
                raise RuntimeError("topography did not converge after %d steps."%TOPO_ITER_MAX)
    v=flow.getVelocity()
    #for d in range(DIM):
         #print("range %d-velocity "%d,inf(v[d]),sup(v[d]))
    #print("Courant = ",inf(dom.getSize()/(length(v)+1e-19)), inf(dom.getSize()**2))
    print("<%s> flow solver completed."%time.asctime())
    n+=1
    t+=dt
    #print("influx= ",integrate(inner(v,dom.getNormal())), sqrt(integrate(length(v)**2,FunctionOnBoundary(dom))), integrate(1., FunctionOnBoundary(dom)))
    print("<%s> Time step %d (t=%e) completed."%(time.asctime(),n,t))

    #====================== setup temperature problem ========================
    heat.setValue(v=v,Q=CHI_REF*flow.getTau()**2/flow.getCurrentEtaEff())
    dt=heat.getSafeTimeStepSize()
    print("<%s> New time step size is %e"%(time.asctime(),dt))
    #=========================== setup topography ============================
    if CREATE_TOPO:
        dt=min(mts.getSafeTimeStepSize()*0.5,dt)
        #print("<%s> New time step size is %e"%(time.asctime(),dt))
    #print("<%s> Start time step %d (t=%e)."%(time.asctime(),n+1,t+dt))
    #
    #  solve temperature:
    #
    T=heat.getSolution(dt)
    print("Temperature range ",inf(T),sup(T))
    print("<%s> temperature update completed."%time.asctime())
    #============================== analysis =================================
    #
    #  .... Nusselt number
    #
    dTdz=grad(T)[DIM-1]
    Nu=1.-integrate(v[DIM-1]*T)/integrate(dTdz)
    eta_bar=integrate(flow.getTau())/integrate(flow.getTau()/flow.getCurrentEtaEff())
    Ra_eff= (t_REF*RHO_0*G*H*(T_1-T_0)*ALPHA_0)/eta_bar
    #print("nusselt number = ",Nu)
    #print("avg. eta = ",eta_bar)
    #print("effective Rayleigh number = ",Ra_eff)
    if CREATE_TOPO:
        topo_level=integrate(topography*top_boundary_mask)/surface_area
        valleys_deep=inf(topography)
        mountains_heigh=sup(topography)
        #print("topography level = ",topo_level)
        #print("valleys deep = ",valleys_deep)
        #print("mountains_heigh = ",mountains_heigh)
        #diagnostics_file.write("%e %e %e %e %e %e %e\n"%(t,Nu, topo_level, valleys_deep, mountains_heigh, eta_bar, Ra_eff))
    #else:
        #diagnostics_file.write("%e %e %e %e\n"%(t,Nu, eta_bar, Ra_eff))
    # ========================================================================
    #
    #    create restart/visualization files:
    #
    if t>=t_vis or n%DN_VIS==0:
        dataMgr.setTime(t)
        t_vis+=DT_VIS
        dataMgr.addData(t=t,n=n,t_vis=t_vis,dt=dt,T=T,v=v,eta=flow.getCurrentEtaEff(),stress=stress,p=p)
        if CREATE_TOPO: dataMgr.addData(topography=topography)
        dataMgr.export()
        print(("<%s> Cycle %d (time %s) exported."%(time.asctime(),dataMgr.getCycle(),t)))

print(("<%s> Calculation finalized after %s seconds."%(time.asctime(),time.time()-t1)))

