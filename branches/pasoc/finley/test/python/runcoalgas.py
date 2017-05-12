#######################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
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
Coal Seam gasL ECLIPSE test case
"""

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"


from esys.escript import *
from esys.weipa import saveVTK
from esys.escript import unitsSI as U
from coalgas import *
import time
from esys.finley import ReadMesh, Rectangle, Brick
from esys.escript.pdetools import Locator

SAVE_VTK=True and False
CONST_G = 9.81 * U.m/U.sec**2
P_0=1.*U.atm

CELL_X=2640*U.ft
CELL_Y=2640*U.ft
CELL_Z=33*U.ft

TOP=2310*U.ft

L_Z=CELL_Z

L_X=CELL_X*210
L_Y=CELL_Y*210
L_X=CELL_X*70
L_Y=CELL_Y*70

N_X=int(L_X/CELL_X)
N_Y=int(L_Y/CELL_Y)
N_Z=int(L_Z/CELL_Z)

OUTPUT_DIR="results"

PERM_F_X = 100 * U.mDarcy
PERM_F_Y = 100 * U.mDarcy
PERM_F_Z = 1e-4 * U.mDarcy

EQUIL = {
    "DATUM_DEPTH" : 2310. * U.ft,
    "DATUM_PRESS" : 1000. * U.psi,
    "GWC_DEPTH" : -1000. * U.ft,
    "GWC_PCOW" : 0. * U.psi
}

TOPS = 2310 * U.ft
PHI_F_0=0.01
SIGMA = 1. /U.ft**2

#DT=[0.1* U.day]*9+[1 * U.day,3* U.day,9* U.day, 17.5*U.day] + [ 30.5*U.day ] *20
DT=[1 * U.day,3* U.day,9* U.day, 17.5*U.day] + [ 30.5*U.day ] *20
DT=[.5 * U.day, .5 * U.day, 0.5*3* U.day, 0.5*3* U.day, 0.5*9* U.day, 0.5*9* U.day, 17.5 *.5 *U.day, 17.5*0.5*U.day] + [ 15.25*U.day ] *40
DT=[.25 * U.day, .25 * U.day, .25 * U.day, .25 * U.day,
     0.25*3* U.day, 0.25*3* U.day, 0.25*3* U.day, 0.25*3* U.day,
     0.25*9* U.day, 0.25*9* U.day, 0.25*9* U.day, 0.25*9* U.day,
     17.5 *.25 *U.day, 17.5*0.25*U.day, 17.5 *.25 *U.day, 17.5*0.25*U.day] + \
    [ 15.25*U.day * 0.5] *80
DT=[.125 * U.day, .125 * U.day, .125 * U.day, .125 * U.day, .125 * U.day, .125 * U.day, .125 * U.day, .125 * U.day,
     0.125*3* U.day, 0.125*3* U.day, 0.125*3* U.day, 0.125*3* U.day, 0.125*3* U.day, 0.125*3* U.day, 0.125*3* U.day, 0.125*3* U.day,
     0.125*9* U.day, 0.125*9* U.day, 0.125*9* U.day, 0.125*9* U.day, 0.125*9* U.day, 0.125*9* U.day, 0.125*9* U.day, 0.125*9* U.day,
     17.5 *.125 *U.day, 17.5*0.125*U.day, 17.5 *.125 *U.day, 17.5*0.125*U.day, 17.5 *.125 *U.day, 17.5*0.125*U.day, 17.5 *.125 *U.day, 17.5*0.125*U.day] + \
    [ 15.25*U.day * 0.25] *160
DT=[.125 * U.day, .125 * U.day, .125 * U.day, .125 * U.day, .125 * U.day, .125 * U.day, .125 * U.day, .125 * U.day]
# DT=[1 * U.day,2* U.day] + [4*U.day ] *200
DT=[.25 * U.day, .25 * U.day, .25 * U.day, .25 * U.day]
#DT=[.125 * U.day, .125 * U.day, .125 * U.day, .125 * U.day, .125 * U.day, .125 * U.day, .125 * U.day, .125 * U.day]
DT=[1./10. * U.day]*10 + [3./10. * U.day]*10 + [ 9./10.* U.day ] *10 + [ 17.5/10 *U.day] * 10 + [ 30.5/10*U.day ] *200

#[0.1 * U.day ] *20

PVTW={ "p_ref" :   1000 * U.psi,
       "B_ref" :  0.997,
       "C" :  3.084E-06/U.psi,
       "mu_ref" : 0.68673 * U.cPoise,
       "C_v" : 0/U.psi
     }

GRAVITY = { "water" : 1.0,
            "gas" : .553 }

ROCK = { "p_ref" :   1000 * U.psi,
         "C" : 3.3E-4 * 1./U.psi }

DIFFCOAL = { "D" : 0.005 * U.ft**2/U.day,
             "f_r": 1.}

LANGMUIR = [
[ 0     * U.psi , 0.00000000 * U.Mscf/U.ft**3],
[ 100   * U.psi , 0.00213886 * U.Mscf/U.ft**3],
[ 200   * U.psi , 0.00383259 * U.Mscf/U.ft**3],
[ 300   * U.psi , 0.00520706 * U.Mscf/U.ft**3],
[ 400   * U.psi , 0.00634474 * U.Mscf/U.ft**3],
[ 500   * U.psi , 0.00730199 * U.Mscf/U.ft**3],
[ 600   * U.psi , 0.00811857 * U.Mscf/U.ft**3],
[ 700   * U.psi , 0.00882336 * U.Mscf/U.ft**3],
[ 800   * U.psi , 0.00943786 * U.Mscf/U.ft**3],
[ 900   * U.psi , 0.00997836 * U.Mscf/U.ft**3],
[ 1000  * U.psi , 0.01045748 * U.Mscf/U.ft**3],
[ 1200  * U.psi , 0.01126912 * U.Mscf/U.ft**3] ]


PVDG = [
[ 14.70 * U.psi ,200.3800 * U.Barrel/U.Mscf , 0.012025 * U.cPoise ] , # psi, rb/Mscf,
[ 20.00 * U.psi ,146.0600 * U.Barrel/U.Mscf , 0.012030 * U.cPoise ] ,
[ 25.00 * U.psi ,116.1461 * U.Barrel/U.Mscf , 0.012034 * U.cPoise ] ,
[ 30.00 * U.psi ,96.3132 * U.Barrel/U.Mscf , 0.012038 * U.cPoise ] ,
[ 35.00 * U.psi ,82.2113 * U.Barrel/U.Mscf , 0.012043 * U.cPoise ] ,
[ 49.33 * U.psi ,57.7891 * U.Barrel/U.Mscf , 0.012055 * U.cPoise ] ,
[ 59.00 * U.psi ,48.0866 * U.Barrel/U.Mscf , 0.012064 * U.cPoise ] ,
[ 69.00 * U.psi ,40.9441 * U.Barrel/U.Mscf , 0.012073 * U.cPoise ] ,
[ 75.00 * U.psi ,37.5839 * U.Barrel/U.Mscf , 0.012078 * U.cPoise ] ,
[ 83.00 * U.psi ,33.8685 * U.Barrel/U.Mscf , 0.012085 * U.cPoise ] ,
[ 90.00 * U.psi ,31.1661 * U.Barrel/U.Mscf , 0.012092 * U.cPoise ] ,
[ 95.00 * U.psi ,29.4827 * U.Barrel/U.Mscf , 0.012097 * U.cPoise ] ,
[ 100.00 * U.psi ,27.9698 * U.Barrel/U.Mscf , 0.012101 * U.cPoise ] ,
[ 105.00 * U.psi ,26.6028 * U.Barrel/U.Mscf , 0.012106 * U.cPoise ] ,
[ 118.60 * U.psi ,23.4749 * U.Barrel/U.Mscf , 0.012119 * U.cPoise ] ,
[ 120.00 * U.psi ,23.1937 * U.Barrel/U.Mscf , 0.012120 * U.cPoise ] ,
[ 140.00 * U.psi ,19.7977 * U.Barrel/U.Mscf , 0.012140 * U.cPoise ] ,
[ 153.23 * U.psi ,18.0443 * U.Barrel/U.Mscf , 0.012153 * U.cPoise ] ,
[ 160.00 * U.psi ,17.2607 * U.Barrel/U.Mscf , 0.012159 * U.cPoise ] ,
[ 170.00 * U.psi ,16.2188 * U.Barrel/U.Mscf , 0.012169 * U.cPoise ] ,
[ 187.86 * U.psi ,14.6373 * U.Barrel/U.Mscf , 0.012188 * U.cPoise ] ,
[ 222.49 * U.psi ,12.3027 * U.Barrel/U.Mscf , 0.012224 * U.cPoise ] ,
[ 257.13 * U.psi ,10.6038 * U.Barrel/U.Mscf , 0.012262 * U.cPoise ] ,
[ 291.76 * U.psi ,9.3134 * U.Barrel/U.Mscf , 0.012301 * U.cPoise ] ,
[ 326.39 * U.psi ,8.3001 * U.Barrel/U.Mscf , 0.012341 * U.cPoise ] ,
[ 361.02 * U.psi ,7.4835 * U.Barrel/U.Mscf , 0.012383 * U.cPoise ] ,
[ 395.66 * U.psi ,6.8114 * U.Barrel/U.Mscf , 0.012425 * U.cPoise ] ,
[ 430.29 * U.psi ,6.2491 * U.Barrel/U.Mscf , 0.012470 * U.cPoise ] ,
[ 464.92 * U.psi ,5.7715 * U.Barrel/U.Mscf , 0.012515 * U.cPoise ] ,
[ 499.55 * U.psi ,5.3610 * U.Barrel/U.Mscf , 0.012562 * U.cPoise ] ,
[ 534.19 * U.psi ,5.0043* U.Barrel/U.Mscf  , 0.012610 * U.cPoise ] ,
[ 568.82 * U.psi ,4.6917 * U.Barrel/U.Mscf , 0.012659 * U.cPoise ] ,
[ 603.45 * U.psi ,4.4154 * U.Barrel/U.Mscf , 0.012710 * U.cPoise ] ,
[ 638.08 * U.psi ,4.1695 * U.Barrel/U.Mscf , 0.012762 * U.cPoise ] ,
[ 672.72 * U.psi ,3.9491 * U.Barrel/U.Mscf , 0.012815 * U.cPoise ] ,
[ 707.35 * U.psi ,3.7507 * U.Barrel/U.Mscf , 0.012869 * U.cPoise ] ,
[ 741.98 * U.psi ,3.5711 * U.Barrel/U.Mscf , 0.012925 * U.cPoise ] ,
[ 776.61 * U.psi ,3.4076 * U.Barrel/U.Mscf , 0.012982 * U.cPoise ] ,
[ 811.25 * U.psi ,3.2583 * U.Barrel/U.Mscf , 0.013041 * U.cPoise ] ,
[ 845.88 * U.psi ,3.1214 * U.Barrel/U.Mscf , 0.013100 * U.cPoise ] ,
[ 880.51 * U.psi ,2.9953 * U.Barrel/U.Mscf , 0.013161 * U.cPoise ] ,
[ 915.14 * U.psi ,2.8790 * U.Barrel/U.Mscf , 0.013223 * U.cPoise ] ,
[ 949.78 * U.psi ,2.7712 * U.Barrel/U.Mscf , 0.013287 * U.cPoise ] ,
[ 984.41 * U.psi ,2.6711 * U.Barrel/U.Mscf , 0.013352 * U.cPoise ] ,
[ 1019.00 * U.psi ,2.5781 * U.Barrel/U.Mscf , 0.013418 * U.cPoise ] ,
[ 1053.70 * U.psi ,2.4909 * U.Barrel/U.Mscf , 0.013486 * U.cPoise ] ,
[ 1088.30 * U.psi ,2.4096 * U.Barrel/U.Mscf , 0.013554 * U.cPoise ] ,
[ 1122.90 * U.psi ,2.3334 * U.Barrel/U.Mscf , 0.013624 * U.cPoise ] ,
[ 1157.60 * U.psi ,2.2616 * U.Barrel/U.Mscf , 0.013696 * U.cPoise ] ,
[ 1192.20 * U.psi ,2.1942 * U.Barrel/U.Mscf , 0.013768 * U.cPoise ] ,
[ 1226.80 * U.psi ,2.1307 * U.Barrel/U.Mscf , 0.013842 * U.cPoise ] ,
[ 1261.50 * U.psi ,2.0705 * U.Barrel/U.Mscf , 0.013917 * U.cPoise ] ,
[ 1296.10 * U.psi ,2.0138 * U.Barrel/U.Mscf , 0.013994 * U.cPoise ] ,
[ 1330.70 * U.psi ,1.9600 * U.Barrel/U.Mscf , 0.014072 * U.cPoise ] ,
[ 1365.40 * U.psi ,1.9089 * U.Barrel/U.Mscf , 0.014151  * U.cPoise ] ]


SGFN = [
[ 0  , 0  , 0 * U.psi],
[ 0.05   , 0  , 0  * U.psi ],
[ 0.1333  , 0.00610   , 0  * U.psi ],
[ 0.2167  , 0.02990   , 0 * U.psi ],
[ 0.3  , 0.0759   , 0  * U.psi],
[ 0.3833  , 0.1471     , 0  * U.psi],
[ 0.46667  , 0.2458     , 0  * U.psi],
[ 0.55  , 0.3739     , 0  * U.psi],
[ 0.6333  , 0.53300    , 0  * U.psi],
[ 0.7167  , 0.7246     , 0  * U.psi],
[ 0.8  , 0.95       , 0  * U.psi] ]
SWFN = [
[ 0.20000  , 0.00000, 0  * U.psi],
[ 0.28330  , 0.03280, 0  * U.psi],
[ 0.36670  , 0.09270, 0  * U.psi],
[ 0.45000  , 0.17030, 0  * U.psi],
[ 0.53330  , 0.26220, 0  * U.psi],
[ 0.61670  , 0.36650, 0  * U.psi],
[ 0.70000  , 0.48170, 0  * U.psi],
[ 0.78330  , 0.60710, 0  * U.psi],
[ 0.86670  , 0.74170, 0  * U.psi],
[ 0.95000  , 0.88500, 0  * U.psi],
[ 1.00000  , 1.00000, 0  * U.psi] ]


wellspecs = {
  'P1' : { "X0" :[ (N_X/2+0.5)*CELL_X,  (N_Y/2+0.5)*CELL_Y, 0.5*CELL_Z],
           "r"  : 0.8333 * U.ft,
           "s"  : 0,
           "Q"  : [0., 2000*U.Barrel/U.day ],
           "BHP" : 75*U.psi,
           "schedule" : [0.*U.yr,  2*U.yr]
         }
}


# print input
print(("<%s> Execution started."%time.asctime()))
DIM=2
domain=Rectangle(N_X, N_Y, l0=L_X, l1=L_Y)

N=1000
for I in wellspecs:
     domain.setTagMap(I, N)
     N+=1
     domain.addDiracPoint(wellspecs[I]["X0"][:DIM], I)
     print(("<%s> Well %s introduced to domain."%(time.asctime(), I)))

#domain=Brick(N_X, N_Y,N_Z,l0=L_X, l1=L_Y,l2=L_Z)

print(("<%s> Domain has been generated."%time.asctime()))

print("length x-direction = %f km"%(sup(domain.getX()[0])/U.km))
print("cell size in x direction = %f m"%(CELL_X/U.m))
print("length y-direction = %f km"%(sup(domain.getX()[1])/U.km))
print("cell size in y direction = %f m"%(CELL_Y/U.m))
print("fracture permeability in x direction= %f mD"%(PERM_F_X/(U.mDarcy)))
print("fracture permeability in y direction= %f mD"%(PERM_F_Y/(U.mDarcy)))
print("fracture permeability in z direction= %f mD"%(PERM_F_Z/(U.mDarcy)))


mkDir(OUTPUT_DIR)

print(("<%s> Mesh set up completed."%time.asctime()))
well_P1=VerticalPeacemanWell('P1', domain, BHP_limit=wellspecs['P1' ]["BHP"],
                                Q=wellspecs['P1']["Q"],
                                r=wellspecs['P1']["r"],
                                X0=[ wellspecs['P1' ]["X0"][0], wellspecs['P1']["X0"][1], wellspecs['P1']["X0"][2]] ,
                                D=[CELL_X, CELL_Y, CELL_Z],
                                perm=[PERM_F_X, PERM_F_Y, PERM_F_Z],
                                schedule=wellspecs['P1']["schedule"],
                                s=wellspecs['P1']["s"])
rho_w = WaterDensity(B_ref=PVTW["B_ref"], p_ref = PVTW["p_ref"], C=PVTW["C"], gravity=GRAVITY["water"])
p_top = EQUIL["DATUM_PRESS"] + P_0
p_bottom=p_top + CONST_G * CELL_Z * rho_w(p_top)

model = PorosityOneHalfModel(domain,
                             phi_f=Porosity(phi_0=PHI_F_0, p_0=(p_bottom +p_top)/2., p_ref=ROCK["p_ref"], C = ROCK["C"]),
                             L_g=InterpolationTable([ l[0] for l in LANGMUIR ], [ l[1] for l in LANGMUIR ] ),
                 perm_f_0=PERM_F_X,
                 perm_f_1=PERM_F_Y,
                 perm_f_2=PERM_F_Z,
                 k_w =InterpolationTable([ l[0] for l in SWFN ], [ l[1] for l in SWFN ], obeyBounds=False ),
                 k_g= InterpolationTable([ l[0] for l in SGFN ], [ l[1] for l in SGFN ], obeyBounds=False ),
                 mu_w = WaterViscosity(mu_ref = PVTW["mu_ref"], p_ref=PVTW["p_ref"], C=PVTW["C_v"]),
                 mu_g = InterpolationTable([ l[0] for l in PVDG ], [ l[2] for l in PVDG ] ),
                 rho_w = rho_w,
                 rho_g=GasDensity( p = [ l[0] for l in PVDG ], B = [ l[1] for l in PVDG ], gravity=GRAVITY["gas"]),
                 sigma=SIGMA,
                 A_mg=DIFFCOAL["D"],
                       f_rg=DIFFCOAL["f_r"],
                 wells=[ well_P1, ], g= CONST_G)
# this needs to be revised:.
model.setInitialState(S_fg=0,  c_mg=None, p_top=p_top, p_bottom=p_bottom)
model.getPDEOptions().setVerbosityOn()
model.getPDEOptions().setSolverMethod(model.getPDEOptions().DIRECT)
model.setIterationControl(iter_max=10, rtol=1.e-4, verbose=True)
print("<%s> Problem set up completed."%time.asctime())
t=0
n_t = 0

p, S_fg, c_mg, BHP, q_gas,q_water =model.getState()

if SAVE_VTK:
    FN=os.path.join(OUTPUT_DIR, "state.%d.vtu"%n_t)
    saveVTK(FN,p=p, S_fg=S_fg, c_mg=c_mg)
    print("<%s> Initial state saved to file %s."%(time.asctime(),FN))
    print(t/U.day, well_P1.locator(p)/U.psi, well_P1.locator(S_fg), well_P1.locator(c_mg)/U.Mscf*U.ft**3)
    print(t/U.day, well_P1.locator(BHP)/U.psi, well_P1.locator(q_gas)/U.Mcf*U.day, well_P1.locator(q_water)/U.Barrel*U.day)


for dt in DT:
    print("<%s>Time step %d, time = %e days started:"%(time.asctime(), n_t+1, (t+dt)/U.day))

    model.update(dt)

    p, S_fg, c_mg, BHP, q_gas,q_water = model.getState()

    if SAVE_VTK:
        FN=os.path.join(OUTPUT_DIR, "state.%d.vtu"%(n_t+1))
        saveVTK(FN,p=p, S_fg=S_fg, c_mg=c_mg)
        print("<%s>State %s saved to file %s."%(time.asctime(),n_t+1,FN))
    print((t+dt)/U.day, well_P1.locator(p)/U.psi, well_P1.locator(S_fg), well_P1.locator(c_mg)/U.Mscf*U.ft**3)
    print((t+dt)/U.day, well_P1.locator(BHP)/U.psi, well_P1.locator(q_gas)/U.Mcf*U.day, well_P1.locator(q_water)/U.Barrel*U.day)

    n_t += 1
    t += dt

