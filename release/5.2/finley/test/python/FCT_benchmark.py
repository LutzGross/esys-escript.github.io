# -*- coding: utf-8 -*-
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

#
#  Flux corrected transport solver benchmark
#  we are moving a Gaussian hill around 
#
#     we solve a* U_{,t} - b *u_{,ii} + c_i u_{,i} + (d_i* u)_{,i}}=0 
#
#               U(0)= U0 * exp ( - |x-x_0(t)|^2/(4*s**2) )
#
#  with a>0, b>=0, s>0
#
#  we set E=b/a v=c/a w=d/a
#
#  the solution is given as   u(x,t)=U0*s^dim/(s**2+E*t)^{dim/2} * exp ( - |x-x_0(t)|^2/(4*(s**2+E*t)) ) 
#
#        with x_0(t) = X0 + (v+w)*t 
#
#
#    the region |x-x_0(t)|^2/(4*(s**2+E*t)) < - log(TAU)  is within the domain for all time 
#
#
#    this holds if
#
#       |x_i-X0_i-(v_i+w_i)*tend | < sqrt(- log(TAU) * 4*(s**2+E*tend))=b0  and
#       |x_i-X0_i | < sqrt(- log(TAU)) * 2*s = b1 implies 0<=x_i<=l_i
#
from math import pi, ceil
from time import time as clock
from esys.finley import Rectangle, Brick
from esys.escript import *
from esys.escript.linearPDEs import LinearSinglePDE, TransportPDE
from esys.weipa import saveVTK

#  
DIM=2
NE_MAX=300000
VERBOSITY=True
TOL=1.e-8
TAU=1e-10
VTK_DIR="output"

#==================
S_MIN=0
TABS = [ 'dx', 'dt', 'peclet', 'error',  'sup', 'integral',  'center', 'width', 'time' ]

a=1.
#==================


mkDir(VTK_DIR)

def uRef(dom,t,E,s,v,x0, onElements=False):
    if onElements:
       x=Function(dom).getX()
    else:
      x=dom.getX()
    X=x0[:dom.getDim()]+v[:dom.getDim()]*t
    u=(s**2/(s**2+E*t))**(dom.getDim()/2.) * exp(-length(x-X)**2/(4*(s**2+E*t))) 
    return u


def getDirection(dim, d="x"):
     k=kronecker(dim)
     if d=="x":
         return k[0]
     elif d=="y":
         return k[1]
     elif d=="z" and dim>2:
         return k[2]
     elif d=="xy":
         return (k[0]+k[1])/sqrt(2.)
     elif d=="yz":
         return (k[1]+k[2])/sqrt(2.)
     elif d=="zx" and dim>2:
         return (k[2]+k[0])/sqrt(2.)
     elif d=="xyz" and dim>2:
         return (k[1]+k[2]+k[0])/sqrt(3.)
     else:
         raise ValueError("Cannot identify direction %s"%d)

def QUALITY(u_h,u_ref):
     u_h_e=interpolate(u_h,u_ref.getFunctionSpace())
     x=u_ref.getFunctionSpace().getX()
     out = {}
     out["error"]=sqrt(integrate((u_ref-u_h_e)**2))/sqrt(integrate(u_ref**2))
     out["sup"]=abs(sup(u_h)-sup(u_ref))/abs(sup(u_ref))
     m0_h=integrate(u_h_e)
     m1_h=integrate(x*u_h_e)/m0_h
     m2_h=integrate(length(x-m1_h)**2*u_h_e)

     m0=integrate(u_ref)
     m1=integrate(x*u_ref)/m0_h
     m2=integrate(length(x-m1)**2*u_ref)

     out["m0"]=abs(m0_h-m0)/abs(m0)
     out["m1"]=length(m1_h-m1)/length(m1)
     out["m2"]=abs(m2_h-m2)/abs(m2)

     return out
 
#================

def XXX(dim,tend,dt, s, h,b,c,d,c_dir="x", d_dir="x", a=1., CN=True):
    """
     dim - sparial dimension 
     s - width of initial profile
     h - mesh size
    """
    v_c=c/a*getDirection(dim,c_dir)
    v_d=d/a*getDirection(dim,d_dir)
    v = (v_c+v_d)
    E=b/a 
    if VERBOSITY: 
           print("="*100)
           print("XX Start test dim  = %d , h=%e, b=%e, c=%e, d=%e, c_dir=%s, d_dir=%s, a=%e, s=%e"%(dim, h,b,c,d,c_dir, d_dir, a, s))
           print("="*100)
           print("initial width s = ",s)
           print("diffusion = ",E)
           print("total velocity = ",v)
           print("tend = ", tend)
           print("tolerance = ",TOL)
           print("number of elements over s =",s/h)
    b0=sqrt(- log(TAU) * 4*(s**2+E*tend))
    b1=sqrt(- log(TAU)) * 2*s 
    X0_0=max(b1,-v[0]*tend + b0)
    X0_1=max(b1,-v[1]*tend + b0)
    l_0=X0_0+max(v[0]*tend + b0 , b1)
    l_1=X0_1+max(v[1]*tend + b0 , b1)
    NE_0=max(int(l_0/h+0.5),1)
    NE_1=max(int(l_1/h+0.5),1)
    if dim == 2:
        if VERBOSITY: print("%d x %d grid over %e  x %e with element size %e."%(NE_0,NE_1,l_0,l_1,h))
        if NE_0*NE_1 > NE_MAX:
           raise ValueError("too many elements %s."%(NE_0*NE_1,))
        dom=Rectangle(n0=NE_0,n1=NE_1,l0=l_0,l1=l_1)
        x0=[X0_0, X0_1]
    else:
       X0_2=max(b1,-v[2]*tend + b0)
       l_2=X0_2+max(v[2]*tend + b0 , b1)
       NE_2=max(int(l_2/h+0.5),1)
       if VERBOSITY: print("%d x %d x %d grid over %e  x %e x %e with element size %e."%(NE_0,NE_1,NE_2,l_0,l_1,l_2,h))
       if NE_0*NE_1*NE_2 > NE_MAX:
          raise ValueError("too many elements %s."%(NE_0*NE_1*NE_2,))
       dom=Brick(n0=NE_0,n1=NE_1, ne2=NE_2, l0=l_0,l1=l_1, l2=l_2)
       x0=[X0_0, X0_1, X0_2]
    if VERBOSITY: 
        print("initial location = ",x0)
    print("XX", interpolate(uRef(dom,0.,E,s,v,x0), FunctionOnBoundary(dom)))
     
    fc_BE=TransportPDE(dom,numEquations=1)
    fc_BE.setValue(M=a, A=-b*kronecker(dom), B=-v_d*a, C=-v_c*a)
    fc_BE.getSolverOptions().setVerbosity(VERBOSITY)
    fc_BE.getSolverOptions().setTolerance(TOL)
    #
    fc_BE.getSolverOptions().setPreconditioner(fc_BE.getSolverOptions().GAUSS_SEIDEL)
    fc_BE.getSolverOptions().setNumSweeps(5)  
    if VERBOSITY: print("Backward Euler Transport created")

    fc_CN=TransportPDE(dom,numEquations=1)
    fc_CN.setValue(M=a, A=-b*kronecker(dom), B=-v_d*a, C=-v_c*a)
    fc_CN.getSolverOptions().setVerbosity(VERBOSITY)
    fc_CN.getSolverOptions().setTolerance(TOL)
   
    #fc_CN.getSolverOptions().setPreconditioner(fc_CN.getSolverOptions().GAUSS_SEIDEL) 
    fc_CN.getSolverOptions().setNumSweeps(2)  
    if VERBOSITY: print("Crank Nicolson Transport created")
    dt_CN=fc_CN.getSafeTimeStepSize()
    if VERBOSITY: print("time step size by Crank Nicolson=",dt_CN)

    U0=uRef(dom,0,E,s,v,x0)
    U0_e=uRef(dom,0,E,s,v,x0,True)
    fc_CN.setInitialSolution(U0)
    fc_BE.setInitialSolution(U0)
    initial_error_L2=sqrt(integrate((U0-U0_e)**2))
    if VERBOSITY:
      print("initial Lsup = ",Lsup(U0), Lsup(U0_e))
      print("initial integral = ",integrate(U0_e))
      print("initial error = ",initial_error_L2)
      print("used time step size =",dt) 
      
    if not CN:
       n=int(ceil(tend/dt))
       if VERBOSITY: 
          print("Solve Backward Euler:")
          print("substeps : ",n)
       t0=clock()
       for i in range(n): u=fc_BE.getSolution(dt)
       t0=clock()-t0
    else:
       if VERBOSITY: print("Solve Crank Nicolson:")
       dt=dt_CN
       t0=clock()
       u=fc_CN.getSolution(tend)
       t0=clock()-t0
    out=QUALITY(u,uRef(dom,tend,E,s,v,x0,True))
    print("XX", interpolate(uRef(dom,tend,E,s,v,x0), FunctionOnBoundary(dom)))
    out['time']=t0
    out['tend']=tend
    out['dt']=dt
    out['dx']=h
    if abs(b)>0:
       out["peclet"]=length(v)*s/b
    else:
        out["peclet"]=9999999.
    # saveVTK("bb.vtu",u0=U0,u_CN=u_CN, uRef=uRef(dom,dt2,E,s,v,X0) )
    return out

# (s, peclet, b, h0)  -> error < 0.01
test_set = ( (0.05, 1., 1., 0.024),  )
test_set = ( (0.05, 100000., 1., 0.024),  )

if False:
    S_MAX=0.5/sqrt(-log(TAU))/2
    s=0.05
    peclet = 1000.
    b=1.
    c=peclet*b/s
    h=0.1/4*1.2/1.25
    dt=6.250000e-10

    print(XXX(DIM,dt,dt,s=s,h=h,b=a*b,c=a*c,d=0,c_dir="x", d_dir="x", CN=True))
    1/0

for tst in test_set:
     s=tst[0]
     peclet=tst[1]
     b=tst[2]
     h0=tst[3]
     c=peclet*b/s

     # find appropraiate tend:
     result=XXX(DIM,1e-99,1.,s=s,h=h0,b=a*b,c=a*c,d=0,c_dir="x", d_dir="x", CN=True)
     tend=result['dt']

     f_test=[ 1 , 2, 4 ]
     f_test=[ 1, 2, 4, 8 ]
     out=""
     tab_name="dt"
     tab_name="tend"
     tab_name="error"
     dt_s=[]
     for f_h in f_test:
        out+="h0/%s "%f_h
        h=h0/f_h
        result=XXX(DIM,tend,tend,s=s,h=h,b=a*b,c=a*c,d=0,c_dir="x", d_dir="x", CN=True)
        out+=", %e"%result[tab_name]
        print("XX",result)
        dt_s.insert(0,result['dt'])
        for i in range(len(f_test)-len(dt_s)): out+=", "
        for dt in dt_s:
            result=XXX(DIM,tend,dt,s=s,h=h,b=a*b,c=a*c,d=0,c_dir="x", d_dir="x", CN=False)
            print("XX",result)
            out+=", %e"%result[tab_name]
        out+="\n"
     header="h\dt , "
     for dt in dt_s: header+=", %e"%dt
     out=header+"\n"+out
     print(out)
