# -*- coding: utf-8 -*-
########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

#
#  Flux corrected transport solver benchmark
#  we are moving a Gaussian hill around 
#
#     we solve a* U_{,t} - b *u_{,ii} + c_i u_{,i} + (d_i* u)_{,i}}=0 
#
#               U(0)= U0 * exp ( - |x-x_0(t)|^2/(4*s) )
#
#  with a>0, b>=0, s>0
#
#  we set E=b/a v=c/a w=d/a
#
#  the solution is given as   u(x,t)=U0*s^{dim/2}/(s+E*t)^{dim/2} * exp ( - |x-x_0(t)|^2/(4*(s+E*t)) ) 
#
#  			      with x_0(t) = X0 + (v+w)*t 
#
#
#    |x-x_0(t)|^2/(4*(s+E*t)) < - log(TAU)   for all time and all x in the
# domain where TAU is a given relative tolerance
#
#
#    this holds if
#
#    	|x_i-X0_i-(v_i+w_i)*tend | < - log(TAU) * 4*(s+E*tend)
#    	|x_i-X0_i | < - log(TAU) * 4*s
#
#        X0_i=min(- log(TAU) * 4*s,-(v_i+w_i)*tend - log(TAU) * 4*(s+E*tend))
#        l_i=X0_i+min((v_i+w_i)*tend - log(TAU) * 4*(s+E*tend)) , - log(TAU) * 4*s) 
#  
DIM=2
NE_MAX=35000
VERBOSITY=True
TOL=1.e-8
TAU=1e-10
VTK_DIR="output"
S_RAISE=1.5
S_MIN=0
S_MAX=0
TABS = [ 'dx', 'dt', 'peclet', 'error',  'sup', 'integral',  'center', 'width', 'time' ]
from math import pi, ceil
from time import time as clock
from esys.finley import Rectangle, Brick
from esys.escript import *
from esys.escript.linearPDEs import LinearSinglePDE, TransportPDE
S_MAX=-0.5/log(TAU)/4


mkDir(VTK_DIR)

def uRef(dom,t,E,s,v,x0, onElements=False):
    if onElements:
       x=Function(dom).getX()
    else:
      x=dom.getX()
    X=x0[:dom.getDim()]+v[:dom.getDim()]*t
    u=(s/(s+E*t))**(dom.getDim()/2.) * exp(-length(x-X)**2/(4*(s+E*t))) 
    return u

NE_test = [ 10, 20, 40, 80, 160 ]
NE_test = [ 10, 20, 40, 80, 160, 320]
NE_test = [ 320 *2 ]
b_test=[ 0., 0.001, 0.1, 1., 10, 100.] 
b_test=[ 1.0 *0] 
c_test=[ 0., 0.001, 0.1, 1., 10, 100. ]
c_test=[ 100.]
d_test=[ 0.0  ]
s_test=[0.0005 ]
# s=[ 0.001, 0.003, 0.01 ]

a=1.
X0=[0.5,0.5,0.5]
step=0.1
s_max=0.01  # don't use an s which is bigger than this as the profile will not fit into the domian 


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
         raise ValueError,"Cannot identify direction %s"%d

#================
def QUALITY(u_h,u_ref):
     dom=u_h.getDomain()
     x=dom.getX()
     out = {}
     out["sup"]=abs(sup(u_h)-sup(u_ref))/abs(sup(u_ref))
     i_h=integrate(u_h)
     i_ref=integrate(u_ref)
     out["integral"]=abs(i_h-i_ref)/abs(i_ref)
     x0h=integrate(x*u_h)/i_h
     x0ref=integrate(x*u_ref)/i_ref
     out["center"]=length(x0h-x0ref)/(1.+step)

     sigma_h=sqrt(integrate(length(x-x0h)**2*u_h))
     sigma_ref=sqrt(integrate(length(x-x0ref)**2*u_ref))
     out["width"]=abs(sigma_h-sigma_ref)/abs(sigma_ref)
     out["error"]=sqrt(integrate((u_ref-u_h)**2))
     
     return out
 


def XXX(dim, s, h,b,c,d,c_dir="x", d_dir="x", a=1.):
    v_c=c/a*getDirection(dim,c_dir)
    v_d=d/a*getDirection(dim,d_dir)
    v = (v_c+v_d)
    E=b/a 
    if abs(E)>0:
       if length(v)>0:
	  tend=min(s/E, s/length(v))
       else:
	  tend=s/E
    else:
	 if length(v)>0:	
	    tend=s/length(v)
	 else:
	    tend=0
    tend*=(1.-S_RAISE)
    if VERBOSITY: 
           print "="*100
           print "Start test dim  = %d , h=%e, b=%e, c=%e, d=%e, c_dir=%s, d_dir=%s, a=%e, s=%e"%(dim, h,b,c,d,c_dir, d_dir, a, s)
           print "="*100
           print "initial width s = ",s
           print "diffusion = ",E
           print "total velocity = ",v
           print "tend = ", tend
           print "tolerance = ",TOL
           print "XX number of elements over s =",s/h
    X0_0=min(- log(TAU) * 4*s,-v[0]*tend - log(TAU) * 4*(s+E*tend))
    X0_1=min(- log(TAU) * 4*s,-v[1]*tend - log(TAU) * 4*(s+E*tend))
    l_0=X0_0+min(v[0]*tend - log(TAU) * 4*(s+E*tend) , - log(TAU) * 4*s)
    l_1=X0_1+min(v[1]*tend - log(TAU) * 4*(s+E*tend) , - log(TAU) * 4*s)
    NE_0=max(int(l_0/h+0.5),1)
    NE_1=max(int(l_1/h+0.5),1)
    if dim == 2:
        if VERBOSITY: print "%d x %d grid over %e  x %e with element size %e."%(NE_0,NE_1,l_0,l_1,h)
        if NE_0*NE_1 > NE_MAX:
	   raise ValueError,"too many elements %s."%(NE_0*NE_1,)
        dom=Rectangle(n0=NE_0,n1=NE_1,l0=l_0,l1=l_1)
        x0=[X0_0, X0_1]
    else:
       X0_2=min(- log(TAU) * 4*s,-v[2]*tend - log(TAU) * 4*(s+E*tend))
       l_2=X0_2+min(v[2]*tend - log(TAU) * 4*(s+E*tend) , - log(TAU) * 4*s)
       NE_2=max(int(l_2/h+0.5),1)
       if VERBOSITY: print "%d x %d x %d grid over %e  x %e x %e with element size %e."%(NE_0,NE_1,NE_2,l_0,l_1,l_2,h)
       if NE_0*NE_1*NE_2 > NE_MAX:
	  raise ValueError,"too many elements %s."%(NE_0*NE_1*NE_2,)
       dom=Brick(n0=NE_0,n1=NE_1, ne2=NE_2, l0=l_0,l1=l_1, l2=l_2)
       x0=[X0_0, X0_1, X0_2]
    if VERBOSITY: 
        print "initial location = ",x0
     
    fc_BE=TransportPDE(dom,numEquations=1,useBackwardEuler=True)
    fc_BE.setValue(M=a, A=-b*kronecker(dom), B=-v_d*a, C=-v_c*a)
    fc_BE.getSolverOptions().setVerbosity(VERBOSITY)
    fc_BE.getSolverOptions().setTolerance(TOL)
    #
    # fc_BE.getSolverOptions().setPreconditioner(fc_BE.getSolverOptions().GAUSS_SEIDEL)
    fc_BE.getSolverOptions().setNumSweeps(1)  
    if VERBOSITY: print "Backward Euler Transport created"

    fc_CN=TransportPDE(dom,numEquations=1,useBackwardEuler=False)
    fc_CN.setValue(M=a, A=-b*kronecker(dom), B=-v_d*a, C=-v_c*a)
    fc_CN.getSolverOptions().setVerbosity(VERBOSITY)
    fc_CN.getSolverOptions().setTolerance(TOL)
   
    #fc_CN.getSolverOptions().setPreconditioner(fc_CN.getSolverOptions().GAUSS_SEIDEL) 
    fc_CN.getSolverOptions().setNumSweeps(2)  
    if VERBOSITY: print "Crank Nicolson Transport created"
    dt_CN=fc_CN.getSafeTimeStepSize()
    if VERBOSITY: print "time step size by Crank Nicolson=",dt_CN
    if tend>0:
       dt=main(dt_CN, tend)
    else:
       dt=1.
    U0=uRef(dom,0,E,s,v,X0)
    U0_e=uRef(dom,0,E,s,v,X0,True)
    fc_CN.setInitialSolution(U0)
    fc_BE.setInitialSolution(U0)
    initial_error_L2=sqrt(integrate((U0-U0_e)**2))
    if VERBOSITY:
      print "used time step size =",dt 
      print "initial integral = ",integrate(U0_e)
      print "initial error = ",initial_error_L2
      
    

    return 
    if False:
	# dt2=min(dt2,20*dt)
	if VERBOSITY: print "Solve Backward Euler:"
        t0=clock()
        u_BE=fc_BE.getSolution(10.*dt)
        u_BE=fc_BE.getSolution(10.*dt)
 	t0=clock()-t0
 	out=QUALITY(u_BE,uRef(dom,dt2,E,s,v,X0))
    #     
    # print integrate(U0), U0    
    # print integrate(u_BE), u_BE
    # print integrate(uRef(dom,dt2,E,s,v,X0)), uRef(dom,dt2,E,s,v,X0)
    # saveVTK("bb_be.vtu",u0=U0,u_BE=u_BE, uRef=uRef(dom,dt2,E,s,v,X0) )
    # 1/0
    # 
    else:
	if VERBOSITY: print "Solve Crank Nicolson:"
	t0=clock()
	u_CN=fc_CN.getSolution(dt)
	t0=clock()-t0
	out=QUALITY(u_CN,uRef(dom,dt,E,s,v,X0))
    out['time']=t0/NN
    out['dt']=dt
    out['dx']=1./NE
    if abs(b)>0:
       out["peclet"]=length(v)/(b*NE)
    else:
        out["peclet"]=9999999.
    # saveVTK("bb.vtu",u0=U0,u_CN=u_CN, uRef=uRef(dom,dt2,E,s,v,X0) )
    return out


for s_fac in [ 1., 0.5, 0.25, 0.125 ] :
   for h_fac in [1., 0.5, 0.25, 0.125 ]:
      XXX(DIM,s=S_MAX*s_fac,h=10*S_MAX*s_fac*h_fac,b=0,c=0,d=0,c_dir="x", d_dir="x")
1/0

txt=""
for t in TABS: txt+=" "+t
for b in b_test: 
  for c in c_test:
    for d in d_test:
      for s in s_test:
	for NE in NE_test:
	    dat=XXX(NE, a,b,c,d,c_dir="x", d_dir="x",s=s)
            txt+="\n" 
            for t in TABS: txt+=" %3.2e"%dat[t]
            
print txt
