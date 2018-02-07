
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
#   AXI-SYMMETRIC NEWTONIAN MODEL ; UPDATED LAGRANGIAN FORMULATION
#
#
#    step 1 rho*(v_star-v) = dt * (sigma'_ij,j-teta3*p,i+f_i)
#    step 2 dp=-dt*B*(v_j,j+teta1*v_star_j,j-dt*teta1*((1-teta3)*p_,jj+teta2*dp_,jj))
#    step 3 rho*(v+-v) = -dt*((1-teta3)*p_,jj+teta2*dp_,jj)
#    step 3b p+=1/2(p+dp+abs(p+dp))
#    step 4 sigma'i+_ij,j=f(v+,p+,...)
#
#
from esys.escript import *
from esys.escript.linearPDEs import LinearSinglePDE, LinearPDESystem
from esys.finley import Rectangle
from esys.weipa import *


nel      =   20
H        =   0.5
L        =   1.0

eta      =   1.0       # shear viscosity
ro       =   1.0
g        =   1.00

alpha_w   =   1.00
alpha    =   1.00*1000000.
Pen=0.
B=100.

nstep    =   3000
dt       =   1.
small    =   EPSILON
w_step=max(int(nstep/50),1)*0+1
toler    =   0.001
teta1    =    0.5
teta2    =    0.5
teta3    =    1  # =0 split A; =1 split B

# create domain:
dom=Rectangle(int(nel*L/min(L,H)),int(nel*H/min(L,H)),order=1, l0=L, l1=H)
x=dom.getX()


momentumStep1=LinearPDESystem(dom) 
momentumStep1.setValue(q=whereZero(x[0])*[1.,0.]+whereZero(x[1])*[0.,1.]) # fix x0=0 and x1=0
face_mask=whereZero(FunctionOnBoundary(dom).getX()[1])

pressureStep2=LinearSinglePDE(dom) 
pressureStep2.setReducedOrderOn() 
pressureStep2.setValue(q=whereZero(x[0]-L)+whereZero(x[1]-H))

momentumStep3=LinearPDESystem(dom)
momentumStep3.setValue(q=whereZero(x[0])*[1.,0.]+whereZero(x[1])*[0.,1.])
#
#   initial values:
#
U=Vector(0.,Solution(dom)) 
p=ro*g*(L-ReducedSolution(dom).getX()[0])*(H-ReducedSolution(dom).getX()[1])/3 
p=ro*g*(H-ReducedSolution(dom).getX()[1])
dev_stress=Tensor(0.,Function(dom))

t=dt
istep=0
while istep < nstep:
    istep=istep+1
    print("time step :",istep," t = ",t)
    r=Function(dom).getX()[0]
    r_b=FunctionOnBoundary(dom).getX()[0]
    print("volume : ",integrate(r))
    #
    #  step 1:
    #
    # calculate normal 
    n_d=dom.getNormal()
    t_d=matrixmult(numpy.array([[0.,-1.],[1.,0]]),n_d)
    sigma_d=(sign(inner(t_d,U))*alpha_w*t_d-n_d)*Pen*clip(inner(n_d,U),0.)
    print("sigma_d =",inf(sigma_d),sup(sigma_d))

    momentumStep1.setValue(D=r*ro*kronecker(dom),
                           Y=r*ro*U+dt*r*[0.,-ro*g], 
                           X=-dt*r*(dev_stress-teta3*p*kronecker(dom)), 
                           y=sigma_d*face_mask*r_b)
    U_star=momentumStep1.getSolution()
    saveVTK("u.vtu",u=U_star,u0=U)
    #
    #  step 2:
    #
    # U2=U+teta1*(U_star-U)
    U2=U+teta1*U_star
    gg2=grad(U2)
    div_U2=gg2[0,0]+gg2[1,1]+U2[0]/r

    grad_p=grad(p)

    pressureStep2.setValue(A=r*dt*B*teta1*teta2/ro*dt*kronecker(dom), 
                           D=r,                            
                           Y=-dt*B*r*div_U2,
                           X=-r*B*dt**2/ro*teta1*(1-teta3)*grad_p)
    dp=pressureStep2.getSolution()
    #
    #  step 3:
    #
    p2=(1-teta3)*p+teta2*dp
    grad_p2=grad(p2)
    momentumStep3.setValue(D=r*ro*kronecker(dom),
                           Y=r*(ro*U_star-dt*teta2*grad_p2))
    U_new=momentumStep3.getSolution()
    #
    #   update:
    #
    p+=dp         
    U=U_new
    print("U:",inf(U),sup(U))
    print("P:",inf(p),sup(p)) 


    p_pos=clip(p,small)
    gg=grad(U) 
    vol=gg[0,0]+gg[1,1]+U[0]/r  
    gamma=sqrt(2*((gg[0,0]-vol/3)**2+(gg[1,1]-vol/3)**2+(U[0]/r-vol/3)**2+(gg[1,0]+gg[0,1])**2/2))
    m=whereNegative(eta*gamma-alpha*p_pos) 
    eta_d=m*eta+(1.-m)*alpha*p_pos/(gamma+small)  
    print("viscosity =",inf(eta_d),sup(eta_d)) 
    dev_stress=eta_d*(symmetric(gg)-2./3.*vol*kronecker(dom))
    #
    # step size control:
    #
    len=inf(dom.getSize())
    dt1=inf(dom.getSize()/(length(U)+small))
    dt2=inf(0.5*ro*(len**2)/eta_d)
    dt=dt1*dt2/(dt1+dt2)
    print("new step size = ",dt)
    #
    #  update geometry
    #
    dom.setX(dom.getX()+U*dt)
    t=t+dt
    if (istep-1)%w_step==0:saveVTK("u.%d.vtu"%((istep-1)/w_step),p=p,eta=eta_d,U=U_star,U_star=U_star,gamma=gamma)
    if istep == 3: 1/0
