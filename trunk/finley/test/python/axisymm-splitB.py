#
#   AXI-SYMMETRIC NEWTONIAN MODEL ; UPDATED LAGRANGIAN FORMULATION
#
#
#    step 1 rho*(v_star-v) = dt * (sigma'_ij,j-teta3*p,i+f_i)
#    step 2 dp=-dt*B*(v_j,j+teta1*v_star_j,j-dt*teta1*((1-teta3)*p_,jj+teta2*dp_,jj))
#    step 3 rho*(v+-v) = -dt*((1-teta3)*p_,jj+teta2*dp_,jj)
#    step 4 p+=1/2(p+dp+abs(p+dp))
#    step 4 sigma'i+_ij,j=f(v+,p+,...)
#
#
from esys.escript import *
from esys.escript.linearPDEs import LinearSinglePDE, LinearPDESystem
from esys.finley import Rectangle
iter     =   0
nstep    =   3000
w_step=max(int(nstep/50),1)*0+1
mstep    =   5
H        =   0.5
eta      =   1.0       # shear viscosity
ro       =   1.0
g        =   1.00
tauy0    =   0.9*ro*g*H/sqrt(3)       #0.11/(3*sqrt(3))
Pen      =   100
vtop     =   0.01
dt       =   10**(-10)
small    =   EPSILON
alpha    =   1.00
nel      =   20
# len      =   sqrt(0.5/nel)
toler    =   0.001
alphaw   =   1.00
L        =   1.0
teta1    =    0.5
teta2    =    0.5
teta3    =    0  # =0 split A; =1 split B
Etau=1000000000.

# create domain:
dom=Rectangle(int(nel*L/min(L,H)),int(nel*H/min(L,H)),order=1, l0=L, l1=H)
x=dom.getX()


momentumStep1=LinearPDESystem(dom) # A momentumStep1
momentumStep1.setValue(q=whereZero(x[0])*[1.,0.]) # +whereZero(x[1])*[1.,1.])
face_mask=whereZero(FunctionOnBoundary(dom).getX()[1])

                      # b   1 U1_2={0.0} 
                      # b   1 U1_1={0,0}
                                                 # b   4 U1_1={0.0}

pressureStep2=LinearSinglePDE(dom) # A [reduced] pressureStep2
pressureStep2.setReducedOrderOn() #
pressureStep2.setValue(q=whereZero(x[0]-L)+whereZero(x[1]-H))
                           # b  free U1={0.0} ?



momentumStep3=LinearPDESystem(dom) # A momentumStep3
momentumStep3.setValue(q=whereZero(x[0])*[1.,0.]) # +whereZero(x[1])*[1.,1.])
                      # b   1 U1_2={0.0} 
                      # b   1 U1_1={0,0}
                                                 # b   4 U1_1={0.0}


# A deformation
# e  V700= [grad] V100
# e  V901= V701+V704 +V101/(X1+small)
# e  V601=sqrt(2*((V701)^2+(V704)^2+(V101/(X1+small))^2+\
# (V702+V703)^2/2))
# e  U1={abs(V901)/(V601+small)}




# A stress
# e V700= [grad] V100
# e  V901= V701+V704 +V101/(X1+small)
# e  V601=sqrt(2*((V701-V901)^2+(V704-V901)^2+(V101/(X1+small)-V901)^2+\
# (V702+V703)^2/2))
# e  V602=((eta*V601<alpha*V401)*eta+(eta*V601>(alpha*V401-small))*(alpha*V401/(V601+small)))
# e V704=2*V602*V704-V401
# b 3 [integrated] {V704}
# 
# A surface
# b 3 [integrated] {1}


U=Vector(0.,Solution(dom)) # V100=0
p=ro*g*(L-ReducedSolution(dom).getX()[0])*(H-ReducedSolution(dom).getX()[1])/3 # V401=ro*g*(L-X1)*(H-X2)/3
t=dt
istep=0
En=0


while istep < nstep:
    istep=istep+1
    print "time step :",istep," t = ",t

    n_d=dom.getNormal()
    t_d=matrixmult(numarray.array([[0.,1.],[-1.,0]]),n_d)

    s=-En*inner(U,n_d)*face_mask
    nStressn=s*wherePositive(s)
    u_boundary=Etau*inner(U,t_d)*face_mask
    m=whereNegative(u_boundary-nStressn)
    nStressTau=nStressn*(1.-m)+u_boundary*m
    print nStressn
    # print nStressn*n_d+nStressTau*t_d
    if istep == 20:
        # print nStressTau
        saveVTK("stress.vtu",s=(nStressn*n_d+nStressTau*t_d))
        1/0
    # print nStressTau

    r=dom.getX()[0]
          # A viscosity
    gg=grad(U) # V700= [grad] V100
    vol=gg[0,0]+gg[1,1]+U[0]/(r+small)    # V901= V701+V704 +V101/(X1+small)
    tau=sqrt(2*((gg[0,0]-vol/3)**2+(gg[1,1]-vol/3)**2+(U[0]/(r+small)-vol/3)**2+(gg[1,0]+gg[0,1])**2/2))
              # V601=sqrt(2*((V701-V901/3)^2+(V704-V901/3)^2+(V101/(X1+small)-V901/3)^2+(V702+V703)^2/2))
    m=whereNegative(eta*tau-alpha*p) # eta*V601<alpha*V401
    eta_d=m*eta+(1.-m)*alpha*p/(tau+small) # 
               # U1={(eta*V601<alpha*V401)*eta+(eta*V601>(alpha*V401-small))*(alpha*V401/(V601+small))}
          # viscosity 1 100 V100 200 V200 300 V300 400 V400
          # solve
    print "	viscosity =",inf(eta_d),sup(eta_d) # show V801
    stress=eta_d*(symmetric(gg)-2./3.*vol*kronecker(dom))
               # solve momentum equationStep1
               # momentumStep1 1 100 V100 200 V200 300 V300 400 V400
    momentumStep1.setValue(D=ro/dt*kronecker(dom), # {ro/dt}U1_i
                           Y=stress[:,0]/(r+small)+[0.,-ro*g],           # -{0.0,ro*g}_i
                           X=-(stress+p*kronecker(dom)), # -D_j{V602}D_j{V100}_i-D_i{V602}D_j{V100}_j+D_i{(2/3)*V602}D_j{V100}_j-D_i{V401}
                           y=-(nStressn*n_d+nStressTau*t_d))
    dU2=momentumStep1.getSolution()
          # solve
    #                     dU2-> V100
    #                     U -> V200
    #                     p -> V500 

    # pressureStep2 1 100 V100 200 V200 300 V300 400 V500    
    
    gg2=grad(dU2)
    div_dU2=gg2[0,0]+gg2[1,1]+dU2[0]/(r+small)
    pressureStep2.setValue(A=r*teta1*teta2*dt**2*kronecker(dom),  # -D_j{teta1*teta2*dt^2}D_jU1
                           D=r*ro,                            # {ro/Pen}U1
                           Y=-r*(ro*dt*vol+teta1*ro*dt*div_dU2))  # -{ro*dt}D_j{V200}_j-{teta1*ro*dt}D_j{V100}_j
    # solve 
    dp=pressureStep2.getSolution()
                   # dp -> V100 
                   # dU2 -> V200
                   # U  -> V300
                   # p -> V600
    p+=dp                  # V601=V601+V101
                           # V701=V101 : dp -> V700
    p=(p+abs(p))/2+small # V601=cut V601
                         # V601=V601+small
    print "	pressure increment :",inf(dp),sup(dp) # show V601
    print "	pressure range :",inf(p),sup(p) # show V601
                # popp
                #   dU2 -> V100
                #   U  -> V200
                #   p -> V500
                #   dp -> V600 

    # momentumStep3 1 100 V100 200 V200 300 V300 400 V600
    A2=Tensor4(0.0,Function(dom)) 
    for i in range(dom.getDim()):
       for j in range(dom.getDim()):
          A2[i,i,j,j]+=1
    momentumStep3.setValue(A=Pen*A2,  # -D_i{Pen}D_jU1_j 
                           D=ro/dt*kronecker(dom), # {ro/dt}U1_i
                           Y=(ro/dt*(dU2+U)-teta2*grad(dp)))
                           # Y=r*(ro/dt*(dU2+U)-teta2*dp*[1.,0]),       # {ro/dt}{V200}_i+{ro/dt}{V100}_i
                           # X=r*teta2*dp**kronecker(dom))            # -D_i{teta2*V401}
    U, U_old=momentumStep3.getSolution(), U
                #   U -> V100
                #   dU2 -> V200
                #   U_old  -> V300
                #   p -> V600
                #   dp -> V700 
           # V400=V300
           # V300=V100
           # popp
           # popp
                #   U  -> V100
                #   U_old  -> V200
                #   p -> V400
                #   dp -> V500 
           # deformation
           # solve
           # show ratio
           # show V101
           # popp

    print "	new U:",inf(U),sup(U) # show v100
    print "	old U",inf(U_old),sup(U_old) # show v200
    print "	relative change:",Lsup(U_old)/Lsup(U) # test=L1 V200
                                                    # test0=L1 V100
                                                    # show test/test0
                     # black
                     # arrow
      
    vmax=Lsup(U)  # vmax1=abs(max V100)
                           # vmax2=abs(min V100)
                           # show vmax1
                           # show vmax2
       
                           # vmax=vmax1
                           # if vmax < vmax2
                           # vmax=vmax2
                           # endif
    print "	max velocity = ",vmax # show vmax
    len=inf(dom.getSize())
    dt1=len/vmax  # Courant condition
    dt2=0.5*ro*(len**2)/eta
    dt=dt1*dt2/(dt1+dt2)
    En=1/dt
    t=t+dt
    print "	time , time step ",t,dt # show t
                                        # show dt
        

    dom.setX(dom.getX()+U*dt) # mapm 500
                            # V500=V500+V100*dt
                            # V600=V500
                            # mapm 500

    print "	volume : ",integrate(r)# vol=Out
                                          # show vol
                #   dU  -> V100
                #   U  -> V200
                #   p -> V400
                #   dp -> V600 
               # clear
               # prim
    # saveVTK("u.%d.xml"%(istep),p=p,eta=eta_d,U=U,dU2=dU2,tau=tau,dp=dp)
    if (istep-1)%w_step==0:saveVTK("u.%d.xml"%((istep-1)/w_step),p=p,eta=eta_d,U=U,dU2=dU2,tau=tau,dp=dp)
