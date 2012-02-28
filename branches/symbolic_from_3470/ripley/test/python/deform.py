# Elastic deformation (see user's guide 1.4, pp. 28)
from esys.escript import *
from esys.finley import Brick

dom=Brick(l0=1., l1=1., l2=1., n0=10, n1=10, n2=10)
x=dom.getX()
q=whereZero(x[0])*[1,0,0]+whereZero(x[1])*[0,1,0]+whereZero(x[2])*[0,0,1]
xc=[0.3, 0.3, 1.]
beta=8.
T_0=1.
params={
    'lam'  : 1.,
    'mu'   : 0.1,
    'alpha': 1.e-6,
    'T_ref': 0.,
    'T'    : T_0*exp(-beta*length(x-xc))
}
lam,mu,alpha,T,T_ref=symbols("lam,mu,alpha,T,T_ref")
u=Symbol('u', (dom.getDim(),), dim=dom.getDim())
gu=grad(u)
I=kronecker(dom)
sigma0=(lam+2./3.*mu)*alpha*(T-T_ref)*I
X=lam*trace(gu)*I + mu*(gu + transpose(gu))- sigma0
pde=NonlinearPDE(dom, u, debug=NonlinearPDE.DEBUG2)
pde.setValue(X=X, q=q)
u0=Data(0., u.getShape(), ContinuousFunction(dom))
ui=pde.getSolution(u=u0, **params)

#=== and the stamdart way:
from esys.escript.linearPDEs import LinearPDE
lam=params['lam']
mu=params['mu']
alpha=params['alpha']
T_ref=params['T_ref']
T=params['T']
A=Tensor4(0., Function(dom))
for i in range(dom.getDim()):
  for j in range(dom.getDim()):
     A[i,i,j,j]+=lam
     A[j,i,j,i]+=mu
     A[j,i,i,j]+=mu

sigma0=(lam+2./3.*mu)*alpha*(T-T_ref)*I
p=LinearPDE(dom)
p.setValue(A=A, X=sigma0, q=q)
u_ref=p.getSolution()

print Lsup(u_ref-ui)

#=== and now as a variational problem
v=VariationalProblem(dom, u=u)
sigma=lam*trace(gu)*I + mu*(gu + transpose(gu))
v.setValue(H=0.5*inner(symmetric(gu),sigma-sigma0), q= q)
ui2=pde.getSolution(u=u0, **params)
print Lsup(u_ref-ui2)