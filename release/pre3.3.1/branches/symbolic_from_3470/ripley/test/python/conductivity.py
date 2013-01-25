  # Elastic deformation (see user's guide 1.4, pp. 28)
from esys.escript import *
from esys.finley import Brick, Rectangle

s_val=2.45
dom=Rectangle(n0=20, n1=10, l0=2, l1=1)

s=Symbol('s', (), dim=dom.getDim())
T = Symbol('T', (), dim=dom.getDim())
log_k = Symbol('log_k', (), dim=dom.getDim())

v = VariationalProblem(dom, u=T,p=log_k)
v.setValue(X=exp(-p)*grad(u), Y=s, h=0.5*(u-0.3)**2)
log_k, T, l = v.getSolution(T=0.,log_K=1, s=s_val)

#sT,S_log_K=v.getSensitivity(s, direction=1, T=T, log_K=log_K, s=2.45)