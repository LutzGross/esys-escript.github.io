# Elastic deformation (see user's guide 1.4, pp. 28)
from esys.escript import *
from esys.finley import Rectangle

dom=Rectangle(l0=1., l1=1.,n0=100, n1=100)
x=dom.getX()
mask=whereZero(x[0])

alpha, u_b, x0=symbols("alpha, u_b, x0")
u=Symbol('u', (), dim=dom.getDim())
pde=NonlinearPDE(dom, u, debug=NonlinearPDE.DEBUG3)

pde.setValue(X=grad(u), Y=alpha*u**2-alpha*x0**2,  y= - whereZero(FunctionOnBoundary(dom).getX()[0]-1), q=mask, r=u_b)


ui=pde.getSolution(u=x[0]*14, alpha=100, u_b=x[0], x0=x[0])
print "error =",Lsup(ui-x[0])
