from esys.escript import *
from esys.weipa import saveVTK, saveSilo
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.finley import ReadGmsh
from esys.escript.pdetools import Locator

print("read in mesh")   
domain=ReadGmsh("simplemesh.msh", 3,  optimize=True )

pde = LinearSinglePDE(domain, isComplex=False)
pde.setSymmetryOn()
x = domain.getX()

pde.setValue(A=kronecker(3), Y=1, q=whereZero(x[0]-inf(x[0])))

options = pde.getSolverOptions()
options.setPackage(SolverOptions.TRILINOS)
options.setSolverMethod(SolverOptions.PCG)
options.setPreconditioner(SolverOptions.AMG)
options.setTrilinosParameter("multigrid algorithm", "sa")
options.setTrilinosParameter("sa: damping factor", 1.3)
options.setTrilinosParameter("max levels", 10)         
options.setTrilinosParameter("coarse: max size", 2000)
options.setTrilinosParameter("coarse: type", "SuperLU")
options.setTrilinosParameter("verbosity", "low") 

print("solve pde")
u=pde.getSolution() 

saveSilo("output/asimple",u=u)
print("simpletest.py finished.")
