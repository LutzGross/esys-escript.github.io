from esys.escript import *
from esys.escript.linearPDEs import Poisson
from esys.ripley import Rectangle 

# Each node is in its own world
sw=SplitWorld(getMPISizeWorld())
buildDomains(sw, Rectangle, 100, 100)

#describe the work we want to do
# In this case we solve a Poisson equation
def task(self, **kwargs):
    v=kwargs['v']
    dom=self.domain
    pde=Poisson(dom)
    x=dom.getX()
    gammaD=whereZero(x[0])+whereZero(x[1])
    pde.setValue(f=v, q=gammaD)
    soln=pde.getSolution()
    soln.dump('soln%d.ncdf'%v)

# Now we add some jobs
for i in range(1,20):
    addJob(sw, FunctionJob, task, v=i)
# Run them
sw.runJobs() 

