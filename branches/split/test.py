from esys.escript import *
from esys.escriptcore.splitworld import *

from esys.ripley import *


sw=SplitWorld(getMPISizeWorld())
buildDomains(sw, Brick, 3, 3, 3 )

class J1(Job):
    def __init__(self, **kwargs):
      super(J1,self).__init__(**kwargs)
      print "Constructed with "+str(self.domain)+" and id="+str(self.jobid)
      self.rounds=self.jobid+2
    def work(self):
      print str(self.jobid)+"is doing work  "+str(self.rounds)
      self.rounds-=1
      return self.rounds<1

class BarrierJob(Job):
    def __init__(self, **kwargs):
      super(BarrierJob,self).__init__(**kwargs)
      self.rounds=self.jobid+2
      
    def work(self):
      print str(self.jobid)+"is doing work  "+str(self.rounds)
      self.rounds-=1
      #print self.domain.OnMasterProcessor()
      print self.domain.onMasterProcessor()
      return self.rounds<1
      
      
class PoissonJob(Job):
    def __init__(self, **kwargs):
      super(PoissonJob, self).__init__(**kwargs)
    
    def work(self):
      print "Starting work"
      x = self.domain.getX()
      print "Here"
      q=x[0]
      print ";;;;;;;;;;"
      whereZero(q)
      print "----"
      gammaD = whereZero(x[0])+whereZero(x[1])
      print "There"
      # define PDE and get its solution u
      mypde = Poisson(domain=self.domain)
      print "Passed constructor"
      mypde.setValue(f=1,q=gammaD)
      print "Beginning solve"
      u = mypde.getSolution()
      print "Solution complete"

#addJob(sw, J1)
#addJob(sw, J1)
#for z in range(1):
#  addJob(sw, J1)
#  addJob(sw, BarrierJob)
addJob(sw, PoissonJob)
  
try:  
    sw.runJobs()
except RuntimeError as e:
    print e