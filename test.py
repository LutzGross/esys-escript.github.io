from esys.escript import *
from esys.escriptcore.splitworld import *

from esys.ripley import *


sw=SplitWorld(getMPISizeWorld()/2)
buildDomains(sw, Brick, 3,3,3)

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
      print self.domain.OnMasterProcessor()
      return self.rounds<1
      
for z in range(3):
#  addJob(sw, J1)
  addJob(sw, BarrierJob)

sw.runJobs()