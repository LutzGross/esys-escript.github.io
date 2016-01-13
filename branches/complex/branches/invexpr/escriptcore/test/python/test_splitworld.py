
##############################################################################
#
# Copyright (c) 2015 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

__copyright__="""Copyright (c) 2015 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Tests to ensure that splitworld features operate correctly
"""

import esys.escriptcore.utestselect as unittest
from esys.escript import *

from esys.escriptcore.splitworld import *
from esys.escript.linearPDEs import Poisson, Helmholtz

class Test_SplitWorld(unittest.TestCase):
  """
  Class to test splitworld functions.
  Requires subclasses to supply self.domainpars which is a list of constructor function followed
  by arguments.
  eg:  if your domain is created with Rectangle(3,4), then your domainpars would be [Rectangle,3,4]
  """
  
  
  class PoissonJob1(Job):
    def __init__(self, **kwargs):
      super(PoissonJob1, self).__init__(**kwargs)
    
    def work(self):
      x = self.domain.getX()
      gammaD = whereZero(x[0])+whereZero(x[1])
      # define PDE and get its solution u
      mypde = Poisson(domain=self.domain)
      mypde.setValue(f=1, q=gammaD)
      u = mypde.getSolution()
      self.exportValue("answer", Lsup(u))
      return True

  class PoissonJob(Job):
    def __init__(self, **kwargs):
      super(Test_SplitWorld.PoissonJob, self).__init__(**kwargs)
    
    def work(self):
      x = self.domain.getX()
      gammaD = whereZero(x[0])+whereZero(x[1])
      # define PDE and get its solution u
      mypde = Poisson(domain=self.domain)
      mypde.setValue(f=self.jobid, q=gammaD)
      u = Lsup(mypde.getSolution())	   # we won't actually export the value to make
      self.exportValue("answer", self.jobid) # testing easier
      return True
      
  class HelmholtzJob(Job):
    def __init__(self, **kwargs):
      super(Test_SplitWorld.HelmholtzJob, self).__init__(**kwargs)
    
    def work(self):
      # define PDE and get its solution u
      mypde = Helmholtz(domain=self.domain)
      mypde.setValue(omega=self.jobid)
      u = mypde.getSolution()
      self.exportValue("hanswer", 2*self.jobid)
      self.exportValue("v", self.jobid)
      return True
      
  class InjectJob(Job):
    """
    Tests jobs taking parameters
    """
    def __init__(self, **kwargs):
      super(Test_SplitWorld.InjectJob, self).__init__(**kwargs)
      self.value=kwargs['val']
      self.name=kwargs['name']
    
    def work(self):
      """
      Make use of values passed to constructor
      """
      self.exportValue(self.name, self.value)
      return True      

  class FactorJob(Job):
    """
    Tests jobs taking parameters
    """
    def __init__(self, **kwargs):
      super(Test_SplitWorld.FactorJob, self).__init__(**kwargs)
      self.divisor=kwargs['fact']
    
    def work(self):
      """
      Make use of values passed to constructor
      """
      z=self.importValue("value")
      if (z%self.divisor==0):
          self.exportValue("boolean", self.divisor)
      return True 
      
      
  class ThrowJob(Job):
    """
    Trigger various faults
    """
    def __init__(self, **kwargs):
      super(Test_SplitWorld.ThrowJob, self).__init__(**kwargs)
      self.faultnum=kwargs['fault']             #leaving this out will test error in constructor
      
    def work(self):
      if self.faultnum==1:
        return "zero"           # non-boolean return value
      if self.faultnum==2:
        z=self.importValue("missing")   # unannounced import
        return True
      if self.faultnum==3:
        self.exportValue("missing",0)   # undeclared
        return True
      if self.faultnum==4:
        self.exportValue("answer","answer")  # type-mismatch in export
      return True

  class DummyJob(Job):
    """
    Trigger various faults
    """
    def __init__(self, **kwargs):
      super(Test_SplitWorld.DummyJob, self).__init__(**kwargs)
      
    def work(self):
      return True      
      
      
  def test_faults(self):
      for x in range(1,5):
        sw=SplitWorld(getMPISizeWorld())
        buildDomains(sw,*self.domainpars)
        addVariable(sw, "answer", makeScalarReducer, "MAX") 
        addJob(sw, Test_SplitWorld.ThrowJob, fault=x)
        self.assertRaises(RuntimeError, sw.runJobs)

  @unittest.skipIf(getMPISizeWorld()>97, "Too many ranks for this test")
  def test_factorjobs(self):
      """
      test importing, multiple phases, max as a flag
      """
      sw=SplitWorld(getMPISizeWorld())
      buildDomains(sw,*self.domainpars)
      addVariable(sw, "value", makeScalarReducer, "MAX")
      addVariable(sw, "boolean", makeScalarReducer, "MAX")
         # first we will load in a value to factorise
         # Don't run this test with 99 or more processes
      addJob(sw, Test_SplitWorld.InjectJob, name='value', val=101)      # Feed it a prime  
      addJob(sw, Test_SplitWorld.InjectJob, name='boolean', val=0)              # so we have a value
      sw.runJobs()
      for x in range(2,getMPISizeWorld()+2):
        addJob(sw, Test_SplitWorld.FactorJob, fact=x)
      sw.runJobs()
      self.assertEquals(sw.getDoubleVariable('boolean'),0)
      sw.clearVariable('value')
      sw.clearVariable('boolean')
      addJob(sw, Test_SplitWorld.InjectJob, name='value', val=101)      # Feed it a prime  
      addJob(sw, Test_SplitWorld.InjectJob, name='boolean', val=0)              # so we have a value
      sw.runJobs()
      sw.clearVariable("value")
      
        # Now test with a value which has a factor
      addJob(sw, Test_SplitWorld.InjectJob, name='value', val=100)       # Feed it a prime  
      addJob(sw, Test_SplitWorld.InjectJob, name='boolean', val=0)               # so we have a value
      sw.runJobs()
      m=0
      for x in range(2,getMPISizeWorld()+2):
        addJob(sw, Test_SplitWorld.FactorJob, fact=x)
        if 100%x==0:
          m=x
      sw.runJobs()
      self.assertEquals(sw.getDoubleVariable('boolean'),m)      
      
  def test_split_simple_solve(self):
    """
    Solve a single equation
    """
    sw=SplitWorld(getMPISizeWorld())
    buildDomains(sw,*self.domainpars)
    addVariable(sw, "answer", makeScalarReducer, "SUM")
    addJob(sw, Test_SplitWorld.PoissonJob)
    sw.runJobs()
    self.assertEquals(sw.getDoubleVariable("answer"),1)
    
  def test_split_simple_solve_multiple(self):
    """
    Solve a number of the same equation in one batch
    """
    sw=SplitWorld(getMPISizeWorld())
    buildDomains(sw,*self.domainpars)
    addVariable(sw, "answer", makeScalarReducer, "SUM")
        # this gives us 1 job per world
    total=0
    jobid=1
    for x in range(0,getMPISizeWorld()):
        addJob(sw, Test_SplitWorld.PoissonJob)
        total+=jobid
        jobid+=1
    sw.runJobs()
    self.assertEquals(sw.getDoubleVariable("answer"), total)
    
  def test_split_simple_and_dummy(self):
    """
    Solve a number of the same equation with some worlds doing dummy Jobs
    """
    sw=SplitWorld(getMPISizeWorld())
    buildDomains(sw,*self.domainpars)
    addVariable(sw, "answer", makeScalarReducer, "SUM")
        # this gives us 1 job per world
    total=0
    mid=getMPISizeWorld()//2
    if getMPISizeWorld()%2==1:
      mid=mid+1
    for x in range(0,mid):
        addJob(sw, Test_SplitWorld.PoissonJob)
        total=total+(x+1)
    for x in range(0,mid):
        addJob(sw, Test_SplitWorld.DummyJob)
    sw.runJobs()
      # expecting this to fail until I work out the answer
    self.assertEqual(sw.getDoubleVariable("answer"), total)
    
  def test_split_simple_and_empty(self):
    """
    Solve a number of the same equation with some worlds doing nothing
    """
    sw=SplitWorld(getMPISizeWorld())
    buildDomains(sw, *self.domainpars)
    addVariable(sw, "answer", makeScalarReducer, "SUM")
        # this gives us at most 1 job per world
    total=0    
    mid=getMPISizeWorld()//2
    if getMPISizeWorld()%2==1:
      mid=mid+1    
    for x in range(0,mid):
        addJob(sw, Test_SplitWorld.PoissonJob)
        total=total+(x+1)
    sw.runJobs()
      # expecting this to fail until I work out the answer
    self.assertEquals(sw.getDoubleVariable("answer"),total)    
    
    
  def test_split_multiple_batches(self):
    """
    Solve a number of the same equation in multiple batches
    """
    sw=SplitWorld(getMPISizeWorld())
    buildDomains(sw,*self.domainpars)
    addVariable(sw, "answer", makeScalarReducer, "SUM")
        # this gives us 1 job per world
    total=0
    sw.runJobs()
    for x in range(0,getMPISizeWorld()):
        addJob(sw, Test_SplitWorld.PoissonJob)
        total=total+x
    sw.runJobs()
    sw.runJobs()
    sw.clearVariable("answer")
    total=0
    for x in range(0,getMPISizeWorld()):
        addJob(sw, Test_SplitWorld.PoissonJob)
        total=total+(x+1+getMPISizeWorld())
    sw.runJobs()
      # expecting this to fail until I work out the answer
    self.assertEquals(sw.getDoubleVariable("answer"),total)    
  
  @unittest.skipIf(getMPISizeWorld()%2!=0, "Test requires even number of processes")
  def test_multiple_equations_size2world(self):
    """
    Test two different types of equations to solve mixed in batches
    We will try various combinations to spread them out over the 
    worlds in different patterns.
    This version attempts this with worlds of size 2
    """
    wc=getMPISizeWorld()//2
    sw=SplitWorld(wc)
    buildDomains(sw, *self.domainpars)
    addVariable(sw, "answer", makeScalarReducer, "SUM")   
    addVariable(sw, "hanswer", makeScalarReducer, "SUM")  
    addVariable(sw, "v", makeScalarReducer, "MAX")
    
    tot=0
    jobid=1
       #first put jobs of the same type close.
    for x in range(0, max(wc//3,1)):
      addJob(sw, Test_SplitWorld.PoissonJob)
      jobid+=1
    for x in range(0, max(wc//3,1)):
      addJob(sw, Test_SplitWorld.HelmholtzJob)
      tot+=2*(jobid)
      jobid+=1
    for x in range(0, max(wc//3,1)):
      addJob(sw, Test_SplitWorld.DummyJob)
      jobid+=1
    sw.runJobs()
    ha=sw.getDoubleVariable("hanswer")
    self.assertEquals(ha, tot)
    sw.clearVariable("answer")
    sw.clearVariable("hanswer")
    sw.clearVariable("v")
    tot=0
      # similar but separated by dummy Jobs
    for x in range(0, max(wc//3,1)):
      addJob(sw, Test_SplitWorld.PoissonJob)
      jobid+=1
    for x in range(0, max(wc//3,1)):
      addJob(sw, Test_SplitWorld.DummyJob)      
      jobid+=1
    for x in range(0, max(wc//3,1)):
      addJob(sw, Test_SplitWorld.HelmholtzJob)
      tot+=2*jobid
      jobid+=1
    sw.runJobs()
    ha=sw.getDoubleVariable("hanswer")
    self.assertEquals(ha, tot)
    sw.clearVariable("answer")
    sw.clearVariable("hanswer")
    sw.clearVariable("v")   
      # mixed
    tot=0
    for x in range(0, max(wc//2,1)):
      addJob(sw, Test_SplitWorld.HelmholtzJob)
      tot+=2*jobid
      addJob(sw, Test_SplitWorld.DummyJob)
      jobid+=2
      addJob(sw, Test_SplitWorld.PoissonJob)
      jobid+=1
    sw.runJobs()
    ha=sw.getDoubleVariable("hanswer")
    self.assertEquals(ha, tot)    

  @unittest.skipIf(getMPISizeWorld()%4!=0, "Test requires number of processes divisible by 4")
  def test_multiple_equations_size4world(self):
    """
    Test two different types of equations to solve mixed in batches
    We will try various combinations to spread them out over the 
    worlds in different patterns.
    This version attempts this with worlds of size 2
    """
    wc=getMPISizeWorld()//4
    sw=SplitWorld(wc)
    buildDomains(sw,*self.domainpars)
    addVariable(sw, "answer", makeScalarReducer, "SUM")   
    addVariable(sw, "hanswer", makeScalarReducer, "SUM")  
    addVariable(sw, "v", makeScalarReducer, "MAX")
    
    jobid=1
    tot=0
       #first put jobs of the same type close.
    for x in range(0, max(wc//2,1)):
      addJob(sw, Test_SplitWorld.PoissonJob)
      jobid+=1
    for x in range(0, max(wc//2,1)):
      addJob(sw, Test_SplitWorld.HelmholtzJob)
      tot+=2*jobid
      jobid+=1      
    for x in range(0, max(wc//2,1)):
      addJob(sw, Test_SplitWorld.DummyJob)
      jobid+=1
    sw.runJobs()
    ha=sw.getDoubleVariable("hanswer")
    self.assertEquals(ha, tot)
    sw.clearVariable("answer")
    sw.clearVariable("hanswer")
    sw.clearVariable("v")
    tot=0
      # similar but separated by dummy Jobs
    for x in range(0, max(wc//2,1)):
      addJob(sw, Test_SplitWorld.PoissonJob)
      jobid+=1      
    for x in range(0, max(wc//2,1)):
      addJob(sw, Test_SplitWorld.DummyJob) 
      jobid+=1     
    for x in range(0, max(wc//2,1)):
      addJob(sw, Test_SplitWorld.HelmholtzJob)
      tot+=2*jobid
      jobid+=1
    sw.runJobs()
    ha=sw.getDoubleVariable("hanswer")
    self.assertEquals(ha, tot)
    sw.clearVariable("answer")
    sw.clearVariable("hanswer")
    sw.clearVariable("v")   
    tot=0
      # mixed
    for x in range(0, max(wc//2,1)):
      addJob(sw, Test_SplitWorld.HelmholtzJob)
      tot+=2*jobid
      jobid+=1
      addJob(sw, Test_SplitWorld.DummyJob)
      addJob(sw, Test_SplitWorld.PoissonJob)
      jobid+=2
    sw.runJobs()
    ha=sw.getDoubleVariable("hanswer")
    self.assertEquals(ha, tot)        
    
    
  def test_multiple_equations_smallworld(self):
    """
    Test two different types of equations to solve mixed in batches
    We will try various combinations to spread them out over the 
    worlds in different patterns
    """
    sw=SplitWorld(getMPISizeWorld())
    buildDomains(sw,*self.domainpars)
    addVariable(sw, "answer", makeScalarReducer, "SUM")   
    addVariable(sw, "hanswer", makeScalarReducer, "SUM")  
    addVariable(sw, "v", makeScalarReducer, "MAX")
    
    tot=0
    jobid=1
       #first put jobs of the same type close together.
    for x in range(0,max(getMPISizeWorld()//3,1)):
      addJob(sw, Test_SplitWorld.PoissonJob)
      jobid+=1
    for x in range(0,max(getMPISizeWorld()//3,1)):
      addJob(sw, Test_SplitWorld.HelmholtzJob)
      tot+=2*jobid
      jobid+=1
    for x in range(0,getMPISizeWorld()//3):
      addJob(sw, Test_SplitWorld.DummyJob)
      jobid+=1
    sw.runJobs()
    ha=sw.getDoubleVariable("hanswer")
    self.assertEquals(ha, tot)
    sw.clearVariable("answer")
    sw.clearVariable("hanswer")
    sw.clearVariable("v")
    tot=0
      # similar but separated by dummy Jobs
    for x in range(0,max(getMPISizeWorld()//3,1)):
      addJob(sw, Test_SplitWorld.PoissonJob)
      jobid+=1
    for x in range(0,getMPISizeWorld()//3):
      addJob(sw, Test_SplitWorld.DummyJob)      
      jobid+=1      
    for x in range(0,max(getMPISizeWorld()//3,1)):
      addJob(sw, Test_SplitWorld.HelmholtzJob)
      tot+=2*jobid
      jobid+=1      
    sw.runJobs()
    ha=sw.getDoubleVariable("hanswer")
    self.assertEquals(ha, tot)
    sw.clearVariable("answer")
    sw.clearVariable("hanswer")
    sw.clearVariable("v")   
    tot=0
      # mixed
    for x in range(0, max(getMPISizeWorld()//2,1)):
      addJob(sw, Test_SplitWorld.HelmholtzJob)
      tot+=jobid*2
      addJob(sw, Test_SplitWorld.DummyJob)
      addJob(sw, Test_SplitWorld.PoissonJob)
      jobid+=3
    sw.runJobs()
    ha=sw.getDoubleVariable("hanswer")
    self.assertEquals(ha, tot)     
