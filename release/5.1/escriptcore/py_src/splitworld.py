##############################################################################
#
# Copyright (c) 2014-2017 by The University of Queensland
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

__copyright__="""Copyright (c) 2014-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

__author__="Joel Fenwick"

"""
This module contains the Python side of the SplitWorld functionality.
"""

import warnings
warnings.simplefilter('default', category=DeprecationWarning)

from . import escriptcpp as escore


class SplitWorld(object):
  """
  Wrapper for the C++ class exposed as __SplitWorld.
  This is a namespace consideration, it allows us to make 
  boost::python::raw_functions into members of a class.
  """
  
  def __init__(self, count):
    """
    :var count: How many equally sized subworlds should our compute resources be partitioned into?
    :var type: `int`
    """
    self.cpp_obj=escore.Internal_SplitWorld(count)
    
  def buildDomains(self, fn, *vec, **kwargs):
    """
    Instruct subworlds how to build the domain.
    :var fn: The function/class to call to create a domain.
    :type fn: `callable`
    The remaining parameters are for the arguments of the function.
    """
    escore.internal_buildDomains(self.cpp_obj, fn, *vec, **kwargs)

    
  def addJob(self, jobctr, *vec, **kwargs):
    """
    Submit a job to be run later on an available subworld.
    :var jobctr: class or function to be called to create a job
    :type jobctr: `callable`
    The remaining parameters are for the arguments of the function.
    """
    escore.internal_addJob(self.cpp_obj, jobctr, *vec, **kwargs)
    
  def addJobPerWorld(self, jobctr, *vec, **kwargs):
    """
    Submit one job per subworld to run later.
    :var jobctr: class or function to be called to create a job
    :type jobctr: `callable`
    The remaining parameters are for the arguments of the function.    
    """
    escore.internal_addJobPerWorld(self.cpp_obj, jobctr,  *vec, **kwargs)
    
    
  def addVariable(self, vname, vartype, *vec, **kwargs):
    """
    Create a variable on all subworlds.
    :var vartype: the type of variable to be created
    :type vartype: `str`
    The remaining parameters are for optional arguments depending on the variable type.
    """
    if vartype=="local":
        escore.internal_addVariable(self.cpp_obj, vname, escore.internal_makeLocalOnly);
    elif vartype=="Data":
        escore.internal_addVariable(self.cpp_obj, vname, escore.internal_makeDataReducer, *vec)
    elif vartype=="float":
        escore.internal_addVariable(self.cpp_obj, vname, escore.internal_makeScalarReducer, *vec)
    else:
        raise ValueError("Unknown variable type (%s)"%str(vartype))

  def runJobs(self):
    """
    Executes pending jobs.
    """
    self.cpp_obj.runJobs()

  def removeVariable(self, name):
    """
    Removes the named variable from all subworlds.
    :var name: 
    :type name: `str`
    """
    self.cpp_obj.removeVariable(name)
    
    
  def clearVariable(self, name):
    """
    Clears the value of the named variable.  The variable itself still exists.
    :var name: variable to clear
    :type name: `str`
    """
    self.cpp_obj.clearVariable(name)
    
  def getVarList(self):
    """
    Returns the names of all declared variables and a boolean for each indicating whether they have values.
    """
    return self.cpp_obj.getVarList()
  
  
  def getVarInfo(self):
    """
    Returns the names of all declared variables and a description of type.
    The details of the output are not fixed and may change without notice
    """
    return self.cpp_obj.getVarInfo()
    
    
  def getFloatVariable(self, vname):
    """
    Return the value of a floating point variable
    """
    return self.cpp_obj.getDoubleVariable(vname)
  
  
  def getLocalObjectVariable(self, vname):
    """
    Return the value of a local object variable - that is, an object (eg tuple) which does not need to
    be reduced/shared between worlds
    """
    return self.cpp_obj.getLocalObjectVariable(vname)
    
  def getSubWorldCount(self):
    """
    Return the number of subworlds in this splitworld
    """
    return self.cpp_obj.getSubWorldCount()

  def getSubWorldID(self):
    """
    Return the id of the subworld which _this_ MPI process belongs to.
    """
    return self.cpp_obj.getSubWorldID()
    
    
  def copyVariable(self, src, dest):
    """
    copy the contents of one splitworld variable into another
    :var src: name of variable to copy from
    :type src: `str`
    :var dest: name of variable to copy to 
    :type dest: `str`
    """
    self.cpp_obj.copyVariable(src, dest)
    
class Job(object):
  """
  Describes a sequence of work to be carried out in a subworld.
  The instances of this class used in the subworlds will
  be constructed by the system.
  The majority of the work done by the Job will be in the 
  *overloaded* work() method.
  To do specific work, this class should be subclassed and the work() 
  (and possibly __init__ methods overloaded).
  The majority of the work done by the job will be in the *overloaded* work() method.
  The work() method should retrieve values from the outside using importValue() and pass values to
  the rest of the system using exportValue().
  The rest of the methods should be considered off limits.
  """

  def __init__(self, *args, **kwargs):
    """
    It ignores all of its parameters, except, it requires the following as keyword arguments
    :var domain: Domain to be used as the basis for all ``Data`` and PDEs in this Job.
    :var jobid: sequence number of this job. The first job has id=1
    :type jobid: Positive ``int``
    """
    self.domain=kwargs["domain"]
    self.jobid=kwargs["jobid"]
    self.wantedvalues=[]                # names of shared values this job wishes to import    
    self.importedvalues={}      # name:values of which this jobs wants to use
    self.exportedvalues={}      # name:values exported by this job
    self.swcount=kwargs["swcount"]      # How many subworlds are there?
    self.swid=kwargs["swid"]    # which subworld are we running in?
    
    
    
  def setImportValue(self, name, v):
    """
    Use to make a value available to the job (ie called from outside the job)
    :var name: label used to identify this import
    :type name: ``str``
    :var v: value to be imported
    :type v: python object
    """
    self.importedvalues[name]=v
        
  def exportValue(self, name, v):
    """
    Make value v available to other Jobs under the label name.
    name must have already been registered with the SplitWorld instance.
    For use inside the work() method.
    :var name: registered label for exported value
    :type name: ``str``
    :var v: value to be imported
    :type v: python object
    """
    if type(name)==type([]):
        for x in name:
          if type(x)!=type(""):
            raise RuntimeError("Variable name must be a string or list of strings- instead got [%s]"%(str(type(x))))
    elif type(name)!=type(""):
      raise RuntimeError("Variable name must be a string or list of strings- instead got %s"%(str(type(name))))
    self.exportedvalues[name]=v
    
  def importValue(self, name):
    """
    For use inside the work() method.
    :var name: label for imported value.
    :type name: ``str``
    """
    if name in self.importedvalues:
        return self.importedvalues[name]
    else:
        raise KeyError("Attempt to import variable \'"+name+"\' which is not available to this job.")
        return None

  def clearExports(self):
    """
    Remove exported values from the map
    """
    self.exportedvalues.clear()
    
  def clearImports(self):
    """
    Remove imported values from their map
    """
    self.importedvalues.clear()
    
  def declareImport(self, name):
    """
    Adds name to the list of imports
    """
    if (not isinstance(name, str)) or len(name)==0:
      raise ValueError("Imports must be identified with non-empty strings")
    if not name in self.wantedvalues:
      self.wantedvalues+=name

  def work(self):
    """
    Need to be overloaded for the job to actually do anthing.
    A return value of True indicates this job thinks it is done.
    A return value of False indicates work still to be done
    """
    raise RuntimeError("work() function not overridden as required")

class FunctionJob(Job):
  """
  Takes a python function (with only self and keyword params) to be called as the work method
  """
  def __init__(self, fn, *args, **kwargs):
    super(FunctionJob, self).__init__(*args, **kwargs)
    self.__fn__ = fn
    if fn is None:
      raise ValueError("Attempt to create a Function Job with no function to run (fn argument missing).")
    self.__calldict__ = kwargs
    if "imports" in kwargs:
      if isinstance(kwargs["imports"], str):
        self.declareImport(kwargs["imports"])
      else:
        for n in kwargs["imports"]:
          self.declareImport(n)

  def work(self):
    self.__fn__(self, **self.__calldict__)
    return True

    
