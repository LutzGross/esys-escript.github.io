##############################################################################
#
# Copyright (c) 2014 by University of Queensland
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

__copyright__="""Copyright (c) 2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__author__="Joel Fenwick"

"""
This module contains the Python side of the SplitWorld functionality.

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

import warnings
warnings.simplefilter('default', category=DeprecationWarning)

from . import escriptcpp as escore

class Job(object):
  """
  Describes a sequence of work to be carried out in a subworld.
  The instances of this class used in the subworlds will
  be constructed by the system.
  The majority of the work done by the Job will be in the 
  *overloaded* work() method
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
    self.wantedvalues=[]		# names of shared values this job wishes to import    
    self.importedvalues={}	# name:values of which this jobs wants to use
    self.exportedvalues={}	# name:values exported by this job
    
  def setImportValue(self, name, v):
    """
    Use to make a value available to the job (ie called from outside the job)
    :var name: label used to identify this import
    :type name: ``str``
    :var v: value to be imported
    :type v: ?
    """
    pass
  
  
  def getExportValue(self, name):
    """
    get value exported by work()  [called from outside the job]
    """
    if name in self.exportedvalues:
	return self.exportedvalues[name]
    else:
	return None
	
  def export(self, name, v):
    """
    Make value v available to other Jobs under the label name.
    name must have already been registered with the SplitWorld instance
    :var name: registered label for exported value
    :type name: ``str``
    """
    self.exportedvalues[name]=v
    
  def getImport(self, name):
    """
    :var name: label for imported value. Returns None if import not present.
    :type name: ``str``
    """
    if name in self.importedvalues:
	return self.importedvalues[name]
    else:
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
    
  def requestImport(self, name):
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
    A return value of True, indicates this job thinks it is done.
    A return value of False indicates work still to be done
    """
    return True