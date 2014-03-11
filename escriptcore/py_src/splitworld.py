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

class Job:
  """
  Describes a sequence of work to be carried out in a subworld.
  The instances of this class used in the subworlds will
  be constructed by the system.
  The majority of the work done by the Job will be in the 
  *overloaded* work() method
  """

  def __init__(self, domain, id):
    """
    :var domain: Domain to be used as the basis for all ``Data`` and PDEs in this Job.
    :var id: sequence number of this job [0, number of jobs-1]
    :type id: Positive ``int``
    """
    self.domain=domain
    self.id=id
    self.incratenames={}
    self.outcratenames={}
    
  # Now does this take a crate / hamper / ....
  # There is the issue of only wanting to put in values for initial set up
  # and so merge is not required
  def setImportValue(self, name, v):
    """
    Use to make a value available to the job (ie called from outside the job)
    :var name: label used to identify this import
    :type name: ``str``
    :var v: value to be imported
    :type v: crate
    """
    pass
  
  
  def getExportValue(self, name):
    """
    get value exported by work()
    """
    pass
  
  
    
  def work(self):
    """
    Need to be overloaded for the job to actually do anthing.
    A return value of True, indicates this job thinks it is done.
    A return value of False indicates work still to be done
    """
    return True