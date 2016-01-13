
##############################################################################
#
# Copyright (c) 2003-2016 by The University of Queensland
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

"""Our most general domain representation. Imports submodules into its namespace
"""

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2016 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"


import esys.escript
import esys.pasowrap    #if you don't import this, you won't be able to see methods not in AbstractSystemmatrix
#from esys.escript import *
from .finleycpp import *
from .factorywrappers import *
from .readers import *

__nodocorecursion=['finleycpp', 'factorywrappers', 'readers']
