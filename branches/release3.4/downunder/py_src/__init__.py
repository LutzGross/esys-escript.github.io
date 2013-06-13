
##############################################################################
#
# Copyright (c) 2003-2013 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

"""Data inversion module built on escript"""

__copyright__="""Copyright (c) 2003-2013 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from .costfunctions import *
from .datasources import *
from .domainbuilder import *
from .forwardmodels import *
from .inversioncostfunctions import *
from .inversions import *
from .mappings import *
from .minimizers import *
from .regularizations import *
from .coordinates import *

import logging
logging.basicConfig(format='%(name)s: %(message)s', level=logging.INFO)

#prevents our doc script from processing these packages since they are already incorporated into this one
__nodocorecursion=['costfunctions', 'datasources', ' domainbuilder', 'forwardmodels', 'inversioncostfunctions',
'inversions', 'mappings', 'minimizers', 'regularizations', 'coordinates']