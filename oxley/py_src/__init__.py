##############################################################################
#
# Copyright (c) 2003-2019 by The University of Queensland
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

"""A domain meshed with uniform rectangles or quadrilaterals
"""


__copyright__="""Copyright (c) 2003-2019 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import esys.escript       # This is just to ensure required libraries are loaded
from .oxleycpp import *

from esys.oxley.RefinementZone import *

__nodocorecursion=['oxleycpp']
