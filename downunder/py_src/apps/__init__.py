
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

"""Our most general domain representation. Imports submodules into its namespace
"""

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"


import esys.escript
from .seismicModels import *
from .MTModels import *
from .magneticModels import *
from .gravityModels import *
from .dcModels import *
try:
    from esys.escriptcore.symboliccore.utils import isSymbol, symbols
except:
	pass
try:
	from esys.escriptcore.symboliccore.pretty import pretty_print, pprint
except:
	pass

__nodocorecursion=['seismicModels', 'MTModels', 'magneticModels', 'gravityModels', 'dcModels']
