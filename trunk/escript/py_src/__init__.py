
##############################################################################
#
# Copyright (c) 2013-2017 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2013-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escriptcore.escriptcpp import *
from esys.escriptcore.util import *
from esys.escriptcore.nonlinearPDE import NonlinearPDE
from esys.escriptcore.datamanager import DataManager
from esys.escriptcore.symbolic import *
from esys.escriptcore.splitworld import *

__all__=[x for x in dir() if not x.startswith('internal_') and not x.startswith('Internal_') and not x.startswith('__') and not str(type(eval(x))).find('module')>=0]


