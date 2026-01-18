
##############################################################################
#
# Copyright (c) 2014-2018 by The University of Queensland
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

"""
Unit test module selection for escript.

This module is imported by all escript unit tests as "unittest". It provides
a central point to switch between the standard ``unittest`` module and
``unittest2`` (if installed) without modifying individual test scripts.
"""


__copyright__="""Copyright (c) 2014-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

#This file will be imported by all escript unit tests as "unittest"
#Replace the line below to switch between unittest and unitest2 (assuming you have it installed)
from unittest import *
