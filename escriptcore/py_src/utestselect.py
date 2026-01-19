
##############################################################################
#
# Copyright (c) 2014-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file in the repository root for contributors and development history
# https://github.com/LutzGross/esys-escript.github.io/blob/master/CREDITS
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
__url__="https://github.com/LutzGross/esys-escript.github.io"

#This file will be imported by all escript unit tests as "unittest"
#Replace the line below to switch between unittest and unitest2 (assuming you have it installed)
from unittest import *
