
##############################################################################
#
# Copyright (c) 2003-2016 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2016 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
test for non-overloaded binary operations

:remark: use see `test_util`
:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Joel Fenwick, joelfenwick@uq.edu.au"

import esys.escriptcore.utestselect as unittest
import numpy
from esys.escript import *
from test_util_base import Test_util_values

class Test_util_binary_new(Test_util_values):
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inner_combined(self):
       opstring='inner(a,b)'
       misccheck=None   # How to work out what the result of type should be
       oraclecheck="numpy.tensordot(refa, refb, axes=getRank(refa))"
       opname="inner"
       noshapemismatch=True
       permitscalarmismatch=False
       self.generate_binary_operation_test_batch_large(opstring, misccheck, oraclecheck, opname, no_shape_mismatch=noshapemismatch, permit_scalar_mismatch=permitscalarmismatch)
           
           
           
           
