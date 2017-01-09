
##############################################################################
#
# Copyright (c) 2016 by The University of Queensland
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

__copyright__="""Copyright (c) 2016 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
test for util operations for unary operations without tagged data

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
import math
import cmath
from esys.escript import *
from test_util_base import Test_util_base, Test_util_values

haveLapack = hasFeature('lapack')

class Test_util_unary_new(Test_util_values):
   """
   test for unary operations. No tagged data are tested.
   """
   def iterateops(self, ops, vals):
       """
       """
       for p in ops:
          o,c,z=p
          for v in vals:
            res=o(v)
            if isinstance(v,complex):
               ref=z(v)
            else:
               ref=c(v)
            self.assertTrue(isinstance(res,type(ref)),"wrong type of result for "+str(o))
            self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result for "+str(o))
            d=Data(v)
            res=o(d)
            self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result for data on "+str(o))

   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def xtest_compare_complex_vs_real_data(self):
      # Compare results of unary ops provided by util and their python equivalents
      # Note that in some cases util calls these python versions so it is not a
      # guarantee that they are correct, we also compare with Data implementation
      # so we reduce the risk unless there is a fault in the system's underlying c library
      # Also note that these calls are only testing scalars
      ops=[(sin,math.sin,cmath.sin), (cos,math.cos,cmath.cos), (tan,math.tan,cmath.tan), (log,math.log,cmath.log), (log10, math.log10, cmath.log10), (Abs, abs, abs),
(acos,math.acos,cmath.acos), (acosh,math.acosh,cmath.acosh), (asin,math.asin,cmath.asin), (asinh, math.asinh,cmath.asinh),
(cosh, math.cosh, cmath.cosh), (exp, math.exp, cmath.exp), (sinh, math.sinh, cmath.sinh), (sqrt, math.sqrt, cmath.sqrt)]
      vals=[1+0j,-1+0j,1j, -1j, math.pi*1j,3+4j]
      self.iterateops(ops,vals)
      ops=[(atan,math.atan,cmath.atan)]
      vals=[1+0j,-1+0j, math.pi*1j,3+4j]
      self.iterateops(ops,vals)
      ops=[(atanh,math.atanh,cmath.atanh)]
      vals=[1j, -1j, math.pi*1j,3+4j]
      self.iterateops(ops,vals)
      # test for zero values for those functions which can take it
      vals=[0j]
      ops=[(sin,math.sin,cmath.sin), (cos,math.cos,cmath.cos), (tan,math.tan,cmath.tan),
           (Abs, abs, abs),
           (acos,math.acos,cmath.acos), (acosh,math.acosh,cmath.acosh),
           (asin,math.asin,cmath.asin), (asinh, math.asinh,cmath.asinh),
           (cosh, math.cosh, cmath.cosh), (exp, math.exp, cmath.exp), 
           (sinh, math.sinh, cmath.sinh), (sqrt, math.sqrt, cmath.sqrt),
           (atan,math.atan,cmath.atan), (atanh,math.atanh,cmath.atanh)]
      self.iterateops(ops,vals)      

   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sin_combined(self):
        supportcplx=True
        opstring="sin(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.sin(ref)"
        opname="sin"
        update1="numpy.sin(r2)"    # The updates are a problem here because this is not a reduction
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1, update2, multisteptag=False)  

   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_cos_combined(self):
        supportcplx=True
        opstring="cos(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.cos(ref)"
        opname="cos"
        update1="numpy.cos(r2)"    # The updates are a problem here because this is not a reduction
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1, update2, multisteptag=False)  
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_tan_combined(self):
        supportcplx=True
        opstring="tan(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.tan(ref)"
        opname="tan"
        update1="numpy.tan(r2)"    # The updates are a problem here because this is not a reduction
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1, update2, multisteptag=False) 
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_log_combined(self):
        supportcplx=True
        opstring="log(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.log(ref)"
        opname="log"
        update1="numpy.log(r2)"    # The updates are a problem here because this is not a reduction
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1, update2, multisteptag=False, 
            input_trans=lambda x: numpy.abs(x) if type(x) is numpy.ndarray and x.dtype.kind=='f' else abs(x) if type(x) is Data and not x.isComplex() else x)
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_log10_combined(self):
        supportcplx=True
        opstring="log10(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.log10(ref)"
        opname="log10"
        update1="numpy.log10(r2)"    # The updates are a problem here because this is not a reduction
        update2=None
        self.generate_operation_test_batch_large(supportcplx, opstring, misccheck, oraclecheck, opname, update1, update2, multisteptag=False,  
            input_trans=lambda x: numpy.abs(x) if type(x) is numpy.ndarray and x.dtype.kind=='f' else abs(x) if type(x) is Data and not x.isComplex() else x)         
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_exp_combined(self):
        supportcplx=True
        opstring="exp(a)"
        misccheck="isinstance(res,type(a))"
        oraclecheck="numpy.exp(ref)"
        opname="exp"
        update1="numpy.exp(r2)"    # The updates are a problem here because this is not a reduction
        update2=None
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1, update2, multisteptag=False)         
