
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
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


__copyright__="""Copyright (c) 2003-2022 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
Tests that escript returns the correct type when adding real and complex Data objects together

"""

__author__="Adam Ellery, a.ellery@uq.edu.au"

import esys.escriptcore.utestselect as unittest
from esys.escript import *

import numpy

from test_util_base import Test_util_base

class Test_addition_types(Test_util_base):
   """
   these overloaded operations still fail!

        - wrong return value of Data binaries (Mantis 0000054) 
   """
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_add_real_complex_function(self):
      x=Scalar(0, Function(self.domain))
      self.assertTrue(x.isComplex()==False,"Real data object was initiated as complex on Function")
      x+=ComplexScalar(5+8j, Function(self.domain))
      self.assertTrue(x.isComplex(),"Real + Complex is not Complex")

   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_add_complex_real_function(self):
      x=ComplexScalar(5+8j, Function(self.domain))
      self.assertTrue(x.isComplex()==True,"Complex data object was initiated as real on Function")
      x+=Scalar(5, Function(self.domain))
      self.assertTrue(x.isComplex(),"Real + Complex is not Complex")

   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_add_real_complex_continuousfunction(self):
      x=Scalar(0, ContinuousFunction(self.domain))
      self.assertTrue(x.isComplex()==False,"Real data object was initiated as complex on Function")
      x+=ComplexScalar(5+8j, ContinuousFunction(self.domain))
      self.assertTrue(x.isComplex(),"Real + Complex is not Complex")

   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_add_complex_real_continuousfunction(self):
      x=ComplexScalar(5+8j, ContinuousFunction(self.domain))
      self.assertTrue(x.isComplex()==True,"Complex data object was initiated as real on Function")
      x+=Scalar(5, ContinuousFunction(self.domain))
      self.assertTrue(x.isComplex(),"Real + Complex is not Complex")
   
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_add_real_complex_reducedfunction(self):
      x=Scalar(0, ReducedFunction(self.domain))
      self.assertTrue(x.isComplex()==False,"Real data object was initiated as complex on Function")
      x+=ComplexScalar(5+8j, ReducedFunction(self.domain))
      self.assertTrue(x.isComplex(),"Real + Complex is not Complex")

   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_add_complex_real_reducedfunction(self):
      x=ComplexScalar(5+8j, ReducedFunction(self.domain))
      self.assertTrue(x.isComplex()==True,"Complex data object was initiated as real on Function")
      x+=Scalar(5, ReducedFunction(self.domain))
      self.assertTrue(x.isComplex(),"Real + Complex is not Complex")
   