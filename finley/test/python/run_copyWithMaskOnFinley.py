
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
import numpy as np
from esys.finley import Brick


class TestCopyWithMaskOnFinley(unittest.TestCase):
   def setUp(self):
      self.domain = Brick(10,10,10)
   def tearDown(self):
      del self.domain
   def test_CopyWithMask_Assertions(self):
      # -- mismatch of dimension
      a=Data(0, (3,4), Solution(self.domain))
      b=Data(0, (3,4), Solution(self.domain))
      m=Data(1, (3,4, 6), Solution(self.domain))
      self.assertRaises(RuntimeError, a.copyWithMask, b,m)
      # -- mismatch of dimension
      a=Data(0, (3,2), Solution(self.domain))
      b=Data(0, (3,4), Solution(self.domain))
      m=Data(1, (3,4), Solution(self.domain))
      self.assertRaises(RuntimeError, a.copyWithMask, b,m)
      # -- mismatch of dimension
      a=Data(0, (3,4), Solution(self.domain))
      b=Data(0, (3,2), Solution(self.domain))
      m=Data(1, (3,4), Solution(self.domain))
      self.assertRaises(RuntimeError, a.copyWithMask, b,m)
      # -- match of FS
      a=Data(0, (3,2), Function(self.domain))
      b=Data(0, (3,2), Function(self.domain))
      m=Data(1, (3,2), Function(self.domain))
      a.copyWithMask(b,m)
      # -- make match of FS
      a=Data(0, (3,2), Function(self.domain))
      b=Data(0, (3,2), Solution(self.domain))
      m=Data(1, (3,2), Solution(self.domain))
      a.copyWithMask(b,m)
      # -- make match of FS
      a=Data(0, (3,2), Function(self.domain))
      b=Data(0, (3,2), Function(self.domain))
      m=Data(1, (3,2), Solution(self.domain))
      a.copyWithMask(b,m)
      # -- make match of FS
      a=Data(0, (3,2), Function(self.domain))
      b=Data(0, (3,2), Solution(self.domain))
      m=Data(1, (3,2), Function(self.domain))
      a.copyWithMask(b,m)
      # -- make match of FS
      a=Data(0, (3,2), Function(self.domain))
      b=Data(0, (3,2), Solution(self.domain))
      m=Data(1, (3,2), Solution(self.domain))
      a.copyWithMask(b,m)
      # -- mismatch of FS
      a=Data(0, (3,2), Solution(self.domain))
      b=Data(0, (3,2), Function(self.domain))
      m=Data(1, (3,2), Solution(self.domain))
      self.assertRaises(RuntimeError, a.copyWithMask, b,m)
      # -- mismatch of FS
      a=Data(0, (3,2), Solution(self.domain))
      b=Data(0, (3,2), Solution(self.domain))
      m=Data(1, (3,2), Function(self.domain))
      self.assertRaises(RuntimeError, a.copyWithMask, b,m)
      # -- mismatch of FS
      a=Data(0, (3,2), Solution(self.domain))
      b=Data(0, (3,2), FunctionOnBoundary(self.domain))
      m=Data(1, (3,2), Function(self.domain))
      self.assertRaises(RuntimeError, a.copyWithMask, b,m)
   def test_CopyWithMask_Rank0_fullMask_singleIn_singleMask(self):
      a = Data(1, (), FunctionOnBoundary(self.domain))
      m = Data(1, (), FunctionOnBoundary(self.domain))
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank0_fullMask_singleIn_taggedMask(self):
      a = Data(1, (), FunctionOnBoundary(self.domain))
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank0_fullMask_singleIn_expandedMask(self):
      a = Data(1, (), FunctionOnBoundary(self.domain))
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank0_fullMask_taggedIn_singleMask(self):
      a = Data(0, (), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', 1)
      a.setTaggedValue('bottom', 1)
      m = Data(1, (), FunctionOnBoundary(self.domain))
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank0_fullMask_taggedIn_taggedMask(self):
      a = Data(0, (), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', 1)
      a.setTaggedValue('bottom', 1)
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank0_fullMask_taggedIn_expandedMask(self):
      a = Data(0, (), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', 1)
      a.setTaggedValue('bottom', 1)
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank0_fullMask_expandedIn_singleMask(self):
      a = length(FunctionOnBoundary(self.domain).getX())
      m = Data(1, (), FunctionOnBoundary(self.domain))
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank0_fullMask_expandedIn_taggedMask(self):
      a = length(FunctionOnBoundary(self.domain).getX())
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank0_fullMask_expandedIn_expandedMask(self):
      a = length(FunctionOnBoundary(self.domain).getX())
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank0_scalarMask_singleIn_singleMask(self):
      a = Data(1, (), FunctionOnBoundary(self.domain))
      m = Data(1, (), FunctionOnBoundary(self.domain))
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank0_scalarMask_singleIn_taggedMask(self):
      a = Data(1, (), FunctionOnBoundary(self.domain))
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank0_scalarMask_singleIn_expandedMask(self):
      a = Data(1, (), FunctionOnBoundary(self.domain))
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank0_scalarMask_taggedIn_singleMask(self):
      a = Data(0, (), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', 1)
      a.setTaggedValue('bottom', 1)
      m = Data(1, (), FunctionOnBoundary(self.domain))
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank0_scalarMask_taggedIn_taggedMask(self):
      a = Data(0, (), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', 1)
      a.setTaggedValue('bottom', 1)
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank0_scalarMask_taggedIn_expandedMask(self):
      a = Data(0, (), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', 1)
      a.setTaggedValue('bottom', 1)
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank0_scalarMask_expandedIn_singleMask(self):
      a = length(FunctionOnBoundary(self.domain).getX())
      m = Data(1, (), FunctionOnBoundary(self.domain))
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank0_scalarMask_expandedIn_taggedMask(self):
      a = length(FunctionOnBoundary(self.domain).getX())
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank0_scalarMask_expandedIn_expandedMask(self):
      a = length(FunctionOnBoundary(self.domain).getX())
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank1_fullMask_singleIn_singleMask(self):
      fac = np.array([0.32123118254552296, 0.6787911399871708])
      mfac = np.array([0, 1])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank1_fullMask_singleIn_taggedMask(self):
      fac = np.array([0.5721394995303497, 0.11192321901591529])
      mfac = np.array([0, 0])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (2,), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank1_fullMask_singleIn_expandedMask(self):
      fac = np.array([0.9050259674049648, 0.7607413037163073])
      mfac = np.array([0, 1])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (2,), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      m.expand()
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank1_fullMask_taggedIn_singleMask(self):
      fac = np.array([0.264029706032198, 0.32519592828550925])
      mfac = np.array([0, 0])
      a = Data(0, (2,), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank1_fullMask_taggedIn_taggedMask(self):
      fac = np.array([0.6045753865040275, 0.36450689865866326])
      mfac = np.array([1, 0])
      a = Data(0, (2,), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (2,), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank1_fullMask_taggedIn_expandedMask(self):
      fac = np.array([0.4686446671772245, 0.9010185568721054])
      mfac = np.array([0, 1])
      a = Data(0, (2,), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (2,), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      m.expand()
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank1_fullMask_expandedIn_singleMask(self):
      fac = np.array([0.12794324841979998, 0.3681562857309649])
      mfac = np.array([1, 1])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank1_fullMask_expandedIn_taggedMask(self):
      fac = np.array([0.026217656258289646, 0.5427135201750596])
      mfac = np.array([0, 1])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (2,), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank1_fullMask_expandedIn_expandedMask(self):
      fac = np.array([0.19836825314429662, 0.451721873307932])
      mfac = np.array([0, 1])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (2,), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      m.expand()
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank1_scalarMask_singleIn_singleMask(self):
      fac = np.array([0.08034071875966065, 0.914178123830289])
      mfac = np.array([0, 0])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank1_scalarMask_singleIn_taggedMask(self):
      fac = np.array([0.06925114510367436, 0.6162839652033565])
      mfac = np.array([1, 0])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank1_scalarMask_singleIn_expandedMask(self):
      fac = np.array([0.8453412713112461, 0.782840219347851])
      mfac = np.array([1, 0])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank1_scalarMask_taggedIn_singleMask(self):
      fac = np.array([0.7371409073480271, 0.3237464282946272])
      mfac = np.array([1, 1])
      a = Data(0, (2,), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank1_scalarMask_taggedIn_taggedMask(self):
      fac = np.array([0.39536773878375375, 0.6277234734805064])
      mfac = np.array([1, 1])
      a = Data(0, (2,), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank1_scalarMask_taggedIn_expandedMask(self):
      fac = np.array([0.7377962188401396, 0.2529523396823189])
      mfac = np.array([0, 1])
      a = Data(0, (2,), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank1_scalarMask_expandedIn_singleMask(self):
      fac = np.array([0.4499749243155703, 0.6188535582977776])
      mfac = np.array([0, 0])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank1_scalarMask_expandedIn_taggedMask(self):
      fac = np.array([0.32088049575350064, 0.7896152490779662])
      mfac = np.array([1, 0])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank1_scalarMask_expandedIn_expandedMask(self):
      fac = np.array([0.7030947249695527, 0.2778685376231357])
      mfac = np.array([0, 0])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (2,), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank2_fullMask_singleIn_singleMask(self):
      fac = np.array([[0.8944149221396697, 0.8354190729458718, 0.9100219038577015], [0.8973141936383877, 0.9049123615957256, 0.4857222888445807], [0.37691982906972665, 0.7889486937138986, 0.6316499753220514], [0.9578654531431999, 0.3247215800551574, 0.7386010398473836]])
      mfac = np.array([[0, 1, 1], [0, 0, 0], [1, 0, 1], [0, 0, 1]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank2_fullMask_singleIn_taggedMask(self):
      fac = np.array([[0.9181320038983173, 0.4910613937104493, 0.4127400459919751], [0.20986032405741983, 0.5408372020393524, 0.28405712964111285], [0.8877594288826252, 0.6398042412385332, 0.9937134307582923], [0.7576747455558622, 0.3778381670780182, 0.6567881025241645]])
      mfac = np.array([[0, 1, 1], [1, 0, 1], [1, 0, 1], [1, 1, 1]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (4, 3), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank2_fullMask_singleIn_expandedMask(self):
      fac = np.array([[0.592763958703754, 0.28246561312733554, 0.7490274852415424], [0.9011195928967389, 0.8258519249829543, 0.7410819601250819], [0.43386726728207714, 0.3236791882106125, 0.12973497533924916], [0.06880947545365335, 0.8047889630552482, 0.1645525219965872]])
      mfac = np.array([[1, 0, 0], [1, 0, 1], [1, 1, 1], [1, 0, 0]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (4, 3), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      m.expand()
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank2_fullMask_taggedIn_singleMask(self):
      fac = np.array([[0.17951857007181538, 0.31463014585890836, 0.4429204901066113], [0.4046307543351202, 0.6852538595885757, 0.40114108952580185], [0.9375199751244503, 0.5693769345940233, 0.12491399137117509], [0.6783873159940194, 0.8415577868503648, 0.9273058188815453]])
      mfac = np.array([[1, 1, 0], [0, 1, 1], [0, 1, 1], [1, 1, 0]])
      a = Data(0, (4, 3), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank2_fullMask_taggedIn_taggedMask(self):
      fac = np.array([[0.09737526034492594, 0.9097236982636414, 0.27478497175244854], [0.043617448646725676, 0.21966754350410944, 0.7646893990048953], [0.3408013829464883, 0.9020112163141054, 0.3893583020844449], [0.10125662781169043, 0.8505381531485663, 0.9942342370758974]])
      mfac = np.array([[0, 0, 0], [1, 0, 1], [0, 0, 0], [0, 0, 1]])
      a = Data(0, (4, 3), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (4, 3), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank2_fullMask_taggedIn_expandedMask(self):
      fac = np.array([[0.4861737244238473, 0.6459831066406707, 0.8423847502435381], [0.6945420268721688, 0.4602299158384736, 0.010015953723026438], [0.5982170547182503, 0.16791393581397396, 0.17896368558223985], [0.0372759480646605, 0.9162007286943107, 0.8915604057450836]])
      mfac = np.array([[0, 1, 1], [0, 0, 0], [1, 1, 0], [1, 1, 0]])
      a = Data(0, (4, 3), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (4, 3), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      m.expand()
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank2_fullMask_expandedIn_singleMask(self):
      fac = np.array([[0.295231946512696, 0.8402287420334332, 0.03497931866501014], [0.18945322743567627, 0.9333393516668791, 0.2385414690498604], [0.1829976728060626, 0.42666424220634136, 0.509432722318471], [0.9336233463252984, 0.1487176584768516, 0.9249270511720536]])
      mfac = np.array([[0, 0, 1], [0, 1, 1], [0, 0, 0], [1, 1, 1]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank2_fullMask_expandedIn_taggedMask(self):
      fac = np.array([[0.5837808904630759, 0.10515337087580123, 0.545958237109611], [0.2991177994833264, 0.5782737523951408, 0.7145251438257507], [0.7586006564325822, 0.3115324527442639, 0.8750278737521898], [0.862827718712728, 0.6428578262967539, 0.8567815400295898]])
      mfac = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 0], [0, 1, 0]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (4, 3), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank2_fullMask_expandedIn_expandedMask(self):
      fac = np.array([[0.03251287198782116, 0.9443144160930977, 0.7107270645347499], [0.4412742884502001, 0.5278598570921272, 0.30241636630481805], [0.8548812805838961, 0.8692150784542156, 0.0766626530775556], [0.4366603612512563, 0.23183306012955862, 0.7988880615828665]])
      mfac = np.array([[1, 1, 0], [1, 1, 0], [1, 1, 1], [0, 0, 0]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (4, 3), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      m.expand()
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank2_scalarMask_singleIn_singleMask(self):
      fac = np.array([[0.8435553832253655, 0.87878960522116, 0.04823138100654756], [0.9538632521509074, 0.22186688163016177, 0.6544306430391585], [0.0377406570876313, 0.6063023253814719, 0.7507323761763305], [0.06472973063860799, 0.01438176753856546, 0.06602853971402933]])
      mfac = np.array([[1, 1, 1], [0, 1, 0], [0, 0, 1], [0, 1, 1]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank2_scalarMask_singleIn_taggedMask(self):
      fac = np.array([[0.4013307827812388, 0.7328465214487374, 0.025113520863358052], [0.5792184956773775, 0.5950307543748828, 0.6120943270541], [0.38388224549044037, 0.9695038607960872, 0.377500094599426], [0.4450618571747347, 0.07345797794795694, 0.5400691702129372]])
      mfac = np.array([[0, 0, 0], [0, 0, 1], [1, 1, 1], [1, 0, 1]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank2_scalarMask_singleIn_expandedMask(self):
      fac = np.array([[0.2628633953650107, 0.49028165320185735, 0.8912667203582183], [0.5234464884117314, 0.0783982761513412, 0.41542111212263855], [0.46051564672654854, 0.9775890788893297, 0.14148091013712438], [0.9099279809257299, 0.5825590215004082, 0.17180752891132334]])
      mfac = np.array([[0, 0, 0], [1, 1, 0], [0, 1, 1], [0, 0, 0]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank2_scalarMask_taggedIn_singleMask(self):
      fac = np.array([[0.892154744361936, 0.5305745964218176, 0.6539107984578421], [0.7540434804409933, 0.5106151061257204, 0.6844638885117594], [0.17214498293403202, 0.5597978999618088, 0.35910370641610556], [0.0015844047438517972, 0.2066162129089043, 0.6553650412226485]])
      mfac = np.array([[1, 0, 1], [0, 1, 1], [1, 0, 0], [1, 1, 1]])
      a = Data(0, (4, 3), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank2_scalarMask_taggedIn_taggedMask(self):
      fac = np.array([[0.5115530974920427, 0.22861425239780908, 0.6199212300605544], [0.5885205947075081, 0.6787350454738361, 0.510215342346114], [0.5166155604369884, 0.11227729694732114, 0.4877924432129529], [0.834934020551455, 0.08409592084694917, 0.1402067467510405]])
      mfac = np.array([[0, 1, 1], [0, 0, 0], [1, 0, 0], [1, 0, 1]])
      a = Data(0, (4, 3), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank2_scalarMask_taggedIn_expandedMask(self):
      fac = np.array([[0.9387408027252425, 0.8436948118420232, 0.12219625889848962], [0.32313813072650643, 0.26371862699894166, 0.014014636962692228], [0.8477918442722557, 0.2523252402161945, 0.8149907194963943], [0.4900467956938588, 0.5709208408849484, 0.9226370172502123]])
      mfac = np.array([[1, 0, 0], [0, 1, 1], [0, 1, 0], [1, 0, 1]])
      a = Data(0, (4, 3), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank2_scalarMask_expandedIn_singleMask(self):
      fac = np.array([[0.024017843365157, 0.394717515384597, 0.9206035566187638], [0.3817657857191755, 0.6598748077286158, 0.41368479777852285], [0.6775501113422915, 0.5063051058506731, 0.7925443852093711], [0.7861269240936445, 0.27262952153903974, 0.373028054922268]])
      mfac = np.array([[0, 0, 1], [1, 1, 1], [1, 1, 1], [0, 0, 0]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank2_scalarMask_expandedIn_taggedMask(self):
      fac = np.array([[0.5645421152367375, 0.05270685122085883, 0.02177909615001916], [0.46060131855342756, 0.45025299850417055, 0.8636907833201748], [0.7739038560337427, 0.6142599682735335, 0.5030612559381432], [0.3739208228946518, 0.4435605112683477, 0.4329048829834494]])
      mfac = np.array([[0, 1, 1], [0, 0, 1], [1, 1, 0], [0, 0, 1]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank2_scalarMask_expandedIn_expandedMask(self):
      fac = np.array([[0.4881377268015522, 0.811665122913442, 0.7307567225098984], [0.17266876185757085, 0.9651942560316845, 0.49363604676713624], [0.18148380244295648, 0.4678593650269771, 0.09317841895452239], [0.5383364802091514, 0.7688383351516204, 0.3648319603275062]])
      mfac = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 0], [0, 0, 1]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (4, 3), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertTrue(m.isExpanded())


   def test_CopyWithMask_Rank3_fullMask_singleIn_singleMask(self):
      fac = np.array([[[0.6564036549921607, 0.8633465631752325], [0.729214858258692, 0.4939764445539492], [0.821219012635921, 0.852285979688934], [0.05418707179190685, 0.9774961555873364], [0.14287958922457877, 0.015455488721642485]], [[0.7642955050777925, 0.32156592743650925], [0.7204523901089209, 0.7337910042397893], [0.14485248190411948, 0.5367126538990786], [0.49576690287736225, 0.9642126565313043], [0.1074501458126651, 0.5458558222226234]]])
      mfac = np.array([[[0, 0], [0, 0], [1, 1], [1, 0], [1, 0]], [[0, 0], [0, 1], [0, 1], [1, 0], [1, 1]]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank3_fullMask_singleIn_taggedMask(self):
      fac = np.array([[[0.5624820066544872, 0.2576590099423097], [0.10378741671064573, 0.9632847968660542], [0.09744534201540855, 0.7789006194684606], [0.5706874095582177, 0.3133382276625707], [0.6711257490865259, 0.4370260380471592]], [[0.6767805240190359, 0.5154722101642899], [0.6969618078004199, 0.5293340302185826], [0.6405241222226302, 0.9937054534576855], [0.6180077508599291, 0.30089500853105766], [0.96257104016361, 0.9597187786253485]]])
      mfac = np.array([[[0, 1], [1, 0], [1, 0], [1, 0], [1, 1]], [[1, 1], [0, 0], [1, 0], [0, 1], [1, 0]]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (2, 5, 2), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank3_fullMask_singleIn_expandedMask(self):
      fac = np.array([[[0.13830835106729078, 0.9829526347164134], [0.46140664658172126, 0.5973647270711099], [0.5986085965240916, 0.04054786860666726], [0.224897422032735, 0.5138409418302913], [0.1633630064285556, 0.08136671108296312]], [[0.30352275294212594, 0.12337397527864236], [0.8067534252942998, 0.32508129591317536], [0.8152202807567452, 0.37493922728512086], [0.500200913950108, 0.4046530952740145], [0.9782834588044695, 0.09743924587589803]]])
      mfac = np.array([[[1, 1], [1, 1], [1, 1], [0, 1], [0, 1]], [[1, 1], [0, 0], [1, 1], [0, 0], [1, 1]]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (2, 5, 2), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      m.expand()
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank3_fullMask_taggedIn_singleMask(self):
      fac = np.array([[[0.13794320791634096, 0.8015193739062891], [0.38367031684160335, 0.6133781487021178], [0.4535674230150366, 0.17857242911733617], [0.8260052573386132, 0.4456545419110418], [0.4407693220752249, 0.9950562321650132]], [[0.5043364036241831, 0.29471368695928524], [0.42535983001757915, 0.2807256159535285], [0.5743203468366933, 0.027197189924305976], [0.03707406920293288, 0.517730726313475], [0.8064644348907217, 0.5742286220194238]]])
      mfac = np.array([[[0, 0], [1, 0], [0, 0], [0, 0], [1, 0]], [[1, 0], [0, 1], [0, 1], [0, 0], [1, 1]]])
      a = Data(0, (2, 5, 2), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank3_fullMask_taggedIn_taggedMask(self):
      fac = np.array([[[0.4337239318956304, 0.24881823615198295], [0.8802761891044987, 0.03063939727778575], [0.6019123427526294, 0.9880333350237059], [0.9755676996956182, 0.056539284314686244], [0.08407770420800487, 0.0999986359040198]], [[0.7794921456494565, 0.8479241738510165], [0.2905658561702277, 0.13117855894924413], [0.06158515149494448, 0.5696238293886696], [0.8362108214120513, 0.16954718253278378], [0.2890317189814722, 0.1805049787174814]]])
      mfac = np.array([[[1, 0], [1, 0], [0, 0], [1, 0], [0, 0]], [[1, 1], [0, 1], [0, 1], [1, 0], [0, 0]]])
      a = Data(0, (2, 5, 2), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (2, 5, 2), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank3_fullMask_taggedIn_expandedMask(self):
      fac = np.array([[[0.4940311678790128, 0.611453273729118], [0.9178790607251128, 0.9046676326614964], [0.6771024981737814, 0.7355737710019532], [0.22489958116697728, 0.6518828638730916], [0.025115087881785292, 0.8117599997888941]], [[0.05420869423044883, 0.2918126328720394], [0.03574986139941516, 0.4975687382232511], [0.13014328722008883, 0.034107538473945875], [0.3385268777920627, 0.36839816728279995], [0.6254929724494126, 0.39309009559218855]]])
      mfac = np.array([[[1, 1], [0, 1], [1, 0], [0, 1], [1, 0]], [[1, 1], [0, 0], [0, 1], [1, 0], [0, 1]]])
      a = Data(0, (2, 5, 2), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (2, 5, 2), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      m.expand()
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank3_fullMask_expandedIn_singleMask(self):
      fac = np.array([[[0.034113571858515046, 0.6385531498211743], [0.6983611462803733, 0.09054322278475213], [0.39232619286410786, 0.9541032418062265], [0.6916141086264136, 0.562782200293913], [0.7113685858607771, 0.3045348173590283]], [[0.8692292962800774, 0.6998838724292705], [0.4113764262227759, 0.952732455898779], [0.740441389831644, 0.8852950984705812], [0.17342912585570625, 0.669126157634618], [0.3280617994953763, 0.8126503455507534]]])
      mfac = np.array([[[1, 0], [0, 1], [0, 1], [0, 0], [1, 0]], [[1, 0], [0, 1], [0, 1], [0, 0], [1, 0]]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank3_fullMask_expandedIn_taggedMask(self):
      fac = np.array([[[0.7960389115973113, 0.9707389007033377], [0.20397594723833057, 0.018976162736547675], [0.33373649594706134, 0.6783179183821797], [0.20229218645341096, 0.2651089035688211], [0.5290979204931201, 0.7962044904668862]], [[0.8566409534904948, 0.8530361791560062], [0.047446113368958875, 0.5766910118450649], [0.6961905391100597, 0.8615325303708437], [0.26270992892068434, 0.9436711077080928], [0.9481570806796188, 0.381150164712039]]])
      mfac = np.array([[[1, 0], [0, 0], [1, 0], [0, 1], [0, 0]], [[1, 0], [0, 0], [0, 1], [0, 1], [1, 0]]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (2, 5, 2), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank3_fullMask_expandedIn_expandedMask(self):
      fac = np.array([[[0.61508091257228, 0.6217656191241002], [0.7835567072545406, 0.5459030996946524], [0.9392035134349033, 0.8678077255812334], [0.3532755929484773, 0.6026530499517566], [0.9947508948941503, 0.9264532308031665]], [[0.9207988463088429, 0.3081355597954878], [0.48246911830996597, 0.21638406362349227], [0.23341847559638773, 0.1671362789375912], [0.8541691786977695, 0.7900010585013276], [0.9165816983704193, 0.000872429378492523]]])
      mfac = np.array([[[1, 1], [0, 1], [0, 1], [1, 0], [0, 1]], [[1, 1], [1, 0], [0, 0], [1, 1], [1, 0]]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (2, 5, 2), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      m.expand()
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank3_scalarMask_singleIn_singleMask(self):
      fac = np.array([[[0.826792681549137, 0.6694134357480092], [0.322273837444097, 0.2594209639606151], [0.48772074522596576, 0.0794658146111602], [0.5855427830972223, 0.8307817801925198], [0.3389321399087044, 0.6912261739536742]], [[0.7450964370051907, 0.15162244207115783], [0.40873193698435883, 0.8494489442477529], [0.9907172677609853, 0.007855811462442297], [0.8813941080565812, 0.2592377726608346], [0.9514363341904258, 0.7990322986777129]]])
      mfac = np.array([[[1, 1], [1, 1], [0, 1], [1, 0], [1, 0]], [[1, 0], [0, 1], [0, 0], [0, 0], [1, 0]]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank3_scalarMask_singleIn_taggedMask(self):
      fac = np.array([[[0.17070860118677023, 0.6340284648733763], [0.42018362620839356, 0.5431950173714032], [0.9950879898326694, 0.9473978926928388], [0.10894235368825245, 0.08696812722317337], [0.9156434803413623, 0.9790245567377176]], [[0.48016479515096255, 0.841306631926787], [0.6546204272816352, 0.745869540265048], [0.6925110305595166, 0.5888591014428405], [0.8159822604179795, 0.7063726338821424], [0.006997394503951093, 0.40436878578484325]]])
      mfac = np.array([[[1, 0], [1, 1], [1, 0], [0, 0], [1, 0]], [[0, 1], [1, 1], [1, 1], [1, 0], [0, 1]]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank3_scalarMask_singleIn_expandedMask(self):
      fac = np.array([[[0.14703198145218233, 0.041501172700075806], [0.8701883205997202, 0.3711491916291281], [0.5233420448436713, 0.39709272508108406], [0.791172252273662, 0.0246044868974723], [0.2878615345028841, 0.4528945287209404]], [[0.2691341771983603, 0.9445363862547972], [0.9157450896211511, 0.5783328129915404], [0.6684125007351824, 0.7024178124484725], [0.007432992183239229, 0.48625862919064566], [0.2353669029784644, 0.07431110702348498]]])
      mfac = np.array([[[0, 0], [0, 1], [0, 1], [1, 1], [1, 0]], [[0, 0], [0, 0], [0, 1], [1, 1], [0, 1]]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank3_scalarMask_taggedIn_singleMask(self):
      fac = np.array([[[0.40894321970093017, 0.04830987534644782], [0.36606258579651685, 0.6607611318721096], [0.9921770324058112, 0.5465426535302496], [0.8466111860377398, 0.18043852460580045], [0.7237201997077487, 0.43825380355541765]], [[0.675522904424004, 0.7110352095747112], [0.577357771394378, 0.7963648987115823], [0.2983544936014152, 0.12758037010906376], [0.07663742124718076, 0.20678279139113964], [0.5119087203231665, 0.19425716812829408]]])
      mfac = np.array([[[1, 0], [1, 0], [1, 0], [0, 1], [1, 1]], [[1, 1], [0, 0], [1, 0], [1, 1], [0, 1]]])
      a = Data(0, (2, 5, 2), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank3_scalarMask_taggedIn_taggedMask(self):
      fac = np.array([[[0.4031658578237126, 0.9946445664099264], [0.07027445127691834, 0.3894973463123723], [0.18575665289622978, 0.24383442505942776], [0.2871408138186984, 0.4903995522012169], [0.7497363291089792, 0.4483175606890121]], [[0.17870417081233847, 0.0563213359755107], [0.8894634631059644, 0.8597144261685833], [0.4137380830075119, 0.4413967459390533], [0.3614307706579817, 0.8295560129338752], [0.5971548082584723, 0.9816547245256385]]])
      mfac = np.array([[[1, 0], [0, 1], [1, 1], [1, 1], [0, 0]], [[0, 0], [1, 0], [0, 1], [1, 1], [1, 0]]])
      a = Data(0, (2, 5, 2), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank3_scalarMask_taggedIn_expandedMask(self):
      fac = np.array([[[0.8000698999069984, 0.6579496117743114], [0.23620701151764245, 0.22826331315583837], [0.802794880974708, 0.20591219390400628], [0.9146775794981367, 0.013443827603828673], [0.5584245591072098, 0.46485723741214136]], [[0.1960329402587191, 0.7210581169862785], [0.287884960571663, 0.3833354034311701], [0.45276311899106125, 0.8465376792024234], [0.6814367215879954, 0.8828082701895011], [0.20439388101438083, 0.24027571526121783]]])
      mfac = np.array([[[0, 1], [0, 1], [1, 1], [0, 1], [1, 0]], [[1, 0], [1, 0], [1, 1], [1, 1], [1, 1]]])
      a = Data(0, (2, 5, 2), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank3_scalarMask_expandedIn_singleMask(self):
      fac = np.array([[[0.1978997064518766, 0.2679993089470952], [0.14909268229865513, 0.08364947138788081], [0.12493779860479326, 0.897446647489098], [0.18244483082353524, 0.6069883884846459], [0.5852608334588592, 0.5434667748904877]], [[0.5176389142496258, 0.39925797292510445], [0.5749092751262712, 0.21864597048456125], [0.29507546853203726, 0.5286847817445772], [0.21473798967081825, 0.33170279403942393], [0.05143641197257476, 0.36208729145385354]]])
      mfac = np.array([[[1, 0], [0, 0], [0, 1], [0, 0], [1, 0]], [[0, 0], [1, 0], [0, 0], [1, 1], [0, 0]]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank3_scalarMask_expandedIn_taggedMask(self):
      fac = np.array([[[0.49732969387762105, 0.03457370408018268], [0.9300207618660574, 0.5428190751392273], [0.6759600424397672, 0.2669438681375854], [0.48452053639833026, 0.8507943886797126], [0.8276509466186609, 0.6696672238119361]], [[0.22159753503769175, 0.7045888030424617], [0.14446293121636422, 0.8555989966786768], [0.16145705033473912, 0.6394242959052708], [0.7420525092410043, 0.15111142415243384], [0.8167058759610161, 0.6340774929327975]]])
      mfac = np.array([[[0, 0], [1, 1], [1, 0], [0, 1], [0, 1]], [[0, 1], [0, 0], [0, 1], [0, 0], [0, 1]]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank3_scalarMask_expandedIn_expandedMask(self):
      fac = np.array([[[0.18004062418642097, 0.9186464091137632], [0.25803657940730296, 0.25775714754925205], [0.3497282291296858, 0.6084347140518583], [0.6644614891783129, 0.8956593888539043], [0.04096909701605411, 0.7016648679416216]], [[0.7687823239194516, 0.9557506374454043], [0.6181184209956104, 0.2212179008748718], [0.5944805878990572, 0.24836816717179755], [0.21031468287113853, 0.12584767066680091], [0.6566029596186534, 0.9440779751740608]]])
      mfac = np.array([[[1, 1], [0, 0], [1, 0], [0, 1], [1, 0]], [[0, 1], [1, 1], [1, 1], [0, 1], [0, 0]]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (2, 5, 2), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertTrue(m.isExpanded())


   def test_CopyWithMask_Rank4_fullMask_singleIn_singleMask(self):
      fac = np.array([[[[0.8835338266670086, 0.12691548757706528, 0.8112969519133564, 0.3804644900416345, 0.47285231037604414, 0.7845741508056099], [0.38814768826016144, 0.8740232609637932, 0.07760311006467979, 0.30627244647295404, 0.9064851303662419, 0.80748598567858], [0.7885929000489532, 0.047966734751030526, 0.5715759612808118, 0.7543076139983727, 0.5723841223081648, 0.521602801378494], [0.09697179913722298, 0.7526146106202443, 0.5758314474078271, 0.9935232821205765, 0.41322712196981826, 0.358863476992171], [0.316248756145588, 0.4918532743784195, 0.03149087947891094, 0.7305260815575416, 0.2979630654370272, 0.5309090780637349]], [[0.12784372677104971, 0.9822803190207813, 0.23998537507770268, 0.48309083591869006, 0.184153233684748, 0.27295156513843355], [0.6477574612507196, 0.7976508130959694, 0.3554682959739608, 0.3310838602164349, 0.31248854297901285, 0.8978236996757646], [0.9248395091734618, 0.4642577042046203, 0.07842236542093406, 0.08368231389442227, 0.10604440291463102, 0.959401614118765], [0.3844230762308489, 0.7661553593345756, 0.9014270686675457, 0.6570939010743703, 0.9314220232291852, 0.3958536884179843], [0.770799369970392, 0.7285477263150595, 0.8365843443608341, 0.9742241828086639, 0.29778967826276403, 0.5247661389932801]]], [[[0.5445456892179997, 0.02239667855530758, 0.07093015945053494, 0.5900927379700044, 0.3700115793921469, 0.3154479493139475], [0.5657090945120083, 0.11555991723679193, 0.37682023979278145, 0.2872215103693584, 0.4341004861861747, 0.4832585145988372], [0.11886721752390661, 0.2266731248241678, 0.35536524878354925, 0.7177705446228131, 0.2229868926705313, 0.7133933647326096], [0.7316082518978863, 0.9281473929852657, 0.9316006235362341, 0.7981766002276496, 0.9348297057431996, 0.5877489314182088], [0.8429311451251411, 0.1575381374679108, 0.3089151474956241, 0.6796709880791347, 0.5604160765250541, 0.7027505281291393]], [[0.7385336436344062, 0.9658737000690376, 0.14686556754030755, 0.1480921443660762, 0.8692688665044322, 0.23851578690021713], [0.8768846278069224, 0.19368287535034967, 0.5354437858995682, 0.38164246915774125, 0.6092191372868754, 0.5562593701639089], [0.6171608658845622, 0.4252067135302582, 0.5288408349036915, 0.8824600628697246, 0.22951712012749081, 0.8885523668850565], [0.9117659138437635, 0.7898751688517968, 0.8470547465708546, 0.9041087273533114, 0.17383889666882202, 0.9689766078972756], [0.6469837055295589, 0.7386113177864063, 0.6187065838772068, 0.9266531517257032, 0.6537328752285891, 0.5779830571763149]]], [[[0.5510569517051453, 0.705598963863549, 0.9292067552748694, 0.5277256300736468, 0.2316087793770153, 0.43508551831107745], [0.35100280093257197, 0.3030911385839543, 0.8821179994513518, 0.6501794813649187, 0.2964436314099853, 0.3934938262302614], [0.40671430956584154, 0.1129228619453424, 0.29119430833130266, 0.20429175308260628, 0.7093236605406713, 0.7623597239227305], [0.34946076385439095, 0.6625181496147041, 0.8757728965207844, 0.7380535806461574, 0.49644798752956054, 0.7875562596723811], [0.12508474604709374, 0.02572426163622643, 0.2380535389424061, 0.7826774694225801, 0.06822132768165245, 0.6706248391359011]], [[0.9768962599743218, 0.6489763980128935, 0.09709606085698941, 0.7501545985244358, 0.3561614997510165, 0.8648455557053425], [0.9194600388892558, 0.29709768557653227, 0.0678527018902777, 0.02580294335166522, 0.19948696170650415, 0.9521646720475428], [0.2495631255358992, 0.7343880917085921, 0.09335079806700763, 0.6164229869557399, 0.7989936249744775, 0.017571585188680072], [0.5967841367310877, 0.8469463412902066, 0.08040951630169624, 0.34303458822118515, 0.2714838187456382, 0.6902019172087734], [0.2065130146460351, 0.6092967956859262, 0.9846843685381581, 0.23768636629281426, 0.24457214027787832, 0.7284012786427921]]]])
      mfac = np.array([[[[0, 1, 1, 0, 0, 1], [0, 0, 0, 1, 1, 1], [1, 0, 1, 0, 0, 0], [1, 0, 1, 0, 0, 0], [1, 1, 1, 1, 0, 0]], [[1, 0, 0, 0, 0, 0], [0, 1, 0, 1, 0, 0], [0, 0, 0, 1, 1, 0], [0, 1, 1, 0, 0, 1], [1, 1, 0, 0, 1, 0]]], [[[1, 0, 1, 1, 1, 1], [0, 1, 0, 0, 0, 1], [0, 1, 0, 0, 1, 0], [1, 0, 0, 1, 0, 0], [0, 1, 0, 1, 0, 1]], [[1, 1, 1, 0, 1, 1], [0, 1, 0, 1, 0, 1], [0, 0, 0, 1, 0, 0], [1, 1, 0, 0, 1, 0], [0, 0, 1, 1, 0, 0]]], [[[1, 1, 1, 1, 0, 0], [0, 1, 0, 1, 1, 0], [1, 1, 1, 0, 1, 1], [0, 1, 0, 1, 0, 1], [0, 1, 0, 0, 1, 1]], [[1, 0, 1, 0, 0, 0], [0, 0, 1, 0, 1, 1], [1, 0, 1, 0, 1, 0], [0, 1, 1, 1, 0, 0], [1, 1, 0, 0, 0, 0]]]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank4_fullMask_singleIn_taggedMask(self):
      fac = np.array([[[[0.28474971287125417, 0.3535341676440239, 0.4774280576565706, 0.4632930495890606, 0.331284252186292, 0.10144634065927083], [0.42700476913377283, 0.4505276590513678, 0.30026380490747706, 0.495468761282964, 0.24419752313126009, 0.29945259494027443], [0.24839709846592695, 0.7688306516572189, 0.22057011438166474, 0.2915285884688682, 0.3551782617471212, 0.015801446334563596], [0.807655344982925, 0.9171300786969879, 0.3974029951403697, 0.7060197354420393, 0.5704320264898329, 0.323118530480213], [0.5183855803881253, 0.4192654657751137, 0.5822150325343173, 0.6177536338777756, 0.5740629323483011, 0.024453997287590123]], [[0.027594679938599875, 0.7780782951062702, 0.49908048724516885, 0.5926774673545643, 0.5296609298009975, 0.5594624798782792], [0.42476284797663855, 0.6306903382111915, 0.19814078850547612, 0.0244428994675292, 0.9859352928588756, 0.524683455888915], [0.9186949895348747, 0.5159336720617776, 0.5428009865976677, 0.6743695655741537, 0.3791971135274891, 0.3270715636415854], [0.7191297912925987, 0.00028850400372371077, 0.4672191412019715, 0.48292410304089106, 0.01470833861127141, 0.43249655467622283], [0.15027590152337222, 0.09018899344987963, 0.6176797811082676, 0.012829583085602447, 0.5870240436444062, 0.8669301057523718]]], [[[0.7734660647575453, 0.5977289068696848, 0.16841656065070343, 0.2875375984550894, 0.06209026367760018, 0.2732504945898999], [0.5793172323463114, 0.4148652583352125, 0.20608491363824666, 0.6423420479039753, 0.14247817088202652, 0.5797840189242851], [0.09121252992944162, 0.9621839016193503, 0.12279959029478826, 0.9638908829085734, 0.8541001164653118, 0.23211520780369377], [0.6908390954165956, 0.5282398240990397, 0.5274066715246831, 0.6538584422414843, 0.09270895474518093, 0.09075670470457542], [0.8404184626794099, 0.32512827676667977, 0.8211496312539296, 0.009585489301437433, 0.652544301582929, 0.6632921610894963]], [[0.9602242981185977, 0.6196102719281482, 0.36161380585720415, 0.12798908178646973, 0.5674796273517891, 0.11519792556499031], [0.9248079948628926, 0.8070354572323817, 0.9183047550312536, 0.6348294266550367, 0.26153826018315574, 0.913197762652338], [0.6885599272886335, 0.6487866250389761, 0.15726693518756252, 0.1675498669126806, 0.13892387996520517, 0.5143769242308446], [0.9308400080291018, 0.1470899618863839, 0.07172285990684879, 0.13914373658802281, 0.5368944735347689, 0.816331448252284], [0.3889730804817272, 0.008233098949246953, 0.23100490223688785, 0.20247025896657278, 0.6621645183343428, 0.7505308596896887]]], [[[0.596180891851228, 0.06625459962956137, 0.1633396663352652, 0.7125990220374819, 0.35496196891187815, 0.9403605099493082], [0.33537375332405495, 0.5772396392955728, 0.3297650901321878, 0.3240885255378827, 0.4537829247713565, 0.4212185295737466], [0.5637235799778947, 0.5076671497095357, 0.017332067398027284, 0.8808153649907975, 0.9526906790716846, 0.42115639674036254], [0.6768837393765627, 0.7467610368881578, 0.05888601997019305, 0.1597838731573058, 0.9762322055506684, 0.606386417972246], [0.8348414621133695, 0.4192892480925702, 0.38862472427701034, 0.35128052068453863, 0.7579335321954391, 0.05844316725785392]], [[0.8215185830656196, 0.5095633626237966, 0.5192186109169252, 0.6936542962493606, 0.092203977525662, 0.05548357237013324], [0.18562983148218226, 0.44389638999338643, 0.937731309113861, 0.6443737905400112, 0.46915054513334753, 0.952469497402146], [0.07231331078440151, 0.8714030712843551, 0.07378328352023844, 0.6528123955287729, 0.6801258544133157, 0.5189953689099828], [0.08383357460901963, 0.5570925819053345, 0.42260156215807, 0.052814970334472466, 0.7715018032477665, 0.12863011600533003], [0.7811183203385312, 0.4680999697218089, 0.3248622129468006, 0.9877537708069157, 0.8796167444426422, 0.34701908671059445]]]])
      mfac = np.array([[[[0, 1, 0, 0, 1, 0], [0, 1, 1, 0, 1, 1], [0, 0, 0, 1, 0, 1], [1, 1, 1, 0, 1, 0], [0, 0, 1, 1, 1, 0]], [[1, 0, 1, 1, 0, 0], [1, 1, 1, 1, 1, 0], [1, 0, 0, 0, 0, 0], [0, 1, 0, 1, 1, 1], [0, 1, 1, 1, 0, 0]]], [[[0, 1, 0, 1, 1, 0], [0, 0, 1, 1, 0, 0], [1, 1, 0, 1, 0, 0], [1, 0, 0, 0, 0, 0], [1, 1, 0, 1, 0, 0]], [[0, 0, 0, 0, 0, 0], [1, 1, 1, 0, 0, 1], [1, 0, 1, 0, 1, 0], [1, 0, 0, 1, 0, 1], [1, 1, 0, 1, 1, 1]]], [[[0, 1, 0, 1, 1, 0], [0, 1, 0, 0, 0, 1], [1, 0, 0, 0, 0, 1], [0, 0, 1, 1, 0, 0], [1, 1, 1, 1, 1, 1]], [[1, 1, 1, 1, 1, 1], [1, 1, 0, 0, 0, 0], [1, 0, 0, 1, 1, 0], [0, 0, 1, 0, 0, 0], [1, 1, 1, 0, 0, 1]]]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank4_fullMask_taggedIn_singleMask(self):
      fac = np.array([[[[0.6284052257540734, 0.8561996202440656, 0.19316271602541124, 0.3372364090659762, 0.8711274636323109, 0.9049123445681732], [0.8407360799120847, 0.5735328358867803, 0.4302472151492772, 0.8850560520473593, 0.6853366429969759, 0.47076667636013325], [0.665056484716661, 0.7836096099628395, 0.22377721857801158, 0.5900419395002415, 0.22709681932109682, 0.20695167692798544], [0.33799027739549226, 0.5935483767183248, 0.188665230269623, 0.7820141205006766, 0.09554446902276859, 0.4651359336933878], [0.02587026795784475, 0.5177338171044166, 0.6301948563599605, 0.15817918561970212, 0.6887274739824443, 0.8093359806556251]], [[0.9413516145972828, 0.7132245437492433, 0.07078646109346254, 0.6142384883280997, 0.7696610514191871, 0.7023925235597231], [0.8841928229908906, 0.5820287194272579, 0.04497578637130928, 0.26214851651963555, 0.23516398099194913, 0.6682659683742007], [0.3464212982617677, 0.02480800803648131, 0.9985302336627655, 0.4759238935599748, 0.2488042617181645, 0.7586105809426085], [0.1312471286046052, 0.16863935778109818, 0.7771170973967171, 0.04239252295703966, 0.8801248885849302, 0.26132206761418164], [0.7066209390173489, 0.9968795095110092, 0.5998992369366923, 0.5762513838105194, 0.9854629125778839, 0.527383557531429]]], [[[0.28885342846254647, 0.860269940419838, 0.16105006090155216, 0.2868090551703735, 0.0032647931379355954, 0.4095847719829665], [0.7737670183346086, 0.5518711776478875, 0.27175318688289396, 0.8927361595449447, 0.34178625535127305, 0.7104130149723888], [0.10012756313759308, 0.22707438674859814, 0.6541773977186218, 0.11093150948978625, 0.5377933234714646, 0.11708956967519579], [0.7034455319936345, 0.6298564468194509, 0.5767520216283096, 0.12679311965723894, 0.6830326661189925, 0.23782095357718647], [0.5871869500449636, 0.5335968651652621, 0.6602875540723813, 0.03245011229139294, 0.43631460262301336, 0.43599042147725986]], [[0.12656290217591715, 0.9797320506266635, 0.10382611543369535, 0.01174922275928214, 0.29134840675850926, 0.9425380841368799], [0.44753819443290277, 0.8407810630692446, 0.06337224588120405, 0.31556732031531054, 0.11122312539641421, 0.3264536896365403], [0.2449496699777277, 0.6592399876581244, 0.49866537034439806, 0.6102539910032074, 0.08837452918266364, 0.05734568739802637], [0.5283229076154792, 0.9980006573631596, 0.5331180579619603, 0.7382789375925712, 0.0762410704670955, 0.8753735978066292], [0.5471646312006674, 0.6005384999771388, 0.6231976130105781, 0.45072535197499, 0.7001022590027505, 0.7336715851277552]]], [[[0.26639647219842655, 0.5684186034505964, 0.5907422672622088, 0.16499271767033097, 0.06684681289696792, 0.9211311898672873], [0.29845014962535743, 0.32348982762959055, 0.7729112331078295, 0.24849678763965577, 0.23128636868056673, 0.9359182826559339], [0.2635493837573125, 0.20057959222961108, 0.8178385816852145, 0.2887977292839575, 0.32288940247951836, 0.8233591148745156], [0.04289339389648017, 0.17510094351542782, 0.731543058529446, 0.6421692653751091, 0.5458810518045735, 0.9965277134823539], [0.7697193003124405, 0.0637102931422372, 0.7414824763465548, 0.20920798445067346, 0.8782772426860077, 0.9047379287829074]], [[0.09820530094516933, 0.23369886312183308, 0.5484404848380423, 0.04227201433920047, 0.909940549162931, 0.5127816533929996], [0.7288170096777803, 0.9597966358561947, 0.8137410042090313, 0.4612371702794579, 0.9556234425193119, 0.5639231984662784], [0.2934446071702205, 0.6143957232788751, 0.3670958041412067, 0.3607810923067074, 0.3831831146234057, 0.9143804285555224], [0.5764217655683613, 0.04287132701143992, 0.2341241339894784, 0.20779834532798924, 0.09948307607809215, 0.8274144957198956], [0.4837516330397762, 0.7350038122846315, 0.7998056830007183, 0.7765827173955392, 0.9684715776169768, 0.5234332022514278]]]])
      mfac = np.array([[[[0, 1, 0, 0, 1, 1], [0, 1, 0, 0, 1, 0], [1, 1, 1, 0, 0, 1], [0, 1, 0, 1, 1, 1], [0, 0, 0, 1, 0, 0]], [[0, 1, 0, 0, 1, 1], [1, 1, 1, 0, 0, 0], [0, 0, 1, 0, 0, 1], [0, 0, 1, 0, 0, 0], [1, 0, 0, 1, 1, 1]]], [[[0, 1, 1, 1, 0, 1], [1, 1, 0, 0, 1, 0], [1, 0, 1, 1, 0, 1], [0, 1, 0, 1, 1, 1], [0, 1, 0, 1, 0, 0]], [[0, 0, 1, 0, 0, 1], [0, 1, 0, 1, 0, 1], [0, 1, 1, 0, 0, 0], [1, 0, 0, 1, 0, 1], [0, 1, 0, 1, 1, 1]]], [[[1, 1, 0, 0, 1, 0], [1, 1, 0, 0, 0, 1], [0, 1, 1, 1, 1, 0], [0, 0, 1, 1, 0, 1], [0, 0, 0, 1, 1, 0]], [[0, 1, 1, 1, 1, 1], [0, 1, 1, 0, 0, 0], [1, 0, 0, 0, 0, 1], [1, 1, 0, 1, 0, 0], [1, 1, 1, 1, 1, 1]]]])
      a = Data(0, (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank4_fullMask_taggedIn_taggedMask(self):
      fac = np.array([[[[0.3510644688440855, 0.5183645515249149, 0.4633600225387745, 0.644816248725324, 0.4749987140124231, 0.00599046605614495], [0.0344181213945054, 0.6472575513716537, 0.12359801304418139, 0.1368897395767028, 0.703842544714489, 0.09790937702412872], [0.31394948599186445, 0.012082752446857192, 0.05659031367614631, 0.438011989057035, 0.9773888839787365, 0.8374424822237174], [0.9813420655514463, 0.8780005496043441, 0.7123356486452734, 0.3406250745326369, 0.3061335174951726, 0.6321088459185739], [0.8009546859127842, 0.7127835758948141, 0.11674617423700062, 0.6002134942645535, 0.22570162306661046, 0.13469524466492488]], [[0.37118765263467646, 0.6984769718696697, 0.8272351203103447, 0.19327107197190496, 0.670712348138681, 0.553152830728181], [0.40927220995079994, 0.2956959165293388, 0.7596818866539011, 0.069039896728179, 0.5477669391964974, 0.26916340062009037], [0.02653479505518752, 0.42659160329954515, 0.10345176405864642, 0.6005756997511319, 0.8492173601860834, 0.7431678822547608], [0.752435072428836, 0.6431993748370979, 0.30534793490715595, 0.3533494820892742, 0.7522554751057295, 0.7867109440371933], [0.3134393974364986, 0.3017419299090146, 0.27675762276223326, 0.5487554594890008, 0.8045769009449489, 0.29125035061893345]]], [[[0.9347834302143185, 0.2205310833088142, 0.34455183915327536, 0.1877528281689138, 0.1633666864673371, 0.039930948814584144], [0.7376739252270252, 0.24970351099579735, 0.06306173625983325, 0.6394857958371383, 0.21733594133500345, 0.7247606189134576], [0.15037073220201114, 0.32040901610034267, 0.971589072835442, 0.6409611165397366, 0.6134939442052417, 0.06422719038455305], [0.7175854180865147, 0.07105256733886145, 0.8134513416332522, 0.4914896158354646, 0.4163065917305234, 0.4810717036582882], [0.2806998477217487, 0.8282001731280783, 0.8512330355909834, 0.8829745498048875, 0.8561777272620271, 0.7195440627198264]], [[0.5998825284998051, 0.2080987575992974, 0.16820984037411513, 0.8688722978230603, 0.11764191104666621, 0.513378623510173], [0.5701378093685312, 0.7980121225546761, 0.2243883091253096, 0.8806146470112469, 0.6244186854013817, 0.24712739578949028], [0.03970938144506475, 0.9237832633855682, 0.4345449876857559, 0.2503427049208118, 0.9128124740284069, 0.4258522214212115], [0.11413733339024335, 0.9839784228852547, 0.10942423314649408, 0.6359720158634314, 0.7416784292057219, 0.425583133426237], [0.10305713160589969, 0.053980287141250916, 0.4295234053883418, 0.2958197764984395, 0.8102583711674611, 0.05329861529732138]]], [[[0.7764405621565171, 0.1019553103930867, 0.12170700495885156, 0.8383172403145929, 0.6728673470633261, 0.09039693787664027], [0.36326283135320037, 0.13173481271916854, 0.3055019360941894, 0.9423224913052334, 0.3104569040620868, 0.4476284441055254], [0.7890304466548387, 0.5997753692280207, 0.8097528987122212, 0.6371857147780207, 0.7309377778487307, 0.8885665727048369], [0.13900635580425458, 0.20204245195181925, 0.3653831706281282, 0.47325339199565875, 0.9651340716042628, 0.02657224625121657], [0.053776114567332556, 0.8950282334034564, 0.3494985303860769, 0.09725481502634825, 0.5933688623717693, 0.9963595069898447]], [[0.2013874490559756, 0.013296618279048489, 0.41086119953941025, 0.1422854589121917, 0.36418669764150935, 0.0477716340076616], [0.5443021489005228, 0.37254487430034666, 0.46969334069511837, 0.14723478203443918, 0.15668490948397573, 0.49914556660780807], [0.058046977863655536, 0.6498733908963004, 0.2893753346195944, 0.9972343445692876, 0.9797959399102546, 0.4907591145280673], [0.4497069230511578, 0.9824795145949844, 0.9821623360625206, 0.49695147290771435, 0.9421030227367252, 0.3254897222228912], [0.15395556929098597, 0.8860267634664146, 0.518236422749814, 0.22556584753528208, 0.2505912470498203, 0.5761709357875672]]]])
      mfac = np.array([[[[1, 0, 1, 0, 1, 1], [0, 1, 0, 0, 1, 0], [1, 1, 1, 1, 0, 1], [0, 0, 0, 1, 1, 0], [0, 0, 1, 0, 1, 0]], [[1, 1, 0, 0, 1, 1], [1, 1, 1, 0, 0, 0], [1, 1, 1, 1, 1, 0], [1, 0, 0, 0, 1, 0], [0, 1, 1, 1, 1, 0]]], [[[1, 1, 0, 1, 1, 1], [1, 0, 1, 0, 1, 1], [0, 0, 0, 1, 0, 0], [1, 0, 0, 0, 1, 0], [0, 0, 1, 1, 0, 1]], [[0, 1, 0, 0, 0, 1], [1, 1, 1, 1, 1, 0], [1, 1, 1, 0, 0, 0], [0, 0, 1, 1, 0, 0], [1, 1, 1, 0, 1, 0]]], [[[0, 1, 1, 1, 0, 0], [0, 1, 1, 1, 1, 1], [0, 1, 0, 1, 1, 1], [1, 0, 0, 1, 0, 1], [1, 0, 1, 1, 1, 1]], [[1, 1, 0, 1, 0, 1], [1, 1, 0, 1, 0, 1], [1, 0, 1, 0, 1, 1], [0, 1, 0, 1, 0, 1], [0, 1, 1, 0, 0, 0]]]])
      a = Data(0, (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank4_scalarMask_singleIn_taggedMask(self):
      fac = np.array([[[[0.24340244461877392, 0.003991636543154264, 0.18416071176898752, 0.36278901168532085, 0.9262578823732283, 0.24880345089298495], [0.9213906359188394, 0.6955465167814712, 0.7295625212576448, 0.3968615857374044, 0.28503833736741424, 0.28568573477420933], [0.9529134028604478, 0.23789011148134664, 0.37389597131335883, 0.37768129809067374, 0.25775114783676434, 0.31329912377889035], [0.4113662071468738, 0.1630457030672493, 0.9138709132123679, 0.16966658550088964, 0.6467781403414309, 0.4961942390786247], [0.47381007027221445, 0.6324934317385829, 0.6542127484097704, 0.8034494444059535, 0.08817026554664875, 0.9590672965188575]], [[0.6686295669153601, 0.7064183681590886, 0.18745674914112798, 0.33350212648382016, 0.6828974782029811, 0.38714026060627293], [0.20620768527278122, 0.6600925097895018, 0.8301141522437175, 0.2390903392550926, 0.3312951249639101, 0.46539785769951614], [0.6116610623763291, 0.2554076178426109, 0.879259191507673, 0.5719973587603542, 0.26596743212484375, 0.023164009069908786], [0.9443058406328223, 0.48461211809887383, 0.5814642859694314, 0.9429112694198492, 0.8373785192720604, 0.25549426698367816], [0.3133092649756002, 0.5649226999611084, 0.6312753881336864, 0.7824811698477031, 0.4668493090789264, 0.940835912145762]]], [[[0.8349404033158234, 0.8794923896958541, 0.7249574790093885, 0.22710232413542597, 0.7953353175497335, 0.33785889259432456], [0.8178231819661628, 0.5916360939504626, 0.7998776354366094, 0.4647071218417502, 0.5211495676885265, 0.5103958942125071], [0.919655352627235, 0.7167777329480843, 0.45536615846379347, 0.8689448584329278, 0.9424275675184963, 0.529741601176882], [0.5346717537960046, 0.18361371843664775, 0.9891817638711139, 0.35182430529300845, 0.020281096895918616, 0.33253255782061564], [0.3894266751116703, 0.41581484251522405, 0.6497192140553895, 0.5302452473182707, 0.19932810611822993, 0.5365856270978093]], [[0.14162717063942953, 0.3680049637302537, 0.11662036768783579, 0.497924471598683, 0.8529349395063882, 0.024381916915324764], [0.4365386260186206, 0.888991159913084, 0.8990111732283905, 0.06645529022166896, 0.23746471290253868, 0.42145705263275346], [0.7953894238129054, 0.28881917685861436, 0.2171492457244001, 0.38145645523788285, 0.4341611198480255, 0.797496648962407], [0.29490089531667296, 0.5613328807810927, 0.3111613910026024, 0.47719572128090004, 0.0048269556400570846, 0.7275357242628578], [0.2071888092736358, 0.2107463612035867, 0.9471960957384761, 0.39058257380984374, 0.13111726477176, 0.6684551726507185]]], [[[0.9550141633398366, 0.2216535244872001, 0.7683399936581771, 0.7408190380620381, 0.10208443296786229, 0.21054891863795955], [0.2808826595907715, 0.15491651835236453, 0.6823731637351869, 0.08683004816384432, 0.882123690346417, 0.7039308501763134], [0.3738303217816191, 0.05246197224873572, 0.49434185454690804, 0.9777014607770798, 0.03518570923608255, 0.2787644387255954], [0.6999566918032673, 0.027372375357986645, 0.3209685148748115, 0.39269707459396586, 0.4303028892646239, 0.48502242502039383], [0.5723048541475331, 0.9218396346466207, 0.23689337904795482, 0.22518086888474564, 0.7491335027331667, 0.0869008343674953]], [[0.26880869526127793, 0.0008066874730127127, 0.09777097537349322, 0.6626360272304919, 0.6715956772728239, 0.46286340401642057], [0.7220247053200817, 0.43079032988851707, 0.6158162292957735, 0.5121846986816639, 0.9832010146500367, 0.27065072956328173], [0.40189979187408387, 0.6715715316675711, 0.009885983517540775, 0.9039884989006461, 0.8073219125494864, 0.7933553862480189], [0.48411251648081843, 0.318592685542217, 0.29162128593656567, 0.26959236745156767, 0.023708763615697737, 0.9682669231262091], [0.370097205523383, 0.9123312281155945, 0.9584284058060969, 0.25170456661630347, 0.02240861090504631, 0.6709429346885786]]]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank4_scalarMask_taggedIn_taggedMask(self):
      fac = np.array([[[[0.23255745282243134, 0.25047967477819555, 0.6108146394241556, 0.13385382302604965, 0.13706718723468414, 0.035870123975077384], [0.6246682909880291, 0.836628524775755, 0.5162237988121678, 0.2263212870966289, 0.3932178434508955, 0.7137301287086041], [0.7900326645006635, 0.6877094372486755, 0.08816645194203776, 0.4387899970329403, 0.5316140968721685, 0.45254925198950136], [0.06812485055764994, 0.7389370886889369, 0.5714257830488181, 0.6047955773425038, 0.27078156831807876, 0.6941573119735807], [0.209989381510422, 0.6333956563974882, 0.4393853249053984, 0.42367337314008047, 0.6748957279189746, 0.32338445580230313]], [[0.2134065392799016, 0.344237642144424, 0.006391732823146423, 0.34032022273305673, 0.8948537900601655, 0.9550600240789674], [0.6676388813796433, 0.5897822335186457, 0.665173408410216, 0.5125134002903305, 0.972581534440156, 0.8803074309213411], [0.6859072537342642, 0.24635080019185462, 0.6864759628868313, 0.9230939191552374, 0.34980486904078456, 0.6234666425070711], [0.19941866735162594, 0.38857567889311384, 0.6804329681480833, 0.8232548806451583, 0.7732278079442935, 0.36739061111911364], [0.664959603547849, 0.354321001742274, 0.4978605827740591, 0.8077836084911515, 0.8860377462245522, 0.3558117956286826]]], [[[0.7022481883456659, 0.47531348642788684, 0.1248278950933227, 0.9954658477193026, 0.4956944965275675, 0.3961830612592562], [0.4452609658645168, 0.191864815188864, 0.6933937596113584, 0.8645633537850476, 0.9983598663949573, 0.7085613415980516], [0.24911687532163918, 0.22078545664164073, 0.8751023042377495, 0.4835181837843944, 0.6275907424582854, 0.646833076432497], [0.9611870465060238, 0.1817683233288383, 0.2531600550840636, 0.6850261894292358, 0.5411832954496947, 0.4927044097007134], [0.9074012181760704, 0.8572672440433726, 0.9477982289219312, 0.45441979253778286, 0.3737391828645419, 0.11409904093648116]], [[0.8380059026984883, 0.388963815716435, 0.31715278996373675, 0.7318626863954599, 0.6701674771302204, 0.5912434291529112], [0.17713082907853328, 0.850518305611716, 0.6038880861452677, 0.6921269193732074, 0.5907501856684764, 0.30281670787792814], [0.9171221669289326, 0.5918122893462844, 0.803345345695726, 0.6615609162170348, 0.5361218657437894, 0.19817134402689784], [0.025833264331518335, 0.5391597929790176, 0.5069982470753921, 0.6946324386296594, 0.8316892306640529, 0.20827400975469235], [0.556647989881616, 0.04869571269138173, 0.54793842219008, 0.9680026989139885, 0.6354172990204311, 0.6437668938278768]]], [[[0.6102421024244974, 0.9752237946192167, 0.14990102113175385, 0.6895695218864366, 0.6549824880011706, 0.45691096911884943], [0.6448913705517718, 0.06201186762968036, 0.7188583178468533, 0.6428477397286124, 0.13196843867694052, 0.2978491557975472], [0.037552529064665996, 0.25335017800285786, 0.762959546507413, 0.9512068299376989, 0.41610961223077725, 0.28259515970236326], [0.876868865490853, 0.9546721456496616, 0.45308978762576624, 0.9210016807509445, 0.3498701613630698, 0.9400827751646211], [0.056991230186118624, 0.15185346244516318, 0.9708614689834713, 0.674758101061766, 0.6609086612573717, 0.4109303912176313]], [[0.3041081599987733, 0.13222381804994032, 0.05960230960866997, 0.6026614014752757, 0.20654611198355444, 0.6510247782715143], [0.7330423411648939, 0.769160239751075, 0.4296736300615197, 0.68245732111301, 0.8357937310534077, 0.9533795456954742], [0.7054475239994354, 0.3295096422500373, 0.9114725768144195, 0.16696638267056285, 0.87747809475601, 0.5468752078600541], [0.6262311830041171, 0.6515771382576069, 0.2601125867736055, 0.5487502389168921, 0.7221886149469052, 0.6283978003189019], [0.512448473907148, 0.8509981872455821, 0.03129729886682475, 0.24202853972259597, 0.9885699268147726, 0.39460010787806177]]]])
      a = Data(0, (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank4_scalarMask_expandedIn_taggedMask(self):
      fac = np.array([[[[0.49268598888585025, 0.3107101466027983, 0.7061659759301252, 0.4820591791469292, 0.07683063390615608, 0.5325488700723157], [0.04420855135858914, 0.005560881873779233, 0.2621433669380535, 0.18504134052153887, 0.08160568376292798, 0.05670023480677733], [0.9043588903176409, 0.02919575399852048, 0.8949611254294384, 0.06157655315550026, 0.6971158737342358, 0.32475618176284204], [0.693985711724827, 0.6806368114121419, 0.12133059650646105, 0.8424244012788182, 0.06906627974177604, 0.5909022262062444], [0.749780948088381, 0.6276044097062183, 0.4194183879183876, 0.5327305589966979, 0.7992340156382242, 0.8437319587286426]], [[0.7986671841237303, 0.7068851820140846, 0.7798276950315332, 0.271750596681249, 0.2340016608189216, 0.5385041982464216], [0.47300429365011565, 0.030977476416714778, 0.9893701739125189, 0.980388798864986, 0.7901958331804076, 0.06099877654457009], [0.18486426595119054, 0.7755597554114215, 0.5303371925091425, 0.22325508652428327, 0.4448756612753437, 0.7366564463964204], [0.564593609538506, 0.9490316554120761, 0.047857980350849094, 0.13300129132175575, 0.2420778330906752, 0.6997330992658103], [0.38658386364074027, 0.3509942043170078, 0.3647893015477538, 0.059665228839896556, 0.26597576119988176, 0.4673654378831553]]], [[[0.9533594966837214, 0.7392134654231676, 0.7659267893682623, 0.5396725792383511, 0.3689249107526186, 0.5231012217727014], [0.10742387207125181, 0.11118824908596836, 0.5615690304371522, 0.6988635405789453, 0.6051048574031144, 0.7289252983498301], [0.8953702625860572, 0.6750002460235109, 0.6843324412730589, 0.575520781395198, 0.8898269527106561, 0.43522983090931466], [0.625507204273611, 0.49782954675026125, 0.6315157861495725, 0.11773027425636073, 0.2608588679097419, 0.5488949765950472], [0.16047842410972868, 0.8233887524412568, 0.6578211756989466, 0.7094676538003308, 0.9174125856873595, 0.17548188804940434]], [[0.3213785659504471, 0.7337994983766531, 0.7738398353364391, 0.8533961953184378, 0.23267404859169627, 0.7799730385017423], [0.38842427759408815, 0.20644454852621452, 0.806740660797546, 0.33520173007420273, 0.5778771764793195, 0.5642534246795946], [0.8201331716393573, 0.6404245132053147, 0.3898093104772955, 0.4106845523916336, 0.20366560906921238, 0.3096562969780384], [0.36474681520156327, 0.8198193877254853, 0.978015743111213, 0.5623609921746932, 0.17242057917686826, 0.9648170870625762], [0.8624604314386543, 0.0910172289057597, 0.891129692059383, 0.14722813652870126, 0.10640804598673248, 0.7617703280418022]]], [[[0.7822286625326228, 0.05831617188148874, 0.6379323829392524, 0.13391418863392768, 0.04741041905042376, 0.7889535355003038], [0.8799983595265761, 0.2803028506985983, 0.5539974564824891, 0.8441759777387376, 0.9419454735304278, 0.9203766165016322], [0.9977230238258449, 0.4156371114418579, 0.013787353150302883, 0.018645988560957116, 0.5519134450149124, 0.07421377879062552], [0.7458372669217389, 0.8055447342613018, 0.3176387414916876, 0.20674721162297582, 0.19665852554929342, 0.3836280863569309], [0.7809033316854367, 0.6985702117999811, 0.74328974899762, 0.2817721743533481, 0.4424980551413057, 0.5460556284628488]], [[0.26996633820979554, 0.7183960911874062, 0.42614903168454943, 0.12764172325688905, 0.6193634799399617, 0.8952610971231936], [0.8012311271226786, 0.1511083297224689, 0.16687997344831917, 0.8842063607612395, 0.6388262962329668, 0.8519129589853756], [0.5533202619319474, 0.7149786370337425, 0.2304269702354962, 0.7244507573686573, 0.3444769707705142, 0.79194046291552], [0.3331704477710953, 0.7546398972635849, 0.8254394038899514, 0.7191601148347769, 0.051136068073240004, 0.40962172698688026], [0.8719645795160615, 0.6867180854696139, 0.36117632258806587, 0.35424351852444336, 0.497175577225738, 0.03275333410739223]]]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())


   def test_CopyWithMask_Rank4_fullMask_singleIn_expandedMask(self):
      fac = np.array([[[[0.8454062762825381, 0.7504779003885301, 0.9217068295029763, 0.4917546558187038, 0.5450124055971526, 0.4112335749528959], [0.7055528972822482, 0.1283096400630186, 0.40465419678109016, 0.3339824379523346, 0.8863830119179582, 0.7638264914191869], [0.3425065449703515, 0.375738251424908, 0.6054308741246542, 0.7615095381822662, 0.09158033127368925, 0.3412946409886265], [0.23123337103600117, 0.8899679327382184, 0.5388553661250899, 0.25001110010864847, 0.10545415138814707, 0.8015134274389525], [0.8476306125719111, 0.599090460728146, 0.2815754272641521, 0.4103959576170386, 0.04002675603912731, 0.07322812124605527]], [[0.8184016897496894, 0.8410157041949969, 0.45709308508882585, 0.40795092103581443, 0.28707894765277286, 0.2586827353643443], [0.4238180832417302, 0.41675402596295297, 0.31113217167710205, 0.22246031290545298, 0.7846721632380986, 0.2132160801135322], [0.46310540784138254, 0.7623906782444497, 0.7517860004173377, 0.49106190984038933, 0.6491604856755372, 0.20780179321012748], [0.7274418813103185, 0.4525374355083256, 0.0497427003146339, 0.4569990226295071, 0.9790529417160692, 0.7480806049061649], [0.8589844993345834, 0.7964583405156614, 0.5861977895993796, 0.7629669635142223, 0.5616863586354357, 0.2216172057711734]]], [[[0.7641016689014933, 0.7288189235675914, 0.30228612201109983, 0.2298441990904967, 0.3845465963307573, 0.3993789889977648], [0.36967933257739105, 0.9268936958764797, 0.1524219663597156, 0.7444132295174904, 0.7402272231128104, 0.36403964139641487], [0.3527165829756569, 0.854641107660538, 0.6646046904700037, 0.22515837222120672, 0.018955490181338797, 0.792022216586335], [0.7158329553351587, 0.6910421252553653, 0.516381394721964, 0.4756475531131029, 0.295099773978709, 0.6147327832733148], [0.7172373978776594, 0.038028025442080726, 0.9314818787918023, 0.2901668799579058, 0.8079482376505619, 0.5955668454707365]], [[0.11586283423582933, 0.7739961358765668, 0.397325133589028, 0.9820758924368529, 0.7678006747091599, 0.4641560498658174], [0.18170833840858147, 0.198335519631488, 0.8837231163463453, 0.3107607119453001, 0.7420480305989258, 0.7272901267918069], [0.9620960145351495, 0.5169927631656737, 0.9280872080651508, 0.17553111840117375, 0.9542894571808797, 0.934348936060397], [0.971460512150884, 0.5966779859792725, 0.26867306282167214, 0.16718496583395104, 0.757014159901295, 0.32054058657591833], [0.03868619736434986, 0.9657989010039688, 0.15580977406961027, 0.6572696429923183, 0.25728064028274655, 0.13613810911617685]]], [[[0.7854451754531789, 0.9709851469821951, 0.24427694065480055, 0.7642242835915412, 0.3469680695456534, 0.8504813556326348], [0.2189113944614155, 0.29936770859158246, 0.06208148929181978, 0.5956356165299023, 0.2699309798836502, 0.636329201444455], [0.9947148474604315, 0.3231243700844687, 0.1902544708530597, 0.8996882841166632, 0.012151534122572749, 0.7663236649200924], [0.10942374084001127, 0.6832535616971418, 0.17467572429160128, 0.8068564243238228, 0.6131264839363639, 0.7740036797366973], [0.8343423888874891, 0.4782613429058685, 0.45662834148274734, 0.027397044214439892, 0.05550137700410773, 0.26196385700249514]], [[0.3298885608708306, 0.3737357735693917, 0.6077566611212749, 0.9086067727112499, 0.40094066505554826, 0.31614230568989543], [0.2773145203152416, 0.6435606374467427, 0.2621174976879833, 0.5972422354548184, 0.9192454052623633, 0.40646304118486176], [0.9731859685571184, 0.7856418383694926, 0.565809443060594, 0.6647270746610757, 0.6501901560361719, 0.06912473917802853], [0.3276836417456621, 0.04705012523358065, 0.9229059474029985, 0.1743121580832918, 0.26113799261982995, 0.5050736884118043], [0.3320215221532071, 0.3063101669254886, 0.4570660728162679, 0.6305076961031194, 0.25566654649249276, 0.34584319410339026]]]])
      mfac = np.array([[[[0, 1, 1, 1, 0, 1], [0, 1, 0, 1, 1, 1], [0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0], [1, 1, 0, 0, 1, 1]], [[0, 0, 1, 1, 1, 0], [0, 0, 1, 1, 1, 1], [1, 0, 1, 0, 1, 1], [0, 0, 0, 0, 1, 1], [0, 1, 1, 0, 0, 0]]], [[[0, 0, 0, 0, 0, 0], [1, 0, 0, 1, 1, 1], [0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 1, 0], [1, 1, 1, 0, 1, 0]], [[1, 1, 1, 0, 0, 1], [0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 1, 0], [0, 1, 0, 1, 0, 1]]], [[[1, 0, 0, 1, 0, 1], [0, 1, 1, 1, 1, 0], [1, 0, 1, 0, 1, 0], [0, 1, 0, 1, 1, 0], [0, 0, 0, 1, 1, 1]], [[1, 1, 1, 0, 0, 1], [0, 0, 0, 1, 0, 1], [1, 0, 0, 1, 1, 0], [0, 1, 1, 1, 1, 1], [0, 0, 0, 1, 1, 1]]]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      m.expand()
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank4_fullMask_taggedIn_expandedMask(self):
      fac = np.array([[[[0.2060257956315361, 0.8013255610901385, 0.41762289281480014, 0.959542720474878, 0.07928351470321748, 0.8349606606057056], [0.6521548264799737, 0.5746515998295404, 0.891779747992215, 0.266526168250091, 0.10447420578161548, 0.21894827651284987], [0.7275127896075808, 0.7129963614164497, 0.9449647212890467, 0.9358049424362349, 0.18551534184074547, 0.9136168177036778], [0.16927567661536913, 0.09100574806043416, 0.2818954219027122, 0.303259625492711, 0.980673410996367, 0.5676841700349545], [0.733154697951832, 0.8529570424418184, 0.0024401421407691526, 0.46834062255202225, 0.48611636977521644, 0.24063372489120427]], [[0.18405245055721842, 0.2190358447381665, 0.4407020056330839, 0.8383326542655447, 0.0040420804100976815, 0.25197575688105633], [0.37154846175978373, 0.15124600705897306, 0.5290200968093539, 0.7524634617093142, 0.22056702599038625, 0.32440228209888256], [0.35214489236654445, 0.08615100385200003, 0.07365690131545821, 0.9859132761806823, 0.39878861080894756, 0.49220177634209994], [0.7429747313840203, 0.9465239492595147, 0.859041530088702, 0.1651602847464274, 0.7191238741715912, 0.9953852607548793], [0.9221563545218003, 0.43179453715319527, 0.0861173643319123, 0.7674514721868696, 0.26408986433456194, 0.7500927492488179]]], [[[0.3316441597866495, 0.3172529422474527, 0.1368382259064722, 0.3558837137140465, 0.6063267204017734, 0.7454645570107197], [0.6113486510034358, 0.39025536706470576, 0.5365660681758005, 0.03367671563921826, 0.05527702326110728, 0.8098274391161415], [0.29540141040671963, 0.6471121839227175, 0.8383818104244506, 0.5166789010078267, 0.21265209449641587, 0.9074507629386533], [0.39313647139928676, 0.5108403822258787, 0.929013483900152, 0.182277765341892, 0.753193241568604, 0.045955912020305245], [0.9675883781465827, 0.504297854983803, 0.6760081665549764, 0.00935531249841326, 0.10476094472004105, 0.5797833924894261]], [[0.5577113760443241, 0.3455199462344659, 0.8229678714766903, 0.6694602393057525, 0.45272427237375634, 0.6795813255898726], [0.1444927467955308, 0.0753148198697946, 0.2528663964081176, 0.938160343533968, 0.26073262943365827, 0.8610127934647804], [0.38333630038671607, 0.41217643297334716, 0.1700727845524853, 0.04686142588961195, 0.13848291526054746, 0.5344784743627528], [0.6358088141386015, 0.8573313063912318, 0.11727594473018366, 0.5952186218005178, 0.13871344323176849, 0.8447590769518593], [0.5484333035159497, 0.08308062571072983, 0.011105466093853655, 0.12506750723530247, 0.9018641312393388, 0.15676405069851373]]], [[[0.7452539859014125, 0.24928533887562354, 0.7708016166290149, 0.7665982130474543, 0.2784527159571365, 0.4703697406718521], [0.5385806936432703, 0.039835242602980436, 0.40027217253198855, 0.03129945807997281, 0.45646704537622906, 0.1748690894683561], [0.9042982801026912, 0.7579829223063139, 0.7544216247889571, 0.0440232338287071, 0.8918990195888908, 0.2626766450744108], [0.3440307752514159, 0.6678718876076561, 0.867332877217701, 0.030542505397874953, 0.20056497264329032, 0.7700944895576698], [0.25351267219807894, 0.23685653554023567, 0.579346193877827, 0.2754279080990303, 0.3008387909599478, 0.9421802410056782]], [[0.6950599323904272, 0.7728767759504737, 0.9522544917004189, 0.7914401982848924, 0.33512842352399397, 0.21662825197443603], [0.5042266623712208, 0.5035375099203054, 0.9929729373688609, 0.6823493470930893, 0.3929623102437775, 0.9114323462282601], [0.0673326917554441, 0.8796549119941219, 0.56049997324928, 0.08220191882133732, 0.4028738663799629, 0.8927377165078402], [0.28463952842944673, 0.13106330426199098, 0.04511006712869414, 0.6119608942509642, 0.6796385338836228, 0.899442952205246], [0.9890093345360628, 0.4129587230616203, 0.2929147356217561, 0.49147833157193754, 0.8717749256631869, 0.1644431984532011]]]])
      mfac = np.array([[[[1, 1, 1, 1, 0, 1], [0, 1, 0, 0, 0, 0], [0, 0, 1, 1, 1, 1], [1, 1, 1, 0, 1, 0], [0, 1, 1, 0, 1, 1]], [[1, 1, 1, 1, 1, 0], [1, 1, 0, 0, 0, 1], [1, 1, 0, 1, 1, 0], [1, 0, 0, 1, 1, 1], [0, 0, 1, 1, 0, 0]]], [[[0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 1, 1], [0, 1, 0, 1, 1, 1], [0, 1, 0, 0, 0, 0], [0, 1, 1, 0, 0, 1]], [[0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 1, 0], [0, 1, 1, 0, 1, 1], [0, 1, 1, 1, 0, 0], [0, 1, 1, 1, 0, 1]]], [[[1, 0, 1, 1, 0, 0], [1, 0, 0, 1, 0, 1], [0, 0, 1, 0, 1, 0], [1, 1, 1, 1, 1, 0], [0, 0, 1, 0, 0, 0]], [[1, 0, 0, 1, 1, 0], [0, 1, 1, 1, 1, 1], [1, 1, 0, 1, 0, 0], [0, 0, 0, 1, 1, 0], [1, 0, 1, 1, 0, 0]]]])
      a = Data(0, (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      m.expand()
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())

   def test_CopyWithMask_Rank4_fullMask_expandedIn_singleMask(self):
      fac = np.array([[[[0.34431933070651066, 0.2956216526076568, 0.16648074609050156, 0.82693345266232, 0.1886598646125427, 0.28248798669383157], [0.13000418544404613, 0.27835264447282093, 0.8326724977173048, 0.9008300878670021, 0.18043784267073715, 0.6645532376751698], [0.8301350435484512, 0.89081172329241, 0.8957816677493445, 0.8909221883146838, 0.4511061612187628, 0.5389312869579176], [0.4842943266380919, 0.9003398131285711, 0.52984753547433, 0.9232697366824476, 0.6302858968083069, 0.36299624013190224], [0.9710944331342275, 0.13022737682584695, 0.9901786503417818, 0.28871390808333564, 0.18181891211924206, 0.8952670724497613]], [[0.20088205898480527, 0.9100534874328549, 0.7238861435006663, 0.9865749831604277, 0.9676978004061866, 0.7320074621573627], [0.1666093245780782, 0.7467880806119441, 0.5140284401615117, 0.5247802356840304, 0.6931285544008728, 0.32517652694615107], [0.014349430732579327, 0.12854029396263167, 0.10330487667189614, 0.5238086417176597, 0.6856420574328915, 0.8482862525509353], [0.042573659181832646, 0.54500943012482, 0.2723874719985864, 0.7097063222767848, 0.6694943432798381, 0.7015590439447432], [0.8553367263109058, 0.18637013543863745, 0.7819928039929638, 0.09163535981872284, 0.8789396119735952, 0.12590238658563135]]], [[[0.04777188854505365, 0.15809590007614105, 0.36532952031172605, 0.7237415749235547, 0.2113666204923943, 0.6014039050610795], [0.9973581451380272, 0.5248896095125268, 0.93574111392229, 0.6628935070512522, 0.7208469111271132, 0.1020490914975285], [0.36582594721082773, 0.7671801020383687, 0.7147984068959402, 0.11190042832996172, 0.501558903430659, 0.8189982388510242], [0.3133191971908146, 0.31702928884006276, 0.8820133325874993, 0.8248549245606946, 0.032325735815006196, 0.14088914144065168], [0.6398154338530787, 0.6920417608248364, 0.549624645022837, 0.5919312264504032, 0.6508280891891102, 0.15953916023654457]], [[0.47777808994088133, 0.7410551023995438, 0.7643760608100465, 0.015327519187263161, 0.10510204528011635, 0.9255865573952102], [0.9905208655853166, 0.9858997608233748, 0.4459878121810972, 0.5999604877487645, 0.606717448039868, 0.716125962482642], [0.7489343287314179, 0.3060955650580902, 0.4493771176639858, 0.46986991181035465, 0.6460768369385653, 0.0675214741028034], [0.7966378018006339, 0.14807654301890083, 0.15388604129111616, 0.7461798702719485, 0.5515179411653017, 0.3555644867666008], [0.0877865788441059, 0.08108382397211966, 0.8587025507927749, 0.4887663058786137, 0.35587034928391403, 0.6757176153813076]]], [[[0.2365197632849022, 0.13786695168882923, 0.1408657131413642, 0.7121331931108396, 0.2711583132337283, 0.5158341308203517], [0.6115519919761471, 0.4106891183818927, 0.9078555490781101, 0.39998221881857765, 0.7382994523479791, 0.5755656929177119], [0.7619844945457624, 0.39280238003884926, 0.10770903095363493, 0.2503795700209497, 0.6924993812479509, 0.784861473422038], [0.9256349702699108, 0.5338721744928653, 0.07884986150989504, 0.7134493701620046, 0.41572448920891203, 0.9391185567720463], [0.9932627479342792, 0.7080603974381019, 0.44334878589987914, 0.34795387783348364, 0.33885004178878164, 0.38616185910504075]], [[0.6629015430817785, 0.19003931561705056, 0.5062322966457088, 0.7776233916867864, 0.18837502693136943, 0.6347628310383449], [0.7627326440357979, 0.41661577173850484, 0.406903256521605, 0.34459651712648864, 0.06195833125826755, 0.24501238016304083], [0.9969507384585272, 0.7352747648214577, 0.9870852965223819, 0.7757263610550733, 0.7196472753496949, 0.12383007445092731], [0.7067450565717149, 0.770581940279722, 0.7185842145053877, 0.38109808855165905, 0.06374210510016398, 0.3045715933956282], [0.8638868536670794, 0.6482162477774626, 0.5079325951310301, 0.9648849496026032, 0.6676677495882813, 0.5064499873597063]]]])
      mfac = np.array([[[[1, 1, 0, 0, 1, 1], [0, 0, 1, 0, 0, 1], [0, 1, 0, 0, 1, 0], [1, 1, 1, 1, 0, 0], [1, 0, 1, 1, 0, 1]], [[0, 0, 0, 1, 0, 1], [0, 1, 0, 1, 1, 1], [1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 1, 0], [1, 1, 0, 1, 1, 1]]], [[[0, 1, 1, 0, 0, 0], [1, 1, 0, 0, 0, 1], [1, 0, 1, 1, 1, 0], [0, 0, 1, 1, 0, 1], [0, 1, 0, 1, 1, 0]], [[1, 1, 0, 0, 0, 1], [1, 0, 0, 1, 1, 1], [1, 0, 1, 1, 1, 1], [0, 0, 0, 0, 0, 1], [1, 1, 0, 1, 0, 0]]], [[[1, 0, 1, 0, 0, 1], [1, 1, 1, 0, 0, 1], [1, 1, 1, 0, 0, 1], [1, 1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0]], [[0, 0, 1, 1, 1, 1], [0, 0, 1, 1, 0, 1], [1, 1, 0, 0, 1, 0], [0, 1, 1, 1, 1, 0], [0, 0, 1, 0, 1, 1]]]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank4_fullMask_expandedIn_taggedMask(self):
      fac = np.array([[[[0.7660908342679651, 0.8878210858929412, 0.5858277894367959, 0.9026208451996311, 0.5383564873901119, 0.49839748842496767], [0.4610225357974309, 0.634277084282585, 0.04436008952271775, 0.7977539198451007, 0.6874023073782206, 0.5055858716889677], [0.8307049865521776, 0.8499865959818322, 0.6503893404864046, 0.9385494077843018, 0.6150432936445988, 0.9718898380989767], [0.5878984644416054, 0.8044399105534293, 0.4950024550212865, 0.5240023033536955, 0.9561570890006335, 0.8513116581475224], [0.24728200165923997, 0.18317431221231484, 0.4537081559221918, 0.6502685538122246, 0.9053939106503347, 0.7263585429098217]], [[0.01646379530877573, 0.4629027197995004, 0.15829442175900965, 0.28199176629354816, 0.24002494785463235, 0.20404490708667422], [0.6549187960050376, 0.8524472417879952, 0.2862403466500534, 0.13849383235215607, 0.6056643588931598, 0.34694423869807267], [0.21822992984774803, 0.9764551975454274, 0.1832225717093422, 0.649108423927957, 0.780144142147636, 0.1771554790424399], [0.017251365173236333, 0.32271597208530844, 0.6622760077122637, 0.38832680850269585, 0.9628469376502358, 0.04515507001753638], [0.30246344673762693, 0.7009239734677206, 0.8512704273009822, 0.7477501756068373, 0.5625962065816328, 0.6074311622202323]]], [[[0.8238435166388135, 0.7303492021995375, 0.6348169116881529, 0.7719995151996101, 0.9147790596445141, 0.28195693351816564], [0.6233223859040478, 0.4190169458925912, 0.8719706385554703, 0.6238752630622026, 0.13807200072010695, 0.8788206255264227], [0.27789631373564505, 0.8659497358967189, 0.8802019459812982, 0.5446680257024953, 0.3413831936685545, 0.5601476573368255], [0.3860331163110172, 0.271104526662592, 0.6324697667383743, 0.12030027184928138, 0.8190884346147314, 0.5117331442160478], [0.10563735336995461, 0.946830103684979, 0.3953225605129612, 0.7644287372843432, 0.7504293034440621, 0.3884006063988289]], [[0.3150890658874459, 0.14820750791231496, 0.824974261215091, 0.7248581216162522, 0.3698954512751256, 0.7217278318928815], [0.8917914593403268, 0.0076845180073500385, 0.9301819371933201, 0.9385525532108595, 0.4195988051340753, 0.8291832615187953], [0.04232795159819569, 0.22139529980683204, 0.5845624865783858, 0.7255562100446149, 0.23313140532716992, 0.45632672401036223], [0.6746102447259757, 0.5788804186856973, 0.4696081559614933, 0.6406768821937994, 0.9610528484088848, 0.45440744271999634], [0.9546558049998668, 0.7850572836595879, 0.5358987076108322, 0.9668393015934319, 0.40488379608481984, 0.5584108505950139]]], [[[0.6987027045565316, 0.8857789907802783, 0.4012509521392982, 0.05481639788686199, 0.3025292219931117, 0.2634171662062711], [0.41112259556310105, 0.7802348354132257, 0.9588005215486888, 0.8729903039347805, 0.7407751224011698, 0.9949051782310238], [0.31494435792909325, 0.9299748226416473, 0.18019802370758742, 0.4789352691036355, 0.18971132432286175, 0.7550208168184991], [0.11731225314591742, 0.6232228041260098, 0.8921567427341464, 0.9248011257300917, 0.36638643231222745, 0.5124849771472058], [0.023961286622143363, 0.21298801380361754, 0.8768121589763805, 0.200885175629277, 0.527683790089323, 0.7093854383916601]], [[0.3464576293267032, 0.6091585594867465, 0.6742059916445082, 0.34277957007462556, 0.769450819219595, 0.9218144296101813], [0.9708071371576861, 0.4035715171362836, 0.08635185173547277, 0.6395282900756231, 0.0746178635619863, 0.6004722722364975], [0.9686441928327798, 0.6568063164358476, 0.5840363338264795, 0.7598460777333612, 0.9823051468131581, 0.5104082030519809], [0.1767818542801588, 0.4494109010606381, 0.3199036241240386, 0.0856212661562945, 0.6031571701939905, 0.6862272494116984], [0.10014399258229623, 0.9711348221527085, 0.580195009492141, 0.9109104898309697, 0.46601667298704486, 0.23033113918077175]]]])
      mfac = np.array([[[[0, 0, 0, 1, 0, 1], [0, 1, 0, 1, 0, 0], [1, 1, 1, 1, 1, 1], [0, 1, 1, 1, 0, 0], [1, 1, 0, 1, 0, 0]], [[0, 0, 0, 1, 1, 0], [1, 1, 1, 0, 1, 0], [0, 0, 1, 0, 0, 0], [1, 1, 1, 1, 0, 1], [1, 1, 1, 0, 0, 0]]], [[[1, 0, 1, 1, 0, 0], [0, 0, 0, 1, 1, 1], [1, 1, 1, 1, 1, 1], [0, 0, 1, 0, 1, 1], [0, 0, 1, 1, 0, 1]], [[1, 1, 1, 0, 1, 0], [0, 1, 0, 1, 0, 1], [1, 1, 0, 0, 0, 1], [1, 1, 1, 1, 1, 1], [0, 1, 0, 0, 0, 1]]], [[[0, 0, 0, 1, 0, 1], [1, 0, 0, 1, 1, 1], [0, 1, 1, 1, 0, 1], [1, 1, 0, 0, 1, 1], [1, 0, 0, 0, 1, 1]], [[1, 0, 1, 0, 0, 1], [0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 1, 1], [0, 1, 1, 1, 0, 1], [0, 1, 1, 0, 1, 0]]]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank4_fullMask_expandedIn_expandedMask(self):
      fac = np.array([[[[0.7495828253741884, 0.09742314702038724, 0.2594119714427673, 0.2756056286260855, 0.49483144022079806, 0.023360420688709915], [0.17391507787669735, 0.14747244711144103, 0.1855691963964886, 0.2712328198005859, 0.8161455389005509, 0.09020045621000983], [0.4503813507608023, 0.6729252631971012, 0.5847745165453715, 0.06664027830829922, 0.863412667971416, 0.5997168888707197], [0.04698861421553224, 0.34290722730819845, 0.7962462108763688, 0.9406386673458905, 0.9823888208553497, 0.8989945190544539], [0.4156925570940745, 0.6733301996760395, 0.15674679055452967, 0.32072219382921874, 0.10841703502716693, 0.2031205996483484]], [[0.678945451206336, 0.937292881343487, 0.4249554274670547, 0.294268928643873, 0.26845212741251356, 0.13009384641101973], [0.4011753448421236, 0.6059856747451037, 0.07015373199531694, 0.6099892215499264, 0.2485736254414157, 0.28857201839228075], [0.6336130014705423, 0.9998002441411232, 0.033527569459249174, 0.07943006199293567, 0.6998806829368511, 0.24791505541712355], [0.9028060741523405, 0.2916338595433138, 0.4493397282752437, 0.7152821006988532, 0.3263280618392064, 0.6340304833829333], [0.08814768470469381, 0.8768440284750081, 0.7338737546905696, 0.4683623648140798, 0.5923970017653253, 0.36581301909464137]]], [[[0.6266153736757945, 0.5470716814830711, 0.6217156442003238, 0.2808096073850209, 0.7726354515103347, 0.2134561766822628], [0.25910535257452305, 0.574691309948229, 0.896456171393585, 0.5489143515415948, 0.5406376653597105, 0.835843881555451], [0.7011422100582267, 0.014883970024378024, 0.4693270403330362, 0.7688067733012917, 0.886144678679621, 0.09330842284326679], [0.262639583296009, 0.9641223169730755, 0.8317795503369214, 0.02910799185680102, 0.054809351465662126, 0.012345817377606982], [0.5154809757358075, 0.82875659423914, 0.19100366409461245, 0.9242423550388797, 0.005056772058989467, 0.6065661430798827]], [[0.06348268228587306, 0.8183562560788075, 0.8196310358730861, 0.5492218768679397, 0.5877741158327879, 0.30338378445738534], [0.8815659931808789, 0.7680181307231783, 0.36735508263681393, 0.06447500379297999, 0.8377562326791684, 0.3028448290083472], [0.5610451409793171, 0.2556518595932168, 0.171562083021768, 0.30352701172093854, 0.30920627324824446, 0.7063398130166118], [0.24756692957711124, 0.030721944405571522, 0.7341311461685025, 0.11318128248709503, 0.4603688430245231, 0.5738680557360636], [0.5186862178960542, 0.8862676287364211, 0.6378121323150525, 0.19049724800358103, 0.30040029210174024, 0.9004094227536699]]], [[[0.4690725850886719, 0.07004675396219673, 0.5802076459427117, 0.11537349143130637, 0.2657433025597855, 0.47953470007981436], [0.11338822814297389, 0.9410716808837645, 0.7372815880159569, 0.9305771396952581, 0.8953606190670959, 0.6240991651118434], [0.8947522729167129, 0.5330724110620295, 0.8732058710326043, 0.5248070379090397, 0.4554741987575972, 0.4199306164363059], [0.35156108135201813, 0.3147885178501917, 0.8004943278570734, 0.1652877185704369, 0.9557968053454043, 0.15927109467383904], [0.4509537015136432, 0.4389817470631874, 0.8437214420925491, 0.10769284787838673, 0.8796550204749072, 0.45521926917161415]], [[0.10283559184517077, 0.7876231990128149, 0.18411818740895347, 0.9785405091213976, 0.515828905496374, 0.8951044707269176], [0.7454485516316658, 0.7259927891055884, 0.6581074293203452, 0.5207775441228788, 0.11685957923196921, 0.8551267592210443], [0.6226982118572767, 0.7610286300731718, 0.47009845792526195, 0.42520502800016624, 0.08671787348083682, 0.36071450760668944], [0.9246121022885089, 0.9344428012289228, 0.3030468080191916, 0.9718204572125643, 0.6623317709338582, 0.6281912061997105], [0.03990475166253882, 0.7402367852298091, 0.6960149045187072, 0.44892592746617943, 0.971780309993067, 0.8557812005944285]]]])
      mfac = np.array([[[[0, 0, 1, 1, 1, 0], [1, 1, 0, 1, 0, 1], [1, 0, 0, 0, 1, 0], [0, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 1]], [[0, 1, 1, 0, 0, 0], [0, 1, 1, 1, 0, 0], [1, 1, 0, 1, 0, 1], [1, 1, 0, 0, 1, 1], [1, 1, 1, 1, 0, 0]]], [[[1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0], [0, 1, 0, 1, 1, 1], [1, 1, 0, 1, 0, 0], [1, 0, 1, 0, 1, 1]], [[1, 1, 0, 0, 1, 0], [1, 1, 1, 1, 1, 1], [0, 1, 0, 0, 1, 0], [1, 1, 0, 0, 1, 0], [0, 1, 0, 0, 1, 0]]], [[[0, 1, 0, 0, 0, 1], [0, 1, 1, 0, 0, 0], [0, 0, 0, 1, 0, 1], [1, 0, 1, 1, 1, 0], [1, 0, 0, 1, 0, 1]], [[0, 1, 0, 1, 1, 1], [1, 0, 1, 0, 1, 1], [1, 1, 1, 1, 0, 1], [1, 1, 0, 0, 1, 1], [0, 1, 1, 0, 0, 0]]]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', mfac)
      m.expand()
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertTrue(m.isExpanded())

   def test_CopyWithMask_Rank4_scalarMask_singleIn_singleMask(self):
      fac = np.array([[[[0.0405591925927844, 0.5003199705000441, 0.3331517944229343, 0.3174331323182764, 0.03897392734172478, 0.04450599759213758], [0.4644576527328501, 0.5838187495070248, 0.7681442101907724, 0.2166707080856899, 0.012651816169763186, 0.20425745175258692], [0.1590487166503738, 0.8162114387941605, 0.21384782726578844, 0.49587480419828645, 0.6944011687677136, 0.9923766491717723], [0.49551427986717156, 0.4804296593932488, 0.012234372810717598, 0.42713665671192014, 0.6281477928074175, 0.05701649284315524], [0.07368822070637482, 0.3052773715098245, 0.5499090488415599, 0.40447517991565685, 0.017643545044307762, 0.8040745404560021]], [[0.575416987233021, 0.617073547162977, 0.08962650880642564, 0.4648423426684908, 0.4857009109670438, 0.2020180165012878], [0.18109886254286833, 0.3529897025334756, 0.8311015941371658, 0.5672425817585529, 0.5124824675453133, 0.5137913948019706], [0.7773384279018636, 0.20733759083640346, 0.4704134706896742, 0.8647167662058954, 0.7139557999155491, 0.7743674090214635], [0.7355495591883454, 0.838049068086811, 0.9465392589916887, 0.046022803719543925, 0.2238773476993311, 0.38045760363598746], [0.29290336034043496, 0.2842574170994131, 0.158329651486815, 0.9458540418173924, 0.33759220130844547, 0.8956995337549604]]], [[[0.5518203935184316, 0.5650407122640412, 0.33323458978608744, 0.1207188041350894, 0.2091844802469066, 0.6805524830387261], [0.7917538841214043, 0.12535825789642174, 0.1985467405382244, 0.08973754897737185, 0.2946820933308397, 0.4357636492989565], [0.884946834153279, 0.18718653582551747, 0.0375771987329967, 0.8531641152422417, 0.6215341634235877, 0.947911941468712], [0.30498177442807095, 0.14741780560212403, 0.6901892775719563, 0.8517613852839027, 0.6132200906063471, 0.0028632094012350784], [0.7157986626694783, 0.6776207693596168, 0.7714135609242434, 0.9037766521188643, 0.9789824123682593, 0.264480230543951]], [[0.41851626496425487, 0.17748310745232565, 0.44375007546577905, 0.8529890627690683, 0.03871647148703494, 0.8568777779592461], [0.1971028260424753, 0.8484852312544274, 0.4564204058059168, 0.6732370127973615, 0.13385923204559103, 0.08929081693469887], [0.6077440386347023, 0.9257466638227365, 0.01614529510456275, 0.11903542443679671, 0.6580283766285335, 0.5902970450121807], [0.8849261602486144, 0.3744479875049712, 0.40862528068322024, 0.9026550284589497, 0.7477784813583536, 0.3165346450594968], [0.9579213213269008, 0.7310233794324796, 0.9510859817621533, 0.03347524425056281, 0.7051700976466232, 0.8325439246715473]]], [[[0.3022629415572945, 0.8854255326047229, 0.7902032250575213, 0.4298297754354018, 0.9596582437713539, 0.6469630068869117], [0.23802884104861977, 0.8816058929405501, 0.13885774734149026, 0.08190947688965688, 0.05416166631730868, 0.14174005623054808], [0.24574175180522406, 0.7300366950187861, 0.2207044376616153, 0.6377992945640776, 0.26774294804588206, 0.1620681814336381], [0.08350215544085471, 0.6960393914824349, 0.5673528010086819, 0.11450414275256637, 0.6906123670419849, 0.6579480067228897], [0.11583517888396155, 0.7785426075100241, 0.7105076584510343, 0.9721537689882377, 0.8488681564779157, 0.5827165309769757]], [[0.6254935430187445, 0.4997360283224658, 0.6546194801245278, 0.973377395561055, 0.08927766007496984, 0.8316170728684793], [0.7780843003567377, 0.8449271479455063, 0.8776282770701859, 0.32162143705090074, 0.2287210510866785, 0.40342787132448044], [0.9686530515139824, 0.2203445988515429, 0.717778564088638, 0.9928419407021238, 0.4940048815045055, 0.03481381534452688], [0.45244731728226284, 0.615909170659405, 0.7351438858284444, 0.03992503713507067, 0.973115462029032, 0.8989412317057036], [0.7526512167289706, 0.7515793823552215, 0.7096438200974552, 0.3490604093435402, 0.3391320904790054, 0.7918313771939032]]]])
      mfac = np.array([[[[1, 0, 0, 0, 1, 1], [1, 0, 0, 1, 0, 0], [0, 1, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1], [0, 1, 1, 0, 1, 0]], [[1, 1, 1, 0, 0, 1], [0, 1, 0, 0, 0, 0], [1, 1, 1, 0, 0, 0], [0, 0, 1, 1, 1, 0], [0, 1, 0, 1, 0, 0]]], [[[1, 0, 1, 1, 0, 0], [1, 1, 1, 1, 1, 0], [0, 1, 1, 1, 0, 1], [0, 0, 1, 0, 0, 1], [1, 0, 1, 1, 1, 0]], [[0, 0, 0, 1, 1, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0], [1, 1, 0, 1, 0, 1], [0, 0, 0, 0, 0, 0]]], [[[1, 1, 1, 0, 0, 1], [0, 0, 0, 0, 1, 1], [1, 1, 0, 1, 1, 1], [1, 0, 1, 1, 1, 0], [1, 0, 0, 1, 0, 1]], [[0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 1], [0, 0, 1, 0, 1, 1], [0, 0, 1, 1, 0, 0], [1, 1, 0, 1, 1, 0]]]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank4_scalarMask_singleIn_expandedMask(self):
      fac = np.array([[[[0.275594313828341, 0.41154210910211675, 0.2454902110701156, 0.5650858745183854, 0.944630060088588, 0.14401713634973268], [0.12571110235628835, 0.9286180794137495, 0.4907172184313835, 0.6692962089758864, 0.7432387091108453, 0.3712128400939114], [0.7935020252889288, 0.9917873010612616, 0.2978947284280614, 0.14696257073666108, 0.03698642471167102, 0.4310241389870657], [0.770867584264277, 0.9821340531126714, 0.45301914531586474, 0.5447935173972609, 0.6476240419270644, 0.7875998417884432], [0.027953668150593636, 0.5106915603167067, 0.9599326837794065, 0.8018021845899082, 0.11059239232536522, 0.4259701310705638]], [[0.6412979082452267, 0.21408626742858416, 0.08582433517452526, 0.424490817559135, 0.054749654554169136, 0.2515177957761273], [0.7136918433155619, 0.4526278563733266, 0.4595656024194673, 0.03965000821563136, 0.5201205972256979, 0.15588305687180481], [0.452610871285684, 0.5879012521311522, 0.08427058423020994, 0.8487009192848999, 0.7753588894647998, 0.47136483140954255], [0.05770223867428237, 0.3593144774640542, 0.7311433238865546, 0.6683141498769898, 0.2235486747342702, 0.4918685937580344], [0.28029562337556346, 0.18757575207693855, 0.22402010398315786, 0.7310847115252193, 0.517613223297527, 0.023514808370909357]]], [[[0.9204673675798855, 0.5263954812590653, 0.6878474588588136, 0.9976471440083647, 0.37534750286749285, 0.7089615854368004], [0.8966946145188766, 0.6580564355524887, 0.7398377621095901, 0.5089387251476682, 0.9409305158150265, 0.25133236992497987], [0.2180733778157773, 0.9094329945396639, 0.5420980994232801, 0.10956763547878967, 0.7431707803698663, 0.31766557765414893], [0.6180322279394116, 0.5512489423091602, 0.2331343016655012, 0.43396057877317495, 0.2876211218634741, 0.24888059448981215], [0.8350835157072292, 0.8764726125276231, 0.791729423406322, 0.00475642831619616, 0.4331182396902582, 0.1792063757420479]], [[0.13317001924048733, 0.0053167693667189875, 0.0030429522741030057, 0.7354014197940056, 0.3829172244228799, 0.5970459539243215], [0.29665600690382277, 0.444145359126738, 0.8730931092944029, 0.003109033576882414, 0.0679810800830628, 0.30771067700448385], [0.5085666973607664, 0.8401196826708746, 0.9561794109712713, 0.23170556789351293, 0.9356602370495886, 0.34963820004147295], [0.630853162100056, 0.2951039116357419, 0.007820674205966882, 0.3000077019142088, 0.22817649935923456, 0.5249585769968387], [0.9670721811956088, 0.7311112876475695, 0.08692729220956164, 0.9938498223342334, 0.04115737009827791, 0.21869983722601394]]], [[[0.30678948048473154, 0.17756200486048423, 0.6248875899300802, 0.06518948880240005, 0.8222803091578031, 0.8157731522975481], [0.5735662669261274, 0.5940019144119213, 0.9565465972172765, 0.8374486515373778, 0.258011325581204, 0.6303464331767112], [0.9361429682058336, 0.8399213827415746, 0.27202142443423427, 0.1982585391599786, 0.39755068695451046, 0.3043496168870523], [0.4322232784256722, 0.9484367824045925, 0.05650794951797955, 0.6527399126877607, 0.13406533057855274, 0.8634218208066415], [0.10533375681812307, 0.5608980063963386, 0.710019833987499, 0.12746833207586816, 0.665726317409711, 0.8919474134998877]], [[0.9042522859493781, 0.29776973671403695, 0.8879386052396433, 0.5324050709108777, 0.6551716462512635, 0.48933425549690024], [0.06583866708934483, 0.051889252480703196, 0.0523462829518071, 0.4625815464957439, 0.7181480769946118, 0.446509819658808], [0.7369263168800787, 0.07577428663640817, 0.290978043122469, 0.6937975957305659, 0.49012199754451835, 0.48835486879906187], [0.7246052650283824, 0.4736386580218115, 0.5392588867574749, 0.6512960529700551, 0.6760977329914386, 0.9148009228393057], [0.29118140820535166, 0.7677023042535771, 0.15505535390999292, 0.09983787523956822, 0.0069171424607699095, 0.7319411890432307]]]])
      a = Data(fac, FunctionOnBoundary(self.domain))
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank4_scalarMask_taggedIn_singleMask(self):
      fac = np.array([[[[0.8098490324664231, 0.9687591386127169, 0.4792014512010625, 0.07895373193069466, 0.4241413546747824, 0.8833750348461686], [0.5673081947281586, 0.3366909443736902, 0.8721032887051752, 0.36383027724206274, 0.208942384345568, 0.3183781661343841], [0.1897079023869228, 0.6705831409510108, 0.900710000511267, 0.9061200075581489, 0.928845508836414, 0.1308649409597189], [0.13713038889174756, 0.029153958838314953, 0.9889093830068906, 0.9915228771470027, 0.08297095723497505, 0.48833814909136086], [0.11391854346467079, 0.820486014771434, 0.9731643402103014, 0.36108416496992835, 0.6427459807106543, 0.659679543945053]], [[0.5021490538471951, 0.6361794758651964, 0.06141009663419794, 0.25028158349923946, 0.014927429031961759, 0.5150863590419135], [0.039494473236582284, 0.36454330083609043, 0.9307551041876155, 0.07168022973829413, 0.03795702195424122, 0.7694104560240711], [0.8975322103823478, 0.15275049197316137, 0.9780327809299368, 0.0421104482883502, 0.7944290515721258, 0.9390737908676204], [0.06500236685476057, 0.09956113181028137, 0.29227175764660074, 0.23844119252219387, 0.7097964392883416, 0.7236519662769395], [0.6040749205725907, 0.3453859524808732, 0.676934558476247, 0.34547454479099926, 0.2582666206297728, 0.601085821095242]]], [[[0.6733132372976107, 0.2172327469450802, 0.030166528093177547, 0.9833700833677914, 0.9049852480554057, 0.4859448583843593], [0.7340695381316727, 0.7381681548164974, 0.678383149236984, 0.10070956577376311, 0.31879799444641765, 0.6845015773674425], [0.7944402204635943, 0.1047055381299894, 0.49295001576629405, 0.21889363684651075, 0.2205261134528056, 0.029982138136158687], [0.7332066715312249, 0.3171999697493657, 0.44283022105530656, 0.8402267048360382, 0.18527535815318463, 0.30112785854878965], [0.7424411952405796, 0.9899159675747673, 0.5993903806887425, 0.39055006545124005, 0.2545493927028267, 0.470617884555954]], [[0.9486366787096607, 0.8285930581216092, 0.5659151406086442, 0.9395343666807842, 0.7467169429832504, 0.275545987273651], [0.8264704751799712, 0.16203678179195147, 0.2846671444117276, 0.4441251547911903, 0.49722038758828213, 0.7089822083399726], [0.9928217042862708, 0.8396988064839795, 0.7786204710557252, 0.5714421294580723, 0.8166399475767092, 0.8833258171742644], [0.4196655290019803, 0.4167071900265501, 0.13611415195088072, 0.5857367466628193, 0.6698977999103918, 0.5074448340573086], [0.6520052564503114, 0.6062996161102265, 0.20676644515843168, 0.3769827256743604, 0.20045575002298888, 0.585647309439655]]], [[[0.8914625657051178, 0.629501199558297, 0.07085597782300979, 0.5704022633795256, 0.9721484694091969, 0.42028789207901807], [0.591065464918871, 0.8748972027083405, 0.9137892729097253, 0.7662269107942323, 0.4819224090675953, 0.45969054683241084], [0.38119478469154666, 0.47708239610294134, 0.11063935795115198, 0.11269431182198986, 0.1812262866981904, 0.45704432374289494], [0.681809579553912, 0.7030631185897487, 0.7604064752906864, 0.26004449341107994, 0.851059243685138, 0.12666932254774554], [0.5449929045701609, 0.006494759317066889, 0.27731295010105206, 0.3196717083775391, 0.741744508727598, 0.27901312746498064]], [[0.9278582503487203, 0.6264315611579352, 0.034023631834472434, 0.2721473316948203, 0.8934520161856371, 0.7116781042074987], [0.19302587533548743, 0.620497762989449, 0.1513786251275232, 0.591860901833134, 0.49400259537205826, 0.9809376076506109], [0.03327557827038219, 0.5273409511419921, 0.4143595276270512, 0.9169949397959852, 0.34185236153587983, 0.4728226703131905], [0.5177748179860026, 0.7399349698249703, 0.8954631687902087, 0.6882671712333893, 0.3117103543723815, 0.8362937565288653], [0.07092277444444639, 0.5020898579141378, 0.04833358464983517, 0.7470162436685314, 0.6258094617192077, 0.8728604680394618]]]])
      mfac = np.array([[[[0, 1, 1, 0, 0, 0], [0, 1, 1, 1, 1, 0], [1, 0, 0, 0, 0, 1], [1, 1, 1, 0, 1, 1], [0, 1, 1, 0, 1, 1]], [[1, 0, 0, 1, 1, 1], [0, 0, 0, 0, 0, 1], [1, 1, 1, 0, 1, 0], [0, 1, 1, 1, 0, 1], [1, 0, 0, 0, 0, 1]]], [[[1, 0, 0, 1, 1, 1], [1, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 1], [1, 1, 0, 1, 1, 1], [1, 1, 0, 1, 0, 0]], [[1, 1, 1, 1, 0, 1], [1, 1, 1, 0, 0, 1], [1, 1, 0, 0, 0, 1], [1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1]]], [[[0, 1, 0, 1, 1, 1], [0, 1, 1, 1, 0, 0], [1, 1, 0, 1, 1, 1], [1, 1, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0]], [[0, 0, 1, 0, 0, 1], [0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 1], [0, 1, 1, 1, 1, 0], [1, 0, 0, 1, 1, 0]]]])
      a = Data(0, (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertFalse(m.isExpanded())

   def test_CopyWithMask_Rank4_scalarMask_taggedIn_expandedMask(self):
      fac = np.array([[[[0.09717446472230873, 0.292313448430216, 0.3522991554956635, 0.4822589930900504, 0.13124907159102384, 0.26663552551295977], [0.25766144726358065, 0.1415810181331827, 0.7497920197629707, 0.6729508676645221, 0.8432175322688631, 0.3112828562305856], [0.5007030245926617, 0.37905110598108105, 0.6110308649166155, 0.7591088555262531, 0.3206986503886966, 0.2775356387583843], [0.8931565500315991, 0.48095119349856374, 0.5906366361254767, 0.3474704204068093, 0.3864807190714318, 0.8779719951555572], [0.4787728729742665, 0.6740001291562715, 0.5023322953063103, 0.9327617815691753, 0.27109167320022887, 0.13822326632376813]], [[0.21994883185283964, 0.3249536746469577, 0.5168781931937787, 0.6581885000163618, 0.005661396707903998, 0.804892775452841], [0.6140939969363799, 0.1926410168362931, 0.7425385347831881, 0.030721005852736183, 0.08793799086650145, 0.7229707506583611], [0.32830876872980586, 0.714764324086105, 0.9145574727568575, 0.2426162348790123, 0.7165356794188242, 0.6028908798398672], [0.8277937888752206, 0.7830956861245211, 0.6302001252738818, 0.9769363607697127, 0.6902757390759418, 0.6560334659969447], [0.9757827439807749, 0.7586666085663586, 0.6693789688703089, 0.8534523968273912, 0.7098842716948698, 0.0006996160446720578]]], [[[0.8816545016773334, 0.8312981048203861, 0.26026428131633284, 0.19764043981123613, 0.34077174774353214, 0.12156728295403874], [0.9006775827942727, 0.06696205837122704, 0.37761100461734165, 0.14547640305123444, 0.7664707767892258, 0.40069070508752835], [0.09590944785992528, 0.478698195408657, 0.16239666498550398, 0.6078017359356626, 0.7700002376997095, 0.508799229892264], [0.9806735121436386, 0.9096624786492297, 0.10142795136234817, 0.025690130232941977, 0.19141204697186565, 0.85218135333202], [0.0456871498623842, 0.18297756319756153, 0.5687160970984109, 0.38892884354597634, 0.6816407539424053, 0.9381893821007532]], [[0.22522742526675976, 0.18283806422721505, 0.3286834800686894, 0.7656949216060334, 0.7790388144350446, 0.32835654753328214], [0.1531875261042357, 0.5626586594393258, 0.7484182596703014, 0.594815327255118, 0.2649530089418344, 0.9688145542770774], [0.12974496377590972, 0.6964295177236213, 0.21575615149539895, 0.3213208666342513, 0.3736779757779216, 0.2089987282742536], [0.6966364270130709, 0.5589395656458586, 0.5126742201297763, 0.649266535805653, 0.48208312666720243, 0.19524971660203527], [0.6801765068910064, 0.6882174248322774, 0.6369372084258119, 0.6642739845700506, 0.7772159044388629, 0.697207683926809]]], [[[0.9215692069710915, 0.13950776617151184, 0.056098686158856714, 0.8320585522181098, 0.7226645113459641, 0.5123763031542093], [0.26304274552352147, 0.2835997422478884, 0.4158428571848407, 0.6056974805045771, 0.08757468071144225, 0.20805049662025688], [0.821830523201849, 0.5480184812236741, 0.6760007562825774, 0.29754009252763647, 0.5381845577080989, 0.15388858345643852], [0.8879619413488697, 0.3153802168991545, 0.7625247788729632, 0.9674494355500671, 0.29327711013272517, 0.1792927601130765], [0.4272455533779207, 0.7794572780861365, 0.5729753271160447, 0.1809271377905498, 0.8993437347289156, 0.6266293497711126]], [[0.8673322592554993, 0.4436374455988331, 0.7962246100718354, 0.7280110475367335, 0.3858286568553817, 0.5166926273975868], [0.6110583887202738, 0.4755463577237843, 0.4747281870250969, 0.2145444910983122, 0.6255268020123184, 0.5263656950012158], [0.5097454707046055, 0.14171517785350518, 0.7862225845468278, 0.9583128996335926, 0.417423623980151, 0.77444413268793], [0.44941056590202744, 0.41290635498668626, 0.22235801269859756, 0.44609871755334485, 0.06933898930286742, 0.8311300926448738], [0.27283703374664536, 0.2095121026995419, 0.8804458214414774, 0.1336219021024151, 0.3849716382169671, 0.8405523887390912]]]])
      a = Data(0, (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      a.setTaggedValue('top', fac)
      a.setTaggedValue('bottom', fac)
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertFalse(a.isExpanded())
      self.assertTrue(m.isExpanded())
   def test_CopyWithMask_Rank4_scalarMask_expandedIn_singleMask(self):
      fac = np.array([[[[0.6862010857597054, 0.366434849229648, 0.008458904447123383, 0.9639067231322768, 0.4393689114112703, 0.14606354832424817], [0.675038912856568, 0.811624724334439, 0.92842036233952, 0.06971114478951745, 0.8835838609044451, 0.7104405318244976], [0.5333203964105102, 0.5480566370691494, 0.4794109605567034, 0.3357144690279771, 0.6200634366249721, 0.40761522461349486], [0.09409305876915552, 0.17543676714057888, 0.05597755048457509, 0.4347897630277179, 0.47850926432541685, 0.21756465934354818], [0.4972382001837645, 0.19137659245557193, 0.48853458040323805, 0.9920899181084524, 0.5010156856694368, 0.7517212447829956]], [[0.8312971500726064, 0.5690951374581835, 0.25184868992849196, 0.2778723578995853, 0.6973475036463278, 0.3480241543295559], [0.5527971044229525, 0.9023813480158557, 0.978117948414839, 0.47700807038435433, 0.7394809856201019, 0.4144907641627189], [0.46875947619118474, 0.8476190156328883, 0.6497526487217814, 0.2354395058125761, 0.36863460339219245, 0.889547572276826], [0.2393153728921089, 0.5685252889812312, 0.614331053316104, 0.4786493095669724, 0.8292618760777034, 0.8736279106266235], [0.6724973846989745, 0.7596684009170975, 0.06600218641440347, 0.23726666728760926, 0.12073443104759407, 0.14665443804518807]]], [[[0.8746533041852003, 0.5518781457756115, 0.14992602511169473, 0.7354492649606098, 0.6249239842031767, 0.8750549184716269], [0.7634586500260844, 0.997686735951208, 0.8540708227551214, 0.3109103213784421, 0.45061340446924547, 0.8099927334851077], [0.2942570015823939, 0.3919510565904065, 0.825031241001668, 0.7684769783794273, 0.019582421591461285, 0.09733035305397875], [0.9124143514868379, 0.06707875780152672, 0.7470987074782445, 0.42447044105956777, 0.7330011866211368, 0.6767390532547443], [0.8621679112362975, 0.9851456402135655, 0.6052729731834932, 0.8610319177162125, 0.1691928651036093, 0.1510590426844456]], [[0.9342595852483552, 0.1512686687374789, 0.4110729055997274, 0.2844602066440479, 0.4466698425244965, 0.9003764724767658], [0.6703656955090221, 0.06925535800389504, 0.760027129311875, 0.297344476131984, 0.5604518776440561, 0.35985352459208375], [0.776164760212593, 0.007011725931318824, 0.9168036233692459, 0.5648285814605376, 0.18837633080367755, 0.14759654957443746], [0.4086514797861296, 0.3757770723699311, 0.05084222548914308, 0.9946854789449722, 0.03296838045694617, 0.20772608306819107], [0.7856356982533486, 0.6611750568310324, 0.9833350026635055, 0.24240717251114707, 0.544781825007803, 0.7318153248475671]]], [[[0.30597108670588524, 0.7914276707513803, 0.7561976309648499, 0.811172362219399, 0.6139344975954201, 0.6508037915691791], [0.1025470755692971, 0.039928758874791415, 0.5182580020403135, 0.8723750080978314, 0.25168224919327686, 0.39746383148910025], [0.3816971074208685, 0.0780092947410368, 0.24926703032994946, 0.022079808236906473, 0.3705149267124743, 0.27149760366373266], [0.061423226005362874, 0.6424199257916028, 0.21612143975495024, 0.010533896880240867, 0.889887108444731, 0.6309376152557631], [0.7789462545756207, 0.06507050059214192, 0.16168642691842328, 0.6319331492114201, 0.06849733311799022, 0.990378977268515]], [[0.7774725417618986, 0.16633334825545987, 0.4042953186669169, 0.7208861917496744, 0.8756443161827009, 0.052718620193694554], [0.915800854599193, 0.14792318881362154, 0.0959685677259775, 0.5488434058978711, 0.30882479235922056, 0.9615998911031416], [0.9552531915388789, 0.7371455848985868, 0.5565490323099859, 0.9483112030951525, 0.8877885487366709, 0.8287131987649984], [0.4753486069259568, 0.41523463854702813, 0.2892557954609556, 0.7955402119955531, 0.451471750895153, 0.6142291694595814], [0.7661222786197202, 0.26408919590232205, 0.7480194378293613, 0.5963644145149376, 0.7733829389052881, 0.7356731687826711]]]])
      mfac = np.array([[[[1, 1, 0, 0, 0, 0], [1, 1, 0, 1, 1, 0], [1, 1, 1, 0, 0, 0], [0, 1, 1, 0, 0, 0], [1, 0, 1, 0, 0, 0]], [[0, 0, 1, 0, 0, 1], [1, 1, 1, 1, 1, 0], [1, 1, 1, 0, 0, 0], [1, 1, 0, 0, 0, 0], [0, 1, 1, 0, 1, 1]]], [[[0, 1, 1, 1, 0, 0], [0, 1, 1, 1, 1, 0], [0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 1], [0, 0, 0, 0, 0, 0]], [[0, 0, 1, 0, 1, 1], [1, 0, 1, 1, 0, 1], [0, 1, 0, 1, 0, 1], [0, 1, 1, 0, 0, 1], [1, 1, 0, 1, 0, 0]]], [[[0, 1, 0, 0, 0, 1], [1, 1, 1, 1, 1, 0], [1, 1, 0, 1, 1, 0], [0, 1, 0, 1, 1, 0], [1, 0, 1, 1, 1, 0]], [[1, 1, 0, 1, 1, 1], [0, 1, 0, 0, 0, 0], [0, 1, 1, 0, 0, 1], [1, 1, 0, 0, 0, 1], [0, 1, 1, 1, 1, 0]]]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(mfac, FunctionOnBoundary(self.domain))
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertFalse(m.isExpanded())
   def test_CopyWithMask_Rank4_scalarMask_expandedIn_expandedMask(self):
      fac = np.array([[[[0.4212048533210767, 0.006262620542562103, 0.44324722268823613, 0.4054736411715598, 0.36018725124329054, 0.6120197039857977], [0.7354802817963213, 0.588490444484958, 0.6076199636820221, 0.4431054531009492, 0.9371858958980814, 0.4923334944097869], [0.9533322730248813, 0.006085710155166457, 0.25935794309672255, 0.9285517491961277, 0.18601794235575075, 0.15797657056785797], [0.8423595828232079, 0.10927605368217241, 0.20643696685884816, 0.46196778726343435, 0.85464606176904, 0.38895146157727767], [0.26272468163731566, 0.6502365083699232, 0.2288226000819802, 0.2811740569759913, 0.5423289500867232, 0.4298134425648379]], [[0.12473767594504981, 0.04698201826712822, 0.9840392204302278, 0.11415160899923049, 0.7510736294977112, 0.019585708529396784], [0.11854077886477121, 0.40573862094455815, 0.7685315150833644, 0.6553313884365054, 0.14390162055492417, 0.7055151119711274], [0.9999506559474575, 0.2704540558844354, 0.5683709948766447, 0.08599341073816791, 0.7116698156636411, 0.317795013713965], [0.9917446737410551, 0.9121707653304977, 0.2878918683203301, 0.6725045918998145, 0.0480339698807043, 0.7002576872455493], [0.323550927896634, 0.9566744620745428, 0.44127320347950705, 0.19286002932755364, 0.6316626891196322, 0.7174697977197761]]], [[[0.20868684570406137, 0.0012005676907441698, 0.02988243056049633, 0.7559488664286255, 0.6014084885790496, 0.08473820269450161], [0.3352945470389497, 0.6136427007622645, 0.596265143220559, 0.62681267150319, 0.37255891465571955, 0.706672842050472], [0.10422305487845451, 0.629781978587141, 0.002695163986669691, 0.9857834045729316, 0.2706667607311113, 0.6292806436597079], [0.9141835970477656, 0.1091105052073158, 0.35061213427043947, 0.39144701224725564, 0.6326695896682373, 0.07993743145873411], [0.8926989726260904, 0.44840200393230745, 0.14147926719457005, 0.46461417939283123, 0.7791407152959933, 0.4837373275283483]], [[0.11539127422310835, 0.5223289006245327, 0.6372483735219494, 0.26905408195300196, 0.6660097124277047, 0.2142013227479148], [0.09604391862498107, 0.7080355287560851, 0.18632699673537523, 0.6459322857074266, 0.3803885362354852, 0.5075304930563064], [0.02639413635173049, 0.5629731650774379, 0.26927109611887645, 0.9661576388514624, 0.9915254909169912, 0.47642927752655173], [0.3721333543997577, 0.20492823524643167, 0.46659569648405996, 0.5962597595554039, 0.4086780352073003, 0.6358185199204799], [0.20578881423360773, 0.30404917582408597, 0.5726845062018174, 0.33051978413902205, 0.9587279179038729, 0.1862602813150318]]], [[[0.6424268155441625, 0.8092805432651543, 0.09679615644440909, 0.39552857383378137, 0.9378054254621545, 0.3280917347303226], [0.7118974031553652, 0.7552933176937816, 0.9568724586371049, 0.488700298939606, 0.17389603033960654, 0.8768648270335301], [0.03442005988739216, 0.6039012243770363, 0.2289933239102715, 0.7281935047705579, 0.28423802107013196, 0.10953197759693611], [0.35217493094438257, 0.7270538906321715, 0.6168836316462364, 0.5278429339198165, 0.057582373487480965, 0.8099706026071666], [0.1590690668595608, 0.8438789203535081, 0.07584547040503997, 0.5880532834890001, 0.625006050449454, 0.03101224785132517]], [[0.8433044091711516, 0.6572965630949376, 0.38846149215857684, 0.9744282115827113, 0.735784741961547, 0.32377927024702813], [0.7353303772099484, 0.4550107297342981, 0.06715455944398863, 0.964551610931829, 0.12663755512553432, 0.9408288723450984], [0.9364787338314408, 0.18149196389391375, 0.38356701572431373, 0.6763479686627537, 0.9991652923243992, 0.3451362747202643], [0.2962657141428474, 0.27375967650954136, 0.25412641180056195, 0.34599084115395184, 0.47255212943077285, 0.6844798741477106], [0.2029794507971744, 0.6343729164343606, 0.562920919173987, 0.07867814597684852, 0.7314288037280526, 0.9110691961726349]]]])
      a = length(FunctionOnBoundary(self.domain).getX()) * fac
      m = Data(0, (), FunctionOnBoundary(self.domain))
      m.setTaggedValue('top', 1)
      m.expand()
      b = Data(1., (3, 2, 5, 6), FunctionOnBoundary(self.domain))
      b_true = b * (1-m) + m * a
      b.copyWithMask(a, m)
      self.assertEqual(b.isExpanded(), m.isExpanded() or a.isExpanded())
      self.assertLess(Lsup(b-b_true), 1e-15 * Lsup(b_true))
      self.assertTrue(a.isExpanded())
      self.assertTrue(m.isExpanded())

if __name__ == '__main__':
   run_tests(__name__, exit_on_failure=True)