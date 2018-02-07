
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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"



import logging
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.downunder.coordinates import *
import esys.escript
from esys.escript import Lsup, sin, cos
import esys.escript.unitsSI as U
import os
import sys
from math import pi
# this is mainly to avoid warning messages
logging.basicConfig(format='%(name)s: %(message)s', level=logging.INFO)

HAS_RIPLEY = True
try:
    from esys.ripley import Brick, Rectangle
except ImportError as e:
    HAS_RIPLEY = False

try:
    WORKDIR=os.environ['DOWNUNDER_WORKDIR']
except KeyError:
    WORKDIR='.'

NE=20 
RTOL=1e-8

class Test_Coordinates(unittest.TestCase):
    def test_CartesianReferenceSystem(self):
         cs=CartesianReferenceSystem()
         
         self.assertEqual(cs.getName(), "CARTESIAN", "wrong name")
         self.assertTrue(str(cs).startswith("CARTESIAN (id"), "wrong str")
         self.assertTrue(cs == cs, "wrong cs == cs")
         self.assertFalse(cs != cs, "wrong cs == cs")
         self.assertFalse(cs == None, "wrong cs == None")
         self.assertTrue(cs != None, "wrong cs == None")

         cs2=CartesianReferenceSystem()
         self.assertTrue(cs == cs2, "wrong cs == cs2")
         self.assertFalse(cs != cs2, "wrong cs == cs2")
         self.assertTrue(cs2 == cs, "wrong cs2 == cs")
         self.assertFalse(cs2 != cs, "wrong cs2 == cs")

    def test_GeodeticReferenceSystem(self):

         cs = GeodeticReferenceSystem(a= 4*U.m, f=0.75, angular_unit=0.1, name="A1")
         self.assertEqual(cs.getName(), "A1", "wrong name")
         self.assertEqual(cs.getAngularUnit(), 0.1, "wrong angular_unit")
         self.assertEqual(cs.getFlattening(), 0.75, "wrong flattening")
         self.assertEqual(cs.getSemiMajorAxis(), 4., "wrong SemiMajorAxis")
         self.assertEqual(cs.getSemiMinorAxis(), 1., "wrong SemiMinorAxis")
         self.assertTrue(cs == cs, "wrong cs == cs")
         self.assertFalse(cs != cs, "wrong cs == cs")
         self.assertFalse(cs == None, "wrong cs == None")
         self.assertTrue(cs != None, "wrong cs == None")

         # lets try a few stupid things:
         self.assertRaises(ValueError, GeodeticReferenceSystem, a= -4)
         self.assertRaises(ValueError, GeodeticReferenceSystem, f=-0.75)
         self.assertRaises(ValueError, GeodeticReferenceSystem, f=1.)
         self.assertRaises(ValueError, GeodeticReferenceSystem, angular_unit=-0.1)
         
         # cs == CartesianReferenceSystem() ?
         cs2=CartesianReferenceSystem()
         self.assertFalse(cs == cs2, "wrong cs == cs2")
         self.assertTrue(cs != cs2, "wrong cs == cs2")
         self.assertFalse(cs2 == cs, "wrong cs2 == cs")
         self.assertTrue(cs2 != cs, "wrong cs2 == cs")         
         
         cs2 = GeodeticReferenceSystem(a= 4*U.m, f=0.75, angular_unit=0.1, name="A2")
         self.assertTrue(cs == cs2, "wrong cs == cs2")
         self.assertFalse(cs != cs2, "wrong cs == cs2")
         self.assertTrue(cs2 == cs, "wrong cs2 == cs")
         self.assertFalse(cs2 != cs, "wrong cs2 == cs")
         
         # different semi majr axis?
         cs2 = GeodeticReferenceSystem(a= 2*U.m, f=0.75, angular_unit=0.1, name="A2")
         self.assertFalse(cs == cs2, "wrong cs == cs2")
         self.assertTrue(cs != cs2, "wrong cs == cs2")
         self.assertFalse(cs2 == cs, "wrong cs2 == cs")
         self.assertTrue(cs2 != cs, "wrong cs2 == cs")
          
         # different flattening
         cs2 = GeodeticReferenceSystem(a= 4*U.m, f=0.7, angular_unit=0.1, name="A2")
         self.assertFalse(cs == cs2, "wrong cs == cs2")
         self.assertTrue(cs != cs2, "wrong cs == cs2")
         self.assertFalse(cs2 == cs, "wrong cs2 == cs")
         self.assertTrue(cs2 != cs, "wrong cs2 == cs")
         
         # different angular_unit
         cs2 = GeodeticReferenceSystem(a= 4*U.m, f=0.75, angular_unit=0.2, name="A2")
         self.assertFalse(cs == cs2, "wrong cs == cs2")
         self.assertTrue(cs != cs2, "wrong cs == cs2")
         self.assertFalse(cs2 == cs, "wrong cs2 == cs")
         self.assertTrue(cs2 != cs, "wrong cs2 == cs")   
         
         sp=SphericalReferenceSystem(R=1*U.km)
         self.assertEqual(sp.getName(), "SPHERE", "wrong name")
         self.assertAlmostEqual(sp.getAngularUnit(), 1*U.DEG, msg="wrong angular_unit")
         self.assertAlmostEqual(sp.getFlattening(), 0., msg="wrong flattening")
         self.assertAlmostEqual(sp.getSemiMajorAxis(), 1000.,msg= "wrong SemiMajorAxis")
         
         sp=SphericalReferenceSystem()
         self.assertEqual(sp.getName(), "SPHERE", msg= "wrong name")
         self.assertAlmostEqual(sp.getAngularUnit(), 1*U.DEG, msg="wrong angular_unit")
         self.assertAlmostEqual(sp.getFlattening(), 0., msg="wrong flattening")
         self.assertAlmostEqual(sp.getSemiMajorAxis(), 6378137, msg="wrong SemiMajorAxis")
           
           
         sp=WGS84ReferenceSystem()
         self.assertEqual(sp.getName(), "WGS84", msg= "wrong name")
         self.assertAlmostEqual(sp.getAngularUnit(), 1*U.DEG, msg="wrong angular_unit")
         self.assertAlmostEqual(sp.getFlattening(), 0.0033528106647474805, msg="wrong flattening")
         self.assertAlmostEqual(sp.getSemiMajorAxis(), 6378137, msg="wrong SemiMajorAxis")
         
         sp=GRS80ReferenceSystem()
         self.assertEqual(sp.getName(), "GRS80", msg= "wrong name")
         self.assertAlmostEqual(sp.getAngularUnit(), 1*U.DEG, msg="wrong angular_unit")
         self.assertAlmostEqual(sp.getFlattening(), 0.003352810681182319, msg="wrong flattening")
         self.assertAlmostEqual(sp.getSemiMajorAxis(), 6378137, msg="wrong SemiMajorAxis")

    @unittest.skipIf(not HAS_RIPLEY, "Ripley domains unavailable")
    def test_CartesianTransformation3D(self):
      
         dom=Brick(NE,NE,NE, l0=10, l1=10, l2=10)
         
         cs=CartesianReferenceSystem()
         tf=cs.createTransformation(dom)
         
         self.assertEqual(tf.getReferenceSystem(),cs , "wrong reference")
         self.assertEqual(tf.getDomain(),dom , "wrong reference")
         self.assertTrue(tf.isCartesian(), "wrong isCartesian check")
         
         
         v=tf.getVolumeFactor()
         self.assertTrue(isinstance(v, esys.escript.Data), "wrong volume factor type")
         self.assertEqual(v.getFunctionSpace(), esys.escript.Function(dom), "wrong volume factor type")
         error=Lsup(v-1.)
         self.assertTrue(error<=RTOL, "volume factor")
         
         s=tf.getScalingFactors()
         self.assertTrue(isinstance(s, esys.escript.Data), "scaling factor type")
         self.assertEqual(s.getShape(), (dom.getDim(),), "scaling factor length")
         self.assertEqual(s.getFunctionSpace(), esys.escript.Function(dom), "wrong 0-th scaling factor type")

         error=Lsup(s[0]-1.)
         self.assertTrue(error<=RTOL, "0-th scaling factor")         
         error=Lsup(s[1]-1.)
         self.assertTrue(error<=RTOL, "1-th scaling factor")  
         error=Lsup(s[2]-1.)
         self.assertTrue(error<=RTOL, "2-th scaling factor")   

    @unittest.skipIf(not HAS_RIPLEY, "Ripley domains unavailable")
    def test_CartesianTransformation2D(self):
      
         dom=Rectangle(NE,NE, l0=10, l1=10)
         
         cs=CartesianReferenceSystem()
         tf=cs.createTransformation(dom)
         
         self.assertEqual(tf.getReferenceSystem(),cs , "wrong reference")
         self.assertEqual(tf.getDomain(),dom , "wrong reference")
         self.assertTrue(tf.isCartesian(), "wrong isCartesian check")
         
         
         v=tf.getVolumeFactor()
         self.assertTrue(isinstance(v, esys.escript.Data), "wrong volume factor type")
         self.assertEqual(v.getFunctionSpace(), esys.escript.Function(dom), "wrong volume factor type")
         error=Lsup(v-1.)
         self.assertTrue(error<=RTOL, "volume factor")
         
         s=tf.getScalingFactors()
         self.assertTrue(isinstance(s, esys.escript.Data), "scaling factor type")
         self.assertEqual(s.getShape(), (dom.getDim(),), "scaling factor length")
         self.assertEqual(s.getFunctionSpace(), esys.escript.Function(dom), "wrong 0-th scaling factor type")
         
         error=Lsup(s[0]-1.)
         self.assertTrue(error<=RTOL, "0-th scaling factor")         
         
         error=Lsup(s[1]-1.)
         self.assertTrue(error<=RTOL, "1-th scaling factor")  

    @unittest.skipIf(not HAS_RIPLEY, "Ripley domains unavailable")
    def test_SphericalTransformation3D(self):
      
         dom=Brick(NE,NE,NE, l0=90, l1=45, l2=10.)
         
         cs=SphericalReferenceSystem()
         tf=cs.createTransformation(dom)
         
         self.assertEqual(tf.getReferenceSystem(),cs , "wrong reference")
         self.assertEqual(tf.getDomain(),dom , "wrong reference")
         self.assertFalse(tf.isCartesian(), "wrong isCartesian check")
         
         R=6378137.0
         x=esys.escript.Function(dom).getX()
         phi=(90.-x[1])/180.*pi
         lam=x[0]/180.*pi
         h=x[2]*1000.
         r=h+R
         
         v=tf.getVolumeFactor()
         self.assertTrue(isinstance(v, esys.escript.Data), "wrong volume factor type")
         self.assertEqual(v.getFunctionSpace(), esys.escript.Function(dom), "wrong volume factor type")
         error=Lsup(v- r**2*sin(phi)*(pi/180.)**2*1000. )
         self.assertTrue(error<=RTOL * R*R*(pi/180.)**2*1000., "volume factor")
         
         s=tf.getScalingFactors()
         self.assertTrue(isinstance(s, esys.escript.Data), "scaling factor type")
         self.assertEqual(s.getShape(), (dom.getDim(),), "scaling factor length")
         self.assertEqual(s.getFunctionSpace(), esys.escript.Function(dom), "wrong 0-th scaling factor type")
         
         error=Lsup(s[1]-1/r/pi*180.)
         self.assertTrue(error<=RTOL/R/pi*180., "0-th scaling factor")         
         
         error=Lsup(s[0]-1/(r*sin(phi))/pi*180.)
         self.assertTrue(error<=RTOL/R/pi*180., "1-th scaling factor")  
         
         error=Lsup(s[2]-1./1000.)
         self.assertTrue(error<=RTOL/1000., "2-th scaling factor")   

    @unittest.skipIf(not HAS_RIPLEY, "Ripley domains unavailable")
    def test_SphericalTransformation2D(self):
      
         dom=Rectangle(NE,NE, l0=45., l1=10.)
         
         cs=SphericalReferenceSystem()
         tf=cs.createTransformation(dom)
         
         self.assertEqual(tf.getReferenceSystem(),cs , "wrong reference")
         self.assertEqual(tf.getDomain(),dom , "wrong reference")
         self.assertFalse(tf.isCartesian(), "wrong isCartesian check")
         
         R=6378137.0
         x=esys.escript.Function(dom).getX()
         phi=(90.-x[0])/180.*pi
         h=x[1]*1000.
         r=h+R
         
         v=tf.getVolumeFactor()
         self.assertTrue(isinstance(v, esys.escript.Data), "wrong volume factor type")
         self.assertEqual(v.getFunctionSpace(), esys.escript.Function(dom), "wrong volume factor type")
         error=Lsup(v-r*pi/180.*1000.)
         self.assertTrue(error<=RTOL*R, "volume factor")
         
         s=tf.getScalingFactors()
         self.assertTrue(isinstance(s, esys.escript.Data), "scaling factor type")
         self.assertEqual(s.getShape(), (dom.getDim(),), "scaling factor length")
         self.assertEqual(s.getFunctionSpace(), esys.escript.Function(dom), "wrong 0-th scaling factor type")

         error=Lsup(s[0]-1./r/pi*180.)
         self.assertTrue(error<=RTOL/R/pi*180., "0-th scaling factor")         
         error=Lsup(s[1]-1./1000)
         self.assertTrue(error<=RTOL/1000., "1-th scaling factor")  
         
         
if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

