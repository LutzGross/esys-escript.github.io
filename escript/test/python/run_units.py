
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import unittest
from esys.escript.unitsSI import *
from esys.escript.util import EPSILON

class UnitsSITestCase(unittest.TestCase):
    TOL=EPSILON*100.
    def testPrefix(self):
       self.failUnless(abs( Yotta * Yocto - 1.) <= self.TOL , " Yotta or Yocto wrong.")
       self.failUnless(abs( Zetta * Zepto - 1.) <= self.TOL , "  Zetta or Zepto  wrong.")
       self.failUnless(abs( Exa * Atto - 1.) <= self.TOL , " Exa or Atto  wrong.")
       self.failUnless(abs( Peta *Femto - 1.) <= self.TOL , " Peta orFemto  wrong.")
       self.failUnless(abs( Tera * Pico - 1.) <= self.TOL , " Tera or Pico  wrong.")
       self.failUnless(abs( Giga * Nano - 1.) <= self.TOL , "  Giga or Nano  wrong.")
       self.failUnless(abs( Mega *Micro - 1.) <= self.TOL , " Mega orMicro  wrong.")
       self.failUnless(abs( Kilo * Milli - 1.) <= self.TOL , "Kilo or Milli   wrong.")
       self.failUnless(abs( Hecto * Centi- 1.) <= self.TOL , " Hecto or Centi  wrong.")
       self.failUnless(abs( Deca * Deci - 1.) <= self.TOL , "Deca or Deci  wrong.")
    def testLength(self):
       self.failUnless(abs( m - 1.) <= self.TOL , " meter wrong.")
       self.failUnless(abs( km - 1000.) <= self.TOL*1000. , " km wrong.")
       self.failUnless(abs( cm - 1./100.) <= self.TOL/100. , " cm wrong.")
       self.failUnless(abs( mm - 1./1000.) <= self.TOL/1000. , " mm wrong.")
    def testTime(self):
       self.failUnless(abs( sec - 1.) <= self.TOL , " second wrong.")
       self.failUnless(abs( minute - 60.) <= self.TOL*60 , " minute wrong.")
       self.failUnless(abs( h - 3600) <= self.TOL*3600., " h wrong.")
       self.failUnless(abs( day - 86400.) <= self.TOL*86400 , " day wrong.")
       self.failUnless(abs( yr - 31556952) <= self.TOL*31556952 , " yr wrong.")
       self.failUnless(abs( Myr - 31556952.e6) <= self.TOL*31556952.e6 , " Myr wrong.")
       self.failUnless(abs( Gyr - 31556952.e9) <= self.TOL*31556952.e9 , " Gyr wrong.")
    def testMass(self):
       self.failUnless(abs( gram - 1./1000.) <= self.TOL/1000. , " gram wrong.")
       self.failUnless(abs( kg - 1.) <= self.TOL , " kg wrong.")
       self.failUnless(abs( ton - 1.e3) <= self.TOL*1.e3 , " ton wrong.")
       self.failUnless(abs( lb - 0.45359237) <= self.TOL , " lb wrong.")
    def testCurrent(self):
       self.failUnless(abs( A- 1.) <= self.TOL , " A wrong.")
    def testFrequency(self):
       self.failUnless(abs( Hz- 1.) <= self.TOL , " Hz wrong.")
    def testForce(self):
       self.failUnless(abs( N- 1.) <= self.TOL , " N wrong.")
    def testPressure(self):
       self.failUnless(abs( Pa- 1.) <= self.TOL , " Pa wrong.")
       self.failUnless(abs( atm- 101325.024) <= self.TOL*101,325.024 , " Pa wrong.")
    def testEnegry(self):
       self.failUnless(abs( J- 1.) <= self.TOL , " J wrong.")
    def testPower(self):
       self.failUnless(abs( W- 1.) <= self.TOL , " W wrong.")
    def testPotential(self):
       self.failUnless(abs( V- 1.) <= self.TOL , " V wrong.")
    def testCharge(self):
       self.failUnless(abs( C- 1.) <= self.TOL , " C wrong.")
    def testResistance(self):
       self.failUnless(abs( Ohm- 1.) <= self.TOL , " Ohm wrong.")
    def testEarth(self):
       R=6367444.65
       self.failUnless(abs( R_Earth_equator-6378137) <= self.TOL*R , " R_Earth_equator wrong.")
       self.failUnless(abs( R_Earth_poles-6356752.3) <= self.TOL*R , " R_Earth_poles wrong.")
       self.failUnless(abs( R_Earth-R) <= self.TOL*R , "R_Earth wrong.")
    def testVlight(self):
       self.failUnless(abs(v_light - 299792458.) <= self.TOL*299792458. , " speed of light is wrong.")

if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(UnitsSITestCase))
    s=unittest.TextTestRunner(verbosity=2).run(suite)
    if not s.wasSuccessful(): sys.exit(1)
