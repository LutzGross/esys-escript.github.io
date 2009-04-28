
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
import sys
from esys.escript.unitsSI import *
from esys.escript.util import EPSILON

class UnitsSITestCase(unittest.TestCase):
    TOL=EPSILON*100.
    def testUnit(self):
       s=Unit("s","something",1.,4.)
       self.failUnless(s.getName() == "s", "wrong name")
       s.setName("t")
       self.failUnless(s.getName() == "t", "wrong reset name")
       self.failUnless(s.getLongName() == "something", "wrong long name")
       s.setLongName("test")
       self.failUnless(s.getLongName() == "test", "wrong reset long name")
       self.failUnless(abs( s(3.)-13.) <= self.TOL*13. , " s(3.) wrong.")
       self.failUnless(abs( 3.*s-13.) <= self.TOL*13. , " 13.s wrong.")
       self.failUnless(abs( 10.*s-41.) <= self.TOL*41. , " 10.s wrong.")
       self.failUnless(abs( 41./s-10.) <= self.TOL*10. , " 10/s wrong.")
       self.failUnless(abs( 13./s-3.) <= self.TOL*3. , " 13/s wrong.")

       x=Unit("x","X",10.,2.)
       y=Unit("y","Y",0.,4.)
       z=Unit("z","Z",0.,0.5)
       p_m=x*s
       self.failUnless(p_m.getName() == "xt", "wrong p_m name")
       self.failUnless(p_m.getLongName() == "X*test", "wrong p_m long name")
       self.failUnless(abs(p_m(0.)-12.) <= self.TOL*12 , " p_m(0) wrong.")
       self.failUnless(abs(p_m(1.)-20.) <= self.TOL*20 , " p_m(1) wrong.")

       p_d=y/z
       self.failUnless(p_d.getName() == "y/z", "wrong p_d name")
       self.failUnless(p_d.getLongName() == "Y/Z", "wrong p_d long name")
       self.failUnless(abs(p_d(1.)-8.) <= self.TOL*8. , " p_d(0) wrong.")
       self.failUnless(abs(p_d(2.)-16.) <= self.TOL*16. , " p_d(1) wrong.")

       p_p=y**3
       self.failUnless(p_p.getName() == "y^3", "p_p mult name")
       self.failUnless(p_p.getLongName() == "Y^3", "wrong p_p long name")
       self.failUnless(abs(p_p(1.)-64.) <= self.TOL*64. , " p_p(0) wrong.")
       self.failUnless(abs(p_p(.5)-32) <= self.TOL*32. , " p_p(1) wrong.")

        
       x=Unit("x","X",0.,2.)
       y=Unit("y","Y",0.,4.)
       z=Unit("z","Z",0.,0.5)

       self.failUnless(((x*x)).getName() == "xx", "wrong name")
       self.failUnless((x*(x*z)).getName() == "xxz", "wrong name")
       self.failUnless((x*(x/z)).getName() == "xx/z", "wrong name")
       self.failUnless((x*(y**3)).getName() == "x y^3", "wrong name")

       self.failUnless(((x*z)*x).getName() == "xzx", "wrong name")
       self.failUnless(((x*z)*(x*z)).getName() == "xzxz", "wrong name")
       self.failUnless(((x*z)*(x/z)).getName() == "xzx/z", "wrong name")
       self.failUnless(((x*z)*(y**3)).getName() == "xz y^3", "wrong name")

       self.failUnless(((x/z)*x).getName() == "x/z x", "wrong name")
       self.failUnless(((x/z)*(x*z)).getName() == "x/z xz", "wrong name")
       self.failUnless(((x/z)*(x/z)).getName() == "x/z x/z", "wrong name")
       self.failUnless(((x/z)*(y**3)).getName() == "x/z y^3", "wrong name")

       self.failUnless(((y**3)*x).getName() == "y^3 x", "wrong name")
       self.failUnless(((y**3)*(x*z)).getName() == "y^3 xz", "wrong name")
       self.failUnless(((y**3)*(x/z)).getName() == "y^3 x/z", "wrong name")
       self.failUnless(((y**3)*(y**3)).getName() == "y^3 y^3", "wrong name")

       self.failUnless((x/x).getName() == "x/x", "wrong name")
       self.failUnless((x/(x*z)).getName() == "x/(xz)", "wrong name")
       self.failUnless((x/(x/z)).getName() == "x/(x/z)", "wrong name")
       self.failUnless((x/(y**3)).getName() == "x/y^3", "wrong name")

       self.failUnless(((x*z)/x).getName() == "xz/x", "wrong name")
       self.failUnless(((x*z)/(x*z)).getName() == "xz/(xz)", "wrong name")
       self.failUnless(((x*z)/(x/z)).getName() == "xz/(x/z)", "wrong name")
       self.failUnless(((x*z)/(y**3)).getName() == "xz/y^3", "wrong name")

       self.failUnless(((x/z)/x).getName() == "x/z/x", "wrong name")
       self.failUnless(((x/z)/(x*z)).getName() == "x/z/(xz)", "wrong name")
       self.failUnless(((x/z)/(x/z)).getName() == "x/z/(x/z)", "wrong name")
       self.failUnless(((x/z)/(y**3)).getName() == "x/z/y^3", "wrong name")

       self.failUnless(((y**3)/x).getName() == "y^3/x", "wrong name")
       self.failUnless(((y**3)/(x*z)).getName() == "y^3/(xz)", "wrong name")
       self.failUnless(((y**3)/(x/z)).getName() == "y^3/(x/z)", "wrong name")
       self.failUnless(((y**3)/(y**3)).getName() == "y^3/y^3", "wrong name")

       self.failUnless((x**2).getName() == "x^2", "wrong name")
       self.failUnless(((x*z)**2).getName() == "(xz)^2", "wrong name")
       self.failUnless(((x/z)**2).getName() == "(x/z)^2", "wrong name")
       self.failUnless(((y**3)**2).getName() == "(y^3)^2", "wrong name")

       self.failUnless(((x*x)).getLongName() == "X*X", "wrong long name")
       self.failUnless((x*(x*z)).getLongName() == "X*X*Z", "wrong long name")
       self.failUnless((x*(x/z)).getLongName() == "X*X/Z", "wrong long name")
       self.failUnless((x*(y**3)).getLongName() == "X*Y^3", "wrong long name")

       self.failUnless(((x*z)*x).getLongName() == "X*Z*X", "wrong long name")
       self.failUnless(((x*z)*(x*z)).getLongName() == "X*Z*X*Z", "wrong long name")
       self.failUnless(((x*z)*(x/z)).getLongName() == "X*Z*X/Z", "wrong long name")
       self.failUnless(((x*z)*(y**3)).getLongName() == "X*Z*Y^3", "wrong long name")

       self.failUnless(((x/z)*x).getLongName() == "X/Z*X", "wrong long name")
       self.failUnless(((x/z)*(x*z)).getLongName() == "X/Z*X*Z", "wrong long name")
       self.failUnless(((x/z)*(x/z)).getLongName() == "X/Z*X/Z", "wrong long name")
       self.failUnless(((x/z)*(y**3)).getLongName() == "X/Z*Y^3", "wrong long name")

       self.failUnless(((y**3)*x).getLongName() == "Y^3*X", "wrong long name")
       self.failUnless(((y**3)*(x*z)).getLongName() == "Y^3*X*Z", "wrong long name")
       self.failUnless(((y**3)*(x/z)).getLongName() == "Y^3*X/Z", "wrong long name")
       self.failUnless(((y**3)*(y**3)).getLongName() == "Y^3*Y^3", "wrong long name")

       self.failUnless((x/x).getLongName() == "X/X", "wrong long name")
       self.failUnless((x/(x*z)).getLongName() == "X/(X*Z)", "wrong long name")
       self.failUnless((x/(x/z)).getLongName() == "X/(X/Z)", "wrong long name")
       self.failUnless((x/(y**3)).getLongName() == "X/Y^3", "wrong long name")

       self.failUnless(((x*z)/x).getLongName() == "X*Z/X", "wrong long name")
       self.failUnless(((x*z)/(x*z)).getLongName() == "X*Z/(X*Z)", "wrong long name")
       self.failUnless(((x*z)/(x/z)).getLongName() == "X*Z/(X/Z)", "wrong long name")
       self.failUnless(((x*z)/(y**3)).getLongName() == "X*Z/Y^3", "wrong long name")

       self.failUnless(((x/z)/x).getLongName() == "X/Z/X", "wrong long name")
       self.failUnless(((x/z)/(x*z)).getLongName() == "X/Z/(X*Z)", "wrong long name")
       self.failUnless(((x/z)/(x/z)).getLongName() == "X/Z/(X/Z)", "wrong long name")
       self.failUnless(((x/z)/(y**3)).getLongName() == "X/Z/Y^3", "wrong long name")

       self.failUnless(((y**3)/x).getLongName() == "Y^3/X", "wrong long name")
       self.failUnless(((y**3)/(x*z)).getLongName() == "Y^3/(X*Z)", "wrong long name")
       self.failUnless(((y**3)/(x/z)).getLongName() == "Y^3/(X/Z)", "wrong long name")
       self.failUnless(((y**3)/(y**3)).getLongName() == "Y^3/Y^3", "wrong long name")

       self.failUnless((x**2).getLongName() == "X^2", "wrong long name")
       self.failUnless(((x*z)**2).getLongName() == "(X*Z)^2", "wrong long name")
       self.failUnless(((x/z)**2).getLongName() == "(X/Z)^2", "wrong long name")
       self.failUnless(((y**3)**2).getLongName() == "(Y^3)^2", "wrong long name")

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
       self.failUnless(abs( atm- 101325.024) <= self.TOL*101325.024 , " Pa wrong.")
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
