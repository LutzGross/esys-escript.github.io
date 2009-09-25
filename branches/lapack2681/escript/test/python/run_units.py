
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
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


       self.failUnless(Yotta.getName() == "Y","Yotta name wrong")
       self.failUnless(Zetta.getName() == "Z","Zetta name wrong")
       self.failUnless(Exa.getName() == "E","Exa name wrong")
       self.failUnless(Peta.getName() == "P","Peta name wrong")
       self.failUnless(Tera.getName() == "T","Tera name wrong")
       self.failUnless(Giga.getName() == "G","Giga name wrong")
       self.failUnless(Mega.getName() == "M","Mega name wrong")
       self.failUnless(Kilo.getName() == "k","Kilo name wrong")
       self.failUnless(Hecto.getName() == "h","Hecto name wrong")
       self.failUnless(Deca.getName() == "da","Deca name wrong")
       self.failUnless(Deci.getName() == "d","Deci name wrong")
       self.failUnless(Centi.getName() == "c","Centi name wrong")
       self.failUnless(Milli.getName() == "m","Milli name wrong")
       self.failUnless(Micro.getName() == "mu","Micro name wrong")
       self.failUnless(Nano.getName() == "n","Nano name wrong")
       self.failUnless(Pico.getName() == "p","Pico name wrong")
       self.failUnless(Femto.getName() == "f","Femto name wrong")
       self.failUnless(Atto.getName() == "a","Atto name wrong")
       self.failUnless(Zepto.getName() == "z","Zepto name wrong")
       self.failUnless(Yocto.getName() == "y","Yocto name wrong")

       self.failUnless(Yotta.getLongName() == "Yotta","wrong long name")
       self.failUnless(Zetta.getLongName() == "Zetta","wrong long name")
       self.failUnless(Exa.getLongName() == "Exa","wrong long name")
       self.failUnless(Peta.getLongName() == "Peta","wrong long name")
       self.failUnless(Tera.getLongName() == "Tera","wrong long name")
       self.failUnless(Giga.getLongName() == "Giga","wrong long name")
       self.failUnless(Mega.getLongName() == "Mega","wrong long name")
       self.failUnless(Kilo.getLongName() == "Kilo","wrong long name")
       self.failUnless(Hecto.getLongName() == "Hecto","wrong long name")
       self.failUnless(Deca.getLongName() == "Deca","wrong long name")
       self.failUnless(Deci.getLongName() == "Deci","wrong long name")
       self.failUnless(Centi.getLongName() == "Centi","wrong long name")
       self.failUnless(Milli.getLongName() == "Milli","wrong long name")
       self.failUnless(Micro.getLongName() == "Micro","wrong long name")
       self.failUnless(Nano.getLongName() == "Nano","wrong long name")
       self.failUnless(Pico.getLongName() == "Pico","wrong long name")
       self.failUnless(Femto.getLongName() == "Femto","wrong long name")
       self.failUnless(Atto.getLongName() == "Atto","wrong long name")
       self.failUnless(Zepto.getLongName() == "Zepto","wrong long name")
       self.failUnless(Yocto.getLongName() == "Yocto","wrong long name")

       self.failUnless(abs(Yotta(1.)/1.e24 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Zetta(1.)/1.e21 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Exa(1.)/1.e18 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Peta(1.)/1.e15 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Tera(1.)/1.e12 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Giga(1.)/1.e9 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Mega(1.)/1.e6 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Kilo(1.)/1.e3 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Hecto(1.)/1.e2 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Deca(1.)/1.e1 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Deci(1.)/1.e-1 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Centi(1.)/1.e-2 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Milli(1.)/1.e-3 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Micro(1.)/1.e-6 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Nano(1.)/1.e-9 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Pico(1.)/1.e-12 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Femto(1.)/1.e-15 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Atto(1.)/1.e-18 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Zepto(1.)/1.e-21 -1.)  <= self.TOL , "wrong b value")
       self.failUnless(abs(Yocto(1.)/1.e-24 -1.)  <= self.TOL , "wrong b value")

       self.failUnless(abs( Yotta(1.) * Yocto(1.) - 1.) <= self.TOL , " Yotta or Yocto wrong.")
       self.failUnless(Yotta(0.) == 0. , " Yotta a value wrong.")
       self.failUnless(Yocto(0.) == 0. , " Yocto a value wrong.")
       self.failUnless(abs( Zetta(1.) * Zepto(1.) - 1.) <= self.TOL , "  Zetta or Zepto  wrong.")
       self.failUnless(Zetta(0.) == 0. , "  Zetta a value  wrong.")
       self.failUnless(Zepto(0.) == 0. , "  Zepto a value  wrong.")
       self.failUnless(abs( Exa(1.) * Atto(1.) - 1.) <= self.TOL , " Exa or Atto  wrong.")
       self.failUnless(Exa(0.) == 0. , " Exa a value  wrong.")
       self.failUnless(Atto(0.) == 0. , "  Atto a value  wrong.")
       self.failUnless(abs( Peta(1.) *Femto(1.) - 1.) <= self.TOL , " Peta or Femto  wrong.")
       self.failUnless(Peta(0.) == 0. , " Peta  a value wrong.")
       self.failUnless(Femto(0.) == 0. , " Femto a value  wrong.")
       self.failUnless(abs( Tera(1.) * Pico(1.) - 1.) <= self.TOL , " Tera or Pico  wrong.")
       self.failUnless(Tera(0.) == 0. , " Tera a value wrong.")
       self.failUnless(Pico(0.) == 0. , " Pico a value  wrong.")
       self.failUnless(abs( Giga(1.) * Nano(1.) - 1.) <= self.TOL , "  Giga or Nano  wrong.")
       self.failUnless(Giga(0.) == 0. , "  Giga a value  wrong.")
       self.failUnless(Nano(0.) == 0. , "  Nano a value  wrong.")
       self.failUnless(abs( Mega(1.) *Micro(1.) - 1.) <= self.TOL , " Mega or Micro  wrong.")
       self.failUnless(Mega(0.) == 0. , "Mega a value wrong.")
       self.failUnless(Micro(0.) == 0. , " Micro a value  wrong.")
       self.failUnless(abs( Kilo(1.) * Milli(1.) - 1.) <= self.TOL , "Kilo or Milli   wrong.")
       self.failUnless( Kilo(0.) == 0. , "Kilo a value wrong.")
       self.failUnless(Milli(0.) == 0. , "Milli a value   wrong.")
       self.failUnless(abs( Hecto(1.) * Centi(1.)- 1.) <= self.TOL , " Hecto or Centi  wrong.")
       self.failUnless(Hecto(0.) == 0., " Hecto a value wrong.")
       self.failUnless(Centi(0.)== 0. , " Centi a value  wrong.")
       self.failUnless(abs( Deca(1.) * Deci(1.) - 1.) <= self.TOL , "Deca or Deci  wrong.")
       self.failUnless(Deca(0.) == 0. , "Deca a value wrong.")
       self.failUnless(Deci(0.) == 0. , "Deci a value  wrong.")
    def testLength(self):
       self.failUnless(m.getName() == "m", "meter name wrong.")
       self.failUnless(m.getLongName() == "meter", "meter long name wrong.")
       self.failUnless(abs( m(1.) - 1.) <= self.TOL , " meter wrong.")
       self.failUnless(abs( km(1.) - 1000.) <= self.TOL*1000. , " km wrong.")
       self.failUnless(abs( cm(1.) - 1./100.) <= self.TOL/100. , " cm wrong.")
       self.failUnless(abs( mm(1.) - 1./1000.) <= self.TOL/1000. , " mm wrong.")
       self.failUnless(m(0.) == 0., " meter wrong.")
       self.failUnless(km(0.) == 0., " km wrong.")
       self.failUnless(cm(0.) == 0. , " cm wrong.")
       self.failUnless(mm(0.) ==0.," mm wrong.")
    def testTime(self):
       self.failUnless(sec.getName() == "sec" , "sec name wrong")
       self.failUnless(sec.getLongName() == "second" , "sec long name wrong")
       self.failUnless(minute.getName() == "min" , "minute name wrong")
       self.failUnless(minute.getLongName() == "minute" , "minute long name wrong")
       self.failUnless(h.getName() == "h" , "hour name wrong")
       self.failUnless(h.getLongName() == "hour" , "hour long name wrong")
       self.failUnless(day.getName() == "d" , "day name wrong")
       self.failUnless(day.getLongName() == "day" , "day long name wrong")
       self.failUnless(yr.getName() == "yr" , "year name wrong")
       self.failUnless(yr.getLongName() == "year" , "year long name wrong")

       self.failUnless(abs( sec(1.) - 1.) <= self.TOL , " second wrong.")
       self.failUnless(abs( minute(1.) - 60.) <= self.TOL*60 , " minute wrong.")
       self.failUnless(abs( h(1.) - 3600) <= self.TOL*3600., " h wrong.")
       self.failUnless(abs( day(1.) - 86400.) <= self.TOL*86400 , " day wrong.")
       self.failUnless(abs( yr(1.) - 31556952) <= self.TOL*31556952 , " yr wrong.")
       self.failUnless(abs( Myr(1.) - 31556952.e6) <= self.TOL*31556952.e6 , " Myr wrong.")
       self.failUnless(abs( Gyr(1.) - 31556952.e9) <= self.TOL*31556952.e9 , " Gyr wrong.")

       self.failUnless(sec(0.) == 0. , " second wrong.")
       self.failUnless(minute(0.) == 0. , " minute wrong.")
       self.failUnless(h(0.) == 0., " h wrong.")
       self.failUnless(day(0.) == 0. , " day wrong.")
       self.failUnless(yr(0.) == 0. , " yr wrong.")
       self.failUnless(Myr(0.) == 0. , " Myr wrong.")
       self.failUnless(Gyr(0.) == 0., " Gyr wrong.")
    def testMass(self):
       self.failUnless(kg.getName() == "kg" , " kg wrong.")
       self.failUnless(lb.getName() == "lb" , " lb wrong.")
       self.failUnless(kg.getLongName() == "kg" , " kg wrong.")
       self.failUnless(lb.getLongName() == "pound" , " lb wrong.")

       self.failUnless(abs( gram(1.) - 1./1000.) <= self.TOL/1000. , " gram wrong.")
       self.failUnless(abs( kg(1.) - 1.) <= self.TOL , " kg wrong.")
       self.failUnless(abs( ton(1.) - 1.e3) <= self.TOL*1.e3 , " ton wrong.")
       self.failUnless(abs( lb(1.) - 0.45359237) <= self.TOL , " lb wrong.")

       self.failUnless(gram(0.) == 0. , " gram wrong.")
       self.failUnless(kg(0.) == 0. , " kg wrong.")
       self.failUnless(ton(0.) == 0. , " ton wrong.")
       self.failUnless(lb(0.) == 0. , " lb wrong.")

    def testCurrent(self):
       self.failUnless(A.getName() == "A" , " A wrong.")
       self.failUnless(A.getLongName() == "Ampere" , " A wrong.")
       self.failUnless(A(0.) == 0. , " A wrong.")
       self.failUnless(abs( A(1.)- 1.) <= self.TOL , " A wrong.")
    def testFrequency(self):
       self.failUnless(Hz(0.) == 0. , " Hz wrong.")
       self.failUnless(abs( Hz(1)- 1.) <= self.TOL , " Hz wrong.")

    def testForce(self):
       self.failUnless(N.getName() == "N" , "N wrong.")
       self.failUnless(N.getLongName() == "Newton" , " N wrong.")
       self.failUnless(N(0.) == 0. , " N wrong.")
       self.failUnless(abs( N(1)- 1.) <= self.TOL , " N wrong.")
    def testPressure(self):
       self.failUnless(Pa.getName() == "Pa" , " Pa wrong.")
       self.failUnless(Pa.getLongName() == "Pascal" , " Pa wrong.")
       self.failUnless(Pa(0.) == 0. , " Pa wrong.")
       self.failUnless(abs( Pa(1)- 1.) <= self.TOL , " Pa wrong.")

       self.failUnless(atm.getName() == "atm" , " atm wrong.")
       self.failUnless(atm.getLongName() == "atmosphere" , " atm wrong.")
       self.failUnless(atm(0.) == 0. , " atm wrong.")
       self.failUnless(abs( atm(1)- 101325.024) <= self.TOL*101325.024 , " atm wrong.")

    def testEnegry(self):
       self.failUnless(J.getName() == "J" , " J wrong.")
       self.failUnless(J.getLongName() == "Joule" , " J wrong.")
       self.failUnless(J(0.) == 0. , " J wrong.")
       self.failUnless(abs( J(1)- 1.) <= self.TOL , " J wrong.")

    def testPower(self):
       self.failUnless(W.getName() == "W" , " W wrong.")
       self.failUnless(W.getLongName() == "Watt" , " W wrong.")
       self.failUnless(W(0.) == 0. , " W wrong.")
       self.failUnless(abs( W(1)- 1.) <= self.TOL , " W wrong.")

    def testPotential(self):
       self.failUnless(V.getName() == "V" , " V wrong.")
       self.failUnless(V.getLongName() == "Volt" , " V wrong.")
       self.failUnless(V(0.) == 0. , " V wrong.")
       self.failUnless(abs( V(1)- 1.) <= self.TOL , " V wrong.")

    def testCharge(self):
       self.failUnless(C.getName() == "C" , " C wrong.")
       self.failUnless(C.getLongName() == "Coulomb" , " C wrong.")
       self.failUnless(C(0.) == 0. , " C wrong.")
       self.failUnless(abs( C(1)- 1.) <= self.TOL , " C wrong.")

    def testResistance(self):
       self.failUnless(Ohm.getName() == "Ohm" , " Ohm wrong.")
       self.failUnless(Ohm.getLongName() == "Ohm" , " Ohm wrong.")
       self.failUnless(Ohm(0.) == 0. , " Ohm wrong.")
       self.failUnless(abs( Ohm(1)- 1.) <= self.TOL , " Ohm wrong.")
    def testTemperature(self):
       self.failUnless(K.getName() == "K" , " Kelvin wrong.")
       self.failUnless(K.getLongName() == "Kelvin" , " Kelvin wrong.")
       self.failUnless(K(0.) == 0. , " Kelvin wrong.")
       self.failUnless(abs( K(1.)- 1.) <= self.TOL , " Kelvin wrong.")

       self.failUnless(Celsius.getName() == "C" , " Celsius wrong.")
       self.failUnless(Celsius.getLongName() == "Celsius" , " Celsius wrong.")
       self.failUnless(abs(Celsius(20.) - 293.15) <= self.TOL*293.15 , " Celsius 20 conversion wrong.")
       self.failUnless(abs(Celsius(-2.)- 271.15) <= self.TOL*271.15 , " Celsius conversion wrong.")
       self.failUnless(abs(Celsius(100.)- 373.15) <= self.TOL*373.15 , " Celsius conversion wrong.")

       self.failUnless(Fahrenheit.getName() == "F" , " Fahrenheit wrong.")
       self.failUnless(Fahrenheit.getLongName() == "Fahrenheit" , " Fahrenheit wrong.")
       self.failUnless(abs( Fahrenheit(0.)- 255.3722222) <= 1.e-8*255.3722222, " Fahrenheit wrong.")
       self.failUnless(abs( Fahrenheit(20.)- 266.4833333) <=  1.e-8*266.4833333 , " Fahrenheit wrong.")
       self.failUnless(abs( Fahrenheit(100.)- 310.9277778 ) <= 1.e-8 * 310.9277778 , " Fahrenheit wrong.")
       self.failUnless(abs( Fahrenheit(-1.)- 254.8166667) <= 1.e-8 *254.8166667, "Fahrenheit wrong.")


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
