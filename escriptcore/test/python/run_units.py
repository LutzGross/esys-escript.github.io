
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


__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
import sys
from esys.escript.unitsSI import *
from esys.escript.util import EPSILON

class Test_UnitsSI(unittest.TestCase):
    TOL=EPSILON*100.
    def testUnit(self):
       s=Unit("s","something",1.,4.)
       self.assertTrue(s.getName() == "s", "wrong name")
       s.setName("t")
       self.assertTrue(s.getName() == "t", "wrong reset name")
       self.assertTrue(s.getLongName() == "something", "wrong long name")
       s.setLongName("test")
       self.assertTrue(s.getLongName() == "test", "wrong reset long name")
       self.assertTrue(abs( s(3.)-13.) <= self.TOL*13. , " s(3.) wrong.")
       self.assertTrue(abs( 3.*s-13.) <= self.TOL*13. , " 13.s wrong.")
       self.assertTrue(abs( 10.*s-41.) <= self.TOL*41. , " 10.s wrong.")
       self.assertTrue(abs( 41./s-10.) <= self.TOL*10. , " 10/s wrong.")
       self.assertTrue(abs( 13./s-3.) <= self.TOL*3. , " 13/s wrong.")

       x=Unit("x","X",10.,2.)
       y=Unit("y","Y",0.,4.)
       z=Unit("z","Z",0.,0.5)
       p_m=x*s
       self.assertTrue(p_m.getName() == "xt", "wrong p_m name")
       self.assertTrue(p_m.getLongName() == "X*test", "wrong p_m long name")
       self.assertTrue(abs(p_m(0.)-12.) <= self.TOL*12 , " p_m(0) wrong.")
       self.assertTrue(abs(p_m(1.)-20.) <= self.TOL*20 , " p_m(1) wrong.")

       p_d=y/z
       self.assertTrue(p_d.getName() == "y/z", "wrong p_d name")
       self.assertTrue(p_d.getLongName() == "Y/Z", "wrong p_d long name")
       self.assertTrue(abs(p_d(1.)-8.) <= self.TOL*8. , " p_d(0) wrong.")
       self.assertTrue(abs(p_d(2.)-16.) <= self.TOL*16. , " p_d(1) wrong.")

       p_p=y**3
       self.assertTrue(p_p.getName() == "y^3", "p_p mult name")
       self.assertTrue(p_p.getLongName() == "Y^3", "wrong p_p long name")
       self.assertTrue(abs(p_p(1.)-64.) <= self.TOL*64. , " p_p(0) wrong.")
       self.assertTrue(abs(p_p(.5)-32) <= self.TOL*32. , " p_p(1) wrong.")

        
       x=Unit("x","X",0.,2.)
       y=Unit("y","Y",0.,4.)
       z=Unit("z","Z",0.,0.5)

       self.assertTrue(((x*x)).getName() == "xx", "wrong name")
       self.assertTrue((x*(x*z)).getName() == "xxz", "wrong name")
       self.assertTrue((x*(x/z)).getName() == "xx/z", "wrong name")
       self.assertTrue((x*(y**3)).getName() == "x y^3", "wrong name")

       self.assertTrue(((x*z)*x).getName() == "xzx", "wrong name")
       self.assertTrue(((x*z)*(x*z)).getName() == "xzxz", "wrong name")
       self.assertTrue(((x*z)*(x/z)).getName() == "xzx/z", "wrong name")
       self.assertTrue(((x*z)*(y**3)).getName() == "xz y^3", "wrong name")

       self.assertTrue(((x/z)*x).getName() == "x/z x", "wrong name")
       self.assertTrue(((x/z)*(x*z)).getName() == "x/z xz", "wrong name")
       self.assertTrue(((x/z)*(x/z)).getName() == "x/z x/z", "wrong name")
       self.assertTrue(((x/z)*(y**3)).getName() == "x/z y^3", "wrong name")

       self.assertTrue(((y**3)*x).getName() == "y^3 x", "wrong name")
       self.assertTrue(((y**3)*(x*z)).getName() == "y^3 xz", "wrong name")
       self.assertTrue(((y**3)*(x/z)).getName() == "y^3 x/z", "wrong name")
       self.assertTrue(((y**3)*(y**3)).getName() == "y^3 y^3", "wrong name")

       self.assertTrue((x/x).getName() == "x/x", "wrong name")
       self.assertTrue((x/(x*z)).getName() == "x/(xz)", "wrong name")
       self.assertTrue((x/(x/z)).getName() == "x/(x/z)", "wrong name")
       self.assertTrue((x/(y**3)).getName() == "x/y^3", "wrong name")

       self.assertTrue(((x*z)/x).getName() == "xz/x", "wrong name")
       self.assertTrue(((x*z)/(x*z)).getName() == "xz/(xz)", "wrong name")
       self.assertTrue(((x*z)/(x/z)).getName() == "xz/(x/z)", "wrong name")
       self.assertTrue(((x*z)/(y**3)).getName() == "xz/y^3", "wrong name")

       self.assertTrue(((x/z)/x).getName() == "x/z/x", "wrong name")
       self.assertTrue(((x/z)/(x*z)).getName() == "x/z/(xz)", "wrong name")
       self.assertTrue(((x/z)/(x/z)).getName() == "x/z/(x/z)", "wrong name")
       self.assertTrue(((x/z)/(y**3)).getName() == "x/z/y^3", "wrong name")

       self.assertTrue(((y**3)/x).getName() == "y^3/x", "wrong name")
       self.assertTrue(((y**3)/(x*z)).getName() == "y^3/(xz)", "wrong name")
       self.assertTrue(((y**3)/(x/z)).getName() == "y^3/(x/z)", "wrong name")
       self.assertTrue(((y**3)/(y**3)).getName() == "y^3/y^3", "wrong name")

       self.assertTrue((x**2).getName() == "x^2", "wrong name")
       self.assertTrue(((x*z)**2).getName() == "(xz)^2", "wrong name")
       self.assertTrue(((x/z)**2).getName() == "(x/z)^2", "wrong name")
       self.assertTrue(((y**3)**2).getName() == "(y^3)^2", "wrong name")

       self.assertTrue(((x*x)).getLongName() == "X*X", "wrong long name")
       self.assertTrue((x*(x*z)).getLongName() == "X*X*Z", "wrong long name")
       self.assertTrue((x*(x/z)).getLongName() == "X*X/Z", "wrong long name")
       self.assertTrue((x*(y**3)).getLongName() == "X*Y^3", "wrong long name")

       self.assertTrue(((x*z)*x).getLongName() == "X*Z*X", "wrong long name")
       self.assertTrue(((x*z)*(x*z)).getLongName() == "X*Z*X*Z", "wrong long name")
       self.assertTrue(((x*z)*(x/z)).getLongName() == "X*Z*X/Z", "wrong long name")
       self.assertTrue(((x*z)*(y**3)).getLongName() == "X*Z*Y^3", "wrong long name")

       self.assertTrue(((x/z)*x).getLongName() == "X/Z*X", "wrong long name")
       self.assertTrue(((x/z)*(x*z)).getLongName() == "X/Z*X*Z", "wrong long name")
       self.assertTrue(((x/z)*(x/z)).getLongName() == "X/Z*X/Z", "wrong long name")
       self.assertTrue(((x/z)*(y**3)).getLongName() == "X/Z*Y^3", "wrong long name")

       self.assertTrue(((y**3)*x).getLongName() == "Y^3*X", "wrong long name")
       self.assertTrue(((y**3)*(x*z)).getLongName() == "Y^3*X*Z", "wrong long name")
       self.assertTrue(((y**3)*(x/z)).getLongName() == "Y^3*X/Z", "wrong long name")
       self.assertTrue(((y**3)*(y**3)).getLongName() == "Y^3*Y^3", "wrong long name")

       self.assertTrue((x/x).getLongName() == "X/X", "wrong long name")
       self.assertTrue((x/(x*z)).getLongName() == "X/(X*Z)", "wrong long name")
       self.assertTrue((x/(x/z)).getLongName() == "X/(X/Z)", "wrong long name")
       self.assertTrue((x/(y**3)).getLongName() == "X/Y^3", "wrong long name")

       self.assertTrue(((x*z)/x).getLongName() == "X*Z/X", "wrong long name")
       self.assertTrue(((x*z)/(x*z)).getLongName() == "X*Z/(X*Z)", "wrong long name")
       self.assertTrue(((x*z)/(x/z)).getLongName() == "X*Z/(X/Z)", "wrong long name")
       self.assertTrue(((x*z)/(y**3)).getLongName() == "X*Z/Y^3", "wrong long name")

       self.assertTrue(((x/z)/x).getLongName() == "X/Z/X", "wrong long name")
       self.assertTrue(((x/z)/(x*z)).getLongName() == "X/Z/(X*Z)", "wrong long name")
       self.assertTrue(((x/z)/(x/z)).getLongName() == "X/Z/(X/Z)", "wrong long name")
       self.assertTrue(((x/z)/(y**3)).getLongName() == "X/Z/Y^3", "wrong long name")

       self.assertTrue(((y**3)/x).getLongName() == "Y^3/X", "wrong long name")
       self.assertTrue(((y**3)/(x*z)).getLongName() == "Y^3/(X*Z)", "wrong long name")
       self.assertTrue(((y**3)/(x/z)).getLongName() == "Y^3/(X/Z)", "wrong long name")
       self.assertTrue(((y**3)/(y**3)).getLongName() == "Y^3/Y^3", "wrong long name")

       self.assertTrue((x**2).getLongName() == "X^2", "wrong long name")
       self.assertTrue(((x*z)**2).getLongName() == "(X*Z)^2", "wrong long name")
       self.assertTrue(((x/z)**2).getLongName() == "(X/Z)^2", "wrong long name")
       self.assertTrue(((y**3)**2).getLongName() == "(Y^3)^2", "wrong long name")

    def testPrefix(self):


       self.assertTrue(Yotta.getName() == "Y","Yotta name wrong")
       self.assertTrue(Zetta.getName() == "Z","Zetta name wrong")
       self.assertTrue(Exa.getName() == "E","Exa name wrong")
       self.assertTrue(Peta.getName() == "P","Peta name wrong")
       self.assertTrue(Tera.getName() == "T","Tera name wrong")
       self.assertTrue(Giga.getName() == "G","Giga name wrong")
       self.assertTrue(Mega.getName() == "M","Mega name wrong")
       self.assertTrue(Kilo.getName() == "k","Kilo name wrong")
       self.assertTrue(Hecto.getName() == "h","Hecto name wrong")
       self.assertTrue(Deca.getName() == "da","Deca name wrong")
       self.assertTrue(Deci.getName() == "d","Deci name wrong")
       self.assertTrue(Centi.getName() == "c","Centi name wrong")
       self.assertTrue(Milli.getName() == "m","Milli name wrong")
       self.assertTrue(Micro.getName() == "mu","Micro name wrong")
       self.assertTrue(Nano.getName() == "n","Nano name wrong")
       self.assertTrue(Pico.getName() == "p","Pico name wrong")
       self.assertTrue(Femto.getName() == "f","Femto name wrong")
       self.assertTrue(Atto.getName() == "a","Atto name wrong")
       self.assertTrue(Zepto.getName() == "z","Zepto name wrong")
       self.assertTrue(Yocto.getName() == "y","Yocto name wrong")

       self.assertTrue(Yotta.getLongName() == "Yotta","wrong long name")
       self.assertTrue(Zetta.getLongName() == "Zetta","wrong long name")
       self.assertTrue(Exa.getLongName() == "Exa","wrong long name")
       self.assertTrue(Peta.getLongName() == "Peta","wrong long name")
       self.assertTrue(Tera.getLongName() == "Tera","wrong long name")
       self.assertTrue(Giga.getLongName() == "Giga","wrong long name")
       self.assertTrue(Mega.getLongName() == "Mega","wrong long name")
       self.assertTrue(Kilo.getLongName() == "Kilo","wrong long name")
       self.assertTrue(Hecto.getLongName() == "Hecto","wrong long name")
       self.assertTrue(Deca.getLongName() == "Deca","wrong long name")
       self.assertTrue(Deci.getLongName() == "Deci","wrong long name")
       self.assertTrue(Centi.getLongName() == "Centi","wrong long name")
       self.assertTrue(Milli.getLongName() == "Milli","wrong long name")
       self.assertTrue(Micro.getLongName() == "Micro","wrong long name")
       self.assertTrue(Nano.getLongName() == "Nano","wrong long name")
       self.assertTrue(Pico.getLongName() == "Pico","wrong long name")
       self.assertTrue(Femto.getLongName() == "Femto","wrong long name")
       self.assertTrue(Atto.getLongName() == "Atto","wrong long name")
       self.assertTrue(Zepto.getLongName() == "Zepto","wrong long name")
       self.assertTrue(Yocto.getLongName() == "Yocto","wrong long name")

       self.assertTrue(abs(Yotta(1.)/1.e24 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Zetta(1.)/1.e21 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Exa(1.)/1.e18 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Peta(1.)/1.e15 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Tera(1.)/1.e12 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Giga(1.)/1.e9 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Mega(1.)/1.e6 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Kilo(1.)/1.e3 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Hecto(1.)/1.e2 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Deca(1.)/1.e1 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Deci(1.)/1.e-1 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Centi(1.)/1.e-2 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Milli(1.)/1.e-3 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Micro(1.)/1.e-6 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Nano(1.)/1.e-9 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Pico(1.)/1.e-12 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Femto(1.)/1.e-15 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Atto(1.)/1.e-18 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Zepto(1.)/1.e-21 -1.)  <= self.TOL , "wrong b value")
       self.assertTrue(abs(Yocto(1.)/1.e-24 -1.)  <= self.TOL , "wrong b value")

       self.assertTrue(abs( Yotta(1.) * Yocto(1.) - 1.) <= self.TOL , " Yotta or Yocto wrong.")
       self.assertTrue(Yotta(0.) == 0. , " Yotta a value wrong.")
       self.assertTrue(Yocto(0.) == 0. , " Yocto a value wrong.")
       self.assertTrue(abs( Zetta(1.) * Zepto(1.) - 1.) <= self.TOL , "  Zetta or Zepto  wrong.")
       self.assertTrue(Zetta(0.) == 0. , "  Zetta a value  wrong.")
       self.assertTrue(Zepto(0.) == 0. , "  Zepto a value  wrong.")
       self.assertTrue(abs( Exa(1.) * Atto(1.) - 1.) <= self.TOL , " Exa or Atto  wrong.")
       self.assertTrue(Exa(0.) == 0. , " Exa a value  wrong.")
       self.assertTrue(Atto(0.) == 0. , "  Atto a value  wrong.")
       self.assertTrue(abs( Peta(1.) *Femto(1.) - 1.) <= self.TOL , " Peta or Femto  wrong.")
       self.assertTrue(Peta(0.) == 0. , " Peta  a value wrong.")
       self.assertTrue(Femto(0.) == 0. , " Femto a value  wrong.")
       self.assertTrue(abs( Tera(1.) * Pico(1.) - 1.) <= self.TOL , " Tera or Pico  wrong.")
       self.assertTrue(Tera(0.) == 0. , " Tera a value wrong.")
       self.assertTrue(Pico(0.) == 0. , " Pico a value  wrong.")
       self.assertTrue(abs( Giga(1.) * Nano(1.) - 1.) <= self.TOL , "  Giga or Nano  wrong.")
       self.assertTrue(Giga(0.) == 0. , "  Giga a value  wrong.")
       self.assertTrue(Nano(0.) == 0. , "  Nano a value  wrong.")
       self.assertTrue(abs( Mega(1.) *Micro(1.) - 1.) <= self.TOL , " Mega or Micro  wrong.")
       self.assertTrue(Mega(0.) == 0. , "Mega a value wrong.")
       self.assertTrue(Micro(0.) == 0. , " Micro a value  wrong.")
       self.assertTrue(abs( Kilo(1.) * Milli(1.) - 1.) <= self.TOL , "Kilo or Milli   wrong.")
       self.assertTrue( Kilo(0.) == 0. , "Kilo a value wrong.")
       self.assertTrue(Milli(0.) == 0. , "Milli a value   wrong.")
       self.assertTrue(abs( Hecto(1.) * Centi(1.)- 1.) <= self.TOL , " Hecto or Centi  wrong.")
       self.assertTrue(Hecto(0.) == 0., " Hecto a value wrong.")
       self.assertTrue(Centi(0.)== 0. , " Centi a value  wrong.")
       self.assertTrue(abs( Deca(1.) * Deci(1.) - 1.) <= self.TOL , "Deca or Deci  wrong.")
       self.assertTrue(Deca(0.) == 0. , "Deca a value wrong.")
       self.assertTrue(Deci(0.) == 0. , "Deci a value  wrong.")
    def testLength(self):
       self.assertTrue(m.getName() == "m", "meter name wrong.")
       self.assertTrue(m.getLongName() == "meter", "meter long name wrong.")
       self.assertTrue(abs( m(1.) - 1.) <= self.TOL , " meter wrong.")
       self.assertTrue(abs( km(1.) - 1000.) <= self.TOL*1000. , " km wrong.")
       self.assertTrue(abs( cm(1.) - 1./100.) <= self.TOL/100. , " cm wrong.")
       self.assertTrue(abs( mm(1.) - 1./1000.) <= self.TOL/1000. , " mm wrong.")
       self.assertTrue(m(0.) == 0., " meter wrong.")
       self.assertTrue(km(0.) == 0., " km wrong.")
       self.assertTrue(cm(0.) == 0. , " cm wrong.")
       self.assertTrue(mm(0.) ==0.," mm wrong.")
    def testTime(self):
       self.assertTrue(sec.getName() == "sec" , "sec name wrong")
       self.assertTrue(sec.getLongName() == "second" , "sec long name wrong")
       self.assertTrue(minute.getName() == "min" , "minute name wrong")
       self.assertTrue(minute.getLongName() == "minute" , "minute long name wrong")
       self.assertTrue(h.getName() == "h" , "hour name wrong")
       self.assertTrue(h.getLongName() == "hour" , "hour long name wrong")
       self.assertTrue(day.getName() == "d" , "day name wrong")
       self.assertTrue(day.getLongName() == "day" , "day long name wrong")
       self.assertTrue(yr.getName() == "yr" , "year name wrong")
       self.assertTrue(yr.getLongName() == "year" , "year long name wrong")

       self.assertTrue(abs( sec(1.) - 1.) <= self.TOL , " second wrong.")
       self.assertTrue(abs( minute(1.) - 60.) <= self.TOL*60 , " minute wrong.")
       self.assertTrue(abs( h(1.) - 3600) <= self.TOL*3600., " h wrong.")
       self.assertTrue(abs( day(1.) - 86400.) <= self.TOL*86400 , " day wrong.")
       self.assertTrue(abs( yr(1.) - 31556952) <= self.TOL*31556952 , " yr wrong.")
       self.assertTrue(abs( Myr(1.) - 31556952.e6) <= self.TOL*31556952.e6 , " Myr wrong.")
       self.assertTrue(abs( Gyr(1.) - 31556952.e9) <= self.TOL*31556952.e9 , " Gyr wrong.")

       self.assertTrue(sec(0.) == 0. , " second wrong.")
       self.assertTrue(minute(0.) == 0. , " minute wrong.")
       self.assertTrue(h(0.) == 0., " h wrong.")
       self.assertTrue(day(0.) == 0. , " day wrong.")
       self.assertTrue(yr(0.) == 0. , " yr wrong.")
       self.assertTrue(Myr(0.) == 0. , " Myr wrong.")
       self.assertTrue(Gyr(0.) == 0., " Gyr wrong.")
    def testMass(self):
       self.assertTrue(kg.getName() == "kg" , " kg wrong.")
       self.assertTrue(lb.getName() == "lb" , " lb wrong.")
       self.assertTrue(kg.getLongName() == "kg" , " kg wrong.")
       self.assertTrue(lb.getLongName() == "pound" , " lb wrong.")

       self.assertTrue(abs( gram(1.) - 1./1000.) <= self.TOL/1000. , " gram wrong.")
       self.assertTrue(abs( kg(1.) - 1.) <= self.TOL , " kg wrong.")
       self.assertTrue(abs( ton(1.) - 1.e3) <= self.TOL*1.e3 , " ton wrong.")
       self.assertTrue(abs( lb(1.) - 0.45359237) <= self.TOL , " lb wrong.")

       self.assertTrue(gram(0.) == 0. , " gram wrong.")
       self.assertTrue(kg(0.) == 0. , " kg wrong.")
       self.assertTrue(ton(0.) == 0. , " ton wrong.")
       self.assertTrue(lb(0.) == 0. , " lb wrong.")

    def testCurrent(self):
       self.assertTrue(A.getName() == "A" , " A wrong.")
       self.assertTrue(A.getLongName() == "Ampere" , " A wrong.")
       self.assertTrue(A(0.) == 0. , " A wrong.")
       self.assertTrue(abs( A(1.)- 1.) <= self.TOL , " A wrong.")
    def testFrequency(self):
       self.assertTrue(Hz(0.) == 0. , " Hz wrong.")
       self.assertTrue(abs( Hz(1)- 1.) <= self.TOL , " Hz wrong.")

    def testForce(self):
       self.assertTrue(N.getName() == "N" , "N wrong.")
       self.assertTrue(N.getLongName() == "Newton" , " N wrong.")
       self.assertTrue(N(0.) == 0. , " N wrong.")
       self.assertTrue(abs( N(1)- 1.) <= self.TOL , " N wrong.")
    def testPressure(self):
       self.assertTrue(Pa.getName() == "Pa" , " Pa wrong.")
       self.assertTrue(Pa.getLongName() == "Pascal" , " Pa wrong.")
       self.assertTrue(Pa(0.) == 0. , " Pa wrong.")
       self.assertTrue(abs( Pa(1)- 1.) <= self.TOL , " Pa wrong.")

       self.assertTrue(atm.getName() == "atm" , " atm wrong.")
       self.assertTrue(atm.getLongName() == "atmosphere" , " atm wrong.")
       self.assertTrue(atm(0.) == 0. , " atm wrong.")
       self.assertTrue(abs( atm(1)- 101325.024) <= self.TOL*101325.024 , " atm wrong.")

    def testEnegry(self):
       self.assertTrue(J.getName() == "J" , " J wrong.")
       self.assertTrue(J.getLongName() == "Joule" , " J wrong.")
       self.assertTrue(J(0.) == 0. , " J wrong.")
       self.assertTrue(abs( J(1)- 1.) <= self.TOL , " J wrong.")

    def testPower(self):
       self.assertTrue(W.getName() == "W" , " W wrong.")
       self.assertTrue(W.getLongName() == "Watt" , " W wrong.")
       self.assertTrue(W(0.) == 0. , " W wrong.")
       self.assertTrue(abs( W(1)- 1.) <= self.TOL , " W wrong.")

    def testPotential(self):
       self.assertTrue(V.getName() == "V" , " V wrong.")
       self.assertTrue(V.getLongName() == "Volt" , " V wrong.")
       self.assertTrue(V(0.) == 0. , " V wrong.")
       self.assertTrue(abs( V(1)- 1.) <= self.TOL , " V wrong.")

    def testCharge(self):
       self.assertTrue(C.getName() == "C" , " C wrong.")
       self.assertTrue(C.getLongName() == "Coulomb" , " C wrong.")
       self.assertTrue(C(0.) == 0. , " C wrong.")
       self.assertTrue(abs( C(1)- 1.) <= self.TOL , " C wrong.")

    def testResistance(self):
       self.assertTrue(Ohm.getName() == "Ohm" , " Ohm wrong.")
       self.assertTrue(Ohm.getLongName() == "Ohm" , " Ohm wrong.")
       self.assertTrue(Ohm(0.) == 0. , " Ohm wrong.")
       self.assertTrue(abs( Ohm(1)- 1.) <= self.TOL , " Ohm wrong.")
    def testTemperature(self):
       self.assertTrue(K.getName() == "K" , " Kelvin wrong.")
       self.assertTrue(K.getLongName() == "Kelvin" , " Kelvin wrong.")
       self.assertTrue(K(0.) == 0. , " Kelvin wrong.")
       self.assertTrue(abs( K(1.)- 1.) <= self.TOL , " Kelvin wrong.")

       self.assertTrue(Celsius.getName() == "C" , " Celsius wrong.")
       self.assertTrue(Celsius.getLongName() == "Celsius" , " Celsius wrong.")
       self.assertTrue(abs(Celsius(20.) - 293.15) <= self.TOL*293.15 , " Celsius 20 conversion wrong.")
       self.assertTrue(abs(Celsius(-2.)- 271.15) <= self.TOL*271.15 , " Celsius conversion wrong.")
       self.assertTrue(abs(Celsius(100.)- 373.15) <= self.TOL*373.15 , " Celsius conversion wrong.")

       self.assertTrue(Fahrenheit.getName() == "F" , " Fahrenheit wrong.")
       self.assertTrue(Fahrenheit.getLongName() == "Fahrenheit" , " Fahrenheit wrong.")
       self.assertTrue(abs( Fahrenheit(0.)- 255.3722222) <= 1.e-8*255.3722222, " Fahrenheit wrong.")
       self.assertTrue(abs( Fahrenheit(20.)- 266.4833333) <=  1.e-8*266.4833333 , " Fahrenheit wrong.")
       self.assertTrue(abs( Fahrenheit(100.)- 310.9277778 ) <= 1.e-8 * 310.9277778 , " Fahrenheit wrong.")
       self.assertTrue(abs( Fahrenheit(-1.)- 254.8166667) <= 1.e-8 *254.8166667, "Fahrenheit wrong.")


    def testEarth(self):
       R=6367444.65
       self.assertTrue(abs( R_Earth_equator-6378137) <= self.TOL*R , " R_Earth_equator wrong.")
       self.assertTrue(abs( R_Earth_poles-6356752.3) <= self.TOL*R , " R_Earth_poles wrong.")
       self.assertTrue(abs( R_Earth-R) <= self.TOL*R , "R_Earth wrong.")
    def testVlight(self):
       self.assertTrue(abs(v_light - 299792458.) <= self.TOL*299792458. , " speed of light is wrong.")

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
