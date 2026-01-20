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

"""
Tools supporting physical units and conversion.

:synopsis: some tools supporting physical units and conversion

:var Yotta: prefix yotta, symbol:   Y
:var Zetta: prefix zetta, symbol: Z
:var Exa: prefix exa, symbol: E
:var Peta: prefix peta, symbol: P
:var Tera: prefix tera, symbol: T
:var Giga: prefix giga, symbol: G
:var Mega: prefix mega, symbol: M
:var Kilo: prefix kilo, symbol: k
:var Hecto: prefix hecto, symbol: h
:var Deca:  prefix deca, symbol: da
:var Deci:  prefix deci, symbol: d
:var Centi: prefix centi, symbol: c
:var Milli: prefix milli, symbol: m
:var Micro: prefix micro, symbol: mu
:var Nano:  prefix nano, symbol: n
:var Pico:  prefix pico, symbol: p
:var Femto: prefix femto, symbol: f
:var Atto: prefix atto, symbol: a
:var Zepto: prefix zepto, symbol: z
:var Yocto: prefix yocto, symbol: y

:var km: unit of kilo meter
:var m: unit of meter
:var cm: unit of centi meter
:var mm: unit of milli meter
:var sec: unit of second
:var msec: unit of milli second
:var minute: unit of minute
:var h: unit of hour
:var hour: unit of hour
:var day: unit of day
:var yr: unit of year
:var Myr: unit of mega year
:var Gyr: unit of giga year
:var gram: unit of gram
:var kg: unit of kilo gram
:var lb: unit of pound
:var ton: metric ton
:var A: unit of Ampere
:var Hz: unit of Hertz (frequency)
:var N: unit of Newton (force)
:var Pa: unit of Pascal (pressure, stress)
:var bar: unit of bar (pressure)
:var atm: unit of atmosphere (pressure)
:var J: unit of Joule (energy, work)
:var W: unit of Watt (power)
:var C: unit of Coulomb (electric charge)
:var V: unit of Volt (electric potential)
:var F: unit of Farad (capacitance)
:var Ohm: unit of Ohm (electric resistance)
:var K: unit of Kelvin (temperature)
:var Mol: unit of Mole (amount of substance)
:var Celsius: unit of Celsius (temperature)
:var Fahrenheit: unit of Fahrenheit (temperature)
:var Poise: unit of Poise (dynamic viscosity)
:var R_Earth_equator: Earth's equatorial radius
:var R_Earth_poles: Earth's polar radius
:var R_Earth: Earth's radius
:var v_light: speed of light
:var pi: value of pi accurate to 10 decimal places
:var Gravitational_Constant: gravitational constant
"""


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"
__author__="Lutz Gross, l.gross@uq.edu.au"

## :file unitsSI.py

from math import pi
class Unit(object):
   """
   a general class to define a physical unit and convert from this unit to an appropriate SI unit.

   `Unit` object have a dual purpose: Firstly physical units can be combined through *,/ and ** to form new physical units or to add prefixes such as
   Milli to m to form mm=Milli*m. Moreover, a given floating point number x (or any other arithmetic object) can be converted from the physical unit to 
   the SI system, eg. 10*mm to create the value for 10mm which is the float number 0.01 in the SI system. In addition, a value in the SI unit can be 
   converted back to the given unit, eg. to express 0.01m in physical units of mm use 0.01/mm which will return 10.
   """
   def __init__(self, name, longname,a ,b ):
       """
       initializes the physical unit
       
       :param name: short name of the physical unit or prefix
       :type name: ``str``
       :param longname: long name of the physical unit or prefix
       :type longname: ``str``
       :param a: absolute value in transformation
       :type a: ``float``
       :param b: slope in transformation
       :type b: ``float``
       """
       self.setName(name)
       self.setLongName(longname)
       self.__a=a
       self.__b=b

   def __str__(self):
       """
       Returns the string representation of the unit (its short name).

       :return: the short name of the unit
       :rtype: ``str``
       """
       return self.getName()

   def getName(self):
       """
       Returns the name of the physical unit
 
       :return: name of the physical unit
       :rtype: ``str``
       """
       return self.__name

   def setName(self, name):
       """
       Sets the name of the physical unit
 
       :param name: new name of the physical unit
       :type name: ``str``
       """
       self.__name=name

   def getLongName(self):
       """
       Returns the long name of the physical unit
 
       :return: name of the physical unit
       :rtype: ``str``
       """
       return self.__longname

   def setLongName(self, name):
       """
       Sets the long name of the physical unit
 
       :param name: new long name of the physical unit
       :type name: ``str``
       """
       self.__longname=name

   def __call__(self,x):
       """
       Converts a value x in the physical unit self to SI.

       :param x: value to convert
       :type x: an arithmetic object
       :return: the value converted to SI units
       :rtype: same type as input
       """
       return self.__b*x+self.__a

   def __mul__(self,other):
       """
       Performs self*other operation for two `Unit` objects

       :param other: an other physical unit
       :type other: `Unit`
       :rtype: `Unit` or ``NotImplemented``
       """
       if isinstance(other, Unit):
          a=self(other(0.))
          b=(self(other(1.))-a)
          if isinstance(other, _PowUnit) or isinstance(self,_DivUnit) or  isinstance(self, _PowUnit):
            return _ProdUnit(self.getName()+" "+other.getName(),self.getLongName()+"*"+other.getLongName(),a ,b)
          else:
            return _ProdUnit(self.getName()+other.getName(),self.getLongName()+"*"+other.getLongName(),a ,b)
       else:
          return NotImplemented

   def __rmul__(self,other):
       """
       Performs other*self operation

       :param other: an other `Unit` or an arithmetic object. if other is a arithmetic object such as ``float`` other is assumed to be given in the physical unit ``self`` and is converted into the corresponding SI unit.
       :type other: `Unit` or
       :rtype: `Unit` of or an arithmetic object
       """
       if isinstance(other, Unit):
          a=other(self(0.))
          b=(other(self(1.))-a)
          if isinstance(other, _PowUnit) or  isinstance(self, _PowUnit) or isinstance(other, _DivUnit):
             return _ProdUnit(other.getName()+" "+self.getName(),other.getLongName()+"*"+self.getLongName(),a ,b)
          else:
             return _ProdUnit(other.getName()+self.getName(),other.getLongName()+"*"+self.getLongName(),a ,b)
       else:
          return self(other)

   #Need this for python2       
   def __div__(self, other):
       return self.__truediv__(other)

   def __truediv__(self,other):
       """
       Performs self*other operation for two `Unit` objects

       :param other: an other physical unit
       :type other: `Unit`
       :rtype: `Unit` or ``NotImplemented``
       """
       if isinstance(other, Unit):
          if abs(self(0.))+abs(other(0.))>0:
              raise ValueError("Division of physical units requires 0 absolute values")
          if  isinstance(other, (_ProdUnit, _DivUnit)):
              # X/(c*d) or X/(c/d) requires brackets:
              return _DivUnit(self.getName()+"/("+other.getName()+")",self.getLongName()+"/("+other.getLongName()+")",0 , self(1.)/other(1.))
          else:
              return _DivUnit(self.getName()+"/"+other.getName(),self.getLongName()+"/"+other.getLongName(),0 , self(1.)/other(1.))
       else:
          return NotImplemented

   #Need this for python2       
   def __rdiv__(self, other):
       return self.__rtruediv__(other)          
          
   def __rtruediv__(self,other):
       """
       Performs other/self operation

       :param other: an other `Unit` or an arithmetic object
       :type other: `Unit` or an arithmetic object
       :rtype: `Unit` or an arithmetic object
       """
       if isinstance(other, Unit):
          if abs(self(0.))+abs(other(0.))>0:
              raise ValueError("Division of physical units requires 0 absolute values")
          if  isinstance(self, (_ProdUnit, _DivUnit)):
              # X/(a*b) or X/(a/b) requires brackets:
              return _DivUnit(other.getName()+"/("+self.getName()+")",other.getLongName()+"/("+self.getLongName()+")",0 , other(1.)/self(1.))
          else:
              return _DivUnit(other.getName()+"/"+self.getName(),other.getLongName()+"/"+self.getLongName(),0 , other(1.)/self(1.))
       else:
          return (other-self(0.))/(self(1.)-self(0.))
   def __pow__(self,other):
       """
       Performs self**other operation

       :param other: an exponent
       :type other: ``int`` or ``float``
       :rtype: `Unit` 
       """
       if isinstance(other, float) or isinstance(other, int):
          if abs(self(0.))>0:
              raise ValueError("Power of physical unit requires 0 absolute values")
          if  isinstance(self, (_ProdUnit, _DivUnit, _PowUnit)):
              return _PowUnit("("+self.getName()+")^%s"%other, "("+self.getLongName()+")^%s"%other, 0., self(1.)**other)
          else:
              return _PowUnit(self.getName()+"^%s"%other, self.getLongName()+"^%s"%other, 0., self(1.)**other)
       else:
          return NotImplemented

class _ProdUnit(Unit):
    pass
class _DivUnit(Unit):
    pass
class _PowUnit(Unit):
    pass

#
#  prefixes:
#
Yotta=Unit("Y","Yotta",0.,1.e24)
Zetta=Unit("Z","Zetta",0.,1.e21)
Exa=Unit("E","Exa",0.,1.e18)
Peta=Unit("P","Peta",0.,1.e15)
Tera=Unit("T","Tera",0.,1.e12)
Giga=Unit("G","Giga",0.,1.e9)
Mega=Unit("M","Mega",0.,1.e6)
Kilo=Unit("k","Kilo",0.,1.e3)
Hecto=Unit("h","Hecto",0.,1.e2)
Deca=Unit("da","Deca",0.,1.e1)
one=Unit("1","1",0.,1.)
Deci=Unit("d","Deci",0.,1.e-1)
Centi=Unit("c","Centi",0.,1.e-2)
Milli=Unit("m","Milli",0.,1.e-3)
Micro=Unit("mu","Micro",0.,1.e-6)
Nano=Unit("n","Nano",0.,1.e-9)
Pico=Unit("p","Pico",0.,1.e-12)
Femto=Unit("f","Femto",0.,1.e-15)
Atto=Unit("a","Atto",0.,1.e-18)
Zepto=Unit("z","Zepto",0.,1.e-21)
Yocto=Unit("y","Yocto",0.,1.e-24)
#
#   length
#
m=Unit("m","meter",0.,1.)
km=Kilo*m
cm=Centi*m
mm=Milli*m
liter=(Deci*m)**3
ft=Unit("ft", "feet", 0., 0.3048 * m)
Barrel=Unit("bbl","barrel", 0., 0.158987294928 * m**3)
Mscf=Unit("MSCF", "thousand standard cubic feet", 0., 28.31685 * m**3)
Mcf=1000. * ft**3
ccm=cm**3
#
#  time
#
sec=Unit("sec","second",0.,1.)
msec=Unit("msec","milli-second",0.,1e-3)
minute=Unit("min","minute",0.,60.)
hour=Unit("h","hour",0.,60.* minute)
h=hour
day=Unit("d","day",0., 24. * h)
yr=Unit("yr","year",0.,365.2425 * day)
year=yr
Myr=Mega*yr
Gyr=Giga*yr
#
#  mass
#
kg=Unit("kg","kg",0.,1.)
gram=Milli*kg
lb=Unit("lb","pound",0.,0.45359237 * kg)
ton=Kilo*kg
#
#   electric current
#
A=Unit("A","Ampere",0.,1.)
#
#   Temperature
#
K=Unit("K","Kelvin",0.,1.)
Celsius=Unit("C","Celsius",273.15,1.)
Fahrenheit=Unit("F","Fahrenheit",459.67*5./9.,5./9.)
#
#  others
#
Mol=Unit("mole","Mole",0.,1.)
Hz=one/sec
N = Unit("N","Newton",0.,1.)
Pa = Unit("Pa","Pascal",0.,1.)
bar=100*Kilo*Pa
atm= Unit("atm","atmosphere",0.,101325.024 * Pa)
psi= Unit("psi","Psi",0.,6894.75728034 * Pa)
J = Unit("J","Joule",0.,1.)
W= Unit("W","Watt",0.,1.)
C=Unit("C","Coulomb",0.,1.)
V = Unit("V","Volt",0.,1.)
F = Unit("F","Farad",0.,1.)
Ohm=Unit("Ohm","Ohm",0.,1.)
RAD=Unit("RAD","rad",0.,1.)
DEG=Unit("Ohm","Ohm",0.,pi/180.)
gal=Unit("gal","galileo",0.,0.01)
mgal=gal*Milli
Tesla=V*sec/(m**2)
T=Tesla
#
#  Derived 
#
Poise= gram/cm/sec
cPoise=Poise*Centi
mPoise=Poise*Milli
Darcy= 9.869233e-13*m**2
mDarcy=Darcy*Milli
#
#  some constants
#
R_Earth_equator=6378.1370*km
R_Earth_poles=6356.7523*km
R_Earth=(R_Earth_equator+R_Earth_poles)/2
v_light=299792458.*m/sec
Mu_0 = 4 * pi * 10**(-7)*V*sec/A/m 
Magnetic_Dipole_Moment_Earth = 8.22 * 10**22*A*m**2 
Gravitational_Constant=6.6742e-11*m**3/(kg*sec**2)


