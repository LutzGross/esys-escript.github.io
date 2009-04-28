
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
__author__="Lutz Gross, l.gross@uq.edu.au"

## @file unitsSI.py

"""
some tools supporting the usage of symbols.

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version

@var Yotta : prefix yotta, symbol:   Y
@var Zetta : prefix zetta, symbol: Z
@var Exa : prefix exa, symbol: E
@var Peta : prefix peta, symbol: P
@var Tera : prefix tera, symbol: T
@var Giga : prefix giga, symbol: G
@var Mega : prefix mega, symbol: M
@var Kilo : prefix kilo, symbol: k
@var Hecto : prefix hecto, symbol: h
@var Deca :  prefix deca, symbol: da
@var Deci :  prefix deci, symbol: d
@var Centi : prefix centi, symbol: c
@var Milli : prefix milli, symbol: m
@var Micro : prefix micro, symbol: mu
@var Nano :  prefix nano, symbol: n
@var Pico :  prefix pico, symbol: p
@var Femto : prefix femto, symbol: f
@var Atto : prefix atto, symbol: a
@var Zepto : prefix zepto, symbol: z
@var Yocto : prefix yocto, symbol: y

@var km : unit of kilo meter
@var m : unit of meter
@var cm : unit of centi meter
@var mm : unit of milli meter
@var sec: unit of second
@var minute : unit of minute
@var h : unit of hour
@var day : unit of day
@var yr : unit of year
@var Myr : unit of mega year
@var Gyr : unit of giga year
@var gram : unit of gram
@var kg : unit of kilo gram
@var lb : unit of pound
@var ton : metric ton
@var A : unit of Ampere
@var Hz: unit of Hertz (frequency)
@var N: unit of Newton (force)
@var Pa: unit of Pascal (pressure, stress)
@var atm: unit of atmosphere (pressure)
@var J: unit of Joule (energy, work)
@var W: unit of Watt (power)
@var C: unit of Coulomb (electric charge)
@var V: unit of Volt (electric potential)
@var F: unit of Farad (Capacitance)
@var Ohm: unit of Ohm (electric resistance)
@var R_Earth_equator: Earth's equatorial radius
@var R_Earth_poles: Earth's polar radius
@var R_Earth: Earth's radius
@var v_light: speed of light
"""
class Unit(object):
   """
   a general class to define a physical unit and a linear converter a+b*x to get transfert to a referene unit system (typically SI)
   """
   def __init__(self, name, longname,a ,b ):
       """
       initializes the unit
       
       @param name: short name of the unit or prefix
       @type name: C{str}
       @param longname: long name of the unit or prefix
       @type longname: C{str}
       @param a: absolute value in transformation
       @type a: C{float}
       @param b: slop in translation
       @type b: C{float}
       """
       self.setName(name)
       self.setLongName(longname)
       self.__a=a
       self.__b=b

   def __str__(self):
       return self.getName()

   def getName(self):
       """
       Returns the name of the unit
 
       @return: name of the unit
       @rtype: C{str}
       """
       return self.__name

   def setName(self, name):
       """
       Sets the name of the unit
 
       @param name: new name of the unit
       @type name: C{str}
       """
       self.__name=name

   def getLongName(self):
       """
       Returns the long name of the unit
 
       @return: name of the unit
       @rtype: C{str}
       """
       return self.__longname

   def setLongName(self, name):
       """
       Sets the long name of the unit
 
       @param name: new long name of the unit
       @type name: C{str}
       """
       self.__longname=name

   def __call__(self,x):
       """
       Converts a value x in the unit self to SI 

       @param x: value to convert
       @type x: an arithmetic object
       """
       return self.__b*x+self.__a

   def __mul__(self,other):
       """
       Performs self*other operation for two L{Unit} objects

       @param other: an other unit
       @type other: L{Unit}
       @rtype: L{Unit} or C{NotImplemented}
       """
       if isinstance(other, Unit):
          a=self(other(0.))
          b=(self(other(1.))-a)
          if isinstance(other, PowUnit) or isinstance(self,DivUnit) or  isinstance(self, PowUnit):
            return ProdUnit(self.getName()+" "+other.getName(),self.getLongName()+"*"+other.getLongName(),a ,b)
          else:
            return ProdUnit(self.getName()+other.getName(),self.getLongName()+"*"+other.getLongName(),a ,b)
       else:
          return NotImplemented

   def __rmul__(self,other):
       """
       Performs other*self operation

       @param other: an other L{Unit} or an arithmetic object
       @type other: L{Unit} 
       @rtype: L{Unit} or an arithmetic object
       """
       if isinstance(other, Unit):
          a=other(self(0.))
          b=(other(self(1.))-a)
          if isinstance(other, PowUnit) or  isinstance(self, PowUnit) or isinstance(other, DivUnit):
             return ProdUnit(other.getName()+" "+self.getName(),other.getLongName()+"*"+self.getLongName(),a ,b)
          else:
             return ProdUnit(other.getName()+self.getName(),other.getLongName()+"*"+self.getLongName(),a ,b)
       else:
          return self(other)

   def __div__(self,other):
       """
       Performs self*other operation for two L{Unit} objects

       @param other: an other unit
       @type other: L{Unit}
       @rtype: L{Unit} or C{NotImplemented}
       """
       if isinstance(other, Unit):
          if abs(self(0.))+abs(other(0.))>0:
              raise ValueError,"Division of units requires 0 absolute values"
          if  isinstance(other, (ProdUnit, DivUnit)):
              # X/(c*d) or X/(c/d) requires brackets:
              return DivUnit(self.getName()+"/("+other.getName()+")",self.getLongName()+"/("+other.getLongName()+")",0 , self(1.)/other(1.))
          else:
              return DivUnit(self.getName()+"/"+other.getName(),self.getLongName()+"/"+other.getLongName(),0 , self(1.)/other(1.))
       else:
          return NotImplemented
   def __rdiv__(self,other):
       """
       Performs other/self operation

       @param other: an other L{Unit} or an arithmetic object
       @type other: L{Unit} or an arithmetic object
       @rtype: L{Unit} or an arithmetic object
       """
       if isinstance(other, Unit):
          if abs(self(0.))+abs(other(0.))>0:
              raise ValueError,"Division of units requires 0 absolute values"
          if  isinstance(self, (ProdUnit, DivUnit)):
              # X/(a*b) or X/(a/b) requires brackets:
              return DivUnit(other.getName()+"/("+self.getName()+")",other.getLongName()+"/("+self.getLongName()+")",0 , other(1.)/self(1.))
          else:
              return DivUnit(other.getName()+"/"+self.getName(),other.getLongName()+"/"+self.getLongName(),0 , other(1.)/self(1.))
       else:
          return (other-self(0.))/(self(1.)-self(0.))
   def __pow__(self,other):
       """
       Performs self**other operation

       @param other: an exponent
       @type other: C{int} or C{float}
       @rtype: L{Unit} 
       """
       if isinstance(other, float) or isinstance(other, int):
          if abs(self(0.))>0:
              raise ValueError,"Power of unit requires 0 absolute values"
          if  isinstance(self, (ProdUnit, DivUnit, PowUnit)):
              return PowUnit("("+self.getName()+")^%s"%other, "("+self.getLongName()+")^%s"%other, 0., self(1.)**other)
          else:
              return PowUnit(self.getName()+"^%s"%other, self.getLongName()+"^%s"%other, 0., self(1.)**other)
       else:
          return NotImplemented

class ProdUnit(Unit):
    pass
class DivUnit(Unit):
    pass
class PowUnit(Unit):
    pass
#
#  prefixes:
#
Yotta=1.e24
Zetta=1.e21
Exa=1.e18
Peta=1.e15
Tera=1.e12
Giga=1.e9
Mega=1.e6
Kilo=1.e3
Hecto=1.e2
Deca=1.e1
Deci=1.e-1
Centi=1.e-2
Milli=1.e-3
Micro=1.e-6
Nano=1.e-9
Pico=1.e-12
Femto=1.e-15
Atto=1.e-18
Zepto=1.e-21
Yocto=1.e-24
#
#   length
#
m=1.
km=Kilo*m
cm=Centi*m
mm=Milli*m
#
#  time
#
sec=1.
minute=60.*sec
h=60.*minute
day=h*24.
yr=day*365.2425 
Myr=Mega*yr
Gyr=Giga*yr
#
#  mass
#
kg=1.
gram=Milli*kg
lb=453.59237*gram
ton=Kilo*kg
#
#   electric current
#
A=1.
#
#  others
#
Hz=1./sec
N = m*kg/sec**2

Pa = N/m**2
atm=101325.024*Pa

J = N*m 
W= J/sec
C=sec*A
V = W/A 
F = C/V
Ohm=V/A
#
#  some constants
#
R_Earth_equator=6378.1370*km
R_Earth_poles=6356.7523*km
R_Earth=(R_Earth_equator+R_Earth_poles)/2
v_light=299792458.*m/sec
