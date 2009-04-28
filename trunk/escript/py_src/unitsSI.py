
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
