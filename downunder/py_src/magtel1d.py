# -*- coding: utf-8 -*-
##############################################################################
#
# Copyright (c) 2015 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2015 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
Some models for flow

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Ralf Schaa, r.schaa@uq.edu.au"

import sys
import numpy
import cmath

class MT_1D(object):
  """
  Calculates the electromagnetic fields in the subsurface for a 1D layered earth.

  Partly based on Fortran code by  Phil Wannamaker in MT2D
  (http://marineemlab.ucsd.edu/Projects/Occam/2DMT/index.html)

  """
  def __init__(self, freq, depths, rho, zcoord):
    """
    DESCRIPTION:
    ------------
    Constructor which initialises the 1D magnetotelluric class:
    (*) check for argument type
    (*) check for valid argument values
    (*) initialises required data lists

    ARGUMENTS:
    ----------
    param freq        :: sounding frequency
    type  freq        :: ``float``
    param depths      :: layer depth interfaces
    type  depths      :: ``list`` (number)
    param rho         :: layer resistivities
    type  rho         :: ``list`` (number)
    param zcoord      :: sample coordinate points
    type  zcoord      :: ``list`` (number)


    DATA ATTRIBUTES:
    ----------------
    self.f  = freq    :: sounding frequency
    self.z  = zcoord  :: sample coordinate points
    self.zl = zl      :: layer depths
    self.dl = dl      :: layer thicknesses
    self.rl = rl      :: layer resistivities
    """

    # ---
    # Check input types:
    # ---

    #make python3 compatible, since long disappeared in python 3
    if sys.version_info[0] == 3:
        long_type = int
    else:
        long_type = long

    if not isinstance(freq, (int,long_type,float) ):
      raise ValueError("Input parameter FREQ must be a number")
    if not isinstance(depths, list) or not all(isinstance(d,(int,long_type,float)) for d in depths):
      raise ValueError("Input parameter DEPTHS must be a list of numbers")
    if not isinstance(rho, list) or not all(isinstance(d,(int,long_type,float)) for d in rho):
      raise ValueError("Input parameter RHO must be a list of numbers")
    if not isinstance(zcoord, list) or not all(isinstance(d,(int,long_type,float)) for d in zcoord):
      raise ValueError("Input parameter ZCOORD must be a list of numbers")


    # ---
    # Check valid input values:
    # ---

    if not freq > 0:
      raise ValueError("Input parameter FREQ must be larger than 0")
    if not all(x>y for x, y in zip(depths, depths[1:])):
      raise ValueError("Input parameter DEPTHS must be all strictly decreasing")
    if not len(depths) > 1:
      raise ValueError("Input parameter DEPTHS must have more than 1 element")
    if not len(rho) == len(depths)-1:
      raise ValueError("Input parameter RHO must be exactly the size of DEPTHS minus 1")
    if not all(x>0 for x in rho):
      raise ValueError("Input parameter RHO must be all positive")
    if not all(x<y for x, y in zip(zcoord, zcoord[1:])):
      raise ValueError("Input parameter ZCOORD must be all strictly increasing")


    # ---
    # Now initialise the required lists for mt1d
    # ---

    # Setup layer thicknesses from interface coordinates.
    dl = []
    for i in range(0,len(depths)-1):
        # Don't include air-layer:
        if rho[i] < 1.0e+10:
            dl.append( abs(depths[i+1] - depths[i]) )

    # Setup list for cumulative layer depths:
    zl = [0] * (len(dl)) ; zl[0] = dl[0]
    if len(dl)-1 >=1:
        for n in range(1,len(dl)):
            zl[n] = zl[n-1] + dl[n]

    # Setup resistivity list without air-layer.
    rl = list(rho)
    if rl[0] > 1.0e+10:
          rl.pop(0)


    # ---
    # initialise all required variables as data attributes
    # ---

    self.f  = freq
    self.z  = zcoord
    self.zl = zl
    self.dl = dl
    self.rl = rl

#__________________________________________________________________________________________________



  def mt1d(self):
    """
    DESCRIPTION:
    -----------
    Public method to calculate the MT-1D EM-fields at sample coordinates.

    USES:
    -----
    self.f  :: sounding frequency
    self.z  :: sample coordinate points
    self.zl :: layer depths
    self.dl :: layer thicknesses
    self.rl :: layer resistivities

    """

    # Compute the transmission & reflection coefficients;
    an, rn = self.__coeff(self.f, self.dl, self.rl)

    # Number of evaluation sample points:
    nz = len(self.z)

    # Initialise output arrays.
    te = numpy.zeros( nz, dtype=complex )
    tm = numpy.zeros( nz, dtype=complex )

    # Calculate the fields at the sample points:
    for i in range( nz ):
        z = self.z[i]
        te[i], tm[i] = self.__field(z, an, rn, self.f, self.zl, self.dl, self.rl)

    #<Note>: return reverse list -> [::-1] so that the first value is the bottom value:
    return te[::-1], tm[::-1]

#__________________________________________________________________________________________________



  def __coeff(self, f, dl, rl):
    """
    DESCRIPTION:
    -----------
    Computes the transmission and reflection coefficients.
    Based on Wannamaker's subroutine 'COAMP' in 'MT2D'

    ARGUMENTS:
    -----------
    f    :: sounding frequency.
    dl   :: layer thicknesses.
    rl   :: layer resistivities.
    """

    # ---
    # Initialise (return) lists for coefficients.
    # ---

    # Number of layers.
    nl = len(rl)

    # Transmission and Reflection coefficients "an" and "rn"
    an = [ complex(0.0,00) ]*(nl)
    rn = [ complex(0.0,00) ]*(nl)



    # ---
    # Constant values.
    # ---

    pi = cmath.pi    # Ratio of circle circumference to it's diameter.
    ra = 1.0e+14     # Resistivity of air.
    mu = 4*pi*1e-7   # Free space permeability.
    w  = 2*pi*f      # Angular frequency.
    wm = w*mu        # Shortcut of product.



    # ---
    # Calculate intrinsic wave numbers <quasi-static>.
    # ---

    # Wave number of air:
    k0 = cmath.sqrt( -1j*wm/ra )

    # Cycle layers and compute wave numbers of other layers:
    k = [None]*nl
    for i in range(nl):
        k[i] = cmath.sqrt( -1j*wm/rl[i] )



    # ---
    # Reflection & transmission coefficients for half-space.
    # ---

    # Half-space case:
    if nl == 1:
        an[0] = 2*k0/(k[0] + k0) # = 1+Ro
        rn[0] = (k0 - k[0])/(k0 + k[0])

        # All done, return the coefficients.
        return an, rn




    # ---
    # Prepare calculations for layers.
    # ---

    # Initialise lists for computed values with complex zeros.
    arg = [ complex(0.0,00) ]*(nl-1)
    exp = [ complex(0.0,00) ]*(nl-1)
    ex2 = [ complex(0.0,00) ]*(nl-1)
    tnh = [ complex(0.0,00) ]*(nl-1)

    # Setup arguments for the exponential for each layer..
    # .. and compute the tanh function and also exp(-2kh).
    for j in range(nl-1):
        arg[j] = 1j*k[j]*dl[j]
        tnh[j]= cmath.tanh(arg[j])
        # Save also exponentials for coefficient calculations later:
        exp[j] = cmath.exp( -arg[j] )
        ex2[j] = exp[j]*exp[j]

    # ---
    # Reflection & transmission coefficients for layers.
    # ---

    #<Note>: Following section is based on the formulae by Wannamaker's code.


    # Initialise recursion with intrinsic impedance of basement.
    zn = wm/k[nl-1]

    # Compute the reflection coefficients for all sub-surface layers..
    # ..start the loop at the basement and cycle up to the first layer:
    for j in range(nl-1,0,-1):
        # Wave impedance of next layer-up:
        zu = wm/k[j-1]
        # Ratio of layer impedances of current-layer and layer-up::
        rn[j] = (zn - zu)/(zn + zu)
        # New apparent impedance for up-layer via Wait's formula:
        zn = zu*(zn + zu*tnh[j-1])/(zu + zn*tnh[j-1])
        # <Note>: "zn" is the surface impedance when finishing the loop.

    # For the first sub-surface layer, we also ..
    # ..have to mind the air-layer at index '0':
    zu = wm/k0 ; rn[0] = (zn - zu)/(zn + zu)


    # Transmission coefficient of first layer takes into account air-layer:
    an[0] = (1+rn[0]) / (1+rn[1]*ex2[0]) # exp[0]*
    #<Note>: Wannamaker does not multiply with exp!

    # And now compute the transmission coefficients for rest of the layers:
    if (nl-1) > 1:
        for n in range(1,nl-1):
            #<Note>: Wannamaker uses num: ~ exp[n-1]!
            num   = (1+rn[n] )*exp[n-1]
            den   = (1+rn[n+1]*ex2[n])
            an[n] = an[n-1]*num/den
    # And mind the basement as well (numerator is 1):
    an[nl-1] = an[nl-2]*exp[nl-2]*(1+rn[nl-1])


    # Return the coefficients.
    return an, rn
#__________________________________________________________________________________________________



  def __field(self, z, an, rn, f, zl, dl, rl):
    """
    DESCRIPTION:
    -----------
    Computes the electric and magnetic field for 1D-MT.
    Based on Wannamaker's subroutine 'ZLFLD' in 'MT2D'

    ARGUMENTS:
    -----------
    z    :: sample coordinate
    an   :: transmission coefficients
    rn   :: reflection coefficients
    f    :: sounding frequency.
    zl   :: layer depths
    dl   :: layer thicknesses
    rl   :: layer resistivities
    """

    # ---------------------------------------------------------------------------------------------
    # Initialisations.
    # ---------------------------------------------------------------------------------------------

    # Number of layers.
    nl = len(rl)

    # Return values.
    ex = complex(0.0, 0.0)
    hy = complex(0.0, 0.0)

    # Constant values.
    pi = cmath.pi    # Ratio of circle circumference to it's diameter.
    ra = 1.0e+14     # Resistivity of air.
    mu = 4*pi*1e-7   # Free space permeability.
    w  = 2*pi*f      # Angular frequency.
    wm = w*mu        # Shortcut of product.



    # ---------------------------------------------------------------------------------------------
    # Calculate intrinsic wave numbers <quasi-static>.
    # ---------------------------------------------------------------------------------------------

    # Free space wave number & amplitude factor of E-field:
    k0 = cmath.sqrt( -1j*wm*1/ra )
    e0 = wm/(2*k0)

    # Cycle layers and compute wave numbers of other layers:
    k = [None]*nl
    for i in range(nl):
        k[i] = cmath.sqrt( -1j*wm/rl[i] )



    # ---------------------------------------------------------------------------------------------
    # Air-layer EM fields.
    # ---------------------------------------------------------------------------------------------

    if z < 0:
        # Compute the argument and fields.
        kz = 1j*k0*z
        ex = e0*cmath.exp(-kz)*(1+rn[0]*cmath.exp(2*kz))
        hy = e0*cmath.exp(-kz)*(1-rn[0]*cmath.exp(2*kz))*(k0/wm)

        # All done, leave.
        return ex, hy




    # ---------------------------------------------------------------------------------------------
    # Uniform half-space EM fields.
    # ---------------------------------------------------------------------------------------------
    if nl == 1:

        # Compute the argument and fields; <Note>:z<0.
        kz = 1j*k[0]*z
        ex = e0*an[0]*cmath.exp(-kz)
        hy = e0*an[0]*cmath.exp(-kz)*(k[0]/wm)

        # All done, leave.
        return ex, hy




    # ---------------------------------------------------------------------------------------------
    # Layered half-space EM fields.
    # ---------------------------------------------------------------------------------------------

    if nl > 1:

        # First get the layer index for the current
        # depth 'z' via cycling over layer depths:
        n = 0
        for i in range(1,nl):
            if z > zl[i-1]:
                n = i

        # This handles the case when 'z' is in the basement layer (no reflection):
        if n == nl-1:

            # Compute the fields.
            kh = 1j*k[n]*(z-zl[n-1])
            ex = e0*an[n]*cmath.exp(-kh)
            hy = e0*an[n]*cmath.exp(-kh)*(k[n]/wm)

            # All done, leave.
            return ex, hy

        else: # Other layers:

            # Compute the fields.
            kh = 1j*k[n]*(z-zl[n])
            kd = 1j*k[n]*dl[n]
            exp = cmath.exp(-kh)
            aex = cmath.exp(-kd)
            ex  = e0*an[n]*(exp + rn[n+1]/exp)*aex
            hy  = e0*an[n]*(exp - rn[n+1]/exp)*aex*(k[n]/wm)

            # All done, leave.
            return ex, hy

    return ex, hy
#__________________________________________________________________________________________________
