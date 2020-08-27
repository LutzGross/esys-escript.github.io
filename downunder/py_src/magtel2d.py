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
2D Magnetotelluric modelling for TE and TM mode.

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Ralf Schaa, r.schaa@uq.edu.au"

import os, sys
import numpy
import math
import cmath
import types
from . import magtel1d         as mt1d
import esys.weipa              as weipa
import esys.escript            as escript
try:
    import esys.finley         as finley
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False

import esys.escript.pdetools   as pdetools
import esys.escript.linearPDEs as pde

class MT_2D(object):

   # class options:
  _debug   = False    #
  _solver = "DEFAULT" #

   # 'private' field:
  __version = 0.1     #

  """
  DESCRIPTION:
  ------------
  solves the scalar 2-D electromagnetic diffusion equation,
  (where 'u' is the electric field E or magnetic field H).

  [1]  -div( k*grad(u) ) + q*u = 0  (+ Boundary Conditions)

  In 2D the equation is solved for the transverse electric
  field (TE mode) or transverse magnetic field (TM mode).
  These fields are parallel to the 2D strike direction.
  Based on the actual mode, the coefficients are given by:

  TE: k = 1/mu   , q = i*w*sigma
  TM: k = 1/sigma, q = i*w*mu

  'mu'    is the vacuum permeability,
  'i'     is the imaginary unit
  'w'     is the angular frequency
  'sigma' is the conductivity

  The EM diffusion equation is complex and is solved as
  a coupled PDE for the real and imaginary parts. The
  coupled PDE is given by the following equations, with
  Er, Ei and Hr, Hi are the real and imaginary components
  of the electric and magnetic field, respectively:

  TE:
  [2] div( grad(Er) ) + w*mu*sigma*Ei = 0
  [3] div( grad(Ei) ) - w*mu*sigma*Er = 0

     the complementary magnetic fields
     are calculated via Faraday's Law:

  [4] Hr =-d/dz(Ei) / (w*mu)
  [5] Hi = d/dz(Er) / (w*mu)


  TM:
  [6] div( rho*grad(Hr) ) + w*mu*Hi = 0
  [7] div( rho*grad(Hi) ) - w*mu*Hr = 0

     (resistivity 'rho' is 1/sigma)
     the complementary electric fields
     are calculated via Ampere's Law:

  [8] Er = d/dz(Hr) * rho
  [9] Ei = d/dz(Hi) * rho


  Based on the ratio of electric to magnetic field
  apparent resistivity and phase is calculated, viz:

  rho_a = (1/w*mu) * [ (Er)^2 + (Ei)^2 ] / [ (Hr)^2 + (Hi)^2 ]
  phase = arctan( [Ei*Hr - Er*Hi] / [Er*Hr + Ei*Hi] )


  Boundary conditions:
  --------------------
  the source term on the right-hand-side of equation [1] is zero,
  i.e. no artificial source is employed but instead the 'source'
  is provided via the boundary conditions of the PDE which are
  given as Dirichlet conditions at all boundaries. To calculate
  the Dirichlet values, a 1D response is calculated at the left
  and right boundary (based on the 1D recursion formula for MT).
  Interpolation from the left to the right sides then provides
  the values at the top and bottom boundary. See module 'mt1d'
  for details of the computation of the 1D response. Once the
  values on the boundaries have been calculated, the values
  inside the domain are solved in this class.
  """

  def __init__(self, domain, mode, freq_def, tags, rho, rho_1d, ifc_1d,
        xstep=100, zstep=100, maps=None, plot=False, limits=None):
    """
    DESCRIPTION:
    ------------
    Constructor which initialises the 2D magnetotelluric class:
    (*) check for argument type
    (*) check for valid argument values
    (*) initialises required values

    ARGUMENTS:
    ----------
    param  domain       :: the 2d mesh domain
    type   domain       :: ``escript data object``

    param  mode         :: TE or TM mode
    type   mode         :: ``string``

    param  freq_def     :: highest/lowest frequency & points per decade
    type   freq_def     :: ``dictionary``

    param  tags         :: the tag names of the regions defined in the mesh
    type   tags         :: ``list``

    param  rho          :: the resistivity values of the regions in the mesh
    type   rho          :: ``list``

    param  rho_1d       :: the resistivity values at the left & right boundary
    type   rho_1d       :: ``dictionary``

    param  ifc_1d       :: the layer interface depths of the left & right boundary
    type   ifc_1d       :: ``dictionary``

    param  xstep        :: user-defined step size for horizontal sample list
    type   xstep        :: ``number``  (optional)

    param  zstep        :: user-defined step size for vertical sample list
    type   zstep        :: ``number``  (optional)

    param  maps         :: list with user-defined  functions which map the resistivity for each region
    type   maps         :: ``list``    (optional)

    param  plot         :: user-defined flag to show a plot of apparent resistivity and phase at each frequency
    type   plot         :: ``boolean`` (optional)



    DATA ATTRIBUTES:
    ----------------
    self.domain         :: escript data object of mesh
    self.X              :: escript data object with all mesh coordinates
    self.mode           :: string with TE or TM mode
    self.xmin           :: float with x-coordinate minimum
    self.xmax           :: float with x-coordinate maximum
    self.zmin           :: float with z-coordinate minimum
    self.zmax           :: float with z-coordinate maximum
    self.zstep          :: number with sample step in vertical direction
    self.xstep          :: number with sample step in horizontal direction
    self.rho            :: list with resistivity values of all regions
    self.rho_1d         :: dictionary with resistivity values at boundaries left/right
    self.ifc_1d         :: dictionary with interface depths at boundaries left/right
    self.plot           :: boolean flag to show plots of apparent resistivity and phase
    self.sigma          :: escript data object with the conductivity model (based on 'rho' and 'maps')
    self.frequencies    :: list of sounding frequencies
    self.boundary_mask  :: Dirichlet mask at boundaries
    """

    if not HAVE_FINLEY:
        raise ImportError("Finley module not available")
    #make python3 compatible, since long disappeared in python 3
    if sys.version_info[0] == 3:
        long_type = int
    else:
        long_type = long
    # ---
    # Checks
    # ---

    # Types:
    if not isinstance(domain, escript.Domain):
      raise ValueError("Input parameter DOMAIN must be an Escript mesh")
    if not isinstance(mode, str):
      raise ValueError("Input parameter MODE must be a string")
    if not isinstance(freq_def, dict) or len(freq_def) != 3:
      raise ValueError("Input parameter FREQ_DEF must be a dictionary with length 3")
    if not isinstance(tags, list) or not all(isinstance(t,str) for t in tags):
      raise ValueError("Input parameter TAGS must be a list of strings")
    if not isinstance(rho, list) or not all(isinstance(d,(int,long_type,float)) for d in rho):
      raise ValueError("Input parameter RHO must be a list of numbers")
    if not isinstance(rho_1d, dict) or len(rho_1d) != 2:
      raise ValueError("Input parameter RHO_1D must be a dictionary with length 2")
    if not isinstance(ifc_1d, dict) or len(ifc_1d) != 2:
      raise ValueError("Input parameter IFC_1D must be a dictionary with length 2")
    if not isinstance(xstep, (int,long_type,float)):
        raise ValueError("Optional input parameter XSTEP must be a number")
    if not isinstance(zstep, (int,long_type,float)):
        raise ValueError("Optional input parameter ZSTEP must be a number")
    if maps is not None:
      if not isinstance(maps, list) or not all(isinstance(m,(types.FunctionType, types.NoneType)) for m in maps):
        raise ValueError("Optional input parameter MAPS must be a list of Functions or Nones")
    if plot is not None:
      if not isinstance(plot, bool):
        raise ValueError("Optional input parameter PLOT must be True or False")

    # Values:
    if mode.upper() != "TE" and mode.upper() != "TM": # Check mode:
      raise ValueError("Input parameter mode must be either 'TE' or 'TM'")
    if not 'high' in freq_def and not 'low' in freq_def and not 'step' in freq_def:
       raise ValueError("Input dictionary FREQ_DEF must have keys 'high', 'low' and 'step' defined" )
    if freq_def['high'] < freq_def['low']:
      raise ValueError("High frequency value is smaller than low frequency value in input parameter FREQ_DEF")
    if freq_def['step'] < 1:
      raise ValueError("Step frequency value is smaller than 1 in input parameter FREQ_DEF")
    if not all(r>0 for r in rho): # Check resistivity values:
      raise ValueError("Input parameter RHO must be all positive")
    if len(rho) != len(tags): # Check resistivity list-length:
      raise ValueError("Input parameter RHO must have the same length as input parameter TAGS")
    if not 'left' in rho_1d and not 'right' in rho_1d:
       raise ValueError("Input dictionary RHO_1D must have keys 'left' and 'right' defined" )
    if not 'left' in ifc_1d and not 'right' in ifc_1d:
      raise ValueError("Input dictionary IFC_1D must have keys 'left' and 'right' defined" )
    if len(ifc_1d['left'])-1 != len(rho_1d['left']) and len(ifc_1d['right'])-1 != len(rho_1d['right']):
      raise ValueError("Lists with values in input dictionary RHO_1D must have length equal to IFC_1D" )
    if xstep < 0.5: # Step size should be non-zero but should not be smaller than 0.5m:
      raise ValueError("Input parameter XSTEP must be at least 0.5" )
    if zstep < 0.5: # Step size should be non-zero but should not be smaller than 0.5m:
      raise ValueError("Input parameter ZSTEP must be at least 0.5" )



    # ---
    # Domain coordinates & tags:
    # ---

    # Extract the model coordinates..
    X = escript.Solution(domain).getX()

    # Get the Min/Max coordinates:
    xmin = escript.inf(X[0])
    xmax = escript.sup(X[0])
    zmin = escript.inf(X[1])
    zmax = escript.sup(X[1])

    # Get the tag names from the mesh file
    mesh_tags = escript.getTagNames(domain)

    if xmin >= xmax or zmin >= zmax: raise ValueError("The mesh limits are not valid (min >= max)" )
    if zmin >= 0                   : raise ValueError("The mesh must be defined with a negative vertical axis" )
    if not set(mesh_tags) == set(tags)       :
        print("user-tags:", tags)
        print("mesh-tags:", mesh_tags)
        raise ValueError("Input parameter TAGS does not match the tags defined in the mesh")



    # ---
    # Define the boundary mask:
    # ---

    boundary_mask = self.__setBoundaryMask(X)


    # ---
    # Calculate list of sounding frequencies:
    # ---

    frequencies = self.__getSoundingFrequencies(freq_def)



    # ---
    # Tag the domain with conductivity maps:
    # ---

    sigma_model = self.__tagDomain(domain, X, tags, rho, maps)

    # Check for valid values
    if  escript.inf(sigma_model) < 0 or escript.sup(sigma_model) < 0:
       raise ValueError("Negative conductivity encountered" )
    if cmath.isnan( escript.inf(sigma_model) ) or \
       cmath.isnan( escript.sup(sigma_model) ) or \
       cmath.isinf( escript.sup(sigma_model) ):
       raise ValueError("The conductivity model contains NaNs or INFs" )



    # ---
    # Projector and Locator objects.
    # ---

    print("Setting up Escript Locator and Projector objects...")

    # Setup a list with sample points along the vertical mesh extent, bottom to top:
    xsample = self.__getSamplePoints(escript.inf(X[0]),escript.sup(X[0]),xstep, constant=0.0)

    # Get the locations of mesh points at the surface via the Locator object
    # operating on the continuous function-space (i.e. nodes) of the domain.
    loc  = pdetools.Locator(escript.ContinuousFunction(domain),xsample )

    # Instantiate the Projector class with smoothing on (fast=False);
    # the Projector is used to calculate the gradient correctly.
    proj = pdetools.Projector(domain, reduce=False, fast=False)




    # ---
    # Print information:
    # ---

    print("")
    print("="*72)
    print("Escript MT2D, version", self.__version)
    print("="*72)
    print("Mesh XMin/XMax       : ", xmin, xmax)
    print("Mesh ZMin/ZMax       : ", zmin, zmax)
    print("Number of Tags       : ", len( tags ))
    print("Mapping              : ", {True: 'Yes', False: 'No'}[maps is not None])
    print("Conductivity Model   : ", sigma_model)
    print("Nr of Frequencies    : ", len( frequencies ))
    print("Start/End/Step (Hz)  : ", freq_def["high"], freq_def["low"], freq_def["step"])
    print("Mode                 : ", mode.upper())
    print("Solver               : ", MT_2D._solver)
    print("Show plots           : ", plot)
    print("="*72)
    print("")

    if self._debug:
      print("Mesh-Info     : ")
      print(domain.print_mesh_info(full=False))



    # ---
    # Set all required variables as data attributes
    # ---

    self.domain         = domain
    self.X              = X
    self.mode           = mode
    self.xmin           = xmin
    self.xmax           = xmax
    self.zmin           = zmin
    self.zmax           = zmax
    self.zstep          = zstep
    self.xstep          = xstep
    self.rho            = rho
    self.rho_1d         = rho_1d
    self.ifc_1d         = ifc_1d
    self.plot           = plot
    self.limits         = limits
    self.sigma          = sigma_model
    self.frequencies    = frequencies
    self.boundary_mask  = boundary_mask
    self.proj           = proj
    self.loc            = loc


#_______________________________________________________________________________


  def __interpolLinear(self,dx,x0,x1,y0,y1):
    """
    DESCRIPTION:
    -----------
    Function for simple 1D interpolation using the line-equation.

    ARGUMENTS:
    ----------
    dx :: interpolation step.
    x0 :: first coordinate point of known value y0.
    x1 :: last coordinate point of known value y1.
    y0 :: known value at first coordinate.
    y1 :: known value at last coordinate.

    RETURNS:
    --------
    y  :: list with interpolated values
    """
    # Initialise return lists.
    y = []

    # Test for long enough interval.
    if abs(x1-x0) <= dx: return y
    # Test for correct abscissae.
    if x0 >= x1: return y

    x = x0
    while x <= x1:
        y.append( y0 + (y1-y0)*(x-x0)/(x1-x0)  )
        x = x + dx

    return y

#_______________________________________________________________________________


  def __getSamplePoints(self, min,max,step,constant=None):
    """
    DESCRIPTION:
    -----------
    Function to setup a list with sample points. If a
    constant value was passed a 2D list is returned
    where the second column is set to the constant.

    ARGUMENTS:
    ----------
    min        :: minimum value.
    max        :: maximum value.
    step       :: step value.
    constant   :: optional constant value for 2nd column.

    RETURNS:
    --------
    sample     :: list with samples.
    """

    # Initialise return list.
    sample = []

    # Cycle with step-size and fill sample list.
    dp = min
    while dp <= max:
        if constant is not None:
            sample.append([dp,constant])
        else:
            sample.append(dp)
        # Increment the step.
        dp = dp + step

    # Return the list:
    return sample
    #___________________________________________________________________________


  def __getSoundingFrequencies(self, frequencies):
    """
    DESCRIPTION:
    -----------
    Defines the sounding frequencies in Hz.

    ARGUMENTS:
    ----------
    frequencies :: dictionary with frequency start/stop/step

    RETURNS:
    --------
    sounding_frequencies  :: list with frequency values
    """
    # Output list with frequencies in Hertz:
    sounding_frequencies = []

    # Period definition (from freq to time):
    tme_1 = 1.0/frequencies["high"]
    tme_n = 1.0/frequencies["low"]

    # Number of points per decade:
    tme_p = frequencies["step"]

    # Number of periods in range:
    nt = int(math.log10(tme_n/tme_1) * tme_p) + 1
    # Fill list with times:
    for n in range(nt):
      # Sounding period in seconds:
      period = tme_1*10**( (n)/float(tme_p))
      # And store as frequency in Hertz:
      sounding_frequencies.append( 1.0/period )

    return sounding_frequencies

#_______________________________________________________________________________


  def __getGradField(self, proj, mt2d_field, wm):
    """
    DESCRIPTION:
    -----------
    Calculates the complementary fields via Faraday's Law (TE-mode)
    or via Ampere's Law (TM-mode). Partial derivative w.r.t. the
    vertical coordinate are taken at the surface for which an Escript
    'Projector' object is used to calculate the gradient.

    ARGUMENTS:
    ----------
    proj       :: escript Projection object
    mt2d_field :: calculated magnetotelluric field
    wm         :: number with actual angular sounding frequency * mu

    RETURNS:
    --------
    mt2d_grad  :: dictionary with computed gradient fields
    """

    # Define the derived gradient fields:
    if self.mode.upper() == 'TE':
      # H = -(dE/dz) / iwm
      grad_real =-proj.getValue( escript.grad(mt2d_field["imag"])/wm )
      grad_imag = proj.getValue( escript.grad(mt2d_field["real"])/wm )
       #<Note the coupled dependency on real/imaginary part>:
    else:
      # E = (dH/dz) / sigma
      grad_real = proj.getValue( escript.grad(mt2d_field["real"])/self.sigma )
      grad_imag = proj.getValue( escript.grad(mt2d_field["imag"])/self.sigma )
      #<'sigma' is an Escript data-object and as such the division
      # will use the tagged sigma values of the associated domains>


    # And return as dictionary for real and imaginary parts:
    mt2d_grad = {"real": grad_real[1], "imag":grad_imag[1] }
    #<Note>: the derivative w.r.t. 'z' is used (i.e. '[1]').

    return mt2d_grad

#_______________________________________________________________________________


  def __tagDomain(self, domain, X, tags, rho, maps):
    """
    DESCRIPTION:
    -----------
    Defines the conductivity model. Conductivities of tagged regions can be mapped
    via user-defined functions passed in 'maps'. If no user-defined functions are
    passed, a constant value is applied as provided in list 'rho' for each region.
    User-defined functions have 3 arguments: x-coordinate, z-coordinate, resistivity.

    ARGUMENTS:
    ----------
    domain  :: escript object of mesh
    X       :: escript object with all coordinates
    tags    :: list with domain tags
    rho     :: list with domain resistivities
    maps    :: list with user-defined resistivity mappings

    RETURNS:
    --------
    sigma   :: escript object of conductivity model

    """
    # Setup the conductivity structure (acts on elements and can be discontinuous).
    sigma = escript.Scalar(0, escript.Function(domain))

    # Setup conductivity domains.
    for i in range( len(tags) ):

      # Default: assign conductivity which is the inverse of resistivity:
      m = 1.0/rho[i]

      # Map a user-defined conductivity distribution if given:
      if maps is not None:
            # Guard against undefined elements:
        if maps[i] is not None:
          # Map the conductivity according to the defined functions:
          m = maps[i]( X[0], X[1], rho[i] )

      # Tag the mesh with the conductivity distributions at each iteration:
      sigma += m * escript.insertTaggedValues(escript.Scalar(0,escript.Function(domain)),**{ tags[i] : 1})


    if self._debug == True:
      sigma.expand()
      mydir = os.getcwd()
      dbgfl = mydir + os.sep + "mt2d_sigma_dbg.silo"
      print("")
      print("DEBUG: writing SILO debug output of conductivity model:")
      print(dbgfl)
      print("")
      weipa.saveSilo(dbgfl, sigma = sigma)


    # All done:
    return sigma

#_______________________________________________________________________________


  def __setBoundaryMask(self, X):
    """
    DESCRIPTION:
    -----------
    Define Dirichlet model boundaries conditions.

    ARGUMENTS:
    ----------
    X :: escript object with all coordinates

    RETURNS:
    --------
    boundary_mask :: escript object with mask values at boundaries

    """
    # Boundaries are defined as masks (1 or 0) for all mesh coordinates;
    # values at the boundary are '1', whereas all other values are '0'.
    mask_l = escript.whereZero( X[0] - escript.inf(X[0]) )
    mask_r = escript.whereZero( X[0] - escript.sup(X[0]) )
    mask_t = escript.whereZero( X[1] - escript.inf(X[1]) )
    mask_b = escript.whereZero( X[1] - escript.sup(X[1]) )

    # Combine the mask for all boundaries:
    boundary_mask = mask_t + mask_b + mask_l + mask_r

    return boundary_mask
    #<Note>: this boundary mask is used later on as PDE coefficient 'q'.

#_______________________________________________________________________________


  def __getBoundaryValues(self, mode, X, rho_1d, ifc_1d, xstep, zstep, frequency):
    """
    DESCRIPTION:
    -----------
    Returns a list with boundary values along each Dirichlet boundary.
    Values at the left and right side of the domain are evaluated at
    sample points and interpolated across the domain. The subroutine
    expects that layers at the right- and left-hand-side are defined.

    ARGUMENTS:
    ----------
    mode      :: string with TE or TM mode
    X         :: escript object with all coordinates
    rho_1d    :: dictionary with resistivities at the left/right boundary
    ifc_1d    :: dictionary with layer interfaces at the left/right boundary
    xstep     :: number with step size for horizontal sample list
    zstep     :: number with step size for vertical sample list
    frequency :: number with actual sounding frequency

    RETURNS:
    --------
    bondary_value :: dictionary with lists of boundary values at sample points
    """

    # ---
    # Sample lists at vertical and horizontal boundaries.
    # ---

    # Horizontal extents:
    xmin = escript.inf(X[0])
    xmax = escript.sup(X[0])

    # Vertical extents:
    zmin = escript.inf(X[1])
    zmax = escript.sup(X[1])

    # Setup a list with sample points along the vertical mesh extent, bottom to top:
    zsamples = self.__getSamplePoints(-zmax,-zmin,zstep)


    # ---
    # Calculate the 1D response at the left- and right-hand-side boundaries
    # ---

    # Instantiate an 'mt1d' object for the left- and right-hand-sides:
    mt1d_left = mt1d.MT_1D( frequency, ifc_1d['left'] , rho_1d['left'] , zsamples )
    mt1d_rght = mt1d.MT_1D( frequency, ifc_1d['right'], rho_1d['right'], zsamples )

    # Compute the 1D field values at the sample nodes:
    te1d_left, tm1d_left  = mt1d_left.mt1d(  )
    te1d_rght, tm1d_rght  = mt1d_rght.mt1d(  )

    # Distinguish TE and TM mode and save 1D values in dictionary:
    if mode.upper() == "TE":
      mt_1d = {"left":te1d_left, "right":te1d_rght}
    else:
      mt_1d = {"left":tm1d_left, "right":tm1d_rght}


    # ---
    # Interpolation across mesh.
    # ---

    # Now setup a 2D-table from left to right at each sampled depth for mesh-interpolation.
    table2d_real = []
    table2d_imag = []

     # 1D-interpolation of values from left to right at different depths 'i':
    for i in range( len(zsamples)):
      table2d_real.append( self.__interpolLinear(xstep, xmin, xmax, mt_1d["left"].real[i], mt_1d["right"].real[i]) )
      table2d_imag.append( self.__interpolLinear(xstep, xmin, xmax, mt_1d["left"].imag[i], mt_1d["right"].imag[i]) )

    # 2D-interpolation to map the values on the mesh coordinates:
    bondary_value_real = escript.interpolateTable( table2d_real, X, (xmin,zmin), (xstep,zstep) )
    bondary_value_imag = escript.interpolateTable( table2d_imag, X, (xmin,zmin), (xstep,zstep) )

    # Return the real and imaginary values as a dictionary:
    boundary_value = {"real":bondary_value_real, "imag":bondary_value_imag}


    return boundary_value

#_______________________________________________________________________________


  def __getAppResPhase(self, mt2d_field, mt2d_grad, wm):
    """
    DESCRIPTION:
    -----------
    Calculates the apparent resistivity and phase.

    ARGUMENTS:
    ----------
    mt2d_field :: dictionary with real/imaginary field values
    mt2d_grad  :: dictionary with real/imaginary gradient values

    RETURNS:
    --------
    apparent resistivity and phase
    """

    # Define the associated modelled fields in readable variables:
    if self.mode.upper() == 'TE':
      # Transverse electric field:
      Er = mt2d_field["real"]
      Ei = mt2d_field["imag"]
      Hr = mt2d_grad["real"]
      Hi = mt2d_grad["imag"]
    else:
      # Transverse magnetic field :
      Hr = mt2d_field["real"]
      Hi = mt2d_field["imag"]
      Er = mt2d_grad["real"]
      Ei = mt2d_grad["imag"]


    # Return apparent Resistivity and Phase:
    arho_2d = ( (Er**2 + Ei**2)/(Hr**2 + Hi**2) ) / wm
    aphi_2d = escript.atan( (Ei*Hr - Er*Hi)/(Er*Hr + Ei*Hi) ) * 180.0/cmath.pi

    return arho_2d, aphi_2d
#_______________________________________________________________________________


  def __showPlot(self, loc, rho_2d, phi_2d, f, **kwargs):
    """
    DESCRIPTION:
    -----------
    Plot of apparent resistivity and phase. Requires matplotlib to be available.

    ARGUMENTS:
    ----------
    loc     :: escript Locator object
    rho_2d  :: list with computed apparent resistivities
    phi_2d  :: list with computed phase values
    f       :: sounding frequency

    RETURNS:
    --------
    Plot in window.

    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Warning: matplotlib not available, plot will not be shown")
        return

    # Abscissas/Ordinates:
    x  = numpy.array( loc.getX() )[:,0]
    y0 = numpy.array( loc.getValue(rho_2d) )
    y1 = numpy.array( loc.getValue(phi_2d) )

    # Plot labels:
    title = 'Escript MT-2D ' + '(' + self.mode.upper() + ')' + ' freq: ' + str(f) + ' Hz'
    ylbl0 = r'Apparent Resistivity $(\Omega\cdot\,m)$'
    ylbl1 = r'Phase $(^{\circ})$'
    xlbl1 = 'Easting (m)'


    # Setup the plot window with app. res. on top and phase on bottom:
    f, ax = plt.subplots(2, figsize=(8,8), sharex=True) # Mind shared axis
    f.subplots_adjust(top=0.9)  # Little extra space for 'suptitle'
    f.suptitle(title)           # This is actually the plot-title

    # Top: apparent resistivity on semi-log plot
    ax[0].plot(x,y0, color='red') # semilogy
    ax[0].grid(b=True, which='both', color='grey',linestyle=':')
    ax[0].set_title( ylbl0 )
    # Plot limits in **kwargs:
    if 'limits' in kwargs:
        ax[0].set_xlim(kwargs["limits"])

    # Bottom: phase on linear plot
    ax[1].plot(x,y1, color='blue')
    ax[1].grid(b=True, which='both', color='grey',linestyle=':')
    ax[1].set_xlabel( xlbl1 )
    ax[1].set_title( ylbl1 )
    # Plot limits in **kwargs:
    if 'limits' in kwargs:
        ax[1].set_xlim(kwargs["limits"])

    plt.show()


#_______________________________________________________________________________


  def __setSolver(self, mode, domain, sigma, boundary_mask, boundary_value, f):
    """
    DESCRIPTION:
    -----------
    Setups the coupled PDE for real and complex part.

    ARGUMENTS:
    ----------
    mode           :: string with TE or TM mode
    domain         :: escript object with mesh domain
    sigma          :: escript object with conductivity model
    boundary_mask  :: escript object with boundary mask
    boundary_value :: dictionary with real/imag boundary values
    f              :: sounding frequency

    RETURNS:
    --------
    mt2d_fields    :: dictionary with solved PDE, magnetotelluric fields real/imag
    """

    # Constants:
    pi  = cmath.pi     # Ratio circle circumference to diameter.
    mu0 = 4*pi*1e-7    # Free space permeability in V.s/(A.m).
    wm  = 2*pi*f*mu0   # Angular frequency times mu0.


    # ---
    # Setup the coupled PDE for real/imaginary parts:
    # ---

    # Initialise the PDE object for two coupled equations (real/imaginary).
    mtpde = pde.LinearPDE(domain, numEquations=2)

    # If set, solve the 2D case using the direct solver:
    if MT_2D._solver.upper() == "DIRECT":
       mtpde.getSolverOptions().setSolverMethod(pde.SolverOptions().DIRECT)
    else:
       mtpde.getSolverOptions().setSolverMethod(pde.SolverOptions().DEFAULT)

    # Now initialise the PDE coefficients 'A' and 'D',
    # as well as the Dirichlet variables 'q' and 'r':
    A = mtpde.createCoefficient("A")
    D = mtpde.createCoefficient("D")
    q = mtpde.createCoefficient("q")
    r = mtpde.createCoefficient("r")

    # Set the appropriate values for the coefficients depending on the mode:
    if mode.upper() == "TE":
        a_val = 1.0
        d_val = wm*sigma
    elif mode.upper() == "TM":
        a_val = 1.0/sigma
        d_val = wm


    # ---
    # Define the PDE parameters, mind boundary conditions.
    # ---


    # Now define the rank-4 coefficient A:
    for i in range(domain.getDim()):
        A[0,i,0,i] = a_val
        A[1,i,1,i] = a_val

    # And define the elements of 'D' which are decomposed into real/imaginary values:
    D[0,0] = 0     ; D[1,0] = d_val
    D[0,1] =-d_val ; D[1,1] = 0


    # Set Dirichlet boundaries and values:
    q[0] = boundary_mask ; r[0] = boundary_value['real']
    q[1] = boundary_mask ; r[1] = boundary_value['imag']

    # ---
    # Solve the PDE
    # ---

    mtpde.setValue(A=A, D=D, q=q, r=r  )
    pde_solution = mtpde.getSolution()

    # And return the real and imaginary parts individually:
    mt2d_fields = {"real":pde_solution[0], "imag":pde_solution[1] }
    #<Note>: the electric field is returned for TE-mode.
    #        the magnetic field is returned for TM-mode.

    return mt2d_fields

#_______________________________________________________________________________


  def pdeSolve(self):
    """
    DESCRIPTION:
    -----------
    Solves the PDE for either the TE or the TM mode.
    (TE mode is the transverse Electric field).
    (TM mode is the transverse Magnetic field).

    ARGUMENTS:
    ----------
    (uses `self`)

    RETURNS:
    --------
    mt2d  :: list with real/imag fields for each sounding frequency
    arho  :: list with apparent resistivities for each sounding frequency
    aphi  :: list with phase values for each sounding frequency

    """


    # ---
    # Constants.
    # ---

    # Pi & vacuum permeability:
    pi = cmath.pi
    mu = 4*pi*1e-7

    # Number of frequencies:
    nfreq = len(self.frequencies)


    # ---
    # Solve the PDE for all frequencies.
    # ---

    # Prepare lists to store the values at each frequency:
    arho = []
    aphi = []
    mt2d = []

    # Cycle over all frequencies:
    print("Solving for frequency: ...")
    for n in range( nfreq ):

      f = self.frequencies[n] # actual frequency (Hz)
      wm  = (2*pi*f)*mu       # angular frequency (rad/s)
      T = 1.0 / f             # sounding period (s)

      print(n+1,":", f, "(Hz)")

      # Calculate 1D Dirichlet boundary values:
      boundary_value = self.__getBoundaryValues(self.mode.upper(), self.X,
            self.rho_1d, self.ifc_1d, self.xstep, self.zstep, f)

      # Solve the 2D-MT PDE:
      fld_2d = self.__setSolver(self.mode.upper(),self.domain, self.sigma,
            self.boundary_mask, boundary_value, f)

      # Calculate the field gradients:
      grd_2d = self.__getGradField(self.proj, fld_2d, wm)

      # Calculate the apparent resistivity and phase:
      rho_2d, phi_2d = self.__getAppResPhase(fld_2d, grd_2d, wm)

      # Save in lists for each frequency:
      mt2d.append( fld_2d )
      arho.append( self.loc.getValue(rho_2d) )
      aphi.append( self.loc.getValue(phi_2d) )

      # Optionally plot the apparent resistivity and phase:
      if self.plot:
          self.__showPlot(self.loc, rho_2d, phi_2d, f, limits=self.limits)


    # ---
    # All done
    # ---

    print("field calculations finished.")
    return mt2d, arho, aphi

#_______________________________________________________________________________







