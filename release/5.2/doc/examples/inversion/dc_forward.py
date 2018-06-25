
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

# -------------------------------------------------------------------------------------------------
# DESCRIPTION:
# ------------
# Script to calculate the DC Resistivity response for a half-sphere
# embedded in a half-space (see example from Rucker et al (2006), Figure 4).
#
# REFERENCES:
# -----------
# "Three-dimensional modelling and inversion of DC resistivity data incorporating
# topography -- I. Modelling", Geophysical Journal International (2006) 166, 495-505
# Carsten Rucker, Thomas Gunther and Klaus Spitzer
# -------------------------------------------------------------------------------------------------

from __future__ import division, print_function

import esys.finley      as finley
import esys.escript     as escript
from esys.downunder     import *
import math


# -------------------------------------------------------------------------------------------------
# Constants
# -------------------------------------------------------------------------------------------------

pi = math.pi

# -------------------------------------------------------------------------------------------------
# Input
# -------------------------------------------------------------------------------------------------

mesh_file = "data/HalfSphere_v1.4.msh"

# Tag volume names and conductivity values (Sm/m) for primary and secondary potential:
tag_p = {"domain" : 1/10.0, "sphere" :  1/10.0} # Primary (homogeneous).
tag_s = {"domain" : 1/10.0, "sphere" :  1/1.0 } # Secondary.

xe_0 = -5.0 # start X-coordinate
numEle =  21  # number of electrodes
a =  0.5 # step size
n=9 # 
midPoint = [xe_0 + (((numEle-1)*a)/2.), 0, 0]
current = 1.0 # (Ampere)
domain = finley.ReadGmsh(mesh_file, 3)
mesh_tags = escript.getTagNames(domain)
directionVector = [1,0]

# -------------------------------------------------------------------------------------------------
# Define the primary and secondary resistivity model.
# -------------------------------------------------------------------------------------------------

#<Note>: the mesh was created with domains, which were all 'tagged' with unique id-s;
# based on the tagged domain id-s, the conductivity values are assigned to the domains.
#<Note>: the primary conductivity is the conductivity in the vicinity of the electrodes;
# the secondary conductivity is defined as the difference between primary and the
# actual input conductivity.

# Setup primary and secondary conductivity objects:
sig_p = escript.Scalar(0,escript.ContinuousFunction(domain))
sig_s = escript.Scalar(0,escript.ContinuousFunction(domain))

# Cycle tag_values dictionary and assign conductivity values;
# this is done for both, primary and secondary, potentials as
# both dictionaries are setup with identical dictionary-keys:
for tag in tag_p:
    # All initially defined tags must be in the tag list.
    # Print an error if it doesn't and exit the program.
    if tag in mesh_tags:
        # Assign value:
        sig_p.setTaggedValue( tag, tag_p[tag] )
        sig_s.setTaggedValue( tag, tag_s[tag] )
    else:
        print("Error: the defined tag is not defined in the mesh: " & tag)
        sys.exit()

# Expand the data objects for output.
sig_p.expand()
sig_s.expand()


schs=SchlumbergerSurvey(domain, sig_p, sig_s, current, a, n, midPoint, directionVector, numEle)
pot=schs.getPotentialAnalytic()
primaryApparentRes=schs.getApparentResistivity(pot[0])
SecondaryApparentRes=schs.getApparentResistivity(pot[1])
totalApparentRes=schs.getApparentResistivity(pot[2])


n=1
print ("Total:\n")
for i in totalApparentRes:
    print ("n = %d:"%n)
    print (i,"\n")
    n=n+1
