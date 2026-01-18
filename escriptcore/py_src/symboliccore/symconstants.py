
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


"""
Symbolic mathematical constants for escript.

This module provides symbolic versions of common mathematical constants
that can be used in symbolic expressions.

:var pi: symbolic representation of pi (3.14159...)
:var e: symbolic representation of Euler's number e (2.71828...)
"""

from .symbol import Symbol
import sympy
pi=Symbol(sympy.pi)
e=Symbol(sympy.E)
