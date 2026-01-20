
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
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
