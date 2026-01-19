
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
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
Symbolic mathematics support for escript.

This package provides symbolic computation capabilities for escript using
sympy as the backend. It allows defining PDEs symbolically and automatic
computation of Jacobians for nonlinear PDE solvers.
"""


__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

from .symbol import *
from . import functions as symfn
from . import symconstants as Symconsts
from .pretty import pretty_print, pprint
from .utils import *
from .evaluator import *

#from .evaluator import *

#__nodocorecursion=['symbol', 'evaluator']
#
#from esys.escriptcore.start import HAVE_SYMBOLS
#if HAVE_SYMBOLS:
    #from . import functions as symfn
    #from . import symconstants as Symconsts
    #from .pretty import pretty_print, pprint
    #from .utils import *

    #prefer escript's implementation of functions such as 'sign' etc.
    #from sympy.utilities.lambdify import MODULES
    #ESCRIPT_NAMESPACE = {}
    #ESCRIPT_DEFAULT = {}
    #ESCRIPT_TRANSLATIONS = {
        #"ln":"log",
    #}

    #if len(MODULES['math'])==3:
        #MODULES['escript']=(ESCRIPT_NAMESPACE, ESCRIPT_TRANSLATIONS,('from esys.escript import *',))
    #else:
        #MODULES['escript']=(ESCRIPT_NAMESPACE, ESCRIPT_DEFAULT, ESCRIPT_TRANSLATIONS,('from esys.escript import *',))

    #del ESCRIPT_NAMESPACE
    #del ESCRIPT_DEFAULT
    #del ESCRIPT_TRANSLATIONS
    #del MODULES


#
# vim: expandtab shiftwidth=4:
