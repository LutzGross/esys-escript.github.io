
##############################################################################
#
# Copyright (c) 2003-2015 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2015 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from .symbol import *
from .evaluator import *

__nodocorecursion=['symbol', 'evaluator']

from esys.escriptcore.start import HAVE_SYMBOLS
if HAVE_SYMBOLS:
    from . import functions as symfn
    from . import symconstants as Symconsts
    from .pretty import pretty_print, pprint
    from .utils import *

    # prefer escript's implementation of functions such as 'sign' etc.
    from sympy.utilities.lambdify import MODULES
    ESCRIPT_NAMESPACE = {}
    ESCRIPT_DEFAULT = {}
    ESCRIPT_TRANSLATIONS = {
        #"ln":"log",
    }

    if len(MODULES['math'])==3:
        MODULES['escript']=(ESCRIPT_NAMESPACE, ESCRIPT_TRANSLATIONS,('from esys.escript import *',))
    else:
        MODULES['escript']=(ESCRIPT_NAMESPACE, ESCRIPT_DEFAULT, ESCRIPT_TRANSLATIONS,('from esys.escript import *',))

    del ESCRIPT_NAMESPACE
    del ESCRIPT_DEFAULT
    del ESCRIPT_TRANSLATIONS
    del MODULES


#
# vim: expandtab shiftwidth=4:
