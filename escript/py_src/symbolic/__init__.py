
##############################################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from symbol import *
from evaluator import *

from esys.escript import HAVE_SYMBOLS
if HAVE_SYMBOLS:
    from . import functions as symfn
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
