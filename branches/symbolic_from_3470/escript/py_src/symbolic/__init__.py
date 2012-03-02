
########################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from symbols import *
from evaluator import *
from pretty import pretty_print, pprint
import functions as symfn

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

#
# vim: expandtab shiftwidth=4:
