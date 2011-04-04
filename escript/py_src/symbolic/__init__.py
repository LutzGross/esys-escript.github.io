
########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from sympy.core import Basic,Symbol,diff,evalf,expand,oo,pi,symbols
from evaluator import *

# prefer escript's implementation of functions such as 'sign' etc.
from sympy.utilities.lambdify import MODULES
ESCRIPT_NAMESPACE = {}
ESCRIPT_TRANSLATIONS = {
    #"ln":"log",
}
MODULES['escript']=(ESCRIPT_NAMESPACE, ESCRIPT_TRANSLATIONS,('from esys.escript import *',))

#
# vim: expandtab shiftwidth=4:
