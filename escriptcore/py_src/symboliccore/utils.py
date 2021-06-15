
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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"
__author__="Cihan Altinay"

"""
:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

import numpy
import sympy
from .symbol import Symbol

def symbols(*names, **kwargs):
    """
    Emulates the behaviour of sympy.symbols.
    """

    shape=kwargs.pop('shape', ())

    s = names[0]
    if not isinstance(s, list):
        import re
        s = re.split('\s|,', s)
    res = []
    for t in s:
        # skip empty strings
        if not t:
            continue
        sym = Symbol(t, shape, **kwargs)
        res.append(sym)
    res = tuple(res)
    if len(res) == 0:   # var('')
        res = None
    elif len(res) == 1: # var('x')
        res = res[0]
                        # otherwise var('a b ...')
    return res


def isSymbol(arg):
    """
    Returns True if the argument ``arg`` is an escript ``Symbol`` or
    ``sympy.Basic`` object, False otherwise.
    """
    return isinstance(arg, Symbol) or isinstance(arg, sympy.Basic)

def removeFsFromGrad(sym):
    """
    Returns ``sym`` with all occurrences grad_n(a,b,c) replaced by grad_n(a,b).
    That is, all functionspace parameters are removed.
    """
    from esys.escript import symfn
    gg=sym.atoms(symfn.grad_n)
    for g in gg:
        if len(g.args)==3:
            r=symfn.grad_n(*g.args[:2])
            sym=sym.subs(g, r)
    return sym
