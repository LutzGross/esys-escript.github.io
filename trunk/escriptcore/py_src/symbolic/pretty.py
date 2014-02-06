
##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
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

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import numpy
import sympy
from sympy.printing.pretty.pretty import PrettyPrinter,prettyForm,pretty_symbol
from .symbol import Symbol

__author__="Cihan Altinay"

class ValueMatrix(object):
    def __init__(self, content):
        self._items=numpy.array(content)
        if self._items.ndim>2:
            raise TypeError("ValueMatrix only supports 1-D and 2-D arrays")
        elif self._items.ndim==1:
            self._items=self._items.reshape((1,)+self._items.shape)
        self.rows,self.cols=self._items.shape

    def __getitem__(self, key):
        return self._items[key]

class EscriptPrettyPrinter(PrettyPrinter):
    """
    """
    def __init__(self, profile=None):
        PrettyPrinter.__init__(self, profile)

    def _print_Symbol(self, e):
        # handle escript symbols
        if isinstance(e, Symbol):
            if e.getRank()<=4:
                return self._print(e.__array__())
            return PrettyPrinter._print_Symbol(self,e)

        # e is a sympy Symbol. Remove any brackets from the name in case e is
        # a component
        n,c=Symbol._symComp(e)
        if len(c)==0:
            return PrettyPrinter._print_Symbol(self,e)
        s=sympy.Symbol(n+'_'.join([str(i) for i in c]))
        return PrettyPrinter._print_Symbol(self, s)

    def _print_ndarray(self, e):
        if e.ndim==0:
            return self._print(e.item())
        elif e.ndim==1:
            m=sympy.Matrix(1,e.shape[0],lambda i,j:e[j])
            return self._print_Matrix(m)
        elif e.ndim==2:
            i,j=e.shape
            m=sympy.Matrix(i,j,lambda i,j:e[i,j])
            return self._print_Matrix(m)
        else: #ndim==3 or 4:
            arr=numpy.empty(e.shape[2:],dtype=object)
            for idx in numpy.ndindex(e.shape[2:]):
                arr[idx]=Symbol(e[idx])
            m=ValueMatrix(arr)
            return self._print_Matrix(m)

    def _print_grad_n(self, e):
        s=prettyForm(*self._print(e.args[0]).parens())
        i=pretty_symbol(",_"+str(e.args[1]))
        return prettyForm(*s.right(i))

def pretty(expr, profile=None, **kargs):
    """
    Returns a string containing the prettified form of expr.

    Supported arguments:
        ``expr``
            the expression to print
        ``wrap_line``
            line wrapping enabled/disabled, should be a boolean value
            (default to True)
        ``use_unicode``
            use unicode characters, such as the Greek letter pi instead of
            the string pi. Values should be boolean or None
        ``full_prec``
            use full precision. Default to "auto"
    """
    from sympy.printing.pretty.pretty import pretty_use_unicode
    if profile is not None:
        profile.update(kargs)
    else:
        profile = kargs
    uflag = pretty_use_unicode(kargs.get("use_unicode", None))
    try:
        pp = EscriptPrettyPrinter(profile)
        return pp.doprint(expr)
    finally:
        pretty_use_unicode(uflag)

def pretty_print(expr, use_unicode=None):
    """
    Prints expr in pretty form.

    pprint is just a shortcut for this function
    """
    print(pretty(expr, use_unicode = use_unicode))

pprint = pretty_print

