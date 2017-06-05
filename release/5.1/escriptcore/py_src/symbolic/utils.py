
##############################################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
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

def combineData(array, shape):
    """
    """

    # array could just be a single value
    if not hasattr(array,'__len__') and shape==():
        return array

    from esys.escript import Data
    n=numpy.array(array) # for indexing

    # find function space if any
    dom=set()
    fs=set()
    for idx in numpy.ndindex(shape):
        if isinstance(n[idx], Data):
            fs.add(n[idx].getFunctionSpace())
            dom.add(n[idx].getDomain())

    if len(dom)>1:
        domain=dom.pop()
        while len(dom)>0:
            if domain!=dom.pop():
                raise ValueError("Mixing of domains not supported")

    if len(fs)>0:
        d=Data(0., shape, fs.pop()) # maybe interpolate instead of using first?
    else:
        d=numpy.zeros(shape)
    for idx in numpy.ndindex(shape):
        #z=numpy.zeros(shape)
        #z[idx]=1.
        #d+=n[idx]*z # much slower!
        if hasattr(n[idx], "ndim") and n[idx].ndim==0:
            d[idx]=float(n[idx])
        else:
            d[idx]=n[idx]
    return d

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

def getTotalDifferential(f, x, order=0):
    """
    This function computes::

        | Df/Dx = del_f/del_x + del_f/del_grad(x)*del_grad(x)/del_x + ...
        |            \   /         \   /
        |              a             b
    
    """

    from esys.escript import util
    res=()
    shape=util.getShape(f)
    if not isSymbol(f):
        res+=(numpy.zeros(shape+x.getShape()),)
        for i in range(order):
            x=x.grad()
            res+=numpy.zeros(shape+x.getShape())

    elif x.getRank()==0:
        f=removeFsFromGrad(f)
        dfdx=f.diff(x)
        dgdx=x.grad().diff(x)
        a=numpy.empty(shape, dtype=object)
        if order>0:
            b=numpy.empty(shape+dgdx.getShape(), dtype=object)

        if len(shape)==0:
            for j in numpy.ndindex(dgdx.getShape()):
                y=dfdx
                z=dgdx[j]
                # expand() and coeff() are very expensive so
                # we set the unwanted factors to zero to extract
                # the one we need
                for jj in numpy.ndindex(dgdx.getShape()):
                    if j==jj: continue
                    y=y.subs(dgdx[jj], 0)
                a=y.subs(z,0) # terms in x and constants
                if order>0:
                    b[j]=y.subs(z,1)-a
        else:
            for i in numpy.ndindex(shape):
                for j in numpy.ndindex(dgdx.getShape()):
                    y=dfdx[i]
                    z=dgdx[j]
                    for jj in numpy.ndindex(dgdx.getShape()):
                        if j==jj: continue
                        y=y.subs(dgdx[jj], 0)
                    a[i]=y.subs(z,0) # terms in x and constants
                    if order>0:
                        b[i+j]=y.subs(z,1)-a[i]
        res+=(Symbol(a, dim=f.getDim(), subs=f.getDataSubstitutions()),)
        if order>0:
            res+=(Symbol(b, dim=f.getDim(), subs=f.getDataSubstitutions()),)

    elif x.getRank()==1:
        f=removeFsFromGrad(f)
        dfdx=f.diff(x)
        dgdx=x.grad().diff(x).transpose(2)
        a=numpy.empty(shape+x.getShape(), dtype=object)
        if order>0:
            b=numpy.empty(shape+x.grad().getShape(), dtype=object)

        if len(shape)==0:
            raise NotImplementedError('f scalar, x vector')
        else:
            for i in numpy.ndindex(shape):
                for k,l in numpy.ndindex(x.grad().getShape()):
                    if dgdx[k,k,l]==0:
                        a[i+(k,)]=0
                        if order>0:
                            b[i+(k,l)]=0
                    else:
                        y=dfdx[i+(k,)]
                        z=dgdx[k,k,l]
                        for kk,ll in numpy.ndindex(x.grad().getShape()):
                            if k==kk and l==ll: continue
                            y=y.subs(dgdx[kk,kk,ll], 0)
                        a[i+(k,)]=y.subs(z,0) # terms in x and constants
                        if order>0:
                            b[i+(k,l)]=y.subs(z,1)-a[i+(k,)]

        res+=(Symbol(a, dim=f.getDim(), subs=f.getDataSubstitutions()),)
        if order>0:
            res+=(Symbol(b, dim=f.getDim(), subs=f.getDataSubstitutions()),)

    if len(res)==1:
        return res[0]
    else:
        return res

