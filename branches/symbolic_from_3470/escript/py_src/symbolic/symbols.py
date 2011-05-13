# -*- coding: utf-8 -*-

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

__author__="Cihan Altinay"


class Symbol(object):
    """
    """

    def __init__(self, *args, **kwargs):
        """
        Initializes a new Symbol object.
        """
        #from esys.escript import Data
        if len(args)==1:
            arg=args[0]
            if isinstance(arg, str):
                self.__arr=numpy.array(sympy.Symbol(arg, **kwargs))
            elif hasattr(arg, "__array__"):
                arr=arg.__array__()
                if len(arr.shape)>4:
                    raise ValueError("Symbol only supports tensors up to order 4")
                self.__arr=arr.copy()
            elif isinstance(arg, list) or isinstance(arg, sympy.Basic):
                self.__arr=numpy.array(arg)
            #elif isinstance(arg, Data):
            #    self.__arr=arg
            else:
                raise TypeError("Unsupported argument type %s"%str(type(arg)))
        elif len(args)==2:
            if not isinstance(args[0], str):
                raise TypeError("First argument must be a string")
            if args[0].find('[')>=0 or args[0].find(']')>=0:
                raise TypeError("Name must not contain '[' or ']'")
            if not isinstance(args[1], tuple):
                raise TypeError("Second argument must be a tuple")
            name=args[0]
            shape=args[1]
            if len(shape)>4:
                raise ValueError("Symbol only supports tensors up to order 4")
            if len(shape)==0:
                self.__arr=numpy.array(sympy.Symbol(name, **kwargs))
            else:
                self.__arr=sympy.symarray(shape, '['+name+']')
        else:
            raise TypeError("Unsupported number of arguments")
        if self.__arr.ndim==0:
            self.name=self.__arr.item()
        else:
            self.name=str(self.__arr.tolist())

    def __repr__(self):
        return str(self.__arr)

    def __str__(self):
        return str(self.__arr)

    def __eq__(self, other):
        if type(self) is not type(other):
            return False
        if self.getRank()!=other.getRank():
            return False
        if self.getShape()!=other.getShape():
            return False
        return (self.__arr==other.__arr).all()

    def __getitem__(self, key):
        return self.__arr[key]

    def __setitem__(self, key, value):
        if isinstance(value,Symbol):
            if value.getRank()==0:
                self.__arr[key]=value
            elif hasattr(self.__arr[key], "shape"):
                if self.__arr[key].shape==value.getShape():
                    self.__arr[key]=value
                else:
                    raise ValueError("Wrong shape of value")
            else:
                raise ValueError("Wrong shape of value")
        elif isinstance(value,sympy.Basic):
            self.__arr[key]=value
        elif hasattr(value, "__array__"):
            self.__arr[key]=map(sympy.sympify,value.flat)
        else:
            self.__arr[key]=sympy.sympify(value)

    def getRank(self):
        return self.__arr.ndim

    def getShape(self):
        return self.__arr.shape

    def atoms(self, *types):
        s=set()
        for el in self.__arr.flat:
            atoms=el.atoms(*types)
            for a in atoms:
                if a.is_Symbol:
                    n,c=Symbol._symComp(a)
                    s.add(sympy.Symbol(n))
                else:
                    s.add(a)
        return s

    def _sympystr_(self, printer):
        return self.lambdarepr()

    def lambdarepr(self):
        from sympy.printing.lambdarepr import lambdarepr
        temp_arr=numpy.empty(self.__arr.shape, dtype=object)
        for idx,el in numpy.ndenumerate(self.__arr):
            atoms=el.atoms(sympy.Symbol)
            # create a dictionary to convert names like x_0_0 to x[0,0]
            symdict={}
            for a in atoms:
                n,c=Symbol._symComp(a)
                if len(c)>0:
                    c=[str(i) for i in c]
                    symstr=n+'['+','.join(c)+']'
                else:
                    symstr=n
                symdict[a.name]=symstr
            s=lambdarepr(el)
            for key in symdict:
                s=s.replace(key, symdict[key])
            temp_arr[idx]=s
        res='combineData(%s,%s)'%(str(temp_arr.tolist()).replace("'",""),str(self.__arr.shape))
        return res

    def diff(self, *symbols, **assumptions):
        symbols=Symbol._symbolgen(*symbols)
        result=Symbol(self.__arr)
        for s in symbols:
            if isinstance(s, Symbol):
                if s.getRank()>0:
                    if s.getShape()!=self.getShape():
                        raise ValueError("Incompatible shapes")
                    a=result.__arr.flat
                    b=s.__arr.flat
                    for idx in range(len(a)):
                        a[idx]=a[idx].diff(b.next())
                else:
                    diff_item=lambda item: getattr(item, 'diff')(s.__arr.item(), **assumptions)
                    result=result.applyfunc(diff_item)

            else:
                diff_item=lambda item: getattr(item, 'diff')(s, **assumptions)
                result=result.applyfunc(diff_item)
        return result

    def swap_axes(self, axis0, axis1):
        return Symbol(numpy.swapaxes(self.__arr, axis0, axis1))

    def tensorProduct(self, other, axis_offset):
        arg0_c=self.__arr.copy()
        sh0=self.__arr.shape
        if isinstance(other, Symbol):
            arg1_c=other.__arr.copy()
            sh1=other.getShape()
        else:
            arg1_c=other.copy()
            sh1=other.shape
        d0,d1,d01=1,1,1
        for i in sh0[:self.__arr.ndim-axis_offset]: d0*=i
        for i in sh1[axis_offset:]: d1*=i
        for i in sh1[:axis_offset]: d01*=i
        arg0_c.resize((d0,d01))
        arg1_c.resize((d01,d1))
        out=numpy.zeros((d0,d1),numpy.object)
        for i0 in range(d0):
            for i1 in range(d1):
                out[i0,i1]=numpy.sum(arg0_c[i0,:]*arg1_c[:,i1])
        out.resize(sh0[:self.__arr.ndim-axis_offset]+sh1[axis_offset:])
        return Symbol(out)

    def transposedTensorProduct(self, other, axis_offset):
        arg0_c=self.__arr.copy()
        sh0=self.__arr.shape
        if isinstance(other, Symbol):
            arg1_c=other.__arr.copy()
            sh1=other.getShape()
        else:
            arg1_c=other.copy()
            sh1=other.shape
        d0,d1,d01=1,1,1
        for i in sh0[axis_offset:]: d0*=i
        for i in sh1[axis_offset:]: d1*=i
        for i in sh1[:axis_offset]: d01*=i
        arg0_c.resize((d01,d0))
        arg1_c.resize((d01,d1))
        out=numpy.zeros((d0,d1),numpy.object)
        for i0 in range(d0):
            for i1 in range(d1):
                out[i0,i1]=numpy.sum(arg0_c[:,i0]*arg1_c[:,i1])
        out.resize(sh0[axis_offset:]+sh1[axis_offset:])
        return Symbol(out)

    def tensorTransposedProduct(self, other, axis_offset):
        arg0_c=self.__arr.copy()
        sh0=self.__arr.shape
        if isinstance(other, Symbol):
            arg1_c=other.__arr.copy()
            sh1=other.getShape()
            r1=other.getRank()
        else:
            arg1_c=other.copy()
            sh1=other.shape
            r1=other.ndim
        d0,d1,d01=1,1,1
        for i in sh0[:self.__arr.ndim-axis_offset]: d0*=i
        for i in sh1[:r1-axis_offset]: d1*=i
        for i in sh1[r1-axis_offset:]: d01*=i
        arg0_c.resize((d0,d01))
        arg1_c.resize((d1,d01))
        out=numpy.zeros((d0,d1),numpy.object)
        for i0 in range(d0):
            for i1 in range(d1):
                out[i0,i1]=numpy.sum(arg0_c[i0,:]*arg1_c[i1,:])
        out.resize(sh0[:self.__arr.ndim-axis_offset]+sh1[:r1-axis_offset])
        return Symbol(out)

    def trace(self, axis_offset):
        sh=self.__arr.shape
        s1=1
        for i in range(axis_offset): s1*=sh[i]
        s2=1
        for i in range(axis_offset+2,len(sh)): s2*=sh[i]
        arr_r=numpy.reshape(self.__arr,(s1,sh[axis_offset],sh[axis_offset],s2))
        out=numpy.zeros([s1,s2],object)
        for i1 in range(s1):
            for i2 in range(s2):
                for j in range(sh[axis_offset]):
                    out[i1,i2]+=arr_r[i1,j,j,i2]
        out.resize(sh[:axis_offset]+sh[axis_offset+2:])
        return Symbol(out)

    def transpose(self, axis_offset):
        if axis_offset is None:
            axis_offset=int(self.__arr.ndim/2)
        axes=range(axis_offset, self.__arr.ndim)+range(0,axis_offset)
        return Symbol(numpy.transpose(self.__arr, axes=axes))

    def applyfunc(self, f):
        assert callable(f)
        if self.__arr.ndim==0:
            out=Symbol(f(self.__arr.item()))
        else:
            out=numpy.empty(self.__arr.shape, dtype=object)
            for idx in numpy.ndindex(self.__arr.shape):
                out[idx]=f(self.__arr[idx])
            out=Symbol(out)
        return out

    def _sympy_(self):
        return self.applyfunc(sympy.sympify)

    @staticmethod
    def _symComp(sym):
        n=sym.name
        a=n.split('[')
        if len(a)!=2:
            return n,()
        a=a[1].split(']')
        if len(a)!=2:
            return n,()
        name=a[0]
        comps=[int(i) for i in a[1].split('_')[1:]]
        return name,tuple(comps)

    @staticmethod
    def _symbolgen(*symbols):
        """
        Generator of all symbols in the argument of diff().
        (cf. sympy.Derivative._symbolgen)

        Example:
        >> ._symbolgen(x, 3, y)
        (x, x, x, y)
        >> ._symbolgen(x, 10**6)
        (x, x, x, x, x, x, x, ...)
        """
        from itertools import repeat
        last_s = symbols[len(symbols)-1]
        if not isinstance(last_s, Symbol):
            last_s=sympy.sympify(last_s)
        for i in xrange(len(symbols)):
            s = symbols[i]
            if not isinstance(s, Symbol):
                s=sympy.sympify(s)
            next_s = None
            if s != last_s:
                next_s = symbols[i+1]
                if not isinstance(next_s, Symbol):
                    next_s=sympy.sympify(next_s)

            if isinstance(s, sympy.Integer):
                continue
            elif isinstance(s, Symbol) or isinstance(s, sympy.Symbol):
                # handle cases like (x, 3)
                if isinstance(next_s, sympy.Integer):
                    # yield (x, x, x)
                    for copy_s in repeat(s,int(next_s)):
                        yield copy_s
                else:
                    yield s
            else:
                yield s

    # unary/binary operations follow

    def __pos__(self):
        return self

    def __neg__(self):
        return Symbol(-self.__arr)

    def __abs__(self):
        return Symbol(abs(self.__arr))

    def __add__(self, other):
        if isinstance(other, Symbol):
            return Symbol(self.__arr+other.__arr)
        return Symbol(self.__arr+other)

    def __radd__(self, other):
        if isinstance(other, Symbol):
            return Symbol(other.__arr+self.__arr)
        return Symbol(other+self.__arr)

    def __sub__(self, other):
        if isinstance(other, Symbol):
            return Symbol(self.__arr-other.__arr)
        return Symbol(self.__arr-other)

    def __rsub__(self, other):
        if isinstance(other, Symbol):
            return Symbol(other.__arr-self.__arr)
        return Symbol(other-self.__arr)

    def __mul__(self, other):
        if isinstance(other, Symbol):
            return Symbol(self.__arr*other.__arr)
        return Symbol(self.__arr*other)

    def __rmul__(self, other):
        if isinstance(other, Symbol):
            return Symbol(other.__arr*self.__arr)
        return Symbol(other*self.__arr)

    def __div__(self, other):
        if isinstance(other, Symbol):
            return Symbol(self.__arr/other.__arr)
        return Symbol(self.__arr/other)

    def __rdiv__(self, other):
        if isinstance(other, Symbol):
            return Symbol(other.__arr/self.__arr)
        return Symbol(other/self.__arr)

    def __pow__(self, other):
        if isinstance(other, Symbol):
            return Symbol(self.__arr**other.__arr)
        return Symbol(self.__arr**other)

    def __rpow__(self, other):
        if isinstance(other, Symbol):
            return Symbol(other.__arr**self.__arr)
        return Symbol(other**self.__arr)


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
        d=Data(0., shape, fs.pop()) #FIXME: interpolate instead of using first?
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

