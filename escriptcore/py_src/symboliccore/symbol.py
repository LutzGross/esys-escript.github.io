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

This module provides the `Symbol` class for representing symbolic mathematical
expressions that can be manipulated algebraically and evaluated with escript
Data objects.

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""


__copyright__ = """Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__ = "https://github.com/LutzGross/esys-escript.github.io"
__author__ = "Cihan Altinay"

from esys.escriptcore.start import HAVE_SYMBOLS
import numpy
from esys.escriptcore.escriptcpp import Data, FunctionSpace
if HAVE_SYMBOLS:
    import sympy
    
__all__ = ['Symbol']
   

class Symbol(object):
    """
    `Symbol` objects are placeholders for a single mathematical symbol, such as
    'x', or for arbitrarily complex mathematical expressions such as
    'c*x**4+alpha*exp(x)-2*sin(beta*x)', where 'alpha', 'beta', 'c', and 'x'
    are also Symbols (the symbolic 'atoms' of the expression).

    With the help of the 'Evaluator' class these symbols and expressions can
    be resolved by substituting numeric values and/or escript `Data` objects
    for the atoms. To facilitate the use of `Data` objects a `Symbol` has a
    shape (and thus a rank) as well as a dimension (see constructor).
    Symbols are useful to perform mathematical simplifications, compute
    derivatives and as coefficients for nonlinear PDEs which can be solved by
    the `NonlinearPDE` class.
    """

    # these are for compatibility with sympy.Symbol. lambdify checks these.
    is_Add = False
    is_Float = False

    def __init__(self, *args, **kwargs):
        """
        Initialises a new `Symbol` object in one of three ways::

            u=Symbol('u')

        returns a scalar symbol by the name 'u'.

            alpha=Symbol('alpha', (4,3))

        returns a rank 2 symbol with the shape (4,3), whose elements are
        named '[alpha]_i_j' (with i=0..3, j=0..2).

            a,b,c=symbols('a,b,c')
            x=Symbol([[a+b,0,0],[0,b-c,0],[0,0,c-a]])

        returns a rank 2 symbol with the shape (3,3) whose elements are
        explicitly specified by numeric values and other symbols/expressions
        within a list or numpy array.

        The dimensionality of the symbol can be specified through the ``dim``
        keyword. All other keywords are passed to the underlying symbolic
        library (currently sympy).

        :param args: initialisation arguments as described above
        :keyword dim: dimensionality of the new Symbol (default: 2)
        :type dim: ``int``
        """
        if not HAVE_SYMBOLS:
            raise RuntimeError("Trying to instantiate a Symbol but sympy not available")
        if 'dim' in kwargs:
            self._dim = kwargs.pop('dim')
        else:
            self._dim = -1 # undefined

        if 'subs' in kwargs:
            self._subs = kwargs.pop('subs')
        else:
            self._subs = {}

        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, str):
                if arg.find('[') >= 0 or arg.find(']') >= 0:
                    raise ValueError("Name must not contain '[' or ']'")
                self._arr = numpy.array(sympy.Symbol(arg, **kwargs))
            elif hasattr(arg, "__array__") or isinstance(arg, list):
                if isinstance(arg, list): arg = numpy.array(arg)
                arr = arg.__array__()
                if len(arr.shape) > 4:
                    raise ValueError("Symbol only supports tensors up to order 4")
                res = numpy.empty(arr.shape, dtype=object)
                for idx in numpy.ndindex(arr.shape):
                    if hasattr(arr[idx], "item"):
                        res[idx] = arr[idx].item()
                    else:
                        res[idx] = arr[idx]
                self._arr = res
                if isinstance(arg, Symbol):
                    self._subs.update(arg._subs)
                    if self._dim == -1:
                        self._dim = arg._dim
            elif isinstance(arg, sympy.Basic):
                self._arr = numpy.array(arg)
            else:
                raise TypeError("Unsupported argument type %s"%str(type(arg)))
        elif len(args) == 2:
            if not isinstance(args[0], str):
                raise TypeError("First argument must be a string")
            if not isinstance(args[1], tuple):
                raise TypeError("Second argument must be a tuple")
            name=args[0]
            shape=args[1]
            if name.find('[')>=0 or name.find(']')>=0:
                raise ValueError("Name must not contain '[' or ']'")
            if len(shape)>4:
                raise ValueError("Symbol only supports tensors up to order 4")
            if len(shape)==0:
                self._arr=numpy.array(sympy.Symbol(name, **kwargs))
            else:
                try:
                    self._arr=sympy.symarray(shape, '['+name+']')
                except TypeError:
                    self._arr=sympy.symarray('['+name+']', shape)
        else:
            raise TypeError("Unsupported number of arguments")
        if self._arr.ndim==0:
            self.name=str(self._arr.item())
        else:
            self.name=str(self._arr.tolist())

    def __repr__(self):
        return str(self._arr)

    def __str__(self):
        return str(self._arr)

    def __eq__(self, other):
        if type(self) is not type(other):
            return False
        if self.getRank()!=other.getRank():
            return False
        if self.getShape()!=other.getShape():
            return False
        return (self._arr==other._arr).all()
        
    def __hash__(self):
        return id(self)

    def __getitem__(self, key):
        """
        Returns an element of this symbol which must have rank >0.
        Unlike item() this method converts sympy objects and numpy arrays into
        escript Symbols in order to facilitate expressions that require
        element access, such as: grad(u)[1]+x

        :param key: (nd-)index of item to be returned
        :return: the requested element
        :rtype: ``Symbol``, ``int``, or ``float``
        """
        res=self._arr[key]
        # replace sympy Symbols/expressions by escript Symbols
        if isinstance(res, sympy.Basic) or isinstance(res, numpy.ndarray):
            res=Symbol(res)
            res._dim=self._dim
        res._subs.update(self._subs)
        return res

    def __setitem__(self, key, value):
        if isinstance(value, Symbol):
            self._subs.update(value._subs)
            if value.getRank()==0:
                self._arr[key]=value.item()
            elif hasattr(self._arr[key], "shape"):
                if self._arr[key].shape==value.getShape():
                    for idx in numpy.ndindex(self._arr[key].shape):
                        self._arr[key][idx]=value[idx].item()
                else:
                    raise ValueError("Wrong shape of value")
            else:
                raise ValueError("Wrong shape of value")
        elif isinstance(value, sympy.Basic):
            self._arr[key]=value
        elif hasattr(value, "__array__"):
            self._arr[key]=map(sympy.sympify,value.flat)
        else:
            self._arr[key]=sympy.sympify(value)

    def __iter__(self):
        return self._arr.__iter__

    def __array__(self, t=None):
        if t:
            return self._arr.astype(t)
        else:
            return self._arr

    def _sympy_(self):
        """
        Converts this Symbol to a sympy expression for compatibility.

        :return: sympy-compatible representation of this symbol
        :rtype: ``sympy.Basic`` or ``numpy.ndarray``
        """
        return self.applyfunc(sympy.sympify)

    def getDim(self):
        """
        Returns the spatial dimensionality of this symbol.

        :return: the symbol's spatial dimensionality, or -1 if undefined
        :rtype: ``int``
        """
        return self._dim

    def getRank(self):
        """
        Returns the rank of this symbol.

        :return: the symbol's rank which is equal to the length of the shape.
        :rtype: ``int``
        """
        return self._arr.ndim

    def getShape(self):
        """
        Returns the shape of this symbol.

        :return: the symbol's shape
        :rtype: ``tuple`` of ``int``
        """
        return self._arr.shape

    def getDataSubstitutions(self):
        """
        Returns a dictionary of symbol names and the escript ``Data`` objects
        they represent within this Symbol.

        :return: the dictionary of substituted ``Data`` objects
        :rtype: ``dict``
        """
        return self._subs

    def item(self, *args):
        """
        Returns an element of this symbol.
        This method behaves like the item() method of numpy.ndarray.
        If this is a scalar Symbol, no arguments are allowed and the only
        element in this Symbol is returned.
        Otherwise, 'args' specifies a flat or nd-index and the element at
        that index is returned.

        :param args: index of item to be returned
        :return: the requested element
        :rtype: ``sympy.Symbol``, ``int``, or ``float``
        """
        return self._arr.item(args)

    def atoms(self, *types):
        """
        Returns the atoms that form the current Symbol.

        By default, only objects that are truly atomic and cannot be divided
        into smaller pieces are returned: symbols, numbers, and number
        symbols like I and pi. It is possible to request atoms of any type,
        however.

        Note that if this symbol contains components such as [x]_i_j then
        only their main symbol 'x' is returned.

        :param types: types to restrict result to
        :return: list of atoms of specified type
        :rtype: ``set``
        """
        for t in types:
            if t == type(self):
                types=types+(type(sympy.Symbol("t")),)
        s=set()
        for el in self._arr.flat:
            if isinstance(el,sympy.Basic):
                atoms=el.atoms(*types)
                for a in atoms:
                    if a.is_Symbol:
                        n,c=Symbol._symComp(a)
                        s.add(sympy.Symbol(n))
                    else:
                        s.add(a)
            elif len(types)==0 or type(el) in types:
                s.add(el)
        return s

    def _sympystr_(self, printer):
        # compatibility with sympy 1.6
        return self._sympystr(printer)

    def _sympystr(self, printer):
        return self.lambdarepr()

    def lambdarepr(self):
        """
        Returns a string representation suitable for use with lambdify.

        This method converts the symbol to a string format that can be
        parsed by sympy's lambdify function for numerical evaluation.

        :return: string representation of this symbol for lambdify
        :rtype: ``str``
        """
        from esys.escript.symbolic.printer import EsysEscriptPrinter
        printer = EsysEscriptPrinter()
        temp_arr=numpy.empty(self.getShape(), dtype=object)
        for idx,el in numpy.ndenumerate(self._arr):
            atoms=el.atoms(sympy.Symbol) if isinstance(el,sympy.Basic) else []
            # create a dictionary to convert names like [x]_0_0 to x[0,0]
            symdict={}
            for a in atoms:
                n,c=Symbol._symComp(a)
                if len(c)>0:
                    c=[str(i) for i in c]
                    symstr=n+'['+','.join(c)+']'
                else:
                    symstr=n
                symdict[a.name]=symstr
            s=printer.doprint(el)
            for key in symdict:
                s=s.replace(key, symdict[key])
            temp_arr[idx]=s
        if self.getRank()==0:
            return temp_arr.item()
        else:
            return 'combineData(%s,%s)'%(str(temp_arr.tolist()).replace("'",""),str(self.getShape()))

    def coeff(self, x, expand=True):
        """
        Returns the coefficient of the term "x" or 0 if there is no "x".

        If "x" is a scalar symbol then "x" is searched in all components of
        this symbol. Otherwise the shapes must match and the coefficients are
        checked component by component.

        Example::
        
            x=Symbol('x', (2,2))
            y=3*x
            print y.coeff(x)
            print y.coeff(x[1,1])

        will print::

            [[3 3]
             [3 3]]

            [[0 0]
             [0 3]]

        :param x: the term whose coefficients are to be found
        :type x: ``Symbol``, ``numpy.ndarray``, `list`
        :return: the coefficient(s) of the term
        :rtype: ``Symbol``
        """
        self._ensureShapeCompatible(x)
        if hasattr(x, '__array__'):
            y=x.__array__()
        else:
            y=numpy.array(x)

        if y.ndim>0:
            result=numpy.zeros(self.getShape(), dtype=object)
            for idx in numpy.ndindex(y.shape):
                if y[idx]!=0:
                    res=self._arr[idx].coeff(y[idx], expand)
                    if res is not None:
                        result[idx]=res
        elif y.item()==0:
            result=numpy.zeros(self.getShape(), dtype=object)
        else:
            coeff_item=lambda item: getattr(item, 'coeff')(y.item(), expand)
            none_to_zero=lambda item: 0 if item is None else item
            result=self.applyfunc(coeff_item)
            result=result.applyfunc(none_to_zero)
        res=Symbol(result, dim=self._dim)
        for i in self._subs: res.subs(i, self._subs[i])
        return res

    def subs(self, old, new):
        """
        Substitutes an expression.
        """
        dataSubs={}
        old._ensureShapeCompatible(new)
        if isinstance(new, Data):
            if old.getShape()==() and new.getShape()!=():
                raise ValueError("Only a scalar Data object can be substituted into a scalar\
                        symbol")
            subs=self._subs.copy()
            name='data'+str(id(new))
            newsym=Symbol(name, new.getShape(), dim=self._dim)
            subs.update({Symbol(name):new})
            result=numpy.empty(self.getShape(), dtype=object)
            for idx in numpy.ndindex(self.getShape()):
                result[idx]=self._arr[idx].subs(old._arr[idx], newsym._arr[idx])
            result=Symbol(result, dim=self._dim, subs=subs)
        elif isinstance(old, Symbol) and old.getRank()>0:
            if isinstance(new, Symbol):
                dataSubs=new.getDataSubstitutions()    
            if hasattr(new, '__array__'):
                new=new.__array__()
            else:
                new=numpy.array(new)

            result=numpy.empty(self.getShape(), dtype=object)
            if new.ndim>0:
                for idx in numpy.ndindex(self.getShape()):
                    result[idx]=self._arr[idx].subs(old._arr[idx], new[idx])
            else: # substitute scalar for non-scalar
                for idx in numpy.ndindex(self.getShape()):
                    result[idx]=self._arr[idx].subs(old._arr[idx], new.item())
            result=Symbol(result, dim=self._dim, subs=self._subs)
        else: # scalar
            if isinstance(new, Symbol):
                dataSubs=new.getDataSubstitutions()    
                new=new.item()
            if isinstance(old, Symbol):
                old=old.item()
            subs_item=lambda item: getattr(item, 'subs')(old, new)
            result=self.applyfunc(subs_item)
        result._subs.update(dataSubs)
        return result

    def diff(self, *symbols, **assumptions):
        """
        Computes the derivative of this symbol with respect to given symbols.

        This method differentiates the expression component-wise with respect
        to the specified symbols. Multiple symbols can be provided for higher
        order derivatives.

        :param symbols: symbols to differentiate with respect to
        :type symbols: `Symbol` or ``sympy.Symbol``
        :param assumptions: additional keyword arguments passed to sympy's diff
        :return: the differentiated symbol
        :rtype: `Symbol`
        """
        symbols=Symbol._symbolgen(*symbols)
        result=Symbol(self._arr, dim=self._dim, subs=self._subs)
        for s in symbols:
            if isinstance(s, Symbol):
                if s.getRank()==0:
                    diff_item=lambda item: getattr(item, 'diff')(s._arr.item(), **assumptions)
                    result=result.applyfunc(diff_item)
                elif s.getRank()==1:
                    dim=s.getShape()[0]
                    out=result._arr.copy().reshape(self.getShape()+(1,)).repeat(dim,axis=self.getRank())
                    for d in range(dim):
                        for idx in numpy.ndindex(self.getShape()):
                            index=idx+(d,)
                            out[index]=out[index].diff(s[d].item(), **assumptions)
                    result=Symbol(out, dim=self._dim, subs=self._subs)
                else:
                    raise ValueError("diff: argument must have rank 0 or 1")
            else:
                diff_item=lambda item: getattr(item, 'diff')(s, **assumptions)
                result=result.applyfunc(diff_item)
        return result

    def grad(self, where=None):
        """
        Returns a symbol which represents the gradient of this symbol.

        :param where: optional function space for the gradient
        :type where: `Symbol` or ``FunctionSpace``
        :return: the gradient of this symbol
        :rtype: `Symbol`
        :raises ValueError: if the symbol has undefined dimensionality
        """
        if self._dim < 0:
            raise ValueError("grad: cannot compute gradient as symbol has undefined dimensionality")
        subs=self._subs
        if isinstance(where, Symbol):
            if where.getRank()>0:
                raise ValueError("grad: 'where' must be a scalar symbol")
            where=where._arr.item()
        elif isinstance(where, FunctionSpace):
            name='fs'+str(id(where))
            fssym=Symbol(name)
            subs=self._subs.copy()
            subs.update({fssym:where})
            where=name

        from .functions import grad_n
        out=self._arr.copy().reshape(self.getShape()+(1,)).repeat(self._dim,axis=self.getRank())
        for d in range(self._dim):
            for idx in numpy.ndindex(self.getShape()):
                index=idx+(d,)
                if where is None:
                    out[index]=grad_n(out[index],d)
                else:
                    out[index]=grad_n(out[index],d,where)
        return Symbol(out, dim=self._dim, subs=subs)

    def inverse(self):
        """
        Returns the symbolic inverse of this matrix symbol.

        Only supports rank 2 symbols (matrices) with square shapes up to 3x3.

        :return: the inverse of this matrix symbol
        :rtype: `Symbol`
        :raises TypeError: if the symbol is not rank 2 or not square
        :raises ZeroDivisionError: if the matrix is symbolically singular
        """
        if not self.getRank()==2:
            raise TypeError("inverse: Only rank 2 supported")
        s=self.getShape()
        if not s[0] == s[1]:
            raise ValueError("inverse: Only square shapes supported")
        out=numpy.zeros(s, object)
        arr=self._arr
        if s[0]==1:
            if arr[0,0].is_zero:
                raise ZeroDivisionError("inverse: Symbol not invertible")
            out[0,0]=1./arr[0,0]
        elif s[0]==2:
            A11=arr[0,0]
            A12=arr[0,1]
            A21=arr[1,0]
            A22=arr[1,1]
            D = A11*A22-A12*A21
            if D.is_zero:
                raise ZeroDivisionError("inverse: Symbol not invertible")
            D=1./D
            out[0,0]= A22*D
            out[1,0]=-A21*D
            out[0,1]=-A12*D
            out[1,1]= A11*D
        elif s[0]==3:
            A11=arr[0,0]
            A21=arr[1,0]
            A31=arr[2,0]
            A12=arr[0,1]
            A22=arr[1,1]
            A32=arr[2,1]
            A13=arr[0,2]
            A23=arr[1,2]
            A33=arr[2,2]
            D = A11*(A22*A33-A23*A32)+ A12*(A31*A23-A21*A33)+A13*(A21*A32-A31*A22)
            if D.is_zero:
                raise ZeroDivisionError("inverse: Symbol not invertible")
            D=1./D
            out[0,0]=(A22*A33-A23*A32)*D
            out[1,0]=(A31*A23-A21*A33)*D
            out[2,0]=(A21*A32-A31*A22)*D
            out[0,1]=(A13*A32-A12*A33)*D
            out[1,1]=(A11*A33-A31*A13)*D
            out[2,1]=(A12*A31-A11*A32)*D
            out[0,2]=(A12*A23-A13*A22)*D
            out[1,2]=(A13*A21-A11*A23)*D
            out[2,2]=(A11*A22-A12*A21)*D
        else:
           raise TypeError("inverse: Only matrix dimensions 1,2,3 are supported")
        return Symbol(out, dim=self._dim, subs=self._subs)

    def swap_axes(self, axis0, axis1):
        """
        Returns a symbol with two axes interchanged.

        :param axis0: first axis
        :type axis0: ``int``
        :param axis1: second axis
        :type axis1: ``int``
        :return: symbol with swapped axes
        :rtype: `Symbol`
        """
        return Symbol(numpy.swapaxes(self._arr, axis0, axis1), dim=self._dim, subs=self._subs)

    def tensorProduct(self, other, axis_offset):
        """
        Computes the general tensor product of this symbol with another.

        :param other: the other operand
        :type other: `Symbol` or ``numpy.ndarray``
        :param axis_offset: number of axes to contract
        :type axis_offset: ``int``
        :return: the tensor product
        :rtype: `Symbol`
        """
        arg0_c=self._arr.copy()
        sh0=self.getShape()
        if isinstance(other, Symbol):
            arg1_c=other._arr.copy()
            sh1=other.getShape()
            dim=other._dim if self._dim < 0 else self._dim
        else:
            arg1_c=other.copy()
            sh1=other.shape
            dim=self._dim
        d0,d1,d01=1,1,1
        for i in sh0[:self._arr.ndim-axis_offset]: d0*=i
        for i in sh1[axis_offset:]: d1*=i
        for i in sh1[:axis_offset]: d01*=i
        arg0_c.resize((d0,d01))
        arg1_c.resize((d01,d1))
        out=numpy.zeros((d0,d1), dtype=object)
        for i0 in range(d0):
            for i1 in range(d1):
                out[i0,i1]=numpy.sum(arg0_c[i0,:]*arg1_c[:,i1])
        out.resize(sh0[:self._arr.ndim-axis_offset]+sh1[axis_offset:])
        subs=self._subs.copy()
        subs.update(other._subs)
        return Symbol(out, dim=dim, subs=subs)

    def transposedTensorProduct(self, other, axis_offset):
        """
        Computes the transposed tensor product of this symbol with another.

        :param other: the other operand
        :type other: `Symbol` or ``numpy.ndarray``
        :param axis_offset: number of axes to contract
        :type axis_offset: ``int``
        :return: the transposed tensor product
        :rtype: `Symbol`
        """
        arg0_c=self._arr.copy()
        sh0=self.getShape()
        if isinstance(other, Symbol):
            arg1_c=other._arr.copy()
            sh1=other.getShape()
            dim=other._dim if self._dim < 0 else self._dim
        else:
            arg1_c=other.copy()
            sh1=other.shape
            dim=self._dim
        d0,d1,d01=1,1,1
        for i in sh0[axis_offset:]: d0*=i
        for i in sh1[axis_offset:]: d1*=i
        for i in sh1[:axis_offset]: d01*=i
        arg0_c.resize((d01,d0))
        arg1_c.resize((d01,d1))
        out=numpy.zeros((d0,d1),object)
        for i0 in range(d0):
            for i1 in range(d1):
                out[i0,i1]=numpy.sum(arg0_c[:,i0]*arg1_c[:,i1])
        out.resize(sh0[axis_offset:]+sh1[axis_offset:])
        subs=self._subs.copy()
        subs.update(other._subs)
        return Symbol(out, dim=dim, subs=subs)

    def tensorTransposedProduct(self, other, axis_offset):
        """
        Computes the tensor transposed product of this symbol with another.

        :param other: the other operand
        :type other: `Symbol` or ``numpy.ndarray``
        :param axis_offset: number of axes to contract
        :type axis_offset: ``int``
        :return: the tensor transposed product
        :rtype: `Symbol`
        """
        arg0_c=self._arr.copy()
        sh0=self.getShape()
        if isinstance(other, Symbol):
            arg1_c=other._arr.copy()
            sh1=other.getShape()
            r1=other.getRank()
            dim=other._dim if self._dim < 0 else self._dim
        else:
            arg1_c=other.copy()
            sh1=other.shape
            r1=other.ndim
            dim=self._dim
        d0,d1,d01=1,1,1
        for i in sh0[:self._arr.ndim-axis_offset]: d0*=i
        for i in sh1[:r1-axis_offset]: d1*=i
        for i in sh1[r1-axis_offset:]: d01*=i
        arg0_c.resize((d0,d01))
        arg1_c.resize((d1,d01))
        out=numpy.zeros((d0,d1),object)
        for i0 in range(d0):
            for i1 in range(d1):
                out[i0,i1]=numpy.sum(arg0_c[i0,:]*arg1_c[i1,:])
        out.resize(sh0[:self._arr.ndim-axis_offset]+sh1[:r1-axis_offset])
        subs=self._subs.copy()
        subs.update(other._subs)
        return Symbol(out, dim=dim, subs=subs)

    def trace(self, axis_offset):
        """
        Returns the trace of this Symbol.
        """
        sh=self.getShape()
        s1=1
        for i in range(axis_offset): s1*=sh[i]
        s2=1
        for i in range(axis_offset+2,len(sh)): s2*=sh[i]
        arr_r=numpy.reshape(self._arr,(s1,sh[axis_offset],sh[axis_offset],s2))
        out=numpy.zeros([s1,s2],object)
        for i1 in range(s1):
            for i2 in range(s2):
                for j in range(sh[axis_offset]):
                    out[i1,i2]+=arr_r[i1,j,j,i2]
        out.resize(sh[:axis_offset]+sh[axis_offset+2:])
        return Symbol(out, dim=self._dim, subs=self._subs)

    def transpose(self, axis_offset):
        """
        Returns the transpose of this Symbol.
        """
        if axis_offset is None:
            axis_offset=int(self._arr.ndim/2)
        axes=list(range(axis_offset, self._arr.ndim))+list(range(0,axis_offset))
        return Symbol(numpy.transpose(self._arr, axes=axes), dim=self._dim, subs=self._subs)

    def applyfunc(self, f, on_type=None):
        """
        Applies the function `f` to all elements (if on_type is None) or to
        all elements of type `on_type`.
        """
        assert callable(f)
        if self._arr.ndim==0:
            if on_type is None or isinstance(self._arr.item(), on_type):
                el=f(self._arr.item())
            else:
                el=self._arr.item()
            if el is not None:
                out=Symbol(el, dim=self._dim, subs=self._subs)
            else:
                return el
        else:
            out=numpy.empty(self.getShape(), dtype=object)
            for idx in numpy.ndindex(self.getShape()):
                if on_type is None or isinstance(self._arr[idx],on_type):
                    out[idx]=f(self._arr[idx])
                else:
                    out[idx]=self._arr[idx]
            out=Symbol(out, dim=self._dim, subs=self._subs)
        return out

    def expand(self):
        """
        Applies the sympy.expand operation on all elements in this symbol
        """
        return self.applyfunc(sympy.expand, sympy.Basic)

    def simplify(self):
        """
        Applies the sympy.simplify operation on all elements in this symbol
        """
        return self.applyfunc(sympy.simplify, sympy.Basic)
    
    def evalf(self):
        """
        Applies the sympy.evalf operation on all elements in this symbol
        """
        evalf_s=lambda item: getattr(item, 'evalf')()
        return self.applyfunc(evalf_s, sympy.Basic)
    # unary/binary operators follow

    def __pos__(self):
        return self

    def __neg__(self):
        return Symbol(-self._arr, dim=self._dim, subs=self._subs)

    def __abs__(self):
        return Symbol(abs(self._arr), dim=self._dim, subs=self._subs)

    def _ensureShapeCompatible(self, other):
        """
        Checks for compatible shapes for binary operations.
        Raises TypeError if not compatible.
        """
        sh0=self.getShape()
        if isinstance(other, Symbol) or isinstance(other, Data):
            sh1=other.getShape()
        elif isinstance(other, numpy.ndarray):
            sh1=other.shape
        elif isinstance(other, list):
            sh1=numpy.array(other).shape
        elif isinstance(other,int) or isinstance(other,float) or isinstance(other,sympy.Basic):
            sh1=()
        else:
            raise TypeError("Unsupported argument type '%s' for operation"%other.__class__.__name__)
        if not sh0==sh1 and not sh0==() and not sh1==():
            raise TypeError("Incompatible shapes for operation")

    def __binaryop(self, op, other):
        """
        Helper for binary operations that checks types, shapes etc.
        """
        self._ensureShapeCompatible(other)
        if isinstance(other, Symbol):
            subs=self._subs.copy()
            subs.update(other._subs)
            dim=other._dim if self._dim < 0 else self._dim
            return Symbol(getattr(self._arr, op)(other._arr), dim=dim, subs=subs)
        if isinstance(other, Data):
            name='data'+str(id(other))
            othersym=Symbol(name, other.getShape(), dim=self._dim)
            subs=self._subs.copy()
            subs.update({Symbol(name):other})
            return Symbol(getattr(self._arr, op)(othersym._arr), dim=self._dim, subs=subs)
        return Symbol(getattr(self._arr, op)(other), dim=self._dim, subs=self._subs)

    def __add__(self, other):
        return self.__binaryop('__add__', other)

    def __radd__(self, other):
        return self.__binaryop('__radd__', other)

    def __sub__(self, other):
        return self.__binaryop('__sub__', other)

    def __rsub__(self, other):
        return self.__binaryop('__rsub__', other)

    def __mul__(self, other):
        return self.__binaryop('__mul__', other)

    def __rmul__(self, other):
        return self.__binaryop('__rmul__', other)

    def __div__(self, other):
        return self.__binaryop('__div__', other)

    def __truediv__(self, other):
        return self.__binaryop('__truediv__', other)
            
    def __rdiv__(self, other):
        return self.__binaryop('__rdiv__', other)
    
    def __rtruediv__(self, other):
        return self.__binaryop('__rtruediv__', other)     
    def __pow__(self, other):
        return self.__binaryop('__pow__', other)

    def __rpow__(self, other):
        return self.__binaryop('__rpow__', other)

    # static methods

    @staticmethod
    def _symComp(sym):
        """
        Extracts the name and component indices from a sympy symbol.

        Parses symbol names like '[x]_0_1' into name 'x' and indices (0, 1).

        :param sym: the sympy symbol to parse
        :type sym: ``sympy.Symbol``
        :return: tuple of (name, component_indices)
        :rtype: ``tuple``
        """
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
        for i in range(len(symbols)):
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

