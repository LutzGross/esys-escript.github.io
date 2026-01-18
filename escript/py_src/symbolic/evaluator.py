
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

"""
Symbolic expression evaluator for escript.

This module provides the `Evaluator` class which compiles symbolic expressions
into callable functions that can be evaluated with escript Data objects.
"""

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"
__author__="Cihan Altinay, Lutz Gross"

__all__ = ['Evaluator']

import esys.escriptcore.util as escript
import esys
from .printer import EsysEscriptPrinter
from .utils import combineData
import math

escript_functions = [ "kronecker", "identity", "zeros", "identityTensor", "identityTensor4", "unitVector", "Lsup", "sup", "inf", "testForZero", "log10", "wherePositive",
                      "whereNegative", "whereNonNegative", "whereNonPositive", "whereZero", "whereNonZero", "Abs", "erf", "sin", "cos", "tan", "asin", "acos", "atan",
                      "atan2", "sinh", "cosh", "tanh", "asinh", "acosh", "atanh", "exp", "sqrt", "log", "sign", "minval", "maxval", "length", "trace", "transpose", "swap_axes",
                      "symmetric", "nonsymmetric", "antisymmetric", "hermitian", "antihermitian", "inverse", "eigenvalues", "mult", "maximum", "minimum", "clip", "cross",
                      "inner", "outer", "matrixmult", "matrix_mult", "tensormult", "tensor_mult", "generalTensorProduct", "transposed_matrix_mult", "transposed_tensor_mult",
                      "generalTransposedTensorProduct", "matrix_transposed_mult", "tensor_transposed_mult",
                      "generalTensorTransposedProduct", "grad", "grad_n", "curl", "integrate", "interpolate", "div", "jump", "L2", "normalize", "deviatoric", "meanValue",
                      "reorderComponents", "positive", "negative", "safeDiv", "condEval", "polarToCart", "phase", "real", "imag", "conjugate" ]


translator = { n :  getattr(escript, n) for n in escript_functions }

class Evaluator(object):
    """
    Evaluates symbolic expressions with escript Data objects.

    The Evaluator compiles symbolic expressions into callable functions that
    can be evaluated efficiently with substituted values.
    """

    def __init__(self, *expressions):
        """
        Creates a symbolic evaluator.

        :param expressions: optional expressions to initialize with
        :type expressions: `Symbol` or ``sympy.Expr``
        """
        self.expressions=[]
        self.symbols=[]
        self.lambdas=[]
        self._subsdict={}
        for ex in expressions:
            self.addExpression(ex)

    def __getstate__(self):
        mydict=self.__dict__.copy()
        # lambda functions cannot be pickled
        mydict['lambdas']=[]
        return mydict

    def __setstate__(self, dict):
        expressions=dict['expressions']
        dict['expressions']=[]
        self.__dict__.update(dict)
        # regenerate lambdas
        for ex in expressions:
            self.addExpression(ex)

    def addExpression(self, expression):
        """
        Adds an expression to this evaluator.

        :param expression: symbolic expression to add
        :type expression: `Symbol` or ``sympy.Expr``
        :return: the modified Evaluator object
        :rtype: `Evaluator`
        :raise TypeError: if argument is not an expression
        """
        import sympy
        from esys import escript
        if not hasattr(expression, "atoms"):
            raise TypeError("Argument is not an expression")
        self.expressions.append(expression)
        if isinstance(expression, sympy.Function):
            sym=set()
            for arg in expression.args:
                if arg is not None:
                    sym.update(arg.atoms(sympy.Symbol))
            sym=tuple(sym)
            self.symbols.append(tuple(sym))
        else:
            sym=tuple(set(expression.atoms(sympy.Symbol)))
            self.symbols.append(sym)
        # Create modules dictionary for lambdify (required for SymPy 1.2+)
        modules = [{'combineData': combineData, 'esys': esys}, translator, math, "numpy"]
        if isinstance(expression, escript.Symbol):
            subs=expression.getDataSubstitutions()
            subs_dict={}
            for s in subs.keys():
                subs_dict[s.name] = subs[s]
            self.subs(**subs_dict)
            self.lambdas.append(sympy.lambdify(sym, expression.lambdarepr(), modules=modules, printer=EsysEscriptPrinter, dummify=False))
        else:
            self.lambdas.append(sympy.lambdify(sym, expression, modules=modules, printer=EsysEscriptPrinter, dummify=False))
        return self

    def subs(self, **args):
        """
        Performs symbol substitution.

        :keyword args: symbol name to value mappings
        :return: the modified Evaluator object
        :rtype: `Evaluator`
        """
        self._subsdict.update(args)
        return self

    def evaluate(self, evalf=False, **args):
        """
        Evaluates all expressions in this evaluator and returns the result.

        :param evalf: if True, call evalf() on Symbol results
        :type evalf: ``bool``
        :keyword args: additional symbol substitutions
        :return: the evaluated expression (single) or tuple of expressions
        :rtype: evaluated result or ``tuple``
        :raise RuntimeError: if not all symbols have values
        """
        from esys import escript
        self.subs(**args)
        res=()
        resEvaled=()
        for i in range(len(self.lambdas)):
            x=self.symbols[i]
            subslist=[self._subsdict[a.name] for a in x if a.name in self._subsdict]
            if len(x)==len(subslist):
                res+=self.lambdas[i](*subslist),
            else:
                raise RuntimeError("Not all symbols have a value")
            if evalf:
                if isinstance(res[i], escript.Symbol):
                    resEvaled+=(res[i].evalf(),)
                else:
                    resEvaled+=(res[i],)
            else:
                resEvaled+=(res[i],)
        if len(res)==1:
            return resEvaled[0]
        else:
            return resEvaled

    def __call__(self, **args):
        """
        Convenience method to substitute and evaluate at once.
        """
        return self.subs(**args).evaluate()

    def __getitem__(self, index):
        """
        Expression accessor.
        """
        return self.expressions[index]

    def __len__(self):
        """
        Returns the number of expressions in this evaluator.
        """
        return len(self.expressions)

    def __iadd__(self, expression):
        """
        Same as addExpression(expression).
        """
        return self.addExpression(expression)

    def __str__(self):
        ret="\n".join([str(e) for e in self.expressions])+"\n"
        for k in sorted(self._subsdict.keys()):
            v=self._subsdict[k]
            if v.__class__.__name__=="Data":
                ret+="%s=<Data object>"%k
            elif v.__class__.__name__=="FunctionSpace":
                ret+="%s=<FunctionSpace object>"%k
            else:
                ret+="%s=%s"%(k,v)
            ret+=", "
        return ret

#
# vim: expandtab shiftwidth=4:
