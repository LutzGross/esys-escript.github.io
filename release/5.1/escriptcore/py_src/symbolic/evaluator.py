
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

__all__ = ['Evaluator']

import distutils.version as duv

"""
Symbolic expression evaluator for escript
"""

class Evaluator(object):
    def __init__(self, *expressions):
        """
        Returns a symbolic evaluator.

        :param expressions: optional expressions to initialise with
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

        :return: the modified Evaluator object
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
            self.symbols.append(tuple(sym))
        else:
            sym=set(expression.atoms(sympy.Symbol))
            self.symbols.append(sym)

        if isinstance(expression, escript.Symbol):
            subs=expression.getDataSubstitutions()
            subs_dict={}
            for s in subs.keys():
                subs_dict[s.name] = subs[s]
            self.subs(**subs_dict)
            if duv.LooseVersion(sympy.__version__) < duv.LooseVersion('0.7.6'):
                self.lambdas.append(sympy.lambdify(sym, expression.lambdarepr(), ["escript","numpy"]))
            else:
                self.lambdas.append(sympy.lambdify(sym, expression.lambdarepr(), ["escript","numpy"],dummify=False))
        else:
            if duv.LooseVersion(sympy.__version__) < duv.LooseVersion('0.7.6'):
                self.lambdas.append(sympy.lambdify(sym, expression, ["escript","numpy"]))
            else:
                self.lambdas.append(sympy.lambdify(sym, expression, ["escript","numpy"], dummify=False))
        return self

    def subs(self, **args):
        """
        Symbol substitution.

        :return: the modified Evaluator object
        """
        self._subsdict.update(args)
        return self

    def evaluate(self, evalf=False, **args):
        """
        Evaluates all expressions in this evaluator and returns the result
        as a tuple.

        :return: the evaluated expressions in the order they were added to
                 this Evaluator.
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
