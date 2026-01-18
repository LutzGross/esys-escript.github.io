"""
SymPy code printer for escript expressions.

This module provides the EsysEscriptPrinter class which converts SymPy
expressions to escript-compatible Python code for use with lambdify.

:note: This module is still under development.
"""

from sympy.core import S
from sympy.core.function import Lambda
from sympy.core.power import Pow
from sympy.printing.pycode import PythonCodePrinter, _known_functions_math, _print_known_const, _print_known_func, _unpack_integral_limits, ArrayPrinter
from sympy.printing.codeprinter import CodePrinter


_not_in_eescript = 'erf erfc factorial gamma loggamma'.split()
_in_eescript = [(k, v) for k, v in _known_functions_math.items() if k not in _not_in_eescript]
_known_functions_eescript = dict(_in_eescript, **{
    'acos': 'acos',
    'acosh': 'acosh',
    'asin': 'asin',
    'asinh': 'asinh',
    'atan': 'atan',
    'atan2': 'atan2',
    'atanh': 'atanh',
    'exp2': 'exp2',
    'sign': 'sign',
    'logaddexp': 'logaddexp',
    'logaddexp2': 'logaddexp2',
    'grad_n': 'grad_n',
    'wherePositive': 'wherePositive',
    'whereNonPositive': 'whereNonPositive',
    'whereNegative': 'whereNegative',
    'whereNonNegative': 'whereNonNegative',
    'whereZero': 'whereZero',
    'whereNonZero': 'whereNonZero',
    'clip': 'clip',
    'maximum': 'maximum',
    'minimum': 'minimum',
    'maxval': 'maxval',
    'minval': 'minval'
})
_known_constants_eescript = {
    'Exp1': 'e',
    'Pi': 'pi',
    'EulerGamma': 'euler_gamma',
    'NaN': 'nan',
    'Infinity': 'PINF',
    'NegativeInfinity': 'NINF'
}

_eescript_known_functions = {k: 'esys.escript.' + v for k, v in _known_functions_eescript.items()}
_eescript_known_constants = {k: 'esys.escript.' + v for k, v in _known_constants_eescript.items()}

class EsysEscriptPrinter(ArrayPrinter, PythonCodePrinter):
    """
    eescript printer which handles vectorized piecewise functions,
    logical operators, etc.
    """

    _module = 'esys.escript'
    _kf = _eescript_known_functions
    _kc = _eescript_known_constants

    def __init__(self, settings=None):
        """
        `settings` is passed to CodePrinter.__init__()
        `module` specifies the array module to use, currently 'eescript', 'CuPy'
        or 'JAX'.
        """
        self.language = "Python with {}".format(self._module)
        self.printmethod = "_{}code".format(self._module)

        self._kf = {**PythonCodePrinter._kf, **self._kf}
        super().__init__(settings=settings)

    def _print(self, expr, **kwargs):
        """Override _print to handle escript Symbol objects"""
        # Check if this is an escript Symbol (not a sympy Symbol)
        if hasattr(expr, 'lambdarepr') and hasattr(expr, 'getShape'):
            # This is an escript Symbol, use its lambdarepr method
            return expr.lambdarepr()
        return super()._print(expr, **kwargs)

    def _print_Symbol(self, expr):
        """Convert symbol names from [x]_0 format to x[0] format"""
        import re
        name = expr.name
        # Check if the symbol name matches the pattern [name]_i_j_...
        match = re.match(r'\[(\w+)\](_\d+)+$', name)
        if match:
            base_name = match.group(1)
            # Extract all the index numbers
            indices = re.findall(r'_(\d+)', name)
            # Convert to Python array notation: x[0,1] or x[0]
            return f"{base_name}[{','.join(indices)}]"
        return super()._print_Symbol(expr)

    def _print_seq(self, seq):
        "General sequence printer: converts to tuple"
        # Print tuples here instead of lists because numba supports
        #     tuples in nopython mode.
        delimiter=', '
        return '({},)'.format(delimiter.join(self._print(item) for item in seq))

    def _print_MatMul(self, expr):
        "Matrix multiplication printer"
        if expr.as_coeff_matrices()[0] is not S.One:
            expr_list = expr.as_coeff_matrices()[1]+[(expr.as_coeff_matrices()[0])]
            return '({})'.format(').dot('.join(self._print(i) for i in expr_list))
        return '({})'.format(').dot('.join(self._print(i) for i in expr.args))

    def _print_MatPow(self, expr):
        "Matrix power printer"
        return '{}({}, {})'.format(self._module_format(self._module + '.linalg.matrix_power'),
            self._print(expr.args[0]), self._print(expr.args[1]))

    def _print_Inverse(self, expr):
        "Matrix inverse printer"
        return '{}({})'.format(self._module_format(self._module + '.linalg.inv'),
            self._print(expr.args[0]))

    def _print_DotProduct(self, expr):
        # DotProduct allows any shape order, but eescript.dot does matrix
        # multiplication, so we have to make sure it gets 1 x n by n x 1.
        arg1, arg2 = expr.args
        if arg1.shape[0] != 1:
            arg1 = arg1.T
        if arg2.shape[1] != 1:
            arg2 = arg2.T

        return "%s(%s, %s)" % (self._module_format(self._module + '.dot'),
                               self._print(arg1),
                               self._print(arg2))

    def _print_MatrixSolve(self, expr):
        return "%s(%s, %s)" % (self._module_format(self._module + '.linalg.solve'),
                               self._print(expr.matrix),
                               self._print(expr.vector))

    def _print_ZeroMatrix(self, expr):
        return '{}({})'.format(self._module_format(self._module + '.zeros'),
            self._print(expr.shape))

    def _print_OneMatrix(self, expr):
        return '{}({})'.format(self._module_format(self._module + '.ones'),
            self._print(expr.shape))

    def _print_FunctionMatrix(self, expr):
        from sympy.abc import i, j
        lamda = expr.lamda
        if not isinstance(lamda, Lambda):
            lamda = Lambda((i, j), lamda(i, j))
        return '{}(lambda {}: {}, {})'.format(self._module_format(self._module + '.fromfunction'),
            ', '.join(self._print(arg) for arg in lamda.args[0]),
            self._print(lamda.args[1]), self._print(expr.shape))

    def _print_HadamardProduct(self, expr):
        func = self._module_format(self._module + '.multiply')
        return ''.join('{}({}, '.format(func, self._print(arg)) \
            for arg in expr.args[:-1]) + "{}{}".format(self._print(expr.args[-1]),
            ')' * (len(expr.args) - 1))

    def _print_KroneckerProduct(self, expr):
        func = self._module_format(self._module + '.kron')
        return ''.join('{}({}, '.format(func, self._print(arg)) \
            for arg in expr.args[:-1]) + "{}{}".format(self._print(expr.args[-1]),
            ')' * (len(expr.args) - 1))

    def _print_Adjoint(self, expr):
        return '{}({}({}))'.format(
            self._module_format(self._module + '.conjugate'),
            self._module_format(self._module + '.transpose'),
            self._print(expr.args[0]))

    def _print_DiagonalOf(self, expr):
        vect = '{}({})'.format(
            self._module_format(self._module + '.diag'),
            self._print(expr.arg))
        return '{}({}, (-1, 1))'.format(
            self._module_format(self._module + '.reshape'), vect)

    def _print_DiagMatrix(self, expr):
        return '{}({})'.format(self._module_format(self._module + '.diagflat'),
            self._print(expr.args[0]))

    def _print_DiagonalMatrix(self, expr):
        return '{}({}, {}({}, {}))'.format(self._module_format(self._module + '.multiply'),
            self._print(expr.arg), self._module_format(self._module + '.eye'),
            self._print(expr.shape[0]), self._print(expr.shape[1]))

    def _print_Piecewise(self, expr):
        "Piecewise function printer"
        from sympy.logic.boolalg import ITE, simplify_logic
        def print_cond(cond):
            """ Problem having an ITE in the cond. """
            if cond.has(ITE):
                return self._print(simplify_logic(cond))
            else:
                return self._print(cond)
        exprs = '[{}]'.format(','.join(self._print(arg.expr) for arg in expr.args))
        conds = '[{}]'.format(','.join(print_cond(arg.cond) for arg in expr.args))
        # If [default_value, True] is a (expr, cond) sequence in a Piecewise object
        #     it will behave the same as passing the 'default' kwarg to select()
        #     *as long as* it is the last element in expr.args.
        # If this is not the case, it may be triggered prematurely.
        return '{}({}, {}, default={})'.format(
            self._module_format(self._module + '.select'), conds, exprs,
            self._print(S.NaN))

    def _print_Relational(self, expr):
        "Relational printer for Equality and Unequality"
        op = {
            '==' :'equal',
            '!=' :'not_equal',
            '<'  :'less',
            '<=' :'less_equal',
            '>'  :'greater',
            '>=' :'greater_equal',
        }
        if expr.rel_op in op:
            lhs = self._print(expr.lhs)
            rhs = self._print(expr.rhs)
            return '{op}({lhs}, {rhs})'.format(op=self._module_format(self._module + '.'+op[expr.rel_op]),
                                               lhs=lhs, rhs=rhs)
        return super()._print_Relational(expr)

    def _print_And(self, expr):
        "Logical And printer"
        # We have to override LambdaPrinter because it uses Python 'and' keyword.
        # If LambdaPrinter didn't define it, we could use StrPrinter's
        # version of the function and add 'logical_and' to eescript_TRANSLATIONS.
        return '{}.reduce(({}))'.format(self._module_format(self._module + '.logical_and'), ','.join(self._print(i) for i in expr.args))

    def _print_Or(self, expr):
        "Logical Or printer"
        # We have to override LambdaPrinter because it uses Python 'or' keyword.
        # If LambdaPrinter didn't define it, we could use StrPrinter's
        # version of the function and add 'logical_or' to eescript_TRANSLATIONS.
        return '{}.reduce(({}))'.format(self._module_format(self._module + '.logical_or'), ','.join(self._print(i) for i in expr.args))

    def _print_Not(self, expr):
        "Logical Not printer"
        # We have to override LambdaPrinter because it uses Python 'not' keyword.
        # If LambdaPrinter didn't define it, we would still have to define our
        #     own because StrPrinter doesn't define it.
        return '{}({})'.format(self._module_format(self._module + '.logical_not'), ','.join(self._print(i) for i in expr.args))

    def _print_Pow(self, expr, rational=False):
        # XXX Workaround for negative integer power error
        if expr.exp.is_integer and expr.exp.is_negative:
            expr = Pow(expr.base, expr.exp.evalf(), evaluate=False)
        return self._hprint_Pow(expr, rational=rational, sqrt=self._module + '.sqrt')

    def _print_Min(self, expr):
        return '{}(({}), axis=0)'.format(self._module_format(self._module + '.amin'), ','.join(self._print(i) for i in expr.args))

    def _print_Max(self, expr):
        return '{}(({}), axis=0)'.format(self._module_format(self._module + '.amax'), ','.join(self._print(i) for i in expr.args))

    def _print_arg(self, expr):
        return "%s(%s)" % (self._module_format(self._module + '.angle'), self._print(expr.args[0]))

    def _print_im(self, expr):
        return "%s(%s)" % (self._module_format(self._module + '.imag'), self._print(expr.args[0]))

    def _print_Mod(self, expr):
        return "%s(%s)" % (self._module_format(self._module + '.mod'), ', '.join(
            map(lambda arg: self._print(arg), expr.args)))

    def _print_re(self, expr):
        return "%s(%s)" % (self._module_format(self._module + '.real'), self._print(expr.args[0]))

    def _print_sinc(self, expr):
        return "%s(%s)" % (self._module_format(self._module + '.sinc'), self._print(expr.args[0]/S.Pi))

    def _print_MatrixBase(self, expr):
        func = self.known_functions.get(expr.__class__.__name__, None)
        if func is None:
            func = self._module_format(self._module + '.array')
        return "%s(%s)" % (func, self._print(expr.tolist()))

    def _print_Identity(self, expr):
        shape = expr.shape
        if all(dim.is_Integer for dim in shape):
            return "%s(%s)" % (self._module_format(self._module + '.eye'), self._print(expr.shape[0]))
        else:
            raise NotImplementedError("Symbolic matrix dimensions are not yet supported for identity matrices")

    def _print_BlockMatrix(self, expr):
        return '{}({})'.format(self._module_format(self._module + '.block'),
                                 self._print(expr.args[0].tolist()))

    def _print_NDimArray(self, expr):
        if len(expr.shape) == 1:
            return self._module + '.array(' + self._print(expr.args[0]) + ')'
        if len(expr.shape) == 2:
            return self._print(expr.tomatrix())
        # Should be possible to extend to more dimensions
        return CodePrinter._print_not_supported(self, expr)

    def _print_grad_n(self, expr):
        return "%s(%s)" % (self._module_format(self._module + '.grad_n'), self._print(expr.args[0]))

    def _print_whereZero(self, expr):
        return "%s(%s)" % (self._module_format(self._module + '.whereZero'), self._print(expr.args[0]))

    _add = "add"
    _einsum = "einsum"
    _transpose = "transpose"
    _ones = "ones"
    _zeros = "zeros"

    _print_lowergamma = CodePrinter._print_not_supported
    _print_uppergamma = CodePrinter._print_not_supported
    _print_fresnelc = CodePrinter._print_not_supported
    _print_fresnels = CodePrinter._print_not_supported

for func in _eescript_known_functions:
    setattr(EsysEscriptPrinter, f'_print_{func}', _print_known_func)

for const in _eescript_known_constants:
    setattr(EsysEscriptPrinter, f'_print_{const}', _print_known_const)
