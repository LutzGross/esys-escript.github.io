# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2013 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

__copyright__="""Copyright (c) 2003-2013 by University of Queensland
http://www.uq.edu.au
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
from .start import HAVE_SYMBOLS
import numpy
from time import time
from . import linearPDEs as lpe
from . import util
from .escriptcpp import Data

if HAVE_SYMBOLS:
    import sympy
    import esys.escriptcore.symbolic as symb

__author__="Cihan Altinay, Lutz Gross"

class iteration_steps_maxReached(Exception):
    """
    Exception thrown if the maximum number of iteration steps is reached.
    """
    pass
class DivergenceDetected(Exception):
    """
    Exception thrown if Newton-Raphson did not converge.
    """
    pass
class InadmissiblePDEOrdering(Exception):
    """
    Exception thrown if the ordering of the PDE equations should be revised.
    """
    pass

def concatenateRow(*args):
    """
    """
    if len(args)==1:
        return args[0]

    if len(args)>2:
        return concatenateRow(args[0], concatenateRow(*args[1:]))

    lshape=util.getShape(args[0])
    rshape=util.getShape(args[1])
    if len(lshape)==0:
        if len(rshape)==0:
            shape=(2,)
            res=numpy.zeros(shape, dtype=object)
            res[0]=args[0]
            res[1]=args[1]
        elif len(rshape)==1:
            shape=(rshape[0]+1,)
            res=numpy.zeros(shape, dtype=object)
            res[0]=args[0]
            res[1:]=args[1]
    elif len(lshape)==1:
        if len(rshape)==1:
            shape=(2,)+lshape
            res=numpy.zeros(shape, dtype=object)
            res[0]=args[0]
            res[1]=args[1]
        else:
            shape=(rshape[0]+1,)+lshape
            res=numpy.zeros(shape, dtype=object)
            res[0]=args[0]
            res[1:]=args[1]
    else:
        if len(rshape)==1:
            shape=(lshape[0]+1,)+lshape[1:]
            res=numpy.zeros(shape, dtype=object)
            res[:lshape[0]]=args[0]
            res[lshape[0]]=args[1]
        else:
            shape=(lshape[0]+rshape[0],)+lshape[1:]
            res=numpy.zeros(shape, dtype=object)
            res[:lshape[0]]=args[0]
            res[lshape[0]:]=args[1]

    subs=args[0].getDataSubstitutions().copy()
    subs.update(args[1].getDataSubstitutions())
    dim=args[1].getDim() if args[0].getDim()<0 else args[0].getDim()
    return symb.Symbol(res, dim=dim, subs=subs)

class NonlinearPDE(object):
    """
    This class is used to define a general nonlinear, steady, second order PDE
    for an unknown function *u* on a given domain defined through a `Domain`
    object.

    For a single PDE having a solution with a single component the nonlinear
    PDE is defined in the following form:

    *-div(X) + Y = 0*

    where *X*,*Y*=f(*u*,*grad(u)*), *div(F)* denotes the divergence of *F* and
    *grad(F)* is the spatial derivative of *F*.

    The coefficients *X* (rank 1) and *Y* (scalar) have to be specified through
    `Symbol` objects.

    The following natural boundary conditions are considered:

    *n[j]*X[j] - y = 0*

    where *n* is the outer normal field. Notice that the coefficient *X*
    is defined in the PDE. The coefficient *y* is a scalar `Symbol`.

    Constraints for the solution prescribe the value of the solution at certain
    locations in the domain. They have the form

    *u=r* where *q>0*

    *r* and *q* are each scalar where *q* is the characteristic function
    defining where the constraint is applied. The constraints override any
    other condition set by the PDE or the boundary condition.

    For a system of PDEs and a solution with several components, *u* is rank
    one, while the PDE coefficient *X* is rank two and *Y* and *y* is rank one.

    The PDE is solved by linearising the coefficients and iteratively solving
    the corresponding linear PDE until the error is smaller than a tolerance
    or a maximum number of iterations is reached.

    Typical usage::

        u = Symbol('u', dim=dom.getDim())
        p = NonlinearPDE(dom, u)
        p.setValue(X=grad(u), Y=5*u)
        v = p.getSolution(u=0.)
    """
    DEBUG0=0 # no debug info
    DEBUG1=1 # info on Newton-Raphson solver printed
    DEBUG2=2 # in addition info from linear solver is printed
    DEBUG3=3 # in addition info on linearization is printed
    DEBUG4=4 # in addition info on LinearPDE handling is printed

    ORDER=0
    def __init__(self, domain, u, debug=DEBUG0):
        """
        Initializes a new nonlinear PDE.

        :param domain: domain of the PDE
        :type domain: `Domain`
        :param u: The symbol for the unknown PDE function u.
        :type u: `Symbol`
        :param debug: level of debug information to be printed
        """
        if not HAVE_SYMBOLS:
            raise RuntimeError("Trying to instantiate a NonlinearPDE but sympy not available")

        self.__COEFFICIENTS = [ "X", "X_reduced", "Y", "Y_reduced", "y", "y_reduced", "y_contact", "y_contact_reduced", "y_dirac"]

        self._r=Data()
        self._set_coeffs={}
        self._unknown=u
        self._debug=debug
        self._rtol=1e-6
        self._atol=0.
        self._iteration_steps_max=100
        self._omega_min=0.0001
        self._quadratic_convergence_limit=0.2
        self._simplified_newton_limit=0.1

        self.dim = domain.getDim()
        if u.getRank()==0:
            numEquations=1
        else:
            numEquations=u.getShape()[0]
        numSolutions=numEquations
        self._lpde=lpe.LinearPDE(domain,numEquations,numSolutions,self._debug > self.DEBUG4 )

    def __str__(self):
        """
        Returns the string representation of the PDE

        :return: a simple representation of the PDE
        :rtype: ``str``
        """
        return "<%s %d>"%(self.__class__.__name__, id(self))

    def getUnknownSymbol(self):
        """
        Returns the symbol of the PDE unknown

        :return: the symbol of the PDE unknown
        :rtype: `Symbol`
        """
        return self._unknown

    def trace1(self, text):
        """
        Prints the text message if the debug level is greater than DEBUG0

        :param text: message to be printed
        :type text: ``string``
        """
        if self._debug > self.DEBUG0:
            print("%s: %s"%(str(self), text))

    def trace3(self, text):
        """
        Prints the text message if the debug level is greater than DEBUG3

        :param text: message to be printed
        :type text: ``string``
        """
        if self._debug > self.DEBUG2:
            print("%s: %s"%(str(self), text))

    def getLinearSolverOptions(self):
        """
        Returns the options of the linear PDE solver class
        """
        return self._lpde.getSolverOptions()

    def getLinearPDE(self):
        """
        Returns the linear PDE used to calculate the Newton-Raphson update

        :rtype: `LinearPDE`
        """
        return self._lpde

    def setOptions(self, **opts):
        """
        Allows setting options for the nonlinear PDE.

        The supported options are:
          ``tolerance``
                error tolerance for the Newton method
          ``iteration_steps_max``
                maximum number of Newton iterations
          ``omega_min``
                minimum relaxation factor
          ``atol``
                solution norms less than ``atol`` are assumed to be ``atol``.
                This can be useful if one of your solutions is expected to
                be zero.
          ``quadratic_convergence_limit``
                if the norm of the Newton-Raphson correction is reduced by
                less than ``quadratic_convergence_limit`` between two iteration
                steps quadratic convergence is assumed.
          ``simplified_newton_limit``
                if the norm of the defect is reduced by less than
                ``simplified_newton_limit`` between two iteration steps and
                quadratic convergence is detected the iteration switches to the
                simplified Newton-Raphson scheme.
        """
        for key in opts:
            if key=='tolerance':
                self._rtol=opts[key]
            elif key=='iteration_steps_max':
                self._iteration_steps_max=opts[key]
            elif key=='omega_min':
                self._omega_min=opts[key]
            elif key=='atol':
                self._atol=opts[key]
            elif key=='quadratic_convergence_limit':
                self._quadratic_convergence_limit=opts[key]
            elif key=='simplified_newton_limit':
                self._simplified_newton_limit=opts[key]
            else:
                raise KeyError("Invalid option '%s'"%key)

    def getSolution(self, **subs):
        """
        Returns the solution of the PDE.

        :param subs: Substitutions for all symbols used in the coefficients
                     including the initial value for the unknown *u*.
        :return: the solution
        :rtype: `Data`
        """
        # get the initial value for the iteration process

        # collect components of unknown in u_syms
        u_syms=[]
        simple_u=False
        for i in numpy.ndindex(self._unknown.getShape()):
            u_syms.append(symb.Symbol(self._unknown[i]).atoms(sympy.Symbol).pop().name)
        if len(set(u_syms))==1: simple_u=True
        
        e=symb.Evaluator(self._unknown)
        for sym in u_syms:
            if not subs.has_key(sym):
                raise KeyError("Initial value for '%s' missing."%sym)
            if not isinstance(subs[sym], Data):
                subs[sym]=Data(subs[sym], self._lpde.getFunctionSpaceForSolution())
            e.subs(**{sym:subs[sym]})
        ui=e.evaluate()

        # modify ui so it meets the constraints:
        q=self._lpde.getCoefficient("q")
        if not q.isEmpty():
            if hasattr(self, "_r"):
                r=self._r
                if symb.isSymbol(r):
                    r=symb.Evaluator(r).evaluate(**subs)
                elif not isinstance(r, Data):
                    r=Data(r, self._lpde.getFunctionSpaceForSolution())
                elif r.isEmpty():
                    r=0
            else:
                r=0
            ui = q * r + (1-q) * ui

        # separate symbolic expressions from other coefficients
        constants={}
        expressions={}
        for n, e in self._set_coeffs.items():
            if symb.isSymbol(e):
                expressions[n]=e
            else:
                constants[n]=e

        # set constant PDE values now
        self._lpde.setValue(**constants)
        self._lpde.getSolverOptions().setAbsoluteTolerance(0.)
        self._lpde.getSolverOptions().setVerbosity(self._debug > self.DEBUG1)
        #=====================================================================
        # perform Newton iterations until error is small enough or
        # maximum number of iterations reached
        n=0
        omega=1. # relaxation factor
        use_simplified_Newton=False
        defect_norm=None
        delta_norm=None
        converged=False

        #subs[u_sym]=ui
        if simple_u:
            subs[u_syms[0]]=ui
        else:
            for i in range(len(u_syms)):
                subs[u_syms[i]]=ui[i]

        while not converged:
            if n > self._iteration_steps_max:
               raise iteration_steps_maxReached("maximum number of iteration steps reached, giving up.")
            self.trace1(5*"="+" iteration step %d "%n + 5*"=")
            # calculate the correction delta_u
            if n == 0:
                self._updateLinearPDE(expressions, subs)
                defect_norm=self._getDefectNorm(self._lpde.getRightHandSide())
                LINTOL=0.1
            else:
                if not use_simplified_Newton:
                    self._updateMatrix(expressions, subs)
                if q_u == None:
                    LINTOL = 0.1 * min(qtol/defect_norm)
                else:
                    LINTOL = 0.1 * max( q_u**2, min(qtol/defect_norm) )
                LINTOL=max(1e-4,min(LINTOL, 0.1))
                #LINTOL=1.e-5
            self._lpde.getSolverOptions().setTolerance(LINTOL)
            self.trace1("PDE is solved with rel. tolerance = %e"%LINTOL)
            delta_u=self._lpde.getSolution()

            #check for reduced defect:
            omega=min(2*omega, 1.) # raise omega
            defect_reduced=False
            ui_old=ui
            while not defect_reduced:
                ui=ui_old-delta_u * omega
                if simple_u:
                    subs[u_syms[0]]=ui
                else:
                    for i in range(len(u_syms)):
                        subs[u_syms[i]]=ui[i]
                self._updateRHS(expressions, subs)
                new_defect_norm=self._getDefectNorm(self._lpde.getRightHandSide())
                defect_reduced=False
                for i in range(len( new_defect_norm)):
                    if new_defect_norm[i] < defect_norm[i]: defect_reduced=True
                
                #print new_defect_norm
                #q_defect=max(self._getSafeRatio(new_defect_norm, defect_norm))
                # if defect_norm==0 and new_defect_norm!=0
                # this will be util.DBLE_MAX
                #self.trace1("Defect reduction = %e with relaxation factor %e."%(q_defect, omega))
                if not defect_reduced:
                    omega*=0.5
                    if omega < self._omega_min:
                        raise DivergenceDetected("Underrelaxtion failed to reduce defect, giving up.")
            self.trace1("Defect reduction with relaxation factor %e."%(omega, ))
            delta_norm, delta_norm_old = self._getSolutionNorm(delta_u) * omega, delta_norm
            defect_norm, defect_norm_old = new_defect_norm, defect_norm
            u_norm=self._getSolutionNorm(ui, atol=self._atol)
            # set the tolerance on equation level:
            qtol=self._getSafeRatio(defect_norm_old * u_norm * self._rtol, delta_norm)
            # if defect_norm_old==0 and defect_norm_old!=0 this will be util.DBLE_MAX
            #    -> the ordering of the equations is not appropriate.
            # if defect_norm_old==0 and defect_norm_old==0  this is zero so
            # convergence can happen for defect_norm==0 only.
            if not max(qtol)<util.DBLE_MAX:
                raise InadmissiblePDEOrdering("Review ordering of PDE equations.")
            # check stopping criteria
            if not delta_norm_old == None:
                q_u=max(self._getSafeRatio(delta_norm, delta_norm_old))
                # if delta_norm_old==0 and delta_norm!=0
                # this will be util.DBLE_MAX
                if q_u <= self._quadratic_convergence_limit and not omega<1. :
                    quadratic_convergence=True
                    self.trace1("Quadratic convergence detected (rate %e <= %e)"%(q_u, self._quadratic_convergence_limit))

                    converged = all( [ defect_norm[i] <= qtol[i] for i in range(len(qtol)) ])

                else:
                    self.trace1("No quadratic convergence detected (rate %e > %e, omega=%e)"%(q_u, self._quadratic_convergence_limit,omega ))
                    quadratic_convergence=False
                    converged=False
            else:
                q_u=None
                converged=False
                quadratic_convergence=False
            if self._debug > self.DEBUG0:
                for i in range(len(u_norm)):
                    self.trace1("Component %s: u: %e, du: %e, defect: %e, qtol: %e"%(i, u_norm[i], delta_norm[i], defect_norm[i], qtol[i]))
                if converged: self.trace1("Iteration has converged.")
            # Can we switch to simplified Newton?
            if quadratic_convergence:
                q_defect=max(self._getSafeRatio(defect_norm, defect_norm_old))
                if q_defect < self._simplified_newton_limit:
                    use_simplified_Newton=True
                    self.trace1("Simplified Newton-Raphson is applied (rate %e < %e)."%(q_defect, self._simplified_newton_limit))
            n+=1
        self.trace1(5*"="+" Newton-Raphson iteration completed after %d steps "%n + 5*"=")
        return ui

    def _getDefectNorm(self, f):
        """
        calculates the norm of the defect ``f``

        :param f: defect vector
        :type f: `Data` of rank 0 or 1.
        :return: component-by-component norm of ``f``
        :rtype: ``numpy.array`` of rank 1
        :raise ValueError: if shape if ``f`` is incorrect.

        """
        # this still needs some work!!!
        out=[]
        s=f.getShape()
        if len(s) == 0:
            out.append(util.Lsup(f))
        elif len(s) == 1:
            for i in range(s[0]):
                out.append(util.Lsup(f[i]))
        else:
            raise ValueError("Illegal shape of defect vector: %s"%s)
        return numpy.array(out)

    def _getSolutionNorm(self, u, atol=0.):
        """
        calculates the norm of the solution ``u``

        :param u: solution vector
        :type u: `Data` of rank 0 or 1.
        :return: component-by-component norm of ``u``
        :rtype: ``numpy.array`` of rank 1
        :raise ValueError: if shape of ``u`` is incorrect.

        """
        out=[]
        s=u.getShape()
        if len(s) == 0:
            out.append(max(util.Lsup(u),atol) )
        elif len(s) == 1:
            for i in range(s[0]):
                out.append(max(util.Lsup(u[i]), atol))
        else:
            raise ValueError("Illegal shape of solution: %s"%s)
        return numpy.array(out)

    def _getSafeRatio(self, a , b):
        """
        Returns the component-by-component ratio of ''a'' and ''b''

        If for a component ``i`` the values ``a[i]`` and ``b[i]`` are both
        equal to zero their ratio is set to zero.
        If ``b[i]`` equals zero but ``a[i]`` is positive the ratio is set to
        `util.DBLE_MAX`.

        :param a: numerator
        :param b: denominator
        :type a: ``numpy.array`` of rank 1 with non-negative entries.
        :type b: ``numpy.array`` of rank 1 with non-negative entries.
        :return: ratio of ``a`` and ``b``
        :rtype: ``numpy.array``
        """
        out=0.
        if a.shape !=b.shape:
            raise ValueError("shapes must match.")
        s=a.shape
        if len(s) != 1:
            raise ValueError("rank one is expected.")
        out=numpy.zeros(s)
        for i in range(s[0]):
            if abs(b[i]) > 0:
                out[i]=a[i]/b[i]
            elif abs(a[i]) > 0:
                out[i] = util.DBLE_MAX
            else:
                out[i] = 0
        return out

    def getNumSolutions(self):
        """
        Returns the number of the solution components
        :rtype: ``int``
        """
        s=self._unknown.getShape()
        if len(s) > 0:
            return s[0]
        else:
            return 1

    def getShapeOfCoefficient(self,name):
        """
        Returns the shape of the coefficient ``name``

        :param name: name of the coefficient enquired
        :type name: ``string``
        :return: the shape of the coefficient ``name``
        :rtype: ``tuple`` of ``int``
        :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
        """
        numSol=self.getNumSolutions()
        dim = self.dim
        if name=="X" or name=="X_reduced":
            if numSol > 1:
                return (numSol,dim)
            else:
                return (dim,)
        elif name=="r" or name == "q" :
            if numSol > 1:
                return (numSol,)
            else:
                return ()
        elif name=="Y" or name=="Y_reduced":
            if numSol > 1:
                return (numSol,)
            else:
                return ()
        elif name=="y" or name=="y_reduced":
            if numSol > 1:
                return (numSol,)
            else:
                return ()
        elif name=="y_contact" or name=="y_contact_reduced":
            if numSol > 1:
                return (numSol,)
            else:
                return ()
        elif name=="y_dirac":
            if numSol > 1:
                return (numSol,)
            else:
                return ()
        else:
            raise lpe.IllegalCoefficient("Attempt to request unknown coefficient %s"%name)

    def createCoefficient(self, name):
        """
        Creates a new coefficient ``name`` as Symbol

        :param name: name of the coefficient requested
        :type name: ``string``
        :return: the value of the coefficient
        :rtype: `Symbol` or `Data` (for name = "q")
        :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
        """
        if name == "q":
            return self._lpde.createCoefficient("q")
        else:
           s=self.getShapeOfCoefficient(name)
           return symb.Symbol(name, s, dim=self.dim)

    def getCoefficient(self, name):
        """
        Returns the value of the coefficient ``name`` as Symbol

        :param name: name of the coefficient requested
        :type name: ``string``
        :return: the value of the coefficient
        :rtype: `Symbol`
        :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
        """
        if self._set_coeffs.has_key(name):
            return self._set_coeffs[name]
        elif name == "r":
             if hasattr(self, "_r"):
                 return self._r
             else:
                 raise lpe.IllegalCoefficient("Attempt to request undefined coefficient %s"%name)
        elif name == "q":
             return self._lpde.getCoefficient("q")
        else:
            raise lpe.IllegalCoefficient("Attempt to request undefined coefficient %s"%name)

    def setValue(self,**coefficients):
        """
        Sets new values to one or more coefficients.

        :param coefficients: new values assigned to coefficients

        :param coefficients: new values assigned to coefficients
        :keyword X: value for coefficient ``X``
        :type X: `Symbol` or any type that can be cast to a `Data` object
        :keyword Y: value for coefficient ``Y``
        :type Y: `Symbol` or any type that can be cast to a `Data` object
        :keyword y: value for coefficient ``y``
        :type y: `Symbol` or any type that can be cast to a `Data` object
        :keyword y_contact: value for coefficient ``y_contact``
        :type y_contact: `Symbol` or any type that can be cast to a `Data` object
        :keyword y_dirac: value for coefficient ``y_dirac``
        :type y_dirac: `Symbol` or any type that can be cast to a `Data` object
        :keyword q: mask for location of constraint
        :type q: any type that can be cast to a `Data` object
        :keyword r: value of solution prescribed by constraint
        :type r: `Symbol` or any type that can be cast to a `Data` object
        :raise IllegalCoefficient: if an unknown coefficient keyword is used
        :raise IllegalCoefficientValue: if a supplied coefficient value has an
                                        invalid shape
        """

        u=self._unknown
        for name,val in coefficients.iteritems():
            shape=util.getShape(val)
            if not shape == self.getShapeOfCoefficient(name):
                raise lpe.IllegalCoefficientValue("%s has shape %s but must have shape %s"%(name, shape, self.getShapeOfCoefficient(name)))
            rank=len(shape)
            if name == "q":
                self._lpde.setValue(q=val)
            elif name == "r":
                self._r=val
            elif name=="X" or name=="X_reduced":
                if rank != u.getRank()+1:
                    raise lpe.IllegalCoefficientValue("%s must have rank %d"%(name,u.getRank()+1))
                T0=time()
                B,A=symb.getTotalDifferential(val, u, 1)
                if name=='X_reduced':
                    self.trace3("Computing A_reduced, B_reduced took %f seconds."%(time()-T0))
                    self._set_coeffs['A_reduced']=A
                    self._set_coeffs['B_reduced']=B
                    self._set_coeffs['X_reduced']=val
                else:
                    self.trace3("Computing A, B took %f seconds."%(time()-T0))
                    self._set_coeffs['A']=A
                    self._set_coeffs['B']=B
                    self._set_coeffs['X']=val
            elif name=="Y" or name=="Y_reduced":
                if rank != u.getRank():
                    raise lpe.IllegalCoefficientValue("%s must have rank %d"%(name,u.getRank()))
                T0=time()
                D,C=symb.getTotalDifferential(val, u, 1)
                if name=='Y_reduced':
                    self.trace3("Computing C_reduced, D_reduced took %f seconds."%(time()-T0))
                    self._set_coeffs['C_reduced']=C
                    self._set_coeffs['D_reduced']=D
                    self._set_coeffs['Y_reduced']=val
                else:
                    self.trace3("Computing C, D took %f seconds."%(time()-T0))
                    self._set_coeffs['C']=C
                    self._set_coeffs['D']=D
                    self._set_coeffs['Y']=val
            elif name in ("y", "y_reduced", "y_contact", "y_contact_reduced", \
                    "y_dirac"):
                y=val
                if rank != u.getRank():
                    raise lpe.IllegalCoefficientValue("%s must have rank %d"%(name,u.getRank()))
                if not hasattr(y, 'diff'):
                    d=numpy.zeros(u.getShape())
                else:
                    d=y.diff(u)
                self._set_coeffs[name]=y
                self._set_coeffs['d'+name[1:]]=d
            else:
                raise lpe.IllegalCoefficient("Attempt to set unknown coefficient %s"%name)

    def getSensitivity(self, f, g=None, **subs):
        """
        Calculates the sensitivity of the solution of an input factor ``f``
        in direction ``g``.

        :param f: the input factor to be investigated. ``f`` may be of rank 0
                  or 1.
        :type f: `Symbol`
        :param g: the direction(s) of change.
                  If not present, it is *g=eye(n)* where ``n`` is the number of
                  components of ``f``.
        :type g: ``list`` or single of ``float``, ``numpy.array`` or `Data`.
        :param subs: Substitutions for all symbols used in the coefficients
                     including unknown *u* and the input factor ``f`` to be
                     investigated
        :return: the sensitivity
        :rtype: `Data` with shape  *u.getShape()+(len(g),)* if *len(g)>1* or
                *u.getShape()* if *len(g)==1*
        """
        s_f=f.getShape()
        if len(s_f) == 0:
            len_f=1
        elif len(s_f) == 1:
            len_f=s_f[0]
        else:
            raise ValueError("rank of input factor must be zero or one.")

        if not g == None:
           if len(s_f) == 0:
               if not isinstance(g, list): g=[g]
           else:
              if isinstance(g, list):
                  if len(g) == 0:
                      raise ValueError("no direction given.")
                  if len(getShape(g[0])) == 0: g=[g] # if g[0] is a scalar we assume that the list g is to be interprested a data object
              else:
                 g=[g]
           # at this point g is a list of directions:
           len_g=len(g)
           for g_i in g:
                if not getShape(g_i) == s_f:
                    raise ValueError("shape of direction (=%s) must match rank of input factor (=%s)"%(getShape(g_i) , s_f) )
        else:
           len_g=len_f

        #*** at this point g is a list of direction or None and len_g is the
        #    number of directions to be investigated.

        # now we make sure that the operator in the lpde is set (it is the same
        # as for the Newton-Raphson scheme)
        # if the solution etc are cached this could be omitted:
        constants={}
        expressions={}
        for n, e in self._set_coeffs.items():
            if n not in self.__COEFFICIENTS:
                if symb.isSymbol(e):
                    expressions[n]=e
                else:
                    constants[n]=e
        self._lpde.setValue(**constants)
        self._updateMatrix(self, expressions, subs)
        #=====================================================================
        self._lpde.getSolverOptions().setAbsoluteTolerance(0.)
        self._lpde.getSolverOptions().setTolerance(self._rtol)
        self._lpde.getSolverOptions().setVerbosity(self._debug > self.DEBUG1)
        #=====================================================================
        #
        #   evaluate the derivatives of X, etc with respect to f:
        #
        ev=symb.Evaluator()
        names=[]
        if hasattr(self, "_r"):
             if symb.isSymbol(self._r):
                 names.append('r')
                 ev.addExpression(self._r.diff(f))
        for n in self._set_coeffs.keys():
            if n in self.__COEFFICIENTS and symb.isSymbol(self._set_coeffs[n]):
                   if n=="X" or n=="X_reduced":
                      T0=time()
                      B,A=symb.getTotalDifferential(self._set_coeffs[n], f, 1)
                      if n=='X_reduced':
                          self.trace3("Computing A_reduced, B_reduced took %f seconds."%(time()-T0))
                          names.append('A_reduced'); ev.addExpression(A)
                          names.append('B_reduced'); ev.addExpression(B)
                      else:
                          self.trace3("Computing A, B took %f seconds."%(time()-T0))
                          names.append('A'); ev.addExpression(A)
                          names.append('B'); ev.addExpression(B)
                   elif n=="Y" or n=="Y_reduced":
                      T0=time()
                      D,C=symb.getTotalDifferential(self._set_coeffs[n], f, 1)
                      if n=='Y_reduced':
                         self.trace3("Computing C_reduced, D_reduced took %f seconds."%(time()-T0))
                         names.append('C_reduced'); ev.addExpression(C)
                         names.append('D_reduced'); ev.addExpression(D)
                      else:
                         self.trace3("Computing C, D took %f seconds."%(time()-T0))
                         names.append('C'); ev.addExpression(C)
                         names.append('D'); ev.addExpression(D)
                   elif n in ("y", "y_reduced", "y_contact", "y_contact_reduced",  "y_dirac"):
                          names.append('d'+name[1:]); ev.addExpression(self._set_coeffs[name].diff(f))
                          relevant_symbols['d'+name[1:]]=self._set_coeffs[name].diff(f)
        res=ev.evaluate()
        if len(names)==1: res=[res]
        self.trace3("RHS expressions evaluated in %f seconds."%(time()-T0))
        for i in range(len(names)):
            self.trace3("util.Lsup(%s)=%s"%(names[i],util.Lsup(res[i])))
        coeffs_f=dict(zip(names,res))
        #

        # now we are ready to calculate the right hand side coefficients into
        # args by multiplication with g and grad(g).
        if len_g >1:
            if self.getNumSolutions() == 1:
                u_g=Data(0., (len_g,), self._lpde.getFunctionSpaceForSolution())
            else:
                u_g=Data(0., (self.getNumSolutions(), len_g), self._lpde.getFunctionSpaceForSolution())

        for i in range(len_g):
              # reset coefficients may be set at previous calls:
              args={}
              for n in self.__COEFFICIENTS: args[n]=Data()
              args['r']=Data()
              if g == None: # g_l=delta_{il} and len_f=len_g
                  for n,v in coeffs_f:
                      name=None
                      if len_f > 1:
                          val=v[:,i]
                      else:
                          val=v
                      if n.startswith("d"):
                           name='y'+n[1:]
                      elif n.startswith("D"):
                          name='Y'+n[1:]
                      elif n.startswith("r"):
                          name='r'+n[1:]
                      if name: args[name]=val
              else:
                    g_i=g[i]
                    for n,v in coeffs_f:
                      name=None
                      if n.startswith("d"):
                          name = 'y'+n[1:]
                          val = self.__mm(v, g_i)
                      elif n.startswith("r"):
                          name= 'r'
                          val = self.__mm(v, g_i)
                      elif n.startswith("D"):
                          name = 'Y'+n[1:]
                          val = self.__mm(v, g_i)
                      elif n.startswith("B") and isinstance(g_i, Data):
                          name = 'Y'+n[1:]
                          val = self.__mm(v, grad(g_i))
                      elif n.startswith("C"):
                          name = 'X'+n[1:]
                          val = matrix_multiply(v, g_i)
                      elif n.startswith("A") and isinstance(g_i, Data):
                          name = 'X'+n[1:]
                          val = self.__mm(v, grad(g_i))
                      if name:
                          if args.has_key(name):
                              args[name]+=val
                          else:
                              args[name]=val
              self._lpde.setValue(**args)
              u_g_i=self._lpde.getSolution()

              if len_g >1:
                 if self.getNumSolutions() == 1:
                     u_g[i]=-u_g_i
                 else:
                     u_g[:,i]=-u_g_i
              else:
                  u_g=-u_g_i

        return u_g

    def __mm(self,  m, v):
        """
        a sligtly crude matrix*matrix multiplication
        m is  A-coefficient, u is vector, v is grad vector:   A_ijkl*v_kl
        m is  B-coefficient, u is vector, v is vector:        B_ijk*v_k
        m is  C-coefficient, u is vector, v is grad vector:   C_ikl*v_kl
        m is  D-coefficient, u is vector, v is vector:        D_ij*v_j
        m is  A-coefficient, u is scalar, v is grad vector:   A_jkl*v_kl
        m is  B-coefficient, u is scalar, v is vector:        B_jk*v_k
        m is  C-coefficient, u is scalar, v is grad vector:   C_kl*v_kl
        m is  D-coefficient, u is scalar, v is vector:        D_j*v_j
        m is  A-coefficient, u is vector, v is grad scalar:   A_ijl*v_l
        m is  B-coefficient, u is vector, v is scalar:        B_ij*v
        m is  C-coefficient, u is vector, v is grad scalar:   C_il*v_l
        m is  D-coefficient, u is vector, v is scalar:        D_i*v
        m is  A-coefficient, u is scalar, v is grad scalar:   A_jl*v_l
        m is  B-coefficient, u is scalar, v is scalar:        B_j*v
        m is  C-coefficient, u is scalar, v is grad scalar:   C_l*v_l
        m is  D-coefficient, u is scalar, v is scalar:        D*v
        """
        s_m=getShape(m)
        s_v=getShape(v)

        # m is  B-coefficient, u is vector, v is scalar:        B_ij*v
        # m is  D-coefficient, u is vector, v is scalar:        D_i*v
        # m is  B-coefficient, u is scalar, v is scalar:        B_j*v
        # m is  D-coefficient, u is scalar, v is scalar:        D*v
        if s_m == () or s_v == ():
            return m*v

        # m is  D-coefficient, u is scalar, v is vector:        D_j*v_j
        # m is  C-coefficient, u is scalar, v is grad scalar:   C_l*v_l
        if len(s_m) == 1:
            return inner(m,v)

        # m is  D-coefficient, u is vector, v is vector:        D_ij*v_j
        # m is  B-coefficient, u is scalar, v is vector:        B_jk*v_k
        # m is  C-coefficient, u is vector, v is grad scalar:   C_il*v_l
        # m is  A-coefficient, u is scalar, v is grad scalar:   A_jl*v_l
        if len(s_m) == 2 and len(s_v) == 1:
            return matrix_mult(m,v)

        # m is  C-coefficient, u is scalar, v is grad vector:   C_kl*v_kl
        if len(s_m) == 2 and len(s_v) == 2:
            return inner(m,v)

        # m is  B-coefficient, u is vector, v is vector:        B_ijk*v_k
        # m is  A-coefficient, u is vector, v is grad scalar:   A_ijl*v_l
        if len(s_m) == 3 and len(s_v) == 1:
            return matrix_mult(m,v)

        # m is  A-coefficient, u is scalar, v is grad vector:   A_jkl*v_kl
        # m is  C-coefficient, u is vector, v is grad vector:   C_ikl*v_kl
        if len(s_m) == 3 and len(s_v) == 2:
            return tensor_mult(m,v)

        # m is  A-coefficient, u is vector, v is grad vector:   A_ijkl*v_kl
        if len(s_m) == 4 and len(s_v) == 2:
            return tensor_mult(m,v)

    def _updateRHS(self, expressions, subs):
        """
        """
        ev=symb.Evaluator()
        names=[]
        for name in expressions:
            if name in self.__COEFFICIENTS:
                ev.addExpression(expressions[name])
                names.append(name)
        if len(names)==0:
            return
        self.trace3("Starting expression evaluation.")
        T0=time()
        ev.subs(**subs)
        res=ev.evaluate()
        if len(names)==1: res=[res]
        self.trace3("RHS expressions evaluated in %f seconds."%(time()-T0))
        for i in range(len(names)):
            self.trace3("util.Lsup(%s)=%s"%(names[i],util.Lsup(res[i])))
        args=dict(zip(names,res))
        # reset coefficients may be set at previous calls:
        for n in self.__COEFFICIENTS:
            if not args.has_key(n): args[n]=Data()
        if not args.has_key('r'): args['r']=Data()
        self._lpde.setValue(**args)

    def _updateMatrix(self, expressions, subs):
        """
        """
        ev=symb.Evaluator()
        names=[]
        for name in expressions:
            if not name in self.__COEFFICIENTS:
                ev.addExpression(expressions[name])
                names.append(name)
        if len(names)==0:
            return
        self.trace3("Starting expression evaluation.")
        T0=time()
        ev.subs(**subs)
        res=ev.evaluate()
        if len(names)==1: res=[res]
        self.trace3("Matrix expressions evaluated in %f seconds."%(time()-T0))
        for i in range(len(names)):
            self.trace3("util.Lsup(%s)=%s"%(names[i],util.Lsup(res[i])))
        self._lpde.setValue(**dict(zip(names,res)))

    def _updateLinearPDE(self, expressions, subs):
        self._updateMatrix(expressions, subs)
        self._updateRHS(expressions, subs)


class VariationalProblem(object):
    """
    This class is used to define a general constraint vartional problem
    for an unknown function *u* and (spatially variable) parameter *p* on a
    given domain defined through a `Domain` object. *u* may be a scalar or a
    vector. *p* which may be a scalar or a vector may not be present.

    The solution *u* and the paremeter *p* are given as the solution of the
    minimization problem:

    *min_{u,p} int(H) + int(h)*

    where int{H} refers to integration over the domain and
    *H*=f(*x*,*u*,*grad(u)*,*p*, *grad(p)*) is a function which may depend on
    the location *x* within the domain and is a function of the solution *u*
    and the parameter *p* and their gradients *grad(u)* and *grad(p)*,
    respectively.
    Similarly, int{H} refers to integration over the boundary of the domain and
    *h=f(*x*,*u*, *p*) is a function which may depend on the location *x*
    within the domain boundary and is a function of the solution *u* and the
    parameter *p*.

    If *p* is present, *u* is the solution of a PDE with coefficients depending
    on the parameter *p*. The PDE defines a constraint for the variational
    problem. It is assumed that, if *p* is present, for any given parameter *p*
    a unique solution *u* of the PDE exists.

    For a single PDE having a solution with a single component the nonlinear
    PDE is defined in the following form:

    *-div(X) + Y = 0*

    where *X*,*Y*=f(*x*,*u*,*grad(u)*, *p*, *grad(p)*), *div(F)* denotes the
    divergence of *F* and *grad(F)* is the spatial derivative of *F*.

    The coefficients *X* (rank 1) and *Y* (scalar) have to be specified through
    `Symbol` objects.

    The following natural boundary conditions are considered:

    *n[j]*X[j] - y = 0*

    where *n* is the outer normal field. Notice that the coefficient *X*
    is defined in the PDE. The coefficient *y* is a scalar `Symbol`.

    Constraints for the solution prescribe the value of the solution at certain
    locations in the domain. They have the form

    *u=r* where *q>0*

    *r* and *q* are each scalar where *q* is the characteristic function
    defining where the constraint is applied. The constraints override any
    other condition set by the PDE or the boundary condition.

    For a system of PDEs and a solution with several components, *u* is rank
    one, while the PDE coefficient *X* is rank two and *Y* and *y* is rank one.

    The PDE is solved by linearising the coefficients and iteratively solving
    the corresponding linear PDE until the error is smaller than a tolerance
    or a maximum number of iterations is reached.

    Typical usage:

        s=Symbol('s', dim=dom.getDim())
        T = Symbol('T', dim=dom.getDim())
        log_k = Symbol('log_k', dim=dom.getDim())
        v = VariationalProblem(dom, u=T, p=log_k)
        v.setValue(X=exp(-log_k)*grad(T), Y=s, h=0.3*(T-0.3)**2)
        T, log_k, l = v.getSolution(T=0., log_k=1, s=2.45)
        sT,S_log_k=v.getSensitivity(s, direction=1, T=T, log_k=log_k, s=2.45)

    """
    DEBUG0=0 # no debug info
    DEBUG1=1 # info on Newton-Raphson solver printed
    DEBUG2=2 # in addition info from linear solver is printed
    DEBUG3=3 # in addition info on linearization is printed
    DEBUG4=4 # in addition info on LinearPDE handling is printed
    def __init__(self, domain, u, p=None, debug=DEBUG0):
        """
        Initializes a new variational problem.

        :param domain: domain of the PDE
        :type domain: `Domain`
        :param u: The symbol for the unknown PDE function u.
        :type u: `Symbol`
        :param p: The symbol for the parameter p.
        :type p: `Symbol`
        :param debug: level of debug information to be printed
        """
        self.__COEFFICIENTS = [ "X", "X_reduced", "Y", "Y_reduced", "y", "y_reduced", "y_contact", "y_contact_reduced", "y_dirac", \
                                "H", "h_reduced", "h", "h_reduced", "h_contact", "h_contact_reduced", "h_dirac" ]

        self._set_coeffs={}
        self._unknown=u
        self._parameter=p
        self._debug=debug
        self.dim = domain.getDim()

        if self._parameter == None:
            U=u
        else:
            self._lagrangean=symb.Symbol("lambda%s"%id(self), self._unknown.getShape(), dim=self.dim)
            U=concatenateRow(self._parameter, self._unknown, self._lagrangean)
        self.__PDE=NonlinearPDE(domain, u=U, debug=debug)

    def __str__(self):
        """
        Returns the string representation of the problem.

        :return: a simple representation of the variational problem
        :rtype: ``str``
        """
        return "<%s %d>"%(self.__class__.__name__, id(self))

    def trace1(self, text):
        """
        Prints the text message if the debug level is greater than DEBUG0

        :param text: message to be printed
        :type text: ``string``
        """
        if self._debug > self.DEBUG0:
            print("%s: %s"%(str(self), text))

    def trace3(self, text):
        """
        Prints the text message if the debug level is greater than DEBUG3

        :param text: message to be printed
        :type text: ``string``
        """
        if self._debug > self.DEBUG2:
            print("%s: %s"%(str(self), text))

    def getNonlinearPDE(self):
        """
        Returns the `NonlinearPDE` used to solve the variational problem

        :return: underlying nonlinear PDE
        :rtype: `NonlinearPDE`
        """
        return self.__PDE

    def getNumSolutions(self):
        """
        Returns the number of solution components
        :rtype: ``int``
        """
        s=self._unknown.getShape()
        if len(s) > 0:
            return s[0]
        else:
            return 1

    def getNumParameters(self):
        """
        Returns the number of parameter components. If no parameter is present
        zero is returned.

        :rtype: ``int``
        """
        if self._parameter == None:
            return 0
        else:
            s=self._parameter.getShape()
            if len(s) > 0:
                return s[0]
            else:
                return 1

    def getShapeOfCoefficient(self, name):
        """
        Returns the shape of the coefficient ``name``

        :param name: name of the coefficient enquired
        :type name: ``string``
        :return: the shape of the coefficient ``name``
        :rtype: ``tuple`` of ``int``
        :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
        """
        numSol=self.getNumSolutions()
        numParams=self.getNumParameters()
        dim = self.dim
        if (name=="X" or name=="X_reduced") and numParams>0:
            if numSol > 1:
                return (numSol,dim)
            else:
                return (dim,)
        elif name=="r" or name == "q" :
            if numSol > 1:
                return (numSol,)
            else:
                return ()
        elif ( name=="rp" or name == "qp") and numParams>0:
            if numParams > 1:
                return (numParams,)
            else:
                return ()
        elif ( name=="Y" or name=="Y_reduced") and numParams>0:
            if numSol > 1:
                return (numSol,)
            else:
                return ()
        elif ( name=="y" or name=="y_reduced") and numParams>0:
            if numSol > 1:
                return (numSol,)
            else:
                return ()
        elif ( name=="y_contact" or name=="y_contact_reduced") and numParams>0:
            if numSol > 1:
                return (numSol,)
            else:
                return ()
        elif name=="y_dirac" and numParams>0:
            if numSol > 1:
                return (numSol,)
            else:
                return ()
        elif name=="H" or name=="H_reduced":
                return ()
        elif name=="h" or name=="h_reduced":
                return ()
        elif name=="h_contact" or name=="h_contact_reduced":
                return ()
        elif name=="h_dirac":
                return ()
        else:
            raise lpe.IllegalCoefficient("Attempt to request unknown coefficient %s"%name)

    def createCoefficient(self, name):
        """
        Creates a new coefficient ``name`` as Symbol

        :param name: name of the coefficient requested
        :type name: ``string``
        :return: the value of the coefficient
        :rtype: `Symbol` or `Data`
        :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
        """
        numSol=self.getNumSolutions()
        numParams=self.getNumParameters()
        if name == "q":
            if numParams > 0:
                return self.getNonlinearPDE().createCoefficient("q")[numParams:numParams+numSol]
            else:
                return self.getNonlinearPDE().createCoefficient("q")
        elif name == "qp":
            if numParams > 0:
                return self.getNonlinearPDE().createCoefficient("q")[:numParams]
            else:
                raise lpe.IllegalCoefficient("Attempt to request coefficient %s"%name)
        else:
           s=self.getShapeOfCoefficient(name)
           return symb.Symbol(name, s, dim=self.dim)

    def getCoefficient(self, name):
        """
        Returns the value of the coefficient ``name`` as Symbol

        :param name: name of the coefficient requested
        :type name: ``string``
        :return: the value of the coefficient
        :rtype: `Symbol`
        :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
        """
        if self._set_coeffs.has_key(name):
            return self._set_coeffs[name]
        else:
            raise lpe.IllegalCoefficient("Attempt to request undefined coefficient %s"%name)

    def __getNonlinearPDECoefficient(self, extension, capson=False, order=0):
        """
        """
        if capson:
            H_key="H"+extension
            X_key="X"+extension
            Y_key="Y"+extension
            Z_key="Z"+extension
        else:
            H_key="h"+extension
            X_key="x"+extension
            Y_key="y"+extension
            Z_key="z"+extension

        Z=0
        if self._set_coeffs.has_key(H_key): Z+=self._set_coeffs[H_key]
        if self._set_coeffs.has_key(X_key): Z+=util.inner(self._set_coeffs[X_key], util.grad(self._lagrangean))
        if self._set_coeffs.has_key(Y_key): Z+=util.inner(self._set_coeffs[Y_key], self._lagrangean)

        if self.getNumParameters() > 0:
            if order == 0:
                Yp=symb.getTotalDifferential(Z, self._parameter, order=0)
                Yu=symb.getTotalDifferential(Z, self._unknown, order=0)
                Yl=symb.getTotalDifferential(Z, self._lagrangean, order=0)
                Y=concatenateRow(Yp, Yl, Yu)  # order different from solution!
            else:
                Yp,Xp=symb.getTotalDifferential(Z, self._parameter, order=1)
                Yu,Xu=symb.getTotalDifferential(Z, self._unknown, order=1)
                Yl,Xl=symb.getTotalDifferential(Z, self._lagrangean, order=1)
                Y=concatenateRow(Yp, Yl, Yu)  # order different from solution!
                X=concatenateRow(Xp, Xl, Xu)  # order different from solution!
        else:
            if order == 0:
                Y=symb.getTotalDifferential(Z, self._unknown, order=0)
            else:
                Y,X=symb.getTotalDifferential(Z, self._unknown, order=1)
        if order == 0:
            return Y
        else:
            return Y,X

    def setValue(self,**coefficients):
        """
        Sets new values to one or more coefficients.

        :keyword H: value for coefficient ``H``
        :type H: `Symbol`
        :keyword h: value for coefficient ``h``
        :type h: `Symbol`
        :keyword h_contact: value for coefficient ``h_contact``
        :type h_contact: `Symbol`
        :keyword h_dirac: value for coefficient ``y_dirac``
        :type h_dirac: `Symbol`
        :keyword X: value for coefficient ``X``
        :type X: `Symbol` or any type that can be cast to a `Data` object
        :keyword Y: value for coefficient ``Y`` = r
        :type Y: `Symbol` or any type that can be cast to a `Data` object
        :keyword y: value for coefficient ``y``
        :type y: `Symbol` or any type that can be cast to a `Data` object
        :keyword y_contact: value for coefficient ``y_contact``
        :type y_contact: `Symbol` or any type that can be cast to a `Data` object
        :keyword y_dirac: value for coefficient ``y_dirac``
        :type y_dirac: `Symbol` or any type that can be cast to a `Data` object
        :keyword q: mask for location of constraint
        :type q: any type that can be casted to a `Data` object
        :keyword r: value of solution prescribed by constraint
        :type r: `Symbol` or any type that can be cast to a `Data` object
        :keyword qp: mask for location of constraint fro parameter
        :type qp: any type that can be casted to a `Data` object
        :keyword rp: value of the parameter prescribed by parameter constraint
        :type rp: `Symbol` or any type that can be cast to a `Data` object

        :raise IllegalCoefficient: if an unknown coefficient keyword is used
        :raise IllegalCoefficientValue: if a supplied coefficient value has an
                                        invalid shape
        """
        numSol=self.getNumSolutions()
        numParams=self.getNumParameters()
        update=[]
        for name,val in coefficients.iteritems():
            shape=util.getShape(val)
            if not shape == self.getShapeOfCoefficient(name):
                raise lpe.IllegalCoefficientValue("%s has shape %s but must have shape %s"%(name, shape, self.getShapeOfCoefficient(name)))
            if name == "q":
                self._q = val
                update.append("q")

            elif name == "qp":
                if numParams<1:
                    raise lpe.IllegalCoefficientValue("Illegal coefficient %s - no parameter present."%name)
                self._qp = val
                update.append("q")

            elif name == "r":
                self._r = val
                update.append("r")

            elif name == "rp":
                if numParams<1:
                    raise lpe.IllegalCoefficientValue("Illegal coefficient %s - no parameter present."%name)
                self._rp = val
                update.append("r")

            elif name=="X":
                if numParams<1:
                    raise lpe.IllegalCoefficientValue("Illegal coefficient %s - no parameter present."%name)
                self._set_coeffs['X']=val
                update.append("Y")

            elif name=="X_reduced":
                if numParams<1:
                    raise lpe.IllegalCoefficientValue("Illegal coefficient %s - no parameter present."%name)
                self._set_coeffs['X_reduced']=val
                update.append("Y_reduced")

            elif name=="Y":
                if numParams<1:
                    raise lpe.IllegalCoefficientValue("Illegal coefficient %s - no parameter present."%name)
                self._set_coeffs['Y']=val
                update.append("Y")

            elif name=="Y_reduced":
                if numParams<1:
                    raise lpe.IllegalCoefficientValue("Illegal coefficient %s - no parameter present."%name)
                self._set_coeffs['Y_reduced']=val
                update.append("Y_reduced")

            elif name=="H":
                self._set_coeffs['H']=val
                update.append("Y")

            elif name=="H_reduced":
                self._set_coeffs['H_reduced']=val
                update.append("Y_reduced")

            elif name=="y":
                if numParams<1:
                    raise lpe.IllegalCoefficientValue("Illegal coefficient %s - no parameter present."%name)
                self._set_coeffs['y']=val
                update.append("y")

            elif name=="y_reduced":
                if numParams<1:
                    raise lpe.IllegalCoefficientValue("Illegal coefficient %s - no parameter present."%name)
                self._set_coeffs['y_reduced']=val
                update.append("y_reduced")

            elif name=="h":
                self._set_coeffs['h']=val
                update.append("y")

            elif name=="h_reduced":
                self._set_coeffs['h_reduced']=val
                update.append("y_reduced")

            elif name=="y_contact":
                if numParams<1:
                    raise lpe.IllegalCoefficientValue("Illegal coefficient %s - no parameter present."%name)
                self._set_coeffs['y_contact']=val
                update.append("y_contact")

            elif name=="y_contact_reduced":
                if numParams<1:
                    raise lpe.IllegalCoefficientValue("Illegal coefficient %s - no parameter present."%name)
                self._set_coeffs['y_contact_reduced']=val
                update.append("y_contact_reduced")

            elif name=="h_contact":
                self._set_coeffs['h_contact']=val
                update.append("y_contact")

            elif name=="h_contact_reduced":
                self._set_coeffs['h_contact_reduced']=val
                update.append("y_contact_reduced")

            elif name=="y_dirac":
                if numParams<1:
                    raise lpe.IllegalCoefficientValue("Illegal coefficient %s - no parameter present."%name)
                self._set_coeffs['y_dirac']=val
                update.append("y_dirac")

            elif name=="h_dirac":
                self._set_coeffs['h_diract']=val
                update.append("y_dirac")
            else:
                raise lpe.IllegalCoefficient("Attempt to set unknown coefficient %s"%name)

        # now we can update the coefficients of the nonlinear PDE:
        coeff2={}
        if "q" in update:
            if numParams>0:
                q=self.getNonlinearPDE().createCoefficient("q")
                if hasattr(self, "_qp"): q[:numParams]=self._qp
                if hasattr(self, "_q"):
                    q[numParams:numParams+numSol]=self._q
                    q[numParams+numSol:]=self._q
            else:
                q=self._q
            coeff2["q"]=q
        if "r" in update:
            if numParams>0:
                if hasattr(self, "_rp"):
                    rp=self._rp
                else:
                    rp=numpy.zeros(self.getShapeOfCoefficient('rp'))
                if hasattr(self, "_r"):
                    r=self._r
                else:
                    r=numpy.zeros(self.getShapeOfCoefficient('r'))
                coeff2["r"]=concatenateRow(rp, r, numpy.zeros(self.getShapeOfCoefficient('r')) )
            else:
                coeff2["r"]=self._r

        if "Y" in update:
            Y,X = self.__getNonlinearPDECoefficient("", capson=True, order=1)
            coeff2["Y"]=Y
            coeff2["X"]=X
        if "y" in update:
            coeff2["y"] = self.__getNonlinearPDECoefficient("", capson=False)
        if "y_contact" in update:
            coeff2["y_contact"] = self.__getNonlinearPDECoefficient("_contact", capson=False)
        if "y_dirac" in update:
            coeff2["y_dirac"] = self.__getNonlinearPDECoefficient("_dirac", capson=False)
        if "Y_reduced" in update:
            Y,X = self.__getNonlinearPDECoefficient("_reduced",capson=True, order=1)
            coeff2["Y_reduced"]=Y
            coeff2["X_reduced"]=X
        if "y_reduced" in update:
            coeff2["y_reduced"]= self.__getNonlinearPDECoefficient("_reduced",capson=False)
        if "y_contact_reduced" in update:
            coeff2["y_contact_reduced"] = self.__getNonlinearPDECoefficient("_contact_reduced",capson=False)

        # and now we can set the PDE coefficients:
        self.getNonlinearPDE().setValue(**coeff2)

    def getSolution(self, **subs):
        """
        Returns the solution of the variational problem.

        :param subs: Substitutions for all symbols used in the coefficients
                     including the initial value for solution *u* and for the
                     parameter *p* (if present)
        :return: parameter, corresponding solution and lagrangean multiplier
        :rtype: tuple of `Data` or single `Data` (if no parameter present)
        """
        numSol=self.getNumSolutions()
        numParams=self.getNumParameters()

        if numParams>0:
            fs=self.getNonlinearPDE().getLinearPDE().getFunctionSpaceForSolution()
            subs[self._lagrangean.name]=Data(0., self._lagrangean.getShape(), fs)
            p_sym=self._parameter.atoms().pop().name
            if not subs.has_key(p_sym):
                raise KeyError("Initial value for '%s' missing."%p_sym)
            pi=subs[p_sym]
            if not isinstance(pi,Data):
                fs=self.getNonlinearPDE().getLinearPDE().getFunctionSpaceForSolution()
                pi=Data(pi, fs)
                subs[p_sym]=pi

        Ui=self.getNonlinearPDE().getSolution(**subs)

        if numParams > 0:
            # return parameter, solution, lagrangean multiplier
            if numParams == 1:
                 p=Ui[0]
            else:
                 p=Ui[:numParams]
            if numSol == 1:
                 u=Ui[numParams]
                 l=Ui[numParams+1]
            else:
                 u=Ui[numParams:numParams+numSol]
                 l= Ui[numParams+numSol:]
            return p,u,l
        else:
            return Ui

