# -*- coding: utf-8 -*-

########################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
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
from time import time
from esys.escript.linearPDEs import LinearPDE, IllegalCoefficient, IllegalCoefficientValue
from esys.escript import util, Data, Symbol, Evaluator

__author__="Cihan Altinay"

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

    *u=u0* where *q>0*

    *u0* and *q* are each scalar where *q* is the characteristic function
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

    def __init__(self,domain,u,debug=False):
        """
        Initializes a new nonlinear PDE.

        :param domain: domain of the PDE
        :type domain: `Domain`
        :param u: The symbol for the unknown PDE function u.
        :type u: `Symbol`
        :param debug: if True debug information is printed
        """
        self.__COEFFICIENTS = [ "X", "X_reduced", "Y", "Y_reduced", "y", "y_reduced", "y_contact", "y_contact_reduced", "y_dirac", "y_dirac_reduced", "q" ]
        self._coeffs={}
        self._u=u
        self._debug=debug
        self._tol=1e-6
        self._iterMax=100
        self.dim = domain.getDim()
        if u.getRank()==0:
            numEquations=1
        else:
            numEquations=u.getShape()[0]
        numSolutions=numEquations
        self._lpde=LinearPDE(domain,numEquations,numSolutions,debug)

    def __str__(self):
        """
        Returns the string representation of the PDE.

        :return: a simple representation of the PDE
        :rtype: ``str``
        """
        return "<NonlinearPDE %d>"%id(self)
        
    def getUnknownSymbol(self):
        """
        Returns the symbol of the PDE unknown

        :return: the symbol of the PDE unknown
        :rtype: `Symbol`
        """
        return self._u
      
    def trace(self, text):
        """
        Prints the text message if debug mode is switched on.

        :param text: message to be printed
        :type text: ``string``
        """
        if self._debug:
            print("%s: %s"%(str(self), text))

    def getLinearSolverOptions(self):
        """
        Returns the options of the linear PDE solver class.
        """
        return self._lpde.getSolverOptions()

    def setOptions(self, **opts):
        """
        Allows setting options for the nonlinear PDE.
        The supported options are:
          'tolerance' - error tolerance for the Newton method
          'iterMax'   - maximum number of Newton iterations
        """
        for key in opts:
            if key=='tolerance':
                self._tol=opts[key]
            elif key=='iterMax':
                self._iterMax=opts[key]
            else:
                raise KeyError("Invalid option '%s'"%key)

    def getSolution(self, **subs):
        """
        Returns the solution of the PDE.

        :param subs: Substitutions for all symbols used in the coefficients
                     including the initial value for *u*.
        :return: the solution
        :rtype: `Data`
        """

        u_sym=self._u.atoms().pop().name
        if not subs.has_key(u_sym):
            raise KeyError("Initial value for '%s' missing."%u_sym)

        ui=subs.pop(u_sym)
        if isinstance(ui,float) or isinstance(ui,int):
            ui=Data(ui, self._lpde.getFunctionSpaceForSolution())

        # separate symbolic expressions from other coefficients
        constants={}
        expressions={}
        for n, e in self._coeffs.items():
            if hasattr(e, "atoms"):
                expressions[n]=e
            else:
                constants[n]=e

        # set constant PDE values now
        self._lpde.setValue(**constants)

        # perform Newton iterations until error is small enough or
        # maximum number of iterations reached
        n=0
        while True:
            n=n+1
            subs[u_sym]=ui
            delta_u=self._newtonStep(expressions, subs)
            ui=ui-delta_u
            self.trace("Error after %d Newton iteration(s): %g"%(n,util.Lsup(delta_u)))
            if util.Lsup(delta_u) < self._tol or n >= self._iterMax:
                break

        self.trace("Final error after %d Newton iteration(s): %g"%(n,util.Lsup(delta_u)))
        return ui

    def getNumSolutions(self):
        """
        Returns the number of the solution components
        :rtype: ``int``
        """
        s=self._u.getShape()
        if len(s) > 0:
          return s[0]
        else:
          return 1

    def getShapeOfCoefficient(self,name):
        """
        Returns the shape of the coefficient ``name``.

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
        elif name=="q" :
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
        else:
            raise IllegalCoefficient("Attempt to request unknown coefficient %s"%name) 

    def createCoefficient(self, name):
        """
        create a new coefficient ``name`` as Symbol

        :param name: name of the coefficient requested
        :type name: ``string``
        :return: the value of the coefficient
        :rtype: `Symbol`
        :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
        """
        s=self.getShapeOfCoefficient(name)
        return Symbol(name, s)

    def getCoefficient(self, name):
        """
        Returns the value of the coefficient ``name`` as Symbol

        :param name: name of the coefficient requested
        :type name: ``string``
        :return: the value of the coefficient
        :rtype: `Symbol`
        :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
        """
        if self._coeffs.has_key(name):
            return self._coeffs[name]
        else:
            raise IllegalCoefficient("Attempt to request undefined coefficient %s"%name)

    def setValue(self,**coefficients):
        """
        Sets new values to one or more coefficients.

        :param coefficients: new values assigned to coefficients
        :keyword X: value for coefficient ``X``
        :type X: `Symbol` or any type that can be cast to a `Data` object
        :keyword Y: value for coefficient ``Y``
        :type Y: `Symbol` or any type that can be cast to a `Data` object
        :keyword y: value for coefficient ``y``
        :type y: `Symbol` or any type that can be cast to a `Data` object
        :keyword y_contact: value for coefficient ``y_contact``
        :type y_contact: `Symbol` or any type that can be cast to a `Data`
                         object
        :keyword q: mask for location of constraints
        :type q: `Symbol` or any type that can be cast to a `Data` object

        :raise IllegalCoefficient: if an unknown coefficient keyword is used
        :raise IllegalCoefficientValue: if a supplied coefficient value has an
                                        invalid shape
        """

        u=self._u
        for name,val in coefficients.iteritems():
            shape=util.getShape(val)
            if not shape == self.getShapeOfCoefficient(name):
                raise IllegalCoefficientValue("%s has shape %s but must have shape %d"%(name, self.getShapeOfCoefficient(name), shape))
            rank=len(shape)
            if name=="X" or name=="X_reduced":
                # DX/Du = del_X/del_u + del_X/del_grad(u)*del_grad(u)/del_u
                #            \   /         \   /
                #              B             A
                
                if rank != u.getRank()+1:
                    raise IllegalCoefficientValue("%s must have rank %d"%(name,u.getRank()+1))
                T0=time()
                if not hasattr(val, 'diff'):
                    A=numpy.zeros(shape+u.grad().getShape())
                    B=numpy.zeros(shape+u.getShape())
                elif u.getRank()==0:
                    X=self._removeFsFromGrad(val)
                    dXdu=X.diff(u)
                    dgdu=u.grad().diff(u)
                    A=numpy.empty(shape+dgdu.getShape(), dtype=object)
                    B=numpy.empty(shape, dtype=object)
                    for i in numpy.ndindex(shape):
                        for j in numpy.ndindex(dgdu.getShape()):
                            y=dXdu[i]
                            x=dgdu[j]
                            # expand() and coeff() are very expensive so
                            # we set the unwanted factors to zero to extract
                            # the one we need
                            for jj in numpy.ndindex(dgdu.getShape()):
                                if j==jj: continue
                                y=y.subs(dgdu[jj], 0)
                            B[i]=y.subs(x,0) # terms in u and constants
                            A[i,j]=y.subs(x,1)-B[i]
                    A=Symbol(A)
                    B=Symbol(B)
                    self.trace("Computing A, B took %f seconds."%(time()-ttt0))
                else:  #u.getRank()==1
                    X=self._removeFsFromGrad(val)
                    dXdu=X.diff(u)
                    dgdu=u.grad().diff(u).transpose(2)
                    I,J=shape
                    K,L=u.grad().getShape()
                    A=numpy.empty((I,J,K,L), dtype=object)
                    B=numpy.empty((I,J,K), dtype=object)
                    for i,j,k,l in numpy.ndindex(I,J,K,L):
                        if dgdu[k,k,l]==0:
                            A[i,j,k,l]=0
                            B[i,j,k]=0
                        else:
                            y=dXdu[i,j,k]
                            x=dgdu[k,k,l]
                            for kk,ll in numpy.ndindex(K,L):
                                if k==kk and l==ll: continue
                                y=y.subs(dgdu[kk,kk,ll], 0)
                            B[i,j,k]=y.subs(x,0) # terms in u and constants
                            A[i,j,k,l]=y.subs(x,1)-B[i,j,k]

                    A=Symbol(A)
                    B=Symbol(B)

                if name=='X_reduced':
                    self.trace("Computing A_reduced, B_reduced took %f seconds."%(time()-T0))
                    self._coeffs['A_reduced']=A
                    self._coeffs['B_reduced']=B
                    self._coeffs['X_reduced']=val
                else:
                    self.trace("Computing A, B took %f seconds."%(time()-T0))
                    self._coeffs['A']=A
                    self._coeffs['B']=B
                    self._coeffs['X']=val
            elif name=="Y" or name=="Y_reduced":
                # DY/Du = del_Y/del_u + del_Y/del_grad(u)*del_grad(u)/del_u
                #            \   /         \   /
                #              D             C
                if rank != u.getRank():
                    raise IllegalCoefficientValue("%s must have rank %d"%(name,u.getRank()))
                T0=time()
                if not hasattr(val, 'diff'):
                    C=numpy.zeros(shape+u.grad().getShape())
                    D=numpy.zeros(shape+u.getShape())
                elif u.getRank()==0:
                    Y=self._removeFsFromGrad(val)
                    dYdu=Y.diff(u)
                    dgdu=u.grad().diff(u)
                    C=numpy.empty(dgdu.getShape(), dtype=object)
                    for j in numpy.ndindex(dgdu.getShape()):
                        y=dYdu
                        x=dgdu[j]
                        # expand() and coeff() are very expensive so
                        # we set the unwanted factors to zero to extract
                        # the one we need
                        for jj in numpy.ndindex(dgdu.getShape()):
                            if j==jj: continue
                            y=y.subs(dgdu[jj], 0)
                        D=y.subs(x,0) # terms in u and constants
                        C[j]=y.subs(x,1)-D
                    C=Symbol(C)
                else:  #u.getRank()==1
                    Y=self._removeFsFromGrad(val)
                    dYdu=Y.diff(u)
                    dgdu=u.grad().diff(u).transpose(2)
                    I,=shape
                    J,K=u.grad().getShape()
                    C=numpy.empty((I,J,K), dtype=object)
                    D=numpy.empty((I,J), dtype=object)
                    for i,j,k in numpy.ndindex(I,J,K):
                        if dgdu[j,j,k]==0:
                            C[i,j,k]=0
                            D[i,j]=0
                        else:
                            y=dYdu[i,j]
                            x=dgdu[j,j,k]
                            for jj,kk in numpy.ndindex(J,K):
                                if j==jj and k==kk: continue
                                y=y.subs(dgdu[jj,jj,kk], 0)
                            D[i,j]=y.subs(x,0) # terms in u and constants
                            C[i,j,k]=y.subs(x,1)-D[i,j]
                    C=Symbol(C)
                    D=Symbol(D)

                if name=='Y_reduced':
                    self.trace("Computing C_reduced, D_reduced took %f seconds."%(time()-T0))
                    self._coeffs['C_reduced']=C
                    self._coeffs['D_reduced']=D
                    self._coeffs['Y_reduced']=val
                else:
                    self.trace("Computing C, D took %f seconds."%(time()-T0))
                    self._coeffs['C']=C
                    self._coeffs['D']=D
                    self._coeffs['Y']=val
            elif name=="y" or name=="y_reduced":
                y=val
                if rank != u.getRank():
                    raise IllegalCoefficientValue("%s must have rank %d"%(name,u.getRank()))
                if not hasattr(y, 'diff'):
                    d=numpy.zeros(u.getShape())
                else:
                    d=y.diff(u)
                if name=='y_reduced':
                    self._coeffs['d_reduced']=d
                    self._coeffs['y_reduced']=y
                else:
                    self._coeffs['d']=d
                    self._coeffs['y']=y
            elif name=="y_contact" or name=="y_contact_reduced":
                y_contact=val
                if rank != u.getRank():
                    raise IllegalCoefficientValue("%s must have rank %d"%(name,u.getRank()))
                if not hasattr(y_contact, 'diff'):
                    d_contact=numpy.zeros(u.getShape())
                else:
                    d_contact=y_contact.diff(u)
                if name=='y_contact_reduced':
                    self._coeffs['d_contact_reduced']=d_contact
                    self._coeffs['y_contact_reduced']=y_contact
                else:
                    self._coeffs['d_contact']=d_contact
                    self._coeffs['y_contact']=y_contact
            elif name=="y_dirac" or name=="y_dirac_reduced":
                y=val
                if rank != u.getRank():
                    raise IllegalCoefficientValue("%s must have rank %d"%(name,u.getRank()))
                if not hasattr(y, 'diff'):
                    d=numpy.zeros(u.getShape())
                else:
                    d=y.diff(u)
                if name=='y_dirac_reduced':
                    self._coeffs['d_dirac_reduced']=d
                    self._coeffs['y_dirac_reduced']=y
                else:
                    self._coeffs['d_dirac']=d
                    self._coeffs['y_dirac']=y
            elif name=="q":
                if rank != u.getRank():
                    raise IllegalCoefficientValue("q must have rank %d"%u.getRank())
                self._coeffs['q']=val
            else:
                raise IllegalCoefficient("Attempt to set unknown coefficient %s"%name)

    def _newtonStep(self, expressions, subs):
        """
        """
        self._updateLinearPDE(expressions, subs)
        delta_u=self._lpde.getSolution()
        return delta_u

    def _updateRHS(self, expressions, subs):
        ev=Evaluator()
        names=[]
        for name in expressions:
            if name in self.__COEFFICIENTS:
                ev.addExpression(expressions[name])
                names.append(name)
        if len(names)==0:
            return
        self.trace("Starting expression evaluation.")
        ttt0=time()
        ev.subs(**subs)
        res=ev.evaluate()
        if len(names)==1: res=[res]
        self.trace("RHS expressions evaluated in %f seconds."%(time()-ttt0))
        for i in range(len(names)):
            self.trace("Lsup(%s)=%s"%(names[i],util.Lsup(res[i])))
        self._lpde.setValue(**dict(zip(names,res)))

    def _updateMatrix(self, expressions, subs):
        ev=Evaluator()
        names=[]
        for name in expressions:
            if not name in self.__COEFFICIENTS:
                ev.addExpression(expressions[name])
                names.append(name)
        if len(names)==0:
            return
        self.trace("Starting expression evaluation.")
        ttt0=time()
        ev.subs(**subs)
        res=ev.evaluate()
        if len(names)==1: res=[res]
        self.trace("Matrix expressions evaluated in %f seconds."%(time()-ttt0))
        for i in range(len(names)):
            self.trace("Lsup(%s)=%s"%(names[i],util.Lsup(res[i])))
        self._lpde.setValue(**dict(zip(names,res)))

    def _updateLinearPDE(self, expressions, subs):
        self._updateMatrix(expressions, subs)
        self._updateRHS(expressions, subs)

    def _removeFsFromGrad(self, sym):
        """
        Returns sym with all occurrences grad_n(a,b,c) replaced by grad_n(a,b).
        That is, all functionspace parameters are removed.
        This method is used in setValue() in order to get substitutions work
        properly.
        """
        from esys.escript import symfn
        gg=sym.atoms(symfn.grad_n)
        for g in gg:
            if len(g.args)==3:
                r=symfn.grad_n(*g.args[:2])
                sym=sym.subs(g, r)
        return sym

