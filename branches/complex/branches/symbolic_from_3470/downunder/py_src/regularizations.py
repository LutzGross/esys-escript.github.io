
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

__all__ = ['Regularization']

import numpy as np
from esys.escript.linearPDEs import LinearSinglePDE
from esys.escript import Data, grad, inner, integrate, kronecker, vol

class Regularization(object):
    """
    """
    def __init__(self, domain, m_ref=0, w0=None, w=None, location_of_set_m=Data(), tol=1e-8):
        self.__domain=domain
        self.__m_ref=m_ref
        self.location_of_set_m=location_of_set_m
        if w0 is None:
            self._w0 = None
        else:
            self._w0=np.asarray(w0) / vol(self.__domain)
        if  w is None:
            self._w = None
        else:
            self._w=np.asarray(w)
            

        self.__projector=LinearSinglePDE(domain)
        self.__projector.getSolverOptions().setTolerance(tol)
        self.__projector.setValue(A=kronecker(domain), q=location_of_set_m, D=0.)

    def getInner(self, f0, f1):
        """
        returns the inner product of two gradients.
        """
        return integrate(inner(grad(f0),grad(f1)))

    def project(self, Y=Data(), X=Data()):
        self.__projector.setValue(Y=Y, X=X)
        return  self.__projector.getSolution()

    def getValue(self, m):
        A=0
        if self._w0 is not None:
            A=(m-self.__m_ref)**2 * self._w0
        if self._w is not None:
            A=inner(self._w, grad(m-self.__m_ref)**2) + A
        return integrate(A)

    def getGradient(self, m):
        if not self._w0 == None:
            Y=2. * (m-self.__m_ref) * self._w0
        else:
            Y=0.
        if not self._w == None:
            X=2. * grad(m-self.__m_ref) * self._w
        else:
            X=np.zeros((self.__domain.getDim(),))

        return Y, X

    def getArguments(self, m):
        return ()

