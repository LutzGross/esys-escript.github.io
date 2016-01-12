##############################################################################
#
# Copyright (c) 2003-2016 by The University of Queensland
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

"""Base classes for forward models"""

from __future__ import division, print_function

__copyright__="""Copyright (c) 2003-2016 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['ForwardModel','ForwardModelWithPotential']

from esys.downunder.coordinates import makeTransformation
from esys.escript.linearPDEs import LinearSinglePDE
from esys.escript.util import *
import numpy as np

class ForwardModel(object):
    """
    An abstract forward model that can be plugged into a cost function.
    Subclasses need to implement `getDefect()`, `getGradient()`, and possibly
    `getArguments()` and 'getCoordinateTransformation'.
    """
    def __init__(self):
        pass

    def getArguments(self, x):
        return ()

    def getCoordinateTransformation(self):
        return None

    def getDefect(self, x, *args):
        raise NotImplementedError

    def getGradient(self, x, *args):
        raise NotImplementedError


class ForwardModelWithPotential(ForwardModel):
    """
    Base class for a forward model using a potential such as magnetic or
    gravity. It defines a cost function:

        defect = 1/2 sum_s integrate( ( weight_i[s] * ( r_i - data_i[s] ) )**2 )

    where s runs over the survey, weight_i are weighting factors, data_i are
    the data, and r_i are the results produced by the forward model.
    It is assumed that the forward model is produced through postprocessing
    of the solution of a potential PDE.
    """
    def __init__(self, domain, w, data,  coordinates=None,
                                 fixPotentialAtBottom=False,
                                 tol=1e-8):
        """
        initializes a new forward model with potential.

        :param domain: domain of the model
        :type domain: `Domain`
        :param w: data weighting factors
        :type w: ``Vector`` or list of ``Vector``
        :param data: data
        :type data: ``Vector`` or list of ``Vector``
        :param coordinates: defines coordinate system to be used
        :type coordinates: `ReferenceSystem` or `SpatialCoordinateTransformation`
        :param fixPotentialAtBottom: if true potential is fixed to zero at the bottom of the domain
                                     in addition to the top.
        :type fixPotentialAtBottom: ``bool``
        :param tol: tolerance of underlying PDE
        :type tol: positive ``float``
        """
        super(ForwardModelWithPotential, self).__init__()
        self.__domain = domain
        self.__trafo = makeTransformation(domain, coordinates)

        try:
            n=len(w)
            m=len(data)
            if not m == n:
                raise ValueError("Length of weight and data must be the same.")
            self.__weight = w
            self.__data = data
        except TypeError:
            self.__weight = [w]
            self.__data = [data]

        BX = boundingBox(domain)
        DIM = domain.getDim()
        x = domain.getX()
        self.__pde=LinearSinglePDE(domain)
        self.__pde.getSolverOptions().setTolerance(tol)
        self.__pde.setSymmetryOn()
        z=x[DIM-1]
        q0=whereZero(z-BX[DIM-1][1])
        if fixPotentialAtBottom: q0+=whereZero(z-BX[DIM-1][0])
        self.__pde.setValue(q=q0)

        self.edge_lengths=np.asarray(boundingBoxEdgeLengths(domain))
        self.diameter=1./sqrt(sum(1./self.edge_lengths**2))

        self.__origweight=[]
        for s in range(len(self.__weight)):
            # save a copy of the original weights in case of rescaling
            self.__origweight.append(1.*self.__weight[s])

        if not self.__trafo.isCartesian():
            fd=1./self.__trafo.getScalingFactors()
            fw=self.__trafo.getScalingFactors()*sqrt(self.__trafo.getVolumeFactor())
            for s in range(len(self.__weight)):
                self.__weight[s] = fw * self.__weight[s]
                self.__data[s]   = fd * self.__data[s]

    def _rescaleWeights(self, scale=1., fetch_factor=1.):
        """
        rescales the weights such that

        *sum_s integrate( ( weight_i[s] *data_i[s]) (weight_j[s]*1/L_j) * L**2 * fetch_factor )=scale*
        """
        if not scale > 0:
             raise ValueError("Value for scale must be positive.")
        A=0
        # copy back original weights before rescaling
        self.__weight=[1.*ow for ow in self.__origweight]

        for s in range(len(self.__weight)):
            A += integrate(abs(inner(self.__weight[s], self.__data[s]) * inner(self.__weight[s], 1/self.edge_lengths) * fetch_factor))
        if A > 0:
            A=sqrt(scale/A)/self.diameter
            if not self.__trafo.isCartesian():
                A*=self.__trafo.getScalingFactors()*sqrt(self.__trafo.getVolumeFactor())
            for s in range(len(self.__weight)):
                self.__weight[s]*=A
        else:
            raise ValueError("Rescaling of weights failed.")

    def getDomain(self):
        """
        Returns the domain of the forward model.

        :rtype: `Domain`
        """
        return self.__domain

    def getCoordinateTransformation(self):
        """
        returns the coordinate transformation being used

        :rtype: ``CoordinateTransformation``
        """
        return self.__trafo

    def getPDE(self):
        """
        Return the underlying PDE.

        :rtype: `LinearPDE`
        """
        return self.__pde

    def _getDefect(self, result):
        """
        Returns the defect value.

        :param result: a result vector
        :type result: `Vector`
        :rtype: ``float``
        """
        A=0.
        for s in range(len(self.__weight)):
            A += integrate( inner(self.__weight[s], self.__data[s]-result)**2 )
        return A/2

    def getDefectGradient(self, result):
        Y=0.
        for s in range(len(self.__weight)):
            Y = inner(self.__weight[s], self.__data[s]-result) * self.__weight[s] + Y
        return Y

    def getSurvey(self, index=None):
        """
        Returns the pair (data_index, weight_index), where data_i is the data
        of survey i, weight_i is the weighting factor for survey i.
        If index is None, all surveys will be returned in a pair of lists.
        """
        if index is None:
            return self.__data, self.__weight
        if index>=len(self.__data):
            raise IndexError("Forward model only has %d surveys"%len(self.__data))
        return self.__data[index], self.__weight[index]

