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

"""Forward model for DC Resistivity"""
from __future__ import division, print_function

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['DcRes']

from .base import ForwardModel
from esys.downunder.coordinates import makeTransformation
from esys.escript import Scalar, DiracDeltaFunctions
from esys.escript.pdetools import Locator
from esys.escript.linearPDEs import LinearPDE
from esys.escript.util import *


class DcRes(ForwardModel):
    """
    Forward Model for DC resistivity, with a given source pair.
    The cost function is defined as:

    :math: defect = 1/2 (sum_s sum_pq w_pqs * ((phi_sp-phi_sq)-v_pqs)**2

    """
    def __init__(self, domain, locator, delphiIn, sampleTags, phiPrimary,
                 sigmaPrimary, w=1., coordinates=None, tol=1e-8,
                 saveMemory=True, b=None):
        """
        setup new forward model
        
        :param domain: the domain of the model
        :type: escript domain
        :param locator: contains locator to the measurement pairs
        :type: `list` of ``Locator``
        :param: delphiIn: this is v_pq, the potential difference for the
                          current source and a set of measurement pairs.
                          A list of measured potential differences is expected.
                          Note this should be the secondary potential only.
        :type delphiIn: tuple
        :param sampleTags: tags of measurement points from which potential
                           differences will be calculated.
        :type sampleTags: list of tuples
        :param phiPrimary: primary potential.
        :type phiPrimary: `Scalar`
        """
        super(DcRes, self).__init__()

        if not isinstance(sampleTags, list):
            raise ValueError("sampleTags must be a list.")
        if not len(sampleTags) == len(delphiIn):
            raise ValueError("sampleTags and delphiIn must have the same length.")
        if not len(sampleTags)>0:
            raise ValueError("sampleTags list is empty.")
        if not isinstance(sampleTags[0], tuple) and not isinstance(sampleTags[0], list):
            raise ValueError("sampleTags must be a list of tuples or a list of lists.")

        if isinstance(w, float) or isinstance(w, int):
            w = [float(w) for z in delphiIn]
            self.__w = w
        else:
            self.__w = w
        if not len(w) == len(delphiIn):
            raise ValueError("Number of confidence factors and number of potential input values don't match.")

        self.__trafo = makeTransformation(domain, coordinates)
        if not self.getCoordinateTransformation().isCartesian():
            raise ValueError("Non-Cartesian Coordinates are not supported yet.")
        if not isinstance(locator, Locator):
            raise ValueError("locator must be an escript locator object")

        self.__domain = domain
        self.__tol = tol
        self.__locator = locator
        self.delphiIn = delphiIn
        self.__sampleTags = sampleTags
        self.__sigmaPrimary = sigmaPrimary
        self.__phiPrimary = phiPrimary
        self.__pde = None
        if not saveMemory:
            self.__pde=self.setUpPDE()

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

    def getPrimaryPotential(self):
        """
        returns the primary potential
        :rtype: `Data`
        """
        return self.__phiPrimary

    def setUpPDE(self):
        """
        Return the underlying PDE.

        :rtype: `LinearPDE`
        """
        if self.__pde is None:
            dom=self.__domain
            x = dom.getX()
            DIM=dom.getDim()
            q=whereZero(x[DIM-1]-inf(x[DIM-1]))
            for i in range(DIM-1):
                xi=x[i]
                q+=whereZero(xi-inf(xi))+whereZero(xi-sup(xi))

            pde=LinearPDE(dom, numEquations=1)
            pde.getSolverOptions().setTolerance(self.__tol)
            pde.setSymmetryOn()
            A=pde.createCoefficient('A')
            X=pde.createCoefficient('X')
            pde.setValue(A=A, X=X, q=q)

        else:
            pde=self.__pde
            pde.resetRightHandSideCoefficients()
        return pde

    def getArguments(self, sigma):
        """
        Returns precomputed values shared by `getDefect()` and `getGradient()`.

        :param sigma: conductivity
        :type sigma: ``Data`` of shape (1,)
        :return: phi
        :rtype: ``Data`` of shape (1,)
        """
        pde = self.setUpPDE()
        X = (self.__sigmaPrimary - sigma) * grad(self.__phiPrimary)
        pde.setValue(A=sigma * kronecker(self.__domain), X=X)
        phi = pde.getSolution()
        loc = self.__locator
        loc_phi = loc.getValue(phi)
        return phi, loc_phi

    def getDefect(self, sigma, phi, loc_phi):
        """
        Returns the defect value.

        :param sigma: a suggestion for conductivity
        :type sigma: ``Data`` of shape (1,)
        :param phi: potential field
        :type phi: ``Data`` of shape (1,)
        :rtype: ``float``
        """
        val=loc_phi
        if not isinstance(val,list):
            tmp=val
            val=[tmp]
        # print "val=",val
        length=len(val)
        # print self.__sampleTags[0]
        if ((self.__sampleTags[0][1]!="-" and (length%2) != 0) or \
                (self.__sampleTags[0][1]!="-" and length/2 != len(self.delphiIn))):
            raise ValueError("length of locator is wrong")

        delphi_calc=[]
        if self.__sampleTags[0][1] != "-":
            for i in range(0,length,2):
                delphi_calc.append(val[i+1]-val[i])
        else:
            for i in range(length):
                delphi_calc.append(val[i])
        A=0
        if (self.__sampleTags[0][1] != "-"):
            for i in range(length//2):
                A += (self.__w[i]*(delphi_calc[i]-self.delphiIn[i])**2)
        else:
            for i in range(length):
                A += (self.__w[i]*(delphi_calc[i]-self.delphiIn[i])**2)
        return  A/2

    def getGradient(self, sigma, phi, loc_phi):
        """
        Returns the gradient of the defect with respect to density.

        :param sigma: a suggestison for conductivity
        :type sigma: ``Data`` of shape (1,)
        :param phi: potential field
        :type phi: ``Data`` of shape (1,)
        """
        val=loc_phi
        if not isinstance(val,list):
            tmp=val
            val=[tmp]
        sampleTags=self.__sampleTags

        jointSamples={}
        for i in range(0,2*len(sampleTags),2): #2*len because sample tags is a list of tuples
            half = i//2
            if sampleTags[half][1]!="-":
                tmp=(val[i+1]-val[i]-self.delphiIn[half])*self.__w[i]
            else:
                tmp=(val[i]-self.delphiIn[i//2]) *self.__w[i]
            sample = sampleTags[half]
            if sample[1]!="-":
                if sample[0] in jointSamples.keys():
                    jointSamples[sample[0]].append((sample[1], -tmp))
                else:
                    jointSamples[sampleTags[half][0]]=[(sample[1],-tmp)]
                if sample[1] in jointSamples.keys():
                    jointSamples[sample[1]].append((sample[0], tmp))
                else:
                    jointSamples[sample[1]]=[(sample[0], tmp)]
            else:
                if sample[0] in jointSamples.keys():
                    jointSamples[sample[0]].append((sample[1], tmp))
                else:
                    jointSamples[sampleTags[half][0]]=[(sample[1],tmp)]

        pde =self.setUpPDE()
        dom=self.__domain
        # conPrimary=self.__sigmaPrimary
        # APrimary = conPrimary * kronecker(dom)

        y_dirac = Scalar(0,DiracDeltaFunctions(dom))
        for i in jointSamples:
            total=0
            for j in jointSamples[i]:
                total+=j[1]
            y_dirac.setTaggedValue(i,total)

        pde.setValue(A=sigma*kronecker(dom), y_dirac=y_dirac)
        u=pde.getSolution()
        retVal=-inner(grad(u),grad(phi+self.__phiPrimary))
        return retVal

