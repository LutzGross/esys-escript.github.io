##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

from __future__ import print_function, division

import esys.escript as es
from esys.escript.linearPDEs import LinearPDE
from esys.escript.pdetools import Locator
from math import pi
from esys.weipa import saveSilo

try:
    xrange
except NameError:
    xrange = range



class DcResistivityForward(object):
    """
    This class allows for the solution of dc resistivity forward problems via
    the calculation of a primary and secondary potential. Conductivity values
    are to be provided for the primary problem which is a homogeneous half space
    of a chosen conductivity and for the secondary problem which typically 
    varies it conductivity spatially across the domain. The primary potential
    acts as a reference point typically based of some know reference conductivity
    however any value will suffice. The primary potential is implemented to 
    avoid the use of dirac delta functions.
    """
    def __init__(self):
        """
        This is a skeleton class for all the other forward modeling classes.
        """
        pass

    def getPotential(self):
        raise NotImplementedError


    def getApparentResistivity(self):
        raise NotImplementedError

    def checkBounds(self):
        X = es.ContinuousFunction(self.domain).getX()
        xDim=[es.inf(X[0]),es.sup(X[0])]
        yDim=[es.inf(X[1]),es.sup(X[1])]
        zDim=[es.inf(X[2]),es.sup(X[2])]
        for i in range(self.numElectrodes):
            if (self.electrodes[i][0] < xDim[0] or self.electrodes[i][0] > xDim[1]
                    or self.electrodes[i][1] < yDim[0] or self.electrodes[i][1] > yDim[1]):
                print (self.electrodes[i])
                raise ValueError("Electrode setup extents past domain dimentions")
    def getElectrodes(self):
        """
        retuns the list of electrodes with locations
        """
        return self.electrodes

class SchlumbergerSurvey(DcResistivityForward):
    """
    Schlumberger survey forward calculation
    """
    def __init__(self, domain, primaryConductivity, secondaryConductivity,
            current, a, n, midPoint, directionVector, numElectrodes):
        super(SchlumbergerSurvey, self).__init__()
        """
        :param domain: Domain of the model
        :type domain: `Domain`
        :param primaryConductivity: preset primary conductivity data object
        :type primaryConductivity: data
        :param secondaryConductivity: preset secondary conductivity data object
        :type secondaryConductivity: data
        :param current: amount of current to be injected at the current electrode
        :type current: float or int
        :param a: the spacing between current and potential electrodes
        :type  a: float or int
        :param n: multiple of spacing between electrodes. typicaly interger
        :type  n: float or int
        :param midPoint: midPoint of the survey, as a list containing x,y coords
        :type a: list
        :param directionVector: two element list specifying the direction the
            survey should extend from the midpoint
        :type a: list
        :param numElectrodes: the number of electrodes to be used in the survey
            must be a multiple of 2 for polepole survey:
        :type numElectrodes: int
        """
        self.domain=domain
        self.primaryConductivity=primaryConductivity
        self.secondaryConductivity=secondaryConductivity
        self.a=a
        self.n=n
        self.current=current
        self.numElectrodes=numElectrodes
        self.delPhiPrimaryList=[]
        self.delPhiSecondaryList=[]
        self.samples=[]
        self.sources=[]
        self.electrodeDict={}
        self.electrodeTags=[]
        if (numElectrodes < 4):
            raise ValueError("numElectrodes must atleast 4 for schlumberger surveys")
        if n > ((numElectrodes-2)//2):
            raise ValueError("specified n does not fit max n = %d"%((numElectrodes-2)//2))
        if len(directionVector) == 2:
            electrodes=[]
            start=[]
            start.append(midPoint[0] - (((numElectrodes-1)*a)/2. * directionVector[0]))
            start.append(midPoint[1] - (((numElectrodes-1)*a)/2. * directionVector[1]))
            for i in range(numElectrodes):
                electrodes.append([start[0]+(directionVector[0]*i*a), start[1]+(directionVector[1]*i*a),0])
                self.electrodeTags.append("e%d"%i)
                self.electrodeDict[self.electrodeTags[i]]=electrodes[i]
        else:
            raise NotImplementedError("2d forward model is not yet implemented please provide a 2 component directionVector for a 3d survey")
        self.electrodes=electrodes
        self.checkBounds()

    def getPotentialNumeric(self):
        """
        Returns 3 list each made up of a number of list containing primary, secondary and total
        potentials diferences. Each of the lists contain a list for each value of n.
        """
        primCon=self.primaryConductivity
        coords=self.domain.getX()
        tol=1e-8
        pde=LinearPDE(self.domain, numEquations=1)
        pde.getSolverOptions().setTolerance(tol)
        pde.setSymmetryOn()
        primaryPde=LinearPDE(self.domain, numEquations=1)
        primaryPde.getSolverOptions().setTolerance(tol)
        primaryPde.setSymmetryOn()
        DIM=self.domain.getDim()
        x=self.domain.getX()
        q=es.whereZero(x[DIM-1]-es.inf(x[DIM-1]))
        for i in xrange(DIM-1):
            xi=x[i]
            q+=es.whereZero(xi-es.inf(xi))+es.whereZero(xi-es.sup(xi))
        A = self.secondaryConductivity * es.kronecker(self.domain)
        APrimary = self.primaryConductivity * es.kronecker(self.domain)
        pde.setValue(A=A,q=q)
        primaryPde.setValue(A=APrimary,q=q)

        delPhiSecondaryList = []
        delPhiPrimaryList = []
        delPhiTotalList = []
        for i in range(1,self.n+1): # 1 to n
            maxR = self.numElectrodes - 1 - (2*i) #max amount of readings that will fit in the survey
            delPhiSecondary = []
            delPhiPrimary = []
            delPhiTotal = []
            for j in range(maxR):
                y_dirac=es.Scalar(0,es.DiracDeltaFunctions(self.domain))
                y_dirac.setTaggedValue(self.electrodeTags[j],self.current)
                y_dirac.setTaggedValue(self.electrodeTags[j + ((2*i) + 1)],-self.current)
                self.sources.append([self.electrodeTags[j], self.electrodeTags[j + ((2*i) + 1)]])
                primaryPde.setValue(y_dirac=y_dirac)
                numericPrimaryPot = primaryPde.getSolution()
                X=(primCon-self.secondaryConductivity) * es.grad(numericPrimaryPot)
                pde.setValue(X=X)
                u=pde.getSolution()
                loc=Locator(self.domain,[self.electrodes[j+i],self.electrodes[j+i+1]])
                self.samples.append([self.electrodeTags[j+i],self.electrodeTags[j+i+1]])
                valPrimary=loc.getValue(numericPrimaryPot)
                valSecondary=loc.getValue(u)
                delPhiPrimary.append(valPrimary[1]-valPrimary[0])
                delPhiSecondary.append(valSecondary[1]-valSecondary[0])
                delPhiTotal.append(delPhiPrimary[j]+delPhiSecondary[j])
            delPhiPrimaryList.append(delPhiPrimary)
            delPhiSecondaryList.append(delPhiSecondary)
            delPhiTotalList.append(delPhiTotal)

        self.delPhiPrimaryList=delPhiPrimaryList
        self.delPhiSecondaryList=delPhiSecondaryList
        self.delPhiTotalList = delPhiTotalList
        return [delPhiPrimaryList, delPhiSecondaryList, delPhiTotalList]



    def getPotentialAnalytic(self):
        """
        Returns 3 list each made up of a number of list containing primary, secondary and total
        potentials diferences. Each of the lists contain a list for each value of n.
        """
        coords=self.domain.getX()
        pde=LinearPDE(self.domain, numEquations=1)
        tol=1e-8
        pde.getSolverOptions().setTolerance(tol)
        pde.setSymmetryOn()
        primCon=self.primaryConductivity
        DIM=self.domain.getDim()
        x=self.domain.getX()
        q=es.whereZero(x[DIM-1]-es.inf(x[DIM-1]))
        for i in xrange(DIM-1):
            xi=x[i]
            q+=es.whereZero(xi-es.inf(xi))+es.whereZero(xi-es.sup(xi))
        A = self.secondaryConductivity * es.kronecker(self.domain)
        pde.setValue(A=A,q=q)


        delPhiSecondaryList = []
        delPhiPrimaryList = []
        delPhiTotalList = []
        for i in range(1,self.n+1): # 1 to n
            maxR = self.numElectrodes - 1 - (2*i) #max amount of readings that will fit in the survey
            delPhiSecondary = []
            delPhiPrimary = []
            delPhiTotal = []
            for j in range(maxR):
                analyticRsOne=es.Data(0,(3,),es.ContinuousFunction(self.domain))
                analyticRsOne[0]=(coords[0]-self.electrodes[j][0])
                analyticRsOne[1]=(coords[1]-self.electrodes[j][1])
                analyticRsOne[2]=(coords[2])
                rsMagOne=(analyticRsOne[0]**2+analyticRsOne[1]**2+analyticRsOne[2]**2)**0.5
                analyticRsTwo=es.Data(0,(3,),es.ContinuousFunction(self.domain))
                analyticRsTwo[0]=(coords[0]-self.electrodes[j + ((2*i) + 1)][0])
                analyticRsTwo[1]=(coords[1]-self.electrodes[j + ((2*i) + 1)][1])
                analyticRsTwo[2]=(coords[2])
                rsMagTwo=(analyticRsTwo[0]**2+analyticRsTwo[1]**2+analyticRsTwo[2]**2)**0.5
                self.sources.append([self.electrodeTags[j], self.electrodeTags[j + ((2*i) + 1)]])
                rsMagOne+=(es.whereZero(rsMagOne)*0.0000001)
                rsMagTwo+=(es.whereZero(rsMagTwo)*0.0000001)
                
                analyticPrimaryPot=(self.current/(2*pi*primCon*rsMagOne))-(self.current/(2*pi*primCon*rsMagTwo))
                analyticRsOnePower=(analyticRsOne[0]**2+analyticRsOne[1]**2+analyticRsOne[2]**2)**1.5
                analyticRsOnePower = analyticRsOnePower+(es.whereZero(analyticRsOnePower)*0.0001)
                analyticRsTwoPower=(analyticRsTwo[0]**2+analyticRsTwo[1]**2+analyticRsTwo[2]**2)**1.5
                analyticRsTwoPower = analyticRsTwoPower+(es.whereZero(analyticRsTwoPower)*0.0001)

                gradAnalyticPrimaryPot = es.Data(0,(3,),es.ContinuousFunction(self.domain))
                gradAnalyticPrimaryPot[0] =(self.current/(2*pi*primCon)) * ((-analyticRsOne[0]/analyticRsOnePower) + (analyticRsTwo[0]/analyticRsTwoPower))
                gradAnalyticPrimaryPot[1] =(self.current/(2*pi*primCon)) * ((-analyticRsOne[1]/analyticRsOnePower) + (analyticRsTwo[1]/analyticRsTwoPower))
                gradAnalyticPrimaryPot[2] =(self.current/(2*pi*primCon)) * ((-analyticRsOne[2]/analyticRsOnePower) + (analyticRsTwo[2]/analyticRsTwoPower))
                X=(primCon-self.secondaryConductivity) * (gradAnalyticPrimaryPot)
                pde.setValue(X=X)
                u=pde.getSolution()
                loc=Locator(self.domain,[self.electrodes[j+i],self.electrodes[j+i+1]])
                self.samples.append([self.electrodeTags[j+i],self.electrodeTags[j+i+1]])
                valPrimary=loc.getValue(analyticPrimaryPot)
                valSecondary=loc.getValue(u)
                delPhiPrimary.append(valPrimary[1]-valPrimary[0])
                delPhiSecondary.append(valSecondary[1]-valSecondary[0])
                delPhiTotal.append(delPhiPrimary[j]+delPhiSecondary[j])
            delPhiPrimaryList.append(delPhiPrimary)
            delPhiSecondaryList.append(delPhiSecondary)
            delPhiTotalList.append(delPhiTotal)

        self.delPhiPrimaryList=delPhiPrimaryList
        self.delPhiSecondaryList=delPhiSecondaryList
        self.delPhiTotalList = delPhiTotalList
        return [delPhiPrimaryList, delPhiSecondaryList, delPhiTotalList]

    def getApparentResistivity(self, delPhiList):
        resistivityList = []
        n=self.n
        if (self.delPhiPrimaryList==[]):
            self.getPotential()

        nCount=1
        for i in delPhiList:
            resistivity=[]
            for j in i:
                resistivity.append((-j/self.current)*(nCount*(nCount+1))*pi*self.a)
            nCount=nCount+1
            resistivityList.append(resistivity)
        return resistivityList

    def getSourcesSamples(self):
        """
        return a list of tuples of sample locations followed by a list of tuples
        of source locations.
        """
        return self.samples, self.sources

    def getElectrodeDict(self):
        """
        retuns the electrode dictionary
        """
        return self.electrodeDict

class WennerSurvey(DcResistivityForward):
    """
    WennerSurvey forward calculation
    """
    def __init__(self, domain, primaryConductivity, secondaryConductivity,
             current, a, midPoint, directionVector, numElectrodes):
        """
        :param domain: domain of the model
        :type domain: `Domain`
        :param primaryConductivity: preset primary conductivity data object
        :type primaryConductivity: ``data``
        :param secondaryConductivity: preset secondary conductivity data object
        :type secondaryConductivity: ``data``
        :param current: amount of current to be injected at the current electrode
        :type current: ``float`` or ``int``
        :param a: the spacing between current and potential electrodes
        :type a: ``float`` or ``int``
        :param midPoint: midPoint of the survey, as a list containing x,y coords
        :type a: ``list``
        :param directionVector: two element list specifying the direction the
            survey should extend from the midpoint
        :type a: ``list``
        :param numElectrodes: the number of electrodes to be used in the survey
            must be a multiple of 2 for polepole survey
        :type numElectrodes: ``int``
        """
        super(WennerSurvey, self).__init__()
        self.domain=domain
        self.primaryConductivity=primaryConductivity
        self.secondaryConductivity=secondaryConductivity
        self.a=a
        self.current=current
        self.numElectrodes=numElectrodes
        self.delPhiSecondary=[]
        self.delPhiPrimary=[]
        if (numElectrodes<4 ):
            raise ValueError("numElectrodes must be atleast 4 for Wenner surveys")
        if len(directionVector) == 2:
            electrodes = []
            start=[]
            start.append(midPoint[0] - (((numElectrodes-1)*a)/2. * directionVector[0]))
            start.append(midPoint[1] - (((numElectrodes-1)*a)/2. * directionVector[1]))
            for i in range(numElectrodes):
                electrodes.append([start[0]+(directionVector[0]*i*a),
                        start[1]+(directionVector[1]*i*a),0])
        else:
            raise NotImplementedError("2d forward model is not yet implemented please provide a 2 component directionVector for a 3d survey")
        self.electrodes=electrodes
        self.checkBounds()

    def getPotential(self):
        """
        returns a list containing 3 lists one for each the primary, secondary
        and total potential.
        """

        primCon=self.primaryConductivity
        coords=self.domain.getX()
        pde=LinearPDE(self.domain, numEquations=1)
        tol=1e-8
        pde.getSolverOptions().setTolerance(tol)
        pde.setSymmetryOn()

        DIM=self.domain.getDim()
        x=self.domain.getX()
        q=es.whereZero(x[DIM-1]-es.inf(x[DIM-1]))
        for i in xrange(DIM-1):
            xi=x[i]
            q+=es.whereZero(xi-es.inf(xi))+es.whereZero(xi-es.sup(xi))
        A = self.secondaryConductivity * es.kronecker(self.domain)
        pde.setValue(A=A,q=q)

        delPhiSecondary = []
        delPhiPrimary = []
        delPhiTotal = []
        if(len(self.electrodes[0])==3):

            for i in range(self.numElectrodes-3):
                analyticRsOne=es.Data(0,(3,),es.ContinuousFunction(self.domain))
                analyticRsOne[0]=(coords[0]-self.electrodes[i][0])
                analyticRsOne[1]=(coords[1]-self.electrodes[i][1])
                analyticRsOne[2]=(coords[2])
                rsMagOne=(analyticRsOne[0]**2+analyticRsOne[1]**2+analyticRsOne[2]**2)**0.5
                analyticRsTwo=es.Data(0,(3,),es.ContinuousFunction(self.domain))
                analyticRsTwo[0]=(coords[0]-self.electrodes[i+3][0])
                analyticRsTwo[1]=(coords[1]-self.electrodes[i+3][1])
                analyticRsTwo[2]=(coords[2])
                rsMagTwo=(analyticRsTwo[0]**2+analyticRsTwo[1]**2+analyticRsTwo[2]**2)**0.5
                rsMagOne+=(es.whereZero(rsMagOne)*0.0000001)
                rsMagTwo+=(es.whereZero(rsMagTwo)*0.0000001)
                analyticPrimaryPot=(self.current/(2*pi*primCon*rsMagOne))-(self.current/(2*pi*primCon*rsMagTwo))

                analyticRsOnePower=(analyticRsOne[0]**2+analyticRsOne[1]**2+analyticRsOne[2]**2)**1.5
                analyticRsOnePower = analyticRsOnePower+(es.whereZero(analyticRsOnePower)*0.0001)
                analyticRsTwoPower=(analyticRsTwo[0]**2+analyticRsTwo[1]**2+analyticRsTwo[2]**2)**1.5
                analyticRsTwoPower = analyticRsTwoPower+(es.whereZero(analyticRsTwoPower)*0.0001)

                gradAnalyticPrimaryPot = es.Data(0,(3,),es.ContinuousFunction(self.domain))
                gradAnalyticPrimaryPot[0] =(self.current/(2*pi*primCon)) \
                        * ((-analyticRsOne[0]/analyticRsOnePower) \
                            + (analyticRsTwo[0]/analyticRsTwoPower))
                gradAnalyticPrimaryPot[1] =(self.current/(2*pi*primCon)) \
                        * ((-analyticRsOne[1]/analyticRsOnePower) \
                            + (analyticRsTwo[1]/analyticRsTwoPower))
                gradAnalyticPrimaryPot[2] =(self.current/(2*pi*primCon)) \
                        * ((-analyticRsOne[2]/analyticRsOnePower)
                            + (analyticRsTwo[2]/analyticRsTwoPower))
                X=(primCon-self.secondaryConductivity) * (gradAnalyticPrimaryPot)
                pde.setValue(X=X)
                u=pde.getSolution()
                loc=Locator(self.domain,[self.electrodes[i+1],self.electrodes[i+2]])
                valPrimary=loc.getValue(analyticPrimaryPot)
                valSecondary=loc.getValue(u)
                delPhiPrimary.append(valPrimary[1]-valPrimary[0])
                delPhiSecondary.append(valSecondary[1]-valSecondary[0])
                delPhiTotal.append(delPhiPrimary[i]+delPhiSecondary[i])
        else:
            raise NotImplementedError("2d forward model is not yet implemented")

        self.delPhiSecondary = delPhiSecondary
        self.delPhiPrimary = delPhiPrimary
        self.delPhiTotal=delPhiTotal
        return [delPhiPrimary, delPhiSecondary, delPhiTotal]

    def getApparentResistivityPrimary(self):
        resistivity = []

        if (self.delPhiPrimary==[]):
            self.getPotential()

        for i in self.delPhiPrimary:
            resistivity.append((-i/self.current)*2*pi*self.a)

        return resistivity


    def getApparentResistivitySecondary(self):
        resistivity = []

        if (self.delPhiSecondary==[]):
            self.getPotential()

        for i in self.delPhiSecondary:
            resistivity.append((-i/self.current)*2*pi*self.a)

        return resistivity

    def getApparentResistivityTotal(self):
            resistivity = []

            if (self.delPhiSecondary==[]):
                self.getPotential()

            for i in range(len(self.delPhiSecondary)):
                resistivity.append((-(self.delPhiSecondary[i]+self.delPhiPrimary[i])/self.current)*2*pi*self.a)

            return resistivity


class DipoleDipoleSurvey(DcResistivityForward):
    """
    DipoleDipoleSurvey forward modeling
    """
    def __init__(self, domain, primaryConductivity, secondaryConductivity,
            current, a, n, midPoint, directionVector, numElectrodes):
        super(DipoleDipoleSurvey, self).__init__()
        """
        :param domain: domain of the model
        :type domain: `Domain`
        :param primaryConductivity: preset primary conductivity data object
        :type primaryConductivity: data
        :param secondaryConductivity: preset secondary conductivity data object
        :type secondaryConductivity: data
        :param current: amount of current to be injected at the current electrode
        :type current: float or int
        :param a: the spacing between current and potential electrodes
        :type  a: float or int
        :param n: multiple of spacing between electrodes. typicaly interger
        :type  n: float or int
        :param midPoint: midPoint of the survey, as a list containing x,y coords
        :type a: list
        :param directionVector: two element list specifying the direction the
            survey should extend from the midpoint
        :type a: list
        :param numElectrodes: the number of electrodes to be used in the survey
            must be a multiple of 2 for polepole survey:
        :type numElectrodes: int
        """
        self.domain=domain
        self.primaryConductivity=primaryConductivity
        self.secondaryConductivity=secondaryConductivity
        self.a=a
        self.n=n
        self.current=current
        self.numElectrodes=numElectrodes
        self.delPhiPrimaryList=[]
        self.delPhiSecondaryList=[]
        if (numElectrodes<4 ):
            raise ValueError("numElectrodes must be atleast 4 for DipoleDipole surveys")
        if n > (numElectrodes-3):
            raise ValueError("specified n does not fit max n = %d"%((numElectrodes-2)//2))
        if len(directionVector) == 2:
            electrodes=[]
            start=[]
            start.append(midPoint[0] - (((numElectrodes-1)*a)/2. * directionVector[0]))
            start.append(midPoint[1] - (((numElectrodes-1)*a)/2. * directionVector[1]))
            for i in range(numElectrodes):
                electrodes.append([start[0]+(directionVector[0]*i*a), start[1]+(directionVector[1]*i*a),0])
        else:
            raise NotImplementedError("2d forward model is not yet implemented please provide a 2 component directionVector for a 3d survey")
        self.electrodes=electrodes
        self.checkBounds()


    def getPotential(self):
        """
        Returns 3 list each made up of a number of list containing primary, secondary and total
        potentials diferences. Each of the lists contain a list for each value of n.
        """
        coords=self.domain.getX()
        pde=LinearPDE(self.domain, numEquations=1)
        tol=1e-8
        pde.getSolverOptions().setTolerance(tol)
        pde.setSymmetryOn()
        primCon=self.primaryConductivity
        DIM=self.domain.getDim()
        x=self.domain.getX()
        q=es.whereZero(x[DIM-1]-es.inf(x[DIM-1]))
        for i in xrange(DIM-1):
            xi=x[i]
            q+=es.whereZero(xi-es.inf(xi))+es.whereZero(xi-es.sup(xi))
        A = self.secondaryConductivity * es.kronecker(self.domain)
        pde.setValue(A=A,q=q)


        delPhiSecondaryList = []
        delPhiPrimaryList = []
        delPhiTotalList = []
        for i in range(1,self.n+1): # 1 to n
            maxR = self.numElectrodes - 2 - (i) #max amount of readings that will fit in the survey
            delPhiSecondary = []
            delPhiPrimary = []
            delPhiTotal = []
            for j in range(maxR):
                analyticRsOne=es.Data(0,(3,),es.ContinuousFunction(self.domain))
                analyticRsOne[0]=(coords[0]-self.electrodes[j][0])
                analyticRsOne[1]=(coords[1]-self.electrodes[j][1])
                analyticRsOne[2]=(coords[2])
                rsMagOne=(analyticRsOne[0]**2+analyticRsOne[1]**2+analyticRsOne[2]**2)**0.5
                analyticRsTwo=es.Data(0,(3,),es.ContinuousFunction(self.domain))
                analyticRsTwo[0]=(coords[0]-self.electrodes[j + 1][0])
                analyticRsTwo[1]=(coords[1]-self.electrodes[j + 1][1])
                analyticRsTwo[2]=(coords[2])
                rsMagTwo=(analyticRsTwo[0]**2+analyticRsTwo[1]**2+analyticRsTwo[2]**2)**0.5
                rsMagOne+=(es.whereZero(rsMagOne)*0.0000001)
                rsMagTwo+=(es.whereZero(rsMagTwo)*0.0000001)
                analyticPrimaryPot=(self.current/(2*pi*primCon*rsMagTwo))-(self.current/(2*pi*primCon*rsMagOne))

                analyticRsOnePower=(analyticRsOne[0]**2+analyticRsOne[1]**2+analyticRsOne[2]**2)**1.5
                analyticRsOnePower = analyticRsOnePower+(es.whereZero(analyticRsOnePower)*0.0001)
                analyticRsTwoPower=(analyticRsTwo[0]**2+analyticRsTwo[1]**2+analyticRsTwo[2]**2)**1.5
                analyticRsTwoPower = analyticRsTwoPower+(es.whereZero(analyticRsTwoPower)*0.0001)

                gradAnalyticPrimaryPot = es.Data(0,(3,),es.ContinuousFunction(self.domain))
                gradAnalyticPrimaryPot[0] =(self.current/(2*pi*primCon)) * ((analyticRsOne[0]/analyticRsOnePower) - (analyticRsTwo[0]/analyticRsTwoPower))
                gradAnalyticPrimaryPot[1] =(self.current/(2*pi*primCon)) * ((analyticRsOne[1]/analyticRsOnePower) - (analyticRsTwo[1]/analyticRsTwoPower))
                gradAnalyticPrimaryPot[2] =(self.current/(2*pi*primCon)) * ((analyticRsOne[2]/analyticRsOnePower) - (analyticRsTwo[2]/analyticRsTwoPower))
                X=(primCon-self.secondaryConductivity) * (gradAnalyticPrimaryPot)
                pde.setValue(X=X)
                u=pde.getSolution()
                loc=Locator(self.domain,[self.electrodes[1+j+i],self.electrodes[j+i+2]])
                valPrimary=loc.getValue(analyticPrimaryPot)
                valSecondary=loc.getValue(u)
                delPhiPrimary.append(valPrimary[1]-valPrimary[0])
                delPhiSecondary.append(valSecondary[1]-valSecondary[0])
                delPhiTotal.append(delPhiPrimary[j]+delPhiSecondary[j])
            delPhiPrimaryList.append(delPhiPrimary)
            delPhiSecondaryList.append(delPhiSecondary)
            delPhiTotalList.append(delPhiTotal)

        self.delPhiPrimaryList=delPhiPrimaryList
        self.delPhiSecondaryList=delPhiSecondaryList
        self.delPhiTotalList = delPhiTotalList
        return [delPhiPrimaryList, delPhiSecondaryList, delPhiTotalList]

    def getApparentResistivityPrimary(self):
        resistivityList = []
        n=self.n
        if (self.delPhiPrimaryList==[]):
            self.getPotential()

        nCount=1
        for i in self.delPhiPrimaryList:
            resistivity=[]
            for j in i:
                resistivity.append((-j/self.current)*(nCount*(nCount+1)*(nCount+2))*pi*self.a)
            nCount=nCount+1
            resistivityList.append(resistivity)
        return resistivityList


    def getApparentResistivitySecondary(self):
        resistivityList = []
        n=self.n
        if (self.delPhiSecondaryList==[]):
            self.getPotential()

        nCount=1
        for i in self.delPhiSecondaryList:
            resistivity=[]
            for j in i:
                resistivity.append((-j/self.current)*(nCount*(nCount+1)*(nCount+2))*pi*self.a)
            nCount=nCount+1
            resistivityList.append(resistivity)
        return resistivityList

    def getApparentResistivityTotal(self):
        resistivityList = []
        n=self.n
        if (self.delPhiSecondaryList==[]):
            self.getPotential()

        nCount=1
        for i in self.delPhiTotalList:
            resistivity=[]
            for j in i:
                resistivity.append((-j/self.current)*(nCount*(nCount+1)*(nCount+2))*pi*self.a)
            nCount=nCount+1
            resistivityList.append(resistivity)
        return resistivityList

class PoleDipoleSurvey(DcResistivityForward):
    """
    Forward model class for a poledipole setup
    """
    def __init__(self, domain, primaryConductivity, secondaryConductivity, current, a, n, midPoint, directionVector, numElectrodes):
        """
        :param domain: domain of the model
        :type domain: `Domain`
        :param primaryConductivity: preset primary conductivity data object
        :type primaryConductivity: data
        :param secondaryConductivity: preset secondary conductivity data object
        :type secondaryConductivity: data
        :param current: amount of current to be injected at the current electrode
        :type current: float or int
        :param a: the spacing between current and potential electrodes
        :type  a: float or int
        :param n: multiple of spacing between electrodes. typicaly interger
        :type  n: float or int
        :param midPoint: midPoint of the survey, as a list containing x,y coords
        :type a: list
        :param directionVector: two element list specifying the direction the
            survey should extend from the midpoint
        :type a: list
        :param numElectrodes: the number of electrodes to be used in the survey
            must be a multiple of 2 for polepole survey:
        :type numElectrodes: int
        """

        super(PoleDipoleSurvey, self).__init__()

        self.domain=domain
        self.primaryConductivity=primaryConductivity
        self.secondaryConductivity=secondaryConductivity
        self.a=a
        self.n=n
        self.current=current
        self.numElectrodes=numElectrodes
        self.delPhiPrimaryList=[]
        self.delPhiSecondaryList=[]
        if ((numElectrodes%4) != 0 ):
            raise ValueError("numElectrodes must be a multiple of 4 for schlumberger surveys")
        if n > (numElectrodes-3):
            raise ValueError("specified n does not fit max n = %d"%((numElectrodes-2)//2))
        if len(directionVector) == 2:
            electrodes=[]
            start=[]
            start.append(midPoint[0] - (((numElectrodes-1)*a)/2. * directionVector[0]))
            start.append(midPoint[1] - (((numElectrodes-1)*a)/2. * directionVector[1]))
            for i in range(numElectrodes):
                electrodes.append([start[0]+(directionVector[0]*i*a), start[1]+(directionVector[1]*i*a),0])
        else:
            raise NotImplementedError("2d forward model is not yet implemented please provide a 2 component directionVector for a 3d survey")
        self.electrodes=electrodes
        self.checkBounds()

    def getPotential(self):
        """
        Returns 3 list each made up of a number of list containing primary, secondary and total
        potentials diferences. Each of the lists contain a list for each value of n.
        """

        primCon=self.primaryConductivity
        coords=self.domain.getX()
        pde=LinearPDE(self.domain, numEquations=1)
        tol=1e-8
        pde.getSolverOptions().setTolerance(tol)
        pde.setSymmetryOn()

        DIM=self.domain.getDim()
        x=self.domain.getX()
        q=es.whereZero(x[DIM-1]-es.inf(x[DIM-1]))
        for i in xrange(DIM-1):
            xi=x[i]
            q+=es.whereZero(xi-es.inf(xi))+es.whereZero(xi-es.sup(xi))
        A = self.secondaryConductivity * es.kronecker(self.domain)
        pde.setValue(A=A,q=q)

        delPhiSecondaryList = []
        delPhiPrimaryList = []
        delPhiTotalList = []
        for i in range(1,self.n+1): # 1 to n
            maxR = self.numElectrodes - 1 - i #max amount of readings that will fit in the survey
            delPhiSecondary = []
            delPhiPrimary = []
            delPhiTotal = []
            for j in range(maxR):
                analyticRs=es.Data(0,(3,),es.ContinuousFunction(self.domain))
                analyticRs[0]=(coords[0]-self.electrodes[j][0])
                analyticRs[1]=(coords[1]-self.electrodes[j][1])
                analyticRs[2]=(coords[2])
                rsMag=(analyticRs[0]**2+analyticRs[1]**2+analyticRs[2]**2)**0.5
                analyticPrimaryPot=(self.current*(1./primCon))/(2*pi*(rsMag+(es.whereZero(rsMag)*0.0000001))) #the magic number 0.0000001 is to avoid devide by 0

                analyticRsPolePower=(analyticRs[0]**2+analyticRs[1]**2+analyticRs[2]**2)**1.5
                analyticRsPolePower = analyticRsPolePower+(es.whereZero(analyticRsPolePower)*0.0000001)
                gradUPrimary = es.Data(0,(3,),es.ContinuousFunction(self.domain))
                gradUPrimary[0] =(self.current/(2*pi*primCon)) * (analyticRs[0]/analyticRsPolePower)
                gradUPrimary[1] =(self.current/(2*pi*primCon)) * (analyticRs[1]/analyticRsPolePower)
                gradUPrimary[2] =(self.current/(2*pi*primCon)) * (analyticRs[2]/analyticRsPolePower)
                gradUPrimary=-gradUPrimary
                X=(primCon-self.secondaryConductivity) * gradUPrimary
                pde.setValue(X=X)
                u=pde.getSolution()
                loc=Locator(self.domain,[self.electrodes[i+j],self.electrodes[i+j+1]])
                valPrimary=loc.getValue(analyticPrimaryPot)
                valSecondary=loc.getValue(u)
                delPhiPrimary.append(valPrimary[1]-valPrimary[0])
                delPhiSecondary.append(valSecondary[1]-valSecondary[0])
                delPhiTotal.append(delPhiPrimary[j]+delPhiSecondary[j])

            delPhiPrimaryList.append(delPhiPrimary)
            delPhiSecondaryList.append(delPhiSecondary)
            delPhiTotalList.append(delPhiTotal)



        self.delPhiPrimaryList=delPhiPrimaryList
        self.delPhiSecondaryList=delPhiSecondaryList
        self.delPhiTotalList = delPhiTotalList

        return [delPhiPrimaryList, delPhiSecondaryList, delPhiTotalList]


    def getApparentResistivityPrimary(self):
        resistivityList = []
        n=self.n
        if (self.delPhiPrimaryList==[]):
            self.getPotential()

        nCount=1
        for i in self.delPhiPrimaryList:
            resistivity=[]
            for j in i:
                resistivity.append((-j/self.current)*(nCount*(nCount+1))*pi*self.a*2)
            nCount=nCount+1
            resistivityList.append(resistivity)
        return resistivityList


    def getApparentResistivitySecondary(self):
        resistivityList = []
        n=self.n
        if (self.delPhiSecondaryList==[]):
            self.getPotential()

        nCount=1
        for i in self.delPhiSecondaryList:
            resistivity=[]
            for j in i:
                resistivity.append((-j/self.current)*(nCount*(nCount+1))*pi*self.a*2)
            nCount=nCount+1
            resistivityList.append(resistivity)
        return resistivityList

    def getApparentResistivityTotal(self):
        resistivityList = []
        n=self.n
        if (self.delPhiSecondaryList==[]):
            self.getPotential()

        nCount=1
        for i in self.delPhiTotalList:
            resistivity=[]
            for j in i:
                resistivity.append((-j/self.current)*(nCount*(nCount+1))*pi*self.a*2)
            nCount=nCount+1
            resistivityList.append(resistivity)
        return resistivityList

class PolePoleSurvey(DcResistivityForward):
    """
    Forward model class for a polepole setup
    """
    def __init__(self, domain, primaryConductivity, secondaryConductivity, current, a, midPoint, directionVector, numElectrodes):
        """
        :param domain: domain of the model
        :type domain: `Domain`
        :param primaryConductivity: preset primary conductivity data object
        :type primaryConductivity: data
        :param secondaryConductivity: preset secondary conductivity data object
        :type secondaryConductivity: data
        :param current: amount of current to be injected at the current electrode
        :type current: float or int
        :param a: the spacing between current and potential electrodes
        :type a: float or int
        :param midPoint: midPoint of the survey, as a list containing x,y coords
        :type a: list
        :param directionVector: two element list specifying the direction the
            survey should extend from the midpoint
        :type a: list
        :param numElectrodes: the number of electrodes to be used in the survey
            must be a multiple of 2 for polepole survey:
        :type numElectrodes: int
        """

        super(PolePoleSurvey, self).__init__()
        self.domain=domain
        self.primaryConductivity=primaryConductivity
        self.secondaryConductivity=secondaryConductivity
        self.a=a
        self.current=current
        self.numElectrodes=numElectrodes
        self.delPhiSecondary=[]
        self.delPhiPrimary=[]
        if ((numElectrodes<2) != 0 ):
            raise ValueError("numElectrodes must be atleast 2 for pole-pole surveys")
        if len(directionVector) == 2:
            electrodes = []
            start=[]
            start.append(midPoint[0] - (((numElectrodes-1)*a)/2. * directionVector[0]))
            start.append(midPoint[1] - (((numElectrodes-1)*a)/2. * directionVector[1]))
            for i in range(numElectrodes):
                electrodes.append([start[0]+(directionVector[0]*i*a), start[1]+(directionVector[1]*i*a),0])
        else:
            raise NotImplementedError("2d forward model is not yet implemented please provide a 2 component directionVector for a 3d survey")
        self.electrodes=electrodes
        self.checkBounds()

    def getPotential(self):
        """
        returns a list containing 3 lists one for each the primary, secondary
        and total potential.
        """


        primCon=self.primaryConductivity
        coords=self.domain.getX()
        pde=LinearPDE(self.domain, numEquations=1)
        tol=1e-8
        pde.getSolverOptions().setTolerance(tol)
        pde.setSymmetryOn()

        DIM=self.domain.getDim()
        x=self.domain.getX()
        q=es.whereZero(x[DIM-1]-es.inf(x[DIM-1]))
        for i in xrange(DIM-1):
            xi=x[i]
            q+=es.whereZero(xi-es.inf(xi))+es.whereZero(xi-es.sup(xi))
        A = self.secondaryConductivity * es.kronecker(self.domain)
        pde.setValue(A=A,q=q)

        delPhiSecondary = []
        delPhiPrimary = []
        delPhiTotal = []
        if(len(self.electrodes[0])==3):

            for i in range(self.numElectrodes-1):
                analyticRs=es.Data(0,(3,),es.ContinuousFunction(self.domain))
                analyticRs[0]=(coords[0]-self.electrodes[i][0])
                analyticRs[1]=(coords[1]-self.electrodes[i][1])
                analyticRs[2]=(coords[2])
                rsMag=(analyticRs[0]**2+analyticRs[1]**2+analyticRs[2]**2)**0.5
                analyticPrimaryPot=(self.current*(1./primCon))/(2*pi*(rsMag+(es.whereZero(rsMag)*0.0000001))) #the magic number 0.0000001 is to avoid devide by 0
                analyticRsPolePower=(analyticRs[0]**2+analyticRs[1]**2+analyticRs[2]**2)**1.5
                analyticRsPolePower = analyticRsPolePower+(es.whereZero(analyticRsPolePower)*0.0000001)
                gradUPrimary = es.Data(0,(3,),es.ContinuousFunction(self.domain))
                gradUPrimary[0] =(self.current/(2*pi*primCon)) * (analyticRs[0]/analyticRsPolePower)
                gradUPrimary[1] =(self.current/(2*pi*primCon)) * (analyticRs[1]/analyticRsPolePower)
                gradUPrimary[2] =(self.current/(2*pi*primCon)) * (analyticRs[2]/analyticRsPolePower)
                gradUPrimary=-gradUPrimary
                X=(primCon-self.secondaryConductivity) * gradUPrimary
                pde.setValue(X=X)
                u=pde.getSolution()
                loc=Locator(self.domain,self.electrodes[i+1])
                delPhiSecondary.append(loc.getValue(u))
                delPhiPrimary.append(loc.getValue(analyticPrimaryPot))
        else:
            raise NotImplementedError("2d forward model is not yet implemented")

        self.delPhiSecondary = delPhiSecondary
        self.delPhiPrimary = delPhiPrimary
        for i in range(len(delPhiPrimary)):
            delPhiTotal.append(delPhiPrimary[i] + delPhiSecondary[i])
        self.delPhiTotal=delPhiTotal
        return [delPhiPrimary, delPhiSecondary, delPhiTotal]

    def getApparentResistivityPrimary(self):
        resistivity = []

        if (self.delPhiPrimary==[]):
            self.getPotential()

        for i in self.delPhiPrimary:
            resistivity.append((i/self.current)*2*pi*self.a)

        return resistivity


    def getApparentResistivitySecondary(self):
        resistivity = []

        if (self.delPhiSecondary==[]):
            self.getPotential()

        for i in self.delPhiSecondary:
            resistivity.append((i/self.current)*2*pi*self.a)

        return resistivity

    def getApparentResistivityTotal(self):
            resistivity = []

            if (self.delPhiSecondary==[]):
                self.getPotential()

            for i in range(len(self.delPhiSecondary)):
                resistivity.append(((self.delPhiSecondary[i]+self.delPhiPrimary[i])/self.current)*2*pi*self.a)

            return resistivity









