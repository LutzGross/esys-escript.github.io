#!/usr/bin/python3
from esys.escript import *
from esys.finley import ReadGmsh
import unittest

class Test_Finley_GMSHReader(unittest.TestCase):
    ELECTRODES = [(-64., 0., 0.), (-32., 0., 0.), ( 0.0, 0., 0.), ( 32., 0., 0.), ( 64., 0., 0.)]
    DIRACTAGS = ["e"+str(st) for st in range(5)]
    VOL=2000000.
    AREA=100000.
    DIRACCOUNT=5.0
    def test_Version22_addDirac(self):
        """
        mesh version 2.2 and addition of Dirac point when reading
        it is known that version 2.2 does not contain "Physical Point"
        
        """
        MSHFILE="meshDirac_22.msh"
        domain=ReadGmsh(MSHFILE, 3, optimize=True, diracPoints = self.ELECTRODES, diracTags=self.DIRACTAGS)
        self.checkDomain(domain)
        self.checkTags(domain)
        

    def test_Version41(self):
        MSHFILE="meshDirac_41.msh"
        domain=ReadGmsh(MSHFILE, 3, optimize=True)
        self.checkDomain(domain)
        self.checkTags(domain)

    def test_Version41_addDirac(self):
        MSHFILE="meshDirac_41.msh"
        domain=ReadGmsh(MSHFILE, 3, optimize=True,  diracPoints = self.ELECTRODES, diracTags=self.DIRACTAGS)
        self.checkDomain(domain)
        self.checkTags(domain)


    def test_Version42(self):
        MSHFILE="meshDirac_42.msh"
        domain=ReadGmsh(MSHFILE, 3, optimize=True)
        self.checkDomain(domain)
        self.checkTags(domain)

    def test_Version42_addDirac(self):
        MSHFILE="meshDirac_42.msh"
        domain=ReadGmsh(MSHFILE, 3, optimize=True,  diracPoints = self.ELECTRODES, diracTags=self.DIRACTAGS)
        self.checkDomain(domain)
        self.checkTags(domain)


    def test_Version44(self):

        MSHFILE="meshDirac_44.msh"
        domain=ReadGmsh(MSHFILE, 3, optimize=True)
        self.checkDomain(domain)
        self.checkTags(domain)


    def test_Version44_addDirac(self):

        MSHFILE="meshDirac_44.msh"
        domain=ReadGmsh(MSHFILE, 3, optimize=True,  diracPoints = self.ELECTRODES, diracTags=self.DIRACTAGS)
        self.checkDomain(domain)
        self.checkTags(domain)

    def checkDomain(self, domain):
        one=Scalar(1., Function(domain))
        V=integrate(one)
        self.assertAlmostEqual(V,  self.VOL, msg="wrong volume")
        
        one=Scalar(1., FunctionOnBoundary(domain))
        A=integrate(one)
        self.assertAlmostEqual(A,  self.AREA, msg="wrong surface area")

        one=Scalar(1., DiracDeltaFunctions(domain))
        D=integrate(one)
        self.assertAlmostEqual(D,  self.DIRACCOUNT, msg="wrong dirac point count")

    def checkTags(self, domain, alternative_dirac_tags=None):
        
        mesh_tags=getTagNames(domain)
        self.assertListEqual(mesh_tags, ['back', 'bottom', 'e0', 'e1', 'e2', 'e3', 'e4', 'front', 'left', 'right', 'top', 'volume'])       
        
        tagsF=Function(domain).getListOfTags() 
        tagsF.sort()
        self.assertListEqual(tagsF, [1], 'wrong volume tags')

        tagsFBC=FunctionOnBoundary(domain).getListOfTags()
        tagsFBC.sort()
        self.assertListEqual(tagsFBC, [2, 3, 4, 5, 6, 7], 'wrong surface tags')

        tagsD=DiracDeltaFunctions(domain).getListOfTags() 
        tagsD.sort()
        # this tag ids could be different when tags are added while reading
        if alternative_dirac_tags:
            self.assertListEqual(tagsD, alternative_dirac_tags, 'wrong Dirac point tags')
        else:    
            self.assertListEqual(tagsD, [8, 9, 10, 11, 12], 'wrong Dirac point tags')
        
        

if __name__ == '__main__':
    unittest.main()
