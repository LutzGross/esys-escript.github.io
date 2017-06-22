
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

#here be nan tests
import esys.escript as es
import esys.escriptcore.utestselect as unittest

class Test_util_NaN_funcs(unittest.TestCase):
    RES_TOL=1.e-7 # RES_TOLerance to compare results
    DIFF_TOL=1.e-7 # RES_TOLerance to derivatives

    def test_replaceNaNExpanded(self):
        dom=self.domain
        scl=es.Scalar(0,es.ContinuousFunction(dom))
        scl.expand()
        sclNaN=scl/0
        self.assertTrue(sclNaN.hasNaN(),"sclNaN should contain NaN but its doesn't")
        sclNaN.replaceNaN(15.0)
        self.assertEqual(es.Lsup(sclNaN), 15.0)
        scl=es.Scalar(0,es.ContinuousFunction(dom))
        scl.expand()
        scl.promote()
        if not es.getEscriptParamInt('AUTOLAZY')==1:            
            sclNaN=scl/0
            self.assertTrue(sclNaN.hasNaN(),"sclNaN should contain NaN but its doesn't")
            sclNaN.replaceNaN(3+4j)
            self.assertEqual(es.Lsup(sclNaN), 5.0)

    
    def test_replaceNaNConstant(self):
        dom=self.domain
        dat = es.Data(10,es.ContinuousFunction(dom))
        dat=(dat*0)/0
        self.assertTrue(dat.hasNaN(),"dat should contain NaN but its doesn't")
        dat.replaceNaN(10)
        self.assertEqual(es.Lsup(dat), 10)
        dat = es.Data(10,es.ContinuousFunction(dom))
        dat.promote()
        dat=(dat*0)/0
        self.assertTrue(dat.hasNaN(),"dat should contain NaN but its doesn't")
        dat.replaceNaN(4+3j)
        self.assertEqual(es.Lsup(dat), 5)

    def test_replaceNaNTagged(self):
        dom=self.domain
        sigma = es.Scalar(0,es.FunctionOnBoundary(dom))
        dat=(sigma*0)/0
        sigma.setTaggedValue(1 , es.Lsup(dat))
        sigma.replaceNaN(10)
        self.assertEqual(es.Lsup(sigma), 10)
        sigma = es.Scalar(0,es.FunctionOnBoundary(dom))
        sigma.promote()
        dat=(sigma*0)/0
        sigma.setTaggedValue(1 , es.Lsup(dat))
        sigma.replaceNaN(3+4j)
        self.assertEqual(es.Lsup(sigma), 5)

    def test_replaceInfExpanded(self):
        dom=self.domain
        scl=es.Scalar(1,es.ContinuousFunction(dom))
        scl.expand()
        sclNaN=scl/0
        self.assertTrue(sclNaN.hasInf(),"sclNaN should contain Inf but its doesn't")
        sclNaN.replaceNaN(17.0) # Make sure it is distinguishing between NaN and Inf
        sclNaN.replaceInf(15.0)
        self.assertEqual(es.Lsup(sclNaN), 15.0)
        scl=es.Scalar(10,es.ContinuousFunction(dom))
        scl.expand()
        scl.promote()
        if not es.getEscriptParamInt('AUTOLAZY')==1:            
            sclNaN=scl/0
            self.assertTrue(sclNaN.hasInf(),"sclNaN should contain Inf but its doesn't")
            sclNaN.replaceInf(3+4j)
            sclNaN.replaceNaN(3+19j)    # Make sure it is distinguishing between NaN and Inf
            self.assertEqual(es.Lsup(sclNaN), 5.0)

    
    def test_replaceInfConstant(self):
        dom=self.domain
        dat = es.Data(10,es.ContinuousFunction(dom))
        dat=(dat*1)/0
        self.assertTrue(dat.hasInf(),"dat should contain Inf but its doesn't")
        dat.replaceNaN(18)
        dat.replaceInf(10)
        self.assertEqual(es.Lsup(dat), 10)
        dat = es.Data(10,es.ContinuousFunction(dom))
        dat.promote()
        dat=(dat*1)/0
        self.assertTrue(dat.hasInf(),"dat should contain Inf but its doesn't")
        dat.replaceInf(4+3j)
        dat.replaceNaN(19+3j)
        self.assertEqual(es.Lsup(dat), 5)

    def test_replaceInfTagged(self):
        dom=self.domain
        sigma = es.Scalar(1,es.FunctionOnBoundary(dom))
        dat=(sigma*1)/0
        sigma.setTaggedValue(1 , es.Lsup(dat))
        sigma.replaceNaN(11)
        sigma.replaceInf(10)
        self.assertEqual(es.Lsup(sigma), 10)
        sigma = es.Scalar(1,es.FunctionOnBoundary(dom))
        dat=(sigma*1)/0
        dat.promote()
        sigma.setTaggedValue(1 , es.Lsup(dat))
        sigma.replaceInf(3+4j)
        sigma.replaceNaN(7+12j)
        self.assertEqual(es.Lsup(sigma), 5)


