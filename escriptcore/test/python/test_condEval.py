
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

__author__="Joel Fenwick, joelfenwick@uq.edu.au"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
import os
import numpy
import sys
from esys.escript import *

class Test_CondEval(unittest.TestCase):
   RES_TOL=1.e-7 # RES_TOLerance to compare results     
   r=getTestDomainFunctionSpace(3,3)
   x=r.getDomain().getX()

   def test_Constant(self):
        d1=Data((1,3),self.x.getFunctionSpace())
        d2=Data((2,4),self.x.getFunctionSpace())
        m1=Data(1,self.x.getFunctionSpace())
        mm1=Data(-1,self.x.getFunctionSpace())
        self.assertTrue(Lsup(condEval(m1,d1,d2)-(1,3))<self.RES_TOL)
        self.assertTrue(Lsup(condEval(mm1,d1,d2)-(2,4))<self.RES_TOL)
        #Now we try lazy
        d1=delay(d1)
        d2=delay(d2)
        m1=delay(m1)
        mm1=delay(mm1)
        self.assertTrue(Lsup(condEval(m1,d1,d2)-(1,3))<self.RES_TOL)
        self.assertTrue(Lsup(condEval(mm1,d1,d2)-(2,4))<self.RES_TOL)

   def test_Tagged(self):
        t1=Data((1,3,5),self.x.getFunctionSpace())
        t2=Data((2,4,6),self.x.getFunctionSpace())
        
        t1.tag()
        t2.tag()
        
        t1.setTaggedValue(1, (9,8,7))
        t1.setTaggedValue(2, (0,0,0))
        
        t2.setTaggedValue(1, (-1,-1,-1))
        t2.setTaggedValue(3, (-2,-2,-2))
        
        mt1=Data(1,self.x.getFunctionSpace())
        mt1.tag()
        mt1.setTaggedValue(1,1)
        mt1.setTaggedValue(4,0)
        mt2=Data(-1,self.x.getFunctionSpace())
        mt2.tag()
        mt2.setTaggedValue(1,0)
        mt2.setTaggedValue(4,-1)
        z=condEval(mt1,t1,t2)
        res=Data((1,3,5), self.x.getFunctionSpace())
        res.setTaggedValue(1,(9,8,7))
        res.setTaggedValue(4,(2,4,6))
        self.assertTrue(Lsup(z-res)<self.RES_TOL)
        z=condEval(mt2,t1,t2)
        res=Data((2,4,6), self.x.getFunctionSpace())
        res.setTaggedValue(1,(-1,-1,-1))
        res.setTaggedValue(4,(2,4,6))
        self.assertTrue(Lsup(z-res)<self.RES_TOL)
        
        #Now we try the same but lazy
        mt1=delay(mt1)
        mt2=delay(mt2)
        t1=delay(t1)
        t2=delay(t2)
        res=Data((1,3,5), self.x.getFunctionSpace())
        res.setTaggedValue(1,(9,8,7))
        res.setTaggedValue(4,(2,4,6))
        res=delay(res)
        
        z=condEval(mt1,t1,t2)
        y=z-res
        Lsup(y)
        self.assertTrue(Lsup(condEval(mt1,t1,t2)-res)<self.RES_TOL)
        res=Data((2,4,6), self.x.getFunctionSpace())
        res.setTaggedValue(1,(-1,-1,-1))
        res.setTaggedValue(4,(2,4,6))
        self.assertTrue(Lsup(condEval(mt2,t1,t2)-res)<self.RES_TOL)
        
   def test_Expanded(self):
        e1=Data((1,3,5),self.x.getFunctionSpace(),True)
        e2=Data((2,4,6),self.x.getFunctionSpace(),False)
        
        me1=Data(1,self.x.getFunctionSpace())*wherePositive(self.x-2)
        me2=Data(1,self.x.getFunctionSpace())*(1-wherePositive(self.x-2))
        
        self.assertTrue(Lsup(condEval(me1, e1,e2)+condEval(me2, e1,e2)-(3,7,11))<self.RES_TOL)
        
        le1=delay(e1)
        le2=delay(e2)
        
        ml1=delay(me1)
        ml2=delay(me2)
        
        z=condEval(ml1, le1,le2)
        
        self.assertTrue(Lsup(condEval(ml1, le1,le2)+condEval(ml2, le1,le2)-(3,7,11))<self.RES_TOL)
        
   def test_Errors(self):
        d1=Data((1,3),self.x.getFunctionSpace())
        d2=Data((2,4),self.x.getFunctionSpace())
        d3=Data((2,4,5),self.x.getFunctionSpace())
        m=Data(1,self.x.getFunctionSpace())
        mS1=Data((1,0),self.x.getFunctionSpace())
        #Non-scalar mask
        self.assertRaises(RuntimeError,condEval, mS1,d1,d2)
        #shape mismatch
        self.assertRaises(RuntimeError,condEval, m, d1, d3)

   def test_promote(self):      
        #This is not an exhaustive test of all possible promotion combinataions
        for v in [False, True]:
                mt1=Data(2,self.x.getFunctionSpace())
                mt1.tag()
                if v:
                    mt1=delay(mt1)
                d1=Data((1,3),self.x.getFunctionSpace())
                d2=Data((67,89), self.x.getFunctionSpace(),True)
                self.assertTrue(Lsup(condEval(mt1,d1,d2)-(1,3))<self.RES_TOL)
                
                me1=Data(2,self.x.getFunctionSpace(),True)
                d1=Data((1,3),self.x.getFunctionSpace())
                if v:
                    d1=delay(d1)
                d2=Data((7,19),self.x.getFunctionSpace())
                d2.tag()
                self.assertTrue(Lsup(condEval(me1,d1,d2)-(1,3))<self.RES_TOL)
                
                mc1=Data(0,self.x.getFunctionSpace())
                d1=Data((1,3),self.x.getFunctionSpace())
                d1.tag()
                d2=Data((7,19),self.x.getFunctionSpace())
                d2.tag()
                if v:
                    d2=delay(d2)
                self.assertTrue(Lsup(condEval(mc1,d1,d2)-(7,19))<self.RES_TOL)  
                
                mt1=Data(2,self.x.getFunctionSpace())
                mt1.tag()
                if v:
                    mt1=delay(mt1)
                d1=Data((1,3),self.x.getFunctionSpace())
                d2=Data((67,89), self.x.getFunctionSpace())
                self.assertTrue(Lsup(condEval(mt1,d1,d2)-(1,3))<self.RES_TOL)
        

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
