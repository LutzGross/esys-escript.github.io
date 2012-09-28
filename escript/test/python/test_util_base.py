
##############################################################################
#
# Copyright (c) 2003-2012 by University of Queensland
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

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
basic tests for util.py

:remark: use see `test_util`

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

import unittest
import numpy
import os
from esys.escript import *
from esys import escript

try:
     ESCRIPT_WORKDIR=os.environ['ESCRIPT_WORKDIR']
except KeyError:
     ESCRIPT_WORKDIR='.'

class Test_util_base(unittest.TestCase):
   """
   basic tests on util.py
   """
   RES_TOL=1.e-7 # RES_TOLerance to compare results
   DIFF_TOL=1.e-7 # RES_TOLerance to derivatices
#=========================================================
#  File writer
#=========================================================
   def __checkContent(self,fn,ref_cont):
        cont=open(fn,'r').readlines()
        self.assertTrue(len(cont)==len(ref_cont),"wrong number of records")
        for i in range(len(cont)):
           self.assertTrue(cont[i].strip()==ref_cont[i],"wrong records %s"%i)
   def test_FileWriter_W(self):
        fn=os.path.join(ESCRIPT_WORKDIR, "filewriter_w.txt")
        self.assertRaises(IOError,FileWriter,fn="",append=False)
        f=FileWriter(fn,append=False)
        self.assertTrue(f.name==fn, "wrong file name.")
        self.assertTrue(f.mode=='w', "wrong mode")
        self.assertTrue(f.newlines==os.linesep, "wrong line seps")
        self.assertTrue(not f.closed,"file shuold not be closed.")
        f.write("line1"+f.newlines)
        f.flush()
        self.__checkContent(fn,["line1"])
        f.writelines(["line2"+f.newlines, "line3"+f.newlines])
        f.close()
        self.assertTrue(f.closed,"file shuold be closed.")
        self.__checkContent(fn,["line1", "line2", "line3"])

   def test_FileWriter_A(self):
        fn=os.path.join(ESCRIPT_WORKDIR, "filewriter_a.txt")
        if getMPIRankWorld()==0: open(fn,'w').write("line1"+os.linesep)
        self.assertRaises(IOError,FileWriter,fn="",append=True)
        f=FileWriter(fn,append=True)
        self.assertTrue(f.name==fn, "wrong file name.")
        self.assertTrue(f.mode=='a', "wrong mode")
        self.assertTrue(f.newlines==os.linesep, "wrong line seps")
        self.assertTrue(not f.closed,"file shuold not be closed.")
        f.write("line2"+f.newlines)
        f.flush()
        self.__checkContent(fn,["line1", "line2"])
        f.writelines(["line3"+f.newlines, "line4"+f.newlines])
        f.close()
        self.assertTrue(f.closed,"file shuold be closed.")
        self.__checkContent(fn,["line1", "line2", "line3", "line4"])

   def test_FileWriter_A_loc(self):
        fn=os.path.join(ESCRIPT_WORKDIR, "filewriter_a_loc.txt")
        if getMPIRankWorld()>0:
            fn2=fn+".%s"%getMPIRankWorld()
        else:
            fn2=fn
        open(fn2,'w').write("line1"+os.linesep)
        self.assertRaises(IOError,FileWriter,fn="",append=True, createLocalFiles=True)
        f=FileWriter(fn,append=True,createLocalFiles=True)
        self.assertTrue(f.name==fn, "wrong file name.")
        self.assertTrue(f.mode=='a', "wrong mode")
        self.assertTrue(f.newlines==os.linesep, "wrong line seps")
        self.assertTrue(not f.closed,"file shuold not be closed.")
        f.write("line2"+f.newlines)
        f.flush()
        self.__checkContent(fn2,["line1", "line2"])
        f.writelines(["line3"+f.newlines, "line4"+f.newlines])
        f.close()
        self.assertTrue(f.closed,"file shuold be closed.")
        self.__checkContent(fn2,["line1", "line2", "line3", "line4"])

   def test_FileWriter_W_loc(self):
        fn=os.path.join(ESCRIPT_WORKDIR, "filewriter_w_loc.txt")
        if getMPIRankWorld()>0:
            fn2=fn+".%s"%getMPIRankWorld()
        else:
            fn2=fn
        self.assertRaises(IOError,FileWriter,fn="",append=True, createLocalFiles=True)
        f=FileWriter(fn,append=False,createLocalFiles=True)
        self.assertTrue(f.name==fn, "wrong file name.")
        self.assertTrue(f.mode=='w', "wrong mode")
        self.assertTrue(f.newlines==os.linesep, "wrong line seps")
        self.assertTrue(not f.closed,"file shuold not be closed.")
        f.write("line1"+f.newlines)
        f.flush()
        self.__checkContent(fn2,["line1"])
        f.writelines(["line2"+f.newlines, "line3"+f.newlines])
        f.close()
        self.assertTrue(f.closed,"file shuold be closed.")
        self.__checkContent(fn2,["line1", "line2", "line3"])
#=========================================================
#  constants
#=========================================================
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_kronecker_1(self):
      val=kronecker(d=1)
      self.assertTrue(isinstance(val,numpy.ndarray),"wrong type of result.")
      self.assertEqual(val[0,0],1.0,"wrong value for (0,0)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_kronecker_2(self):
      val=kronecker(d=2)
      self.assertTrue(isinstance(val,numpy.ndarray),"wrong type of result.")
      self.assertEqual(val[0,0],1.0,"wrong value for (0,0)")
      self.assertEqual(val[0,1],0.0,"wrong value for (0,1)")
      self.assertEqual(val[1,0],0.0,"wrong value for (1,0)")
      self.assertEqual(val[1,1],1.0,"wrong value for (1,1)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_kronecker_3(self):
      val=kronecker(d=3)
      self.assertTrue(isinstance(val,numpy.ndarray),"wrong type of result.")
      self.assertEqual(val[0,0],1.0,"wrong value for (0,0)")
      self.assertEqual(val[0,1],0.0,"wrong value for (0,1)")
      self.assertEqual(val[0,2],0.0,"wrong value for (0,2)")
      self.assertEqual(val[1,0],0.0,"wrong value for (1,0)")
      self.assertEqual(val[1,1],1.0,"wrong value for (1,1)")
      self.assertEqual(val[1,2],0.0,"wrong value for (1,2)")
      self.assertEqual(val[2,0],0.0,"wrong value for (2,0)")
      self.assertEqual(val[2,1],0.0,"wrong value for (2,1)")
      self.assertEqual(val[2,2],1.0,"wrong value for (2,2)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_kronecker_domain(self):
      val=kronecker(d=self.domain)
      self.assertTrue(isinstance(val,numpy.ndarray),"wrong type of result.")
      self.assertEqual(val.shape,(self.domain.getDim(),self.domain.getDim()),"wrong shape.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_kronecker_functionspace(self):
      val=kronecker(d=self.functionspace)
      self.assertTrue(isinstance(val,escript.Data),"wrong type of result.")
      self.assertEqual(val.getShape(),(self.functionspace.getDim(),self.functionspace.getDim()),"wrong shape.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor_1(self):
      val=identityTensor(d=1)
      self.assertTrue(isinstance(val,numpy.ndarray),"wrong type of result.")
      self.assertEqual(val[0,0],1.0,"wrong value for (0,0)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor_2(self):
      val=identityTensor(d=2)
      self.assertTrue(isinstance(val,numpy.ndarray),"wrong type of result.")
      self.assertEqual(val[0,0],1.0,"wrong value for (0,0)")
      self.assertEqual(val[0,1],0.0,"wrong value for (0,1)")
      self.assertEqual(val[1,0],0.0,"wrong value for (1,0)")
      self.assertEqual(val[1,1],1.0,"wrong value for (1,1)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor_3(self):
      val=identityTensor(d=3)
      self.assertTrue(isinstance(val,numpy.ndarray),"wrong type of result.")
      self.assertEqual(val[0,0],1.0,"wrong value for (0,0)")
      self.assertEqual(val[0,1],0.0,"wrong value for (0,1)")
      self.assertEqual(val[0,2],0.0,"wrong value for (0,2)")
      self.assertEqual(val[1,0],0.0,"wrong value for (1,0)")
      self.assertEqual(val[1,1],1.0,"wrong value for (1,1)")
      self.assertEqual(val[1,2],0.0,"wrong value for (1,2)")
      self.assertEqual(val[2,0],0.0,"wrong value for (2,0)")
      self.assertEqual(val[2,1],0.0,"wrong value for (2,1)")
      self.assertEqual(val[2,2],1.0,"wrong value for (2,2)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor_domain(self):
      val=identityTensor(d=self.domain)
      self.assertTrue(isinstance(val,numpy.ndarray),"wrong type of result.")
      self.assertEqual(val.shape,(self.domain.getDim(),self.domain.getDim()),"wrong shape.")
      if self.domain.getDim()==2:
         self.assertEqual(val[0,0],1.0,"wrong value for (0,0)")
         self.assertEqual(val[0,1],0.0,"wrong value for (0,1)")
         self.assertEqual(val[1,0],0.0,"wrong value for (1,0)")
         self.assertEqual(val[1,1],1.0,"wrong value for (1,1)")
      if self.domain.getDim()==3:
         self.assertEqual(val[0,0],1.0,"wrong value for (0,0)")
         self.assertEqual(val[0,1],0.0,"wrong value for (0,1)")
         self.assertEqual(val[0,2],0.0,"wrong value for (0,2)")
         self.assertEqual(val[1,0],0.0,"wrong value for (1,0)")
         self.assertEqual(val[1,1],1.0,"wrong value for (1,1)")
         self.assertEqual(val[1,2],0.0,"wrong value for (1,2)")
         self.assertEqual(val[2,0],0.0,"wrong value for (2,0)")
         self.assertEqual(val[2,1],0.0,"wrong value for (2,1)")
         self.assertEqual(val[2,2],1.0,"wrong value for (2,2)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor_functionspace(self):
      val=identityTensor(d=self.functionspace)
      self.assertTrue(isinstance(val,escript.Data),"wrong type of result.")
      self.assertEqual(val.getShape(),(self.functionspace.getDim(),self.functionspace.getDim()),"wrong shape.")
      if self.domain.getDim()==2:
         self.assertTrue(Lsup(val[0,0]-1.0)<=self.RES_TOL,"wrong value for (0,0)")
         self.assertTrue(Lsup(val[0,1]-0.0)<=self.RES_TOL,"wrong value for (0,1)")
         self.assertTrue(Lsup(val[1,0]-0.0)<=self.RES_TOL,"wrong value for (1,0)")
         self.assertTrue(Lsup(val[1,1]-1.0)<=self.RES_TOL,"wrong value for (1,1)")
      if self.domain.getDim()==3:
         self.assertTrue(Lsup(val[0,0]-1.0)<=self.RES_TOL,"wrong value for (0,0)")
         self.assertTrue(Lsup(val[0,1]-0.0)<=self.RES_TOL,"wrong value for (0,1)")
         self.assertTrue(Lsup(val[0,2]-0.0)<=self.RES_TOL,"wrong value for (0,2)")
         self.assertTrue(Lsup(val[1,0]-0.0)<=self.RES_TOL,"wrong value for (1,0)")
         self.assertTrue(Lsup(val[1,1]-1.0)<=self.RES_TOL,"wrong value for (1,1)")
         self.assertTrue(Lsup(val[1,2]-0.0)<=self.RES_TOL,"wrong value for (1,2)")
         self.assertTrue(Lsup(val[2,0]-0.0)<=self.RES_TOL,"wrong value for (2,0)")
         self.assertTrue(Lsup(val[2,1]-0.0)<=self.RES_TOL,"wrong value for (2,1)")
         self.assertTrue(Lsup(val[2,2]-1.0)<=self.RES_TOL,"wrong value for (2,2)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor4_1(self):
      val=identityTensor4(d=1)
      self.assertTrue(isinstance(val,numpy.ndarray),"wrong type of result.")
      self.assertEqual(val[0,0,0,0],1.0,"wrong value for (0,0,0,0)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor4_2(self):
      val=identityTensor4(d=2)
      self.assertTrue(isinstance(val,numpy.ndarray),"wrong type of result.")
      self.assertEqual(val[0,0,0,0],1.0,"wrong value for (0,0,0,0)")
      self.assertEqual(val[0,0,0,1],0.0,"wrong value for (0,0,0,1)")
      self.assertEqual(val[0,0,1,0],0.0,"wrong value for (0,0,1,0)")
      self.assertEqual(val[0,0,1,1],0.0,"wrong value for (0,0,1,1)")
      self.assertEqual(val[0,1,0,0],0.0,"wrong value for (0,1,0,0)")
      self.assertEqual(val[0,1,0,1],1.0,"wrong value for (0,1,0,1)")
      self.assertEqual(val[0,1,1,0],0.0,"wrong value for (0,1,1,0)")
      self.assertEqual(val[0,1,1,1],0.0,"wrong value for (0,1,1,1)")
      self.assertEqual(val[1,0,0,0],0.0,"wrong value for (1,0,0,0)")
      self.assertEqual(val[1,0,0,1],0.0,"wrong value for (1,0,0,1)")
      self.assertEqual(val[1,0,1,0],1.0,"wrong value for (1,0,1,0)")
      self.assertEqual(val[1,0,1,1],0.0,"wrong value for (1,0,1,1)")
      self.assertEqual(val[1,1,0,0],0.0,"wrong value for (1,1,0,0)")
      self.assertEqual(val[1,1,0,1],0.0,"wrong value for (1,1,0,1)")
      self.assertEqual(val[1,1,1,0],0.0,"wrong value for (1,1,1,0)")
      self.assertEqual(val[1,1,1,1],1.0,"wrong value for (1,1,1,1)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor4_3(self):
      val=identityTensor4(d=3)
      self.assertTrue(isinstance(val,numpy.ndarray),"wrong type of result.")
      self.assertEqual(val[0,0,0,0],1.0,"wrong value for (0,0,0,0)")
      self.assertEqual(val[0,0,0,1],0.0,"wrong value for (0,0,0,1)")
      self.assertEqual(val[0,0,0,2],0.0,"wrong value for (0,0,0,2)")
      self.assertEqual(val[0,0,1,0],0.0,"wrong value for (0,0,1,0)")
      self.assertEqual(val[0,0,1,1],0.0,"wrong value for (0,0,1,1)")
      self.assertEqual(val[0,0,1,2],0.0,"wrong value for (0,0,1,2)")
      self.assertEqual(val[0,0,2,0],0.0,"wrong value for (0,0,2,0)")
      self.assertEqual(val[0,0,2,1],0.0,"wrong value for (0,0,2,1)")
      self.assertEqual(val[0,0,2,2],0.0,"wrong value for (0,0,2,2)")
      self.assertEqual(val[0,1,0,0],0.0,"wrong value for (0,1,0,0)")
      self.assertEqual(val[0,1,0,1],1.0,"wrong value for (0,1,0,1)")
      self.assertEqual(val[0,1,0,2],0.0,"wrong value for (0,1,0,2)")
      self.assertEqual(val[0,1,1,0],0.0,"wrong value for (0,1,1,0)")
      self.assertEqual(val[0,1,1,1],0.0,"wrong value for (0,1,1,1)")
      self.assertEqual(val[0,1,1,2],0.0,"wrong value for (0,1,1,2)")
      self.assertEqual(val[0,1,2,0],0.0,"wrong value for (0,1,2,0)")
      self.assertEqual(val[0,1,2,1],0.0,"wrong value for (0,1,2,1)")
      self.assertEqual(val[0,1,2,2],0.0,"wrong value for (0,1,2,2)")
      self.assertEqual(val[0,2,0,0],0.0,"wrong value for (0,2,0,0)")
      self.assertEqual(val[0,2,0,1],0.0,"wrong value for (0,2,0,1)")
      self.assertEqual(val[0,2,0,2],1.0,"wrong value for (0,2,0,2)")
      self.assertEqual(val[0,2,1,0],0.0,"wrong value for (0,2,1,0)")
      self.assertEqual(val[0,2,1,1],0.0,"wrong value for (0,2,1,1)")
      self.assertEqual(val[0,2,1,2],0.0,"wrong value for (0,2,1,2)")
      self.assertEqual(val[0,2,2,0],0.0,"wrong value for (0,2,2,0)")
      self.assertEqual(val[0,2,2,1],0.0,"wrong value for (0,2,2,1)")
      self.assertEqual(val[0,2,2,2],0.0,"wrong value for (0,2,2,2)")
      self.assertEqual(val[1,0,0,0],0.0,"wrong value for (1,0,0,0)")
      self.assertEqual(val[1,0,0,1],0.0,"wrong value for (1,0,0,1)")
      self.assertEqual(val[1,0,0,2],0.0,"wrong value for (1,0,0,2)")
      self.assertEqual(val[1,0,1,0],1.0,"wrong value for (1,0,1,0)")
      self.assertEqual(val[1,0,1,1],0.0,"wrong value for (1,0,1,1)")
      self.assertEqual(val[1,0,1,2],0.0,"wrong value for (1,0,1,2)")
      self.assertEqual(val[1,0,2,0],0.0,"wrong value for (1,0,2,0)")
      self.assertEqual(val[1,0,2,1],0.0,"wrong value for (1,0,2,1)")
      self.assertEqual(val[1,0,2,2],0.0,"wrong value for (1,0,2,2)")
      self.assertEqual(val[1,1,0,0],0.0,"wrong value for (1,1,0,0)")
      self.assertEqual(val[1,1,0,1],0.0,"wrong value for (1,1,0,1)")
      self.assertEqual(val[1,1,0,2],0.0,"wrong value for (1,1,0,2)")
      self.assertEqual(val[1,1,1,0],0.0,"wrong value for (1,1,1,0)")
      self.assertEqual(val[1,1,1,1],1.0,"wrong value for (1,1,1,1)")
      self.assertEqual(val[1,1,1,2],0.0,"wrong value for (1,1,1,2)")
      self.assertEqual(val[1,1,2,0],0.0,"wrong value for (1,1,2,0)")
      self.assertEqual(val[1,1,2,1],0.0,"wrong value for (1,1,2,1)")
      self.assertEqual(val[1,1,2,2],0.0,"wrong value for (1,1,2,2)")
      self.assertEqual(val[1,2,0,0],0.0,"wrong value for (1,2,0,0)")
      self.assertEqual(val[1,2,0,1],0.0,"wrong value for (1,2,0,1)")
      self.assertEqual(val[1,2,0,2],0.0,"wrong value for (1,2,0,2)")
      self.assertEqual(val[1,2,1,0],0.0,"wrong value for (1,2,1,0)")
      self.assertEqual(val[1,2,1,1],0.0,"wrong value for (1,2,1,1)")
      self.assertEqual(val[1,2,1,2],1.0,"wrong value for (1,2,1,2)")
      self.assertEqual(val[1,2,2,0],0.0,"wrong value for (1,2,2,0)")
      self.assertEqual(val[1,2,2,1],0.0,"wrong value for (1,2,2,1)")
      self.assertEqual(val[1,2,2,2],0.0,"wrong value for (1,2,2,2)")
      self.assertEqual(val[2,0,0,0],0.0,"wrong value for (2,0,0,0)")
      self.assertEqual(val[2,0,0,1],0.0,"wrong value for (2,0,0,1)")
      self.assertEqual(val[2,0,0,2],0.0,"wrong value for (2,0,0,2)")
      self.assertEqual(val[2,0,1,0],0.0,"wrong value for (2,0,1,0)")
      self.assertEqual(val[2,0,1,1],0.0,"wrong value for (2,0,1,1)")
      self.assertEqual(val[2,0,1,2],0.0,"wrong value for (2,0,1,2)")
      self.assertEqual(val[2,0,2,0],1.0,"wrong value for (2,0,2,0)")
      self.assertEqual(val[2,0,2,1],0.0,"wrong value for (2,0,2,1)")
      self.assertEqual(val[2,0,2,2],0.0,"wrong value for (2,0,2,2)")
      self.assertEqual(val[2,1,0,0],0.0,"wrong value for (2,1,0,0)")
      self.assertEqual(val[2,1,0,1],0.0,"wrong value for (2,1,0,1)")
      self.assertEqual(val[2,1,0,2],0.0,"wrong value for (2,1,0,2)")
      self.assertEqual(val[2,1,1,0],0.0,"wrong value for (2,1,1,0)")
      self.assertEqual(val[2,1,1,1],0.0,"wrong value for (2,1,1,1)")
      self.assertEqual(val[2,1,1,2],0.0,"wrong value for (2,1,1,2)")
      self.assertEqual(val[2,1,2,0],0.0,"wrong value for (2,1,2,0)")
      self.assertEqual(val[2,1,2,1],1.0,"wrong value for (2,1,2,1)")
      self.assertEqual(val[2,1,2,2],0.0,"wrong value for (2,1,2,2)")
      self.assertEqual(val[2,2,0,0],0.0,"wrong value for (2,2,0,0)")
      self.assertEqual(val[2,2,0,1],0.0,"wrong value for (2,2,0,1)")
      self.assertEqual(val[2,2,0,2],0.0,"wrong value for (2,2,0,2)")
      self.assertEqual(val[2,2,1,0],0.0,"wrong value for (2,2,1,0)")
      self.assertEqual(val[2,2,1,1],0.0,"wrong value for (2,2,1,1)")
      self.assertEqual(val[2,2,1,2],0.0,"wrong value for (2,2,1,2)")
      self.assertEqual(val[2,2,2,0],0.0,"wrong value for (2,2,2,0)")
      self.assertEqual(val[2,2,2,1],0.0,"wrong value for (2,2,2,1)")
      self.assertEqual(val[2,2,2,2],1.0,"wrong value for (2,2,2,2)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor4_domain(self):
      val=identityTensor4(d=self.domain)
      self.assertTrue(isinstance(val,numpy.ndarray),"wrong type of result.")
      self.assertEqual(val.shape,(self.domain.getDim(),self.domain.getDim(),self.domain.getDim(),self.domain.getDim()),"wrong shape.")
      if self.domain.getDim()==2:
         self.assertEqual(val[0,0,0,0],1.0,"wrong value for (0,0,0,0)")
         self.assertEqual(val[0,0,0,1],0.0,"wrong value for (0,0,0,1)")
         self.assertEqual(val[0,0,1,0],0.0,"wrong value for (0,0,1,0)")
         self.assertEqual(val[0,0,1,1],0.0,"wrong value for (0,0,1,1)")
         self.assertEqual(val[0,1,0,0],0.0,"wrong value for (0,1,0,0)")
         self.assertEqual(val[0,1,0,1],1.0,"wrong value for (0,1,0,1)")
         self.assertEqual(val[0,1,1,0],0.0,"wrong value for (0,1,1,0)")
         self.assertEqual(val[0,1,1,1],0.0,"wrong value for (0,1,1,1)")
         self.assertEqual(val[1,0,0,0],0.0,"wrong value for (1,0,0,0)")
         self.assertEqual(val[1,0,0,1],0.0,"wrong value for (1,0,0,1)")
         self.assertEqual(val[1,0,1,0],1.0,"wrong value for (1,0,1,0)")
         self.assertEqual(val[1,0,1,1],0.0,"wrong value for (1,0,1,1)")
         self.assertEqual(val[1,1,0,0],0.0,"wrong value for (1,1,0,0)")
         self.assertEqual(val[1,1,0,1],0.0,"wrong value for (1,1,0,1)")
         self.assertEqual(val[1,1,1,0],0.0,"wrong value for (1,1,1,0)")
         self.assertEqual(val[1,1,1,1],1.0,"wrong value for (1,1,1,1)")
      if self.domain.getDim()==3:
         self.assertEqual(val[0,0,0,0],1.0,"wrong value for (0,0,0,0)")
         self.assertEqual(val[0,0,0,1],0.0,"wrong value for (0,0,0,1)")
         self.assertEqual(val[0,0,0,2],0.0,"wrong value for (0,0,0,2)")
         self.assertEqual(val[0,0,1,0],0.0,"wrong value for (0,0,1,0)")
         self.assertEqual(val[0,0,1,1],0.0,"wrong value for (0,0,1,1)")
         self.assertEqual(val[0,0,1,2],0.0,"wrong value for (0,0,1,2)")
         self.assertEqual(val[0,0,2,0],0.0,"wrong value for (0,0,2,0)")
         self.assertEqual(val[0,0,2,1],0.0,"wrong value for (0,0,2,1)")
         self.assertEqual(val[0,0,2,2],0.0,"wrong value for (0,0,2,2)")
         self.assertEqual(val[0,1,0,0],0.0,"wrong value for (0,1,0,0)")
         self.assertEqual(val[0,1,0,1],1.0,"wrong value for (0,1,0,1)")
         self.assertEqual(val[0,1,0,2],0.0,"wrong value for (0,1,0,2)")
         self.assertEqual(val[0,1,1,0],0.0,"wrong value for (0,1,1,0)")
         self.assertEqual(val[0,1,1,1],0.0,"wrong value for (0,1,1,1)")
         self.assertEqual(val[0,1,1,2],0.0,"wrong value for (0,1,1,2)")
         self.assertEqual(val[0,1,2,0],0.0,"wrong value for (0,1,2,0)")
         self.assertEqual(val[0,1,2,1],0.0,"wrong value for (0,1,2,1)")
         self.assertEqual(val[0,1,2,2],0.0,"wrong value for (0,1,2,2)")
         self.assertEqual(val[0,2,0,0],0.0,"wrong value for (0,2,0,0)")
         self.assertEqual(val[0,2,0,1],0.0,"wrong value for (0,2,0,1)")
         self.assertEqual(val[0,2,0,2],1.0,"wrong value for (0,2,0,2)")
         self.assertEqual(val[0,2,1,0],0.0,"wrong value for (0,2,1,0)")
         self.assertEqual(val[0,2,1,1],0.0,"wrong value for (0,2,1,1)")
         self.assertEqual(val[0,2,1,2],0.0,"wrong value for (0,2,1,2)")
         self.assertEqual(val[0,2,2,0],0.0,"wrong value for (0,2,2,0)")
         self.assertEqual(val[0,2,2,1],0.0,"wrong value for (0,2,2,1)")
         self.assertEqual(val[0,2,2,2],0.0,"wrong value for (0,2,2,2)")
         self.assertEqual(val[1,0,0,0],0.0,"wrong value for (1,0,0,0)")
         self.assertEqual(val[1,0,0,1],0.0,"wrong value for (1,0,0,1)")
         self.assertEqual(val[1,0,0,2],0.0,"wrong value for (1,0,0,2)")
         self.assertEqual(val[1,0,1,0],1.0,"wrong value for (1,0,1,0)")
         self.assertEqual(val[1,0,1,1],0.0,"wrong value for (1,0,1,1)")
         self.assertEqual(val[1,0,1,2],0.0,"wrong value for (1,0,1,2)")
         self.assertEqual(val[1,0,2,0],0.0,"wrong value for (1,0,2,0)")
         self.assertEqual(val[1,0,2,1],0.0,"wrong value for (1,0,2,1)")
         self.assertEqual(val[1,0,2,2],0.0,"wrong value for (1,0,2,2)")
         self.assertEqual(val[1,1,0,0],0.0,"wrong value for (1,1,0,0)")
         self.assertEqual(val[1,1,0,1],0.0,"wrong value for (1,1,0,1)")
         self.assertEqual(val[1,1,0,2],0.0,"wrong value for (1,1,0,2)")
         self.assertEqual(val[1,1,1,0],0.0,"wrong value for (1,1,1,0)")
         self.assertEqual(val[1,1,1,1],1.0,"wrong value for (1,1,1,1)")
         self.assertEqual(val[1,1,1,2],0.0,"wrong value for (1,1,1,2)")
         self.assertEqual(val[1,1,2,0],0.0,"wrong value for (1,1,2,0)")
         self.assertEqual(val[1,1,2,1],0.0,"wrong value for (1,1,2,1)")
         self.assertEqual(val[1,1,2,2],0.0,"wrong value for (1,1,2,2)")
         self.assertEqual(val[1,2,0,0],0.0,"wrong value for (1,2,0,0)")
         self.assertEqual(val[1,2,0,1],0.0,"wrong value for (1,2,0,1)")
         self.assertEqual(val[1,2,0,2],0.0,"wrong value for (1,2,0,2)")
         self.assertEqual(val[1,2,1,0],0.0,"wrong value for (1,2,1,0)")
         self.assertEqual(val[1,2,1,1],0.0,"wrong value for (1,2,1,1)")
         self.assertEqual(val[1,2,1,2],1.0,"wrong value for (1,2,1,2)")
         self.assertEqual(val[1,2,2,0],0.0,"wrong value for (1,2,2,0)")
         self.assertEqual(val[1,2,2,1],0.0,"wrong value for (1,2,2,1)")
         self.assertEqual(val[1,2,2,2],0.0,"wrong value for (1,2,2,2)")
         self.assertEqual(val[2,0,0,0],0.0,"wrong value for (2,0,0,0)")
         self.assertEqual(val[2,0,0,1],0.0,"wrong value for (2,0,0,1)")
         self.assertEqual(val[2,0,0,2],0.0,"wrong value for (2,0,0,2)")
         self.assertEqual(val[2,0,1,0],0.0,"wrong value for (2,0,1,0)")
         self.assertEqual(val[2,0,1,1],0.0,"wrong value for (2,0,1,1)")
         self.assertEqual(val[2,0,1,2],0.0,"wrong value for (2,0,1,2)")
         self.assertEqual(val[2,0,2,0],1.0,"wrong value for (2,0,2,0)")
         self.assertEqual(val[2,0,2,1],0.0,"wrong value for (2,0,2,1)")
         self.assertEqual(val[2,0,2,2],0.0,"wrong value for (2,0,2,2)")
         self.assertEqual(val[2,1,0,0],0.0,"wrong value for (2,1,0,0)")
         self.assertEqual(val[2,1,0,1],0.0,"wrong value for (2,1,0,1)")
         self.assertEqual(val[2,1,0,2],0.0,"wrong value for (2,1,0,2)")
         self.assertEqual(val[2,1,1,0],0.0,"wrong value for (2,1,1,0)")
         self.assertEqual(val[2,1,1,1],0.0,"wrong value for (2,1,1,1)")
         self.assertEqual(val[2,1,1,2],0.0,"wrong value for (2,1,1,2)")
         self.assertEqual(val[2,1,2,0],0.0,"wrong value for (2,1,2,0)")
         self.assertEqual(val[2,1,2,1],1.0,"wrong value for (2,1,2,1)")
         self.assertEqual(val[2,1,2,2],0.0,"wrong value for (2,1,2,2)")
         self.assertEqual(val[2,2,0,0],0.0,"wrong value for (2,2,0,0)")
         self.assertEqual(val[2,2,0,1],0.0,"wrong value for (2,2,0,1)")
         self.assertEqual(val[2,2,0,2],0.0,"wrong value for (2,2,0,2)")
         self.assertEqual(val[2,2,1,0],0.0,"wrong value for (2,2,1,0)")
         self.assertEqual(val[2,2,1,1],0.0,"wrong value for (2,2,1,1)")
         self.assertEqual(val[2,2,1,2],0.0,"wrong value for (2,2,1,2)")
         self.assertEqual(val[2,2,2,0],0.0,"wrong value for (2,2,2,0)")
         self.assertEqual(val[2,2,2,1],0.0,"wrong value for (2,2,2,1)")
         self.assertEqual(val[2,2,2,2],1.0,"wrong value for (2,2,2,2)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor4_functionspace(self):
      val=identityTensor4(d=self.functionspace)
      self.assertTrue(isinstance(val,escript.Data),"wrong type of result.")
      self.assertEqual(val.getShape(),(self.functionspace.getDim(),self.functionspace.getDim(),self.functionspace.getDim(),self.functionspace.getDim()),"wrong shape.")
      if self.domain.getDim()==2:
         self.assertTrue(Lsup(val[0,0,0,0]-1.0)<=self.RES_TOL,"wrong value for (0,0,0,0)")
         self.assertTrue(Lsup(val[0,0,0,1]-0.0)<=self.RES_TOL,"wrong value for (0,0,0,1)")
         self.assertTrue(Lsup(val[0,0,1,0]-0.0)<=self.RES_TOL,"wrong value for (0,0,1,0)")
         self.assertTrue(Lsup(val[0,0,1,1]-0.0)<=self.RES_TOL,"wrong value for (0,0,1,1)")
         self.assertTrue(Lsup(val[0,1,0,0]-0.0)<=self.RES_TOL,"wrong value for (0,1,0,0)")
         self.assertTrue(Lsup(val[0,1,0,1]-1.0)<=self.RES_TOL,"wrong value for (0,1,0,1)")
         self.assertTrue(Lsup(val[0,1,1,0]-0.0)<=self.RES_TOL,"wrong value for (0,1,1,0)")
         self.assertTrue(Lsup(val[0,1,1,1]-0.0)<=self.RES_TOL,"wrong value for (0,1,1,1)")
         self.assertTrue(Lsup(val[1,0,0,0]-0.0)<=self.RES_TOL,"wrong value for (1,0,0,0)")
         self.assertTrue(Lsup(val[1,0,0,1]-0.0)<=self.RES_TOL,"wrong value for (1,0,0,1)")
         self.assertTrue(Lsup(val[1,0,1,0]-1.0)<=self.RES_TOL,"wrong value for (1,0,1,0)")
         self.assertTrue(Lsup(val[1,0,1,1]-0.0)<=self.RES_TOL,"wrong value for (1,0,1,1)")
         self.assertTrue(Lsup(val[1,1,0,0]-0.0)<=self.RES_TOL,"wrong value for (1,1,0,0)")
         self.assertTrue(Lsup(val[1,1,0,1]-0.0)<=self.RES_TOL,"wrong value for (1,1,0,1)")
         self.assertTrue(Lsup(val[1,1,1,0]-0.0)<=self.RES_TOL,"wrong value for (1,1,1,0)")
         self.assertTrue(Lsup(val[1,1,1,1]-1.0)<=self.RES_TOL,"wrong value for (1,1,1,1)")
      if self.domain.getDim()==3:
         self.assertTrue(Lsup(val[0,0,0,0]-1.0)<=self.RES_TOL,"wrong value for (0,0,0,0)")
         self.assertTrue(Lsup(val[0,0,0,1]-0.0)<=self.RES_TOL,"wrong value for (0,0,0,1)")
         self.assertTrue(Lsup(val[0,0,0,2]-0.0)<=self.RES_TOL,"wrong value for (0,0,0,2)")
         self.assertTrue(Lsup(val[0,0,1,0]-0.0)<=self.RES_TOL,"wrong value for (0,0,1,0)")
         self.assertTrue(Lsup(val[0,0,1,1]-0.0)<=self.RES_TOL,"wrong value for (0,0,1,1)")
         self.assertTrue(Lsup(val[0,0,1,2]-0.0)<=self.RES_TOL,"wrong value for (0,0,1,2)")
         self.assertTrue(Lsup(val[0,0,2,0]-0.0)<=self.RES_TOL,"wrong value for (0,0,2,0)")
         self.assertTrue(Lsup(val[0,0,2,1]-0.0)<=self.RES_TOL,"wrong value for (0,0,2,1)")
         self.assertTrue(Lsup(val[0,0,2,2]-0.0)<=self.RES_TOL,"wrong value for (0,0,2,2)")
         self.assertTrue(Lsup(val[0,1,0,0]-0.0)<=self.RES_TOL,"wrong value for (0,1,0,0)")
         self.assertTrue(Lsup(val[0,1,0,1]-1.0)<=self.RES_TOL,"wrong value for (0,1,0,1)")
         self.assertTrue(Lsup(val[0,1,0,2]-0.0)<=self.RES_TOL,"wrong value for (0,1,0,2)")
         self.assertTrue(Lsup(val[0,1,1,0]-0.0)<=self.RES_TOL,"wrong value for (0,1,1,0)")
         self.assertTrue(Lsup(val[0,1,1,1]-0.0)<=self.RES_TOL,"wrong value for (0,1,1,1)")
         self.assertTrue(Lsup(val[0,1,1,2]-0.0)<=self.RES_TOL,"wrong value for (0,1,1,2)")
         self.assertTrue(Lsup(val[0,1,2,0]-0.0)<=self.RES_TOL,"wrong value for (0,1,2,0)")
         self.assertTrue(Lsup(val[0,1,2,1]-0.0)<=self.RES_TOL,"wrong value for (0,1,2,1)")
         self.assertTrue(Lsup(val[0,1,2,2]-0.0)<=self.RES_TOL,"wrong value for (0,1,2,2)")
         self.assertTrue(Lsup(val[0,2,0,0]-0.0)<=self.RES_TOL,"wrong value for (0,2,0,0)")
         self.assertTrue(Lsup(val[0,2,0,1]-0.0)<=self.RES_TOL,"wrong value for (0,2,0,1)")
         self.assertTrue(Lsup(val[0,2,0,2]-1.0)<=self.RES_TOL,"wrong value for (0,2,0,2)")
         self.assertTrue(Lsup(val[0,2,1,0]-0.0)<=self.RES_TOL,"wrong value for (0,2,1,0)")
         self.assertTrue(Lsup(val[0,2,1,1]-0.0)<=self.RES_TOL,"wrong value for (0,2,1,1)")
         self.assertTrue(Lsup(val[0,2,1,2]-0.0)<=self.RES_TOL,"wrong value for (0,2,1,2)")
         self.assertTrue(Lsup(val[0,2,2,0]-0.0)<=self.RES_TOL,"wrong value for (0,2,2,0)")
         self.assertTrue(Lsup(val[0,2,2,1]-0.0)<=self.RES_TOL,"wrong value for (0,2,2,1)")
         self.assertTrue(Lsup(val[0,2,2,2]-0.0)<=self.RES_TOL,"wrong value for (0,2,2,2)")
         self.assertTrue(Lsup(val[1,0,0,0]-0.0)<=self.RES_TOL,"wrong value for (1,0,0,0)")
         self.assertTrue(Lsup(val[1,0,0,1]-0.0)<=self.RES_TOL,"wrong value for (1,0,0,1)")
         self.assertTrue(Lsup(val[1,0,0,2]-0.0)<=self.RES_TOL,"wrong value for (1,0,0,2)")
         self.assertTrue(Lsup(val[1,0,1,0]-1.0)<=self.RES_TOL,"wrong value for (1,0,1,0)")
         self.assertTrue(Lsup(val[1,0,1,1]-0.0)<=self.RES_TOL,"wrong value for (1,0,1,1)")
         self.assertTrue(Lsup(val[1,0,1,2]-0.0)<=self.RES_TOL,"wrong value for (1,0,1,2)")
         self.assertTrue(Lsup(val[1,0,2,0]-0.0)<=self.RES_TOL,"wrong value for (1,0,2,0)")
         self.assertTrue(Lsup(val[1,0,2,1]-0.0)<=self.RES_TOL,"wrong value for (1,0,2,1)")
         self.assertTrue(Lsup(val[1,0,2,2]-0.0)<=self.RES_TOL,"wrong value for (1,0,2,2)")
         self.assertTrue(Lsup(val[1,1,0,0]-0.0)<=self.RES_TOL,"wrong value for (1,1,0,0)")
         self.assertTrue(Lsup(val[1,1,0,1]-0.0)<=self.RES_TOL,"wrong value for (1,1,0,1)")
         self.assertTrue(Lsup(val[1,1,0,2]-0.0)<=self.RES_TOL,"wrong value for (1,1,0,2)")
         self.assertTrue(Lsup(val[1,1,1,0]-0.0)<=self.RES_TOL,"wrong value for (1,1,1,0)")
         self.assertTrue(Lsup(val[1,1,1,1]-1.0)<=self.RES_TOL,"wrong value for (1,1,1,1)")
         self.assertTrue(Lsup(val[1,1,1,2]-0.0)<=self.RES_TOL,"wrong value for (1,1,1,2)")
         self.assertTrue(Lsup(val[1,1,2,0]-0.0)<=self.RES_TOL,"wrong value for (1,1,2,0)")
         self.assertTrue(Lsup(val[1,1,2,1]-0.0)<=self.RES_TOL,"wrong value for (1,1,2,1)")
         self.assertTrue(Lsup(val[1,1,2,2]-0.0)<=self.RES_TOL,"wrong value for (1,1,2,2)")
         self.assertTrue(Lsup(val[1,2,0,0]-0.0)<=self.RES_TOL,"wrong value for (1,2,0,0)")
         self.assertTrue(Lsup(val[1,2,0,1]-0.0)<=self.RES_TOL,"wrong value for (1,2,0,1)")
         self.assertTrue(Lsup(val[1,2,0,2]-0.0)<=self.RES_TOL,"wrong value for (1,2,0,2)")
         self.assertTrue(Lsup(val[1,2,1,0]-0.0)<=self.RES_TOL,"wrong value for (1,2,1,0)")
         self.assertTrue(Lsup(val[1,2,1,1]-0.0)<=self.RES_TOL,"wrong value for (1,2,1,1)")
         self.assertTrue(Lsup(val[1,2,1,2]-1.0)<=self.RES_TOL,"wrong value for (1,2,1,2)")
         self.assertTrue(Lsup(val[1,2,2,0]-0.0)<=self.RES_TOL,"wrong value for (1,2,2,0)")
         self.assertTrue(Lsup(val[1,2,2,1]-0.0)<=self.RES_TOL,"wrong value for (1,2,2,1)")
         self.assertTrue(Lsup(val[1,2,2,2]-0.0)<=self.RES_TOL,"wrong value for (1,2,2,2)")
         self.assertTrue(Lsup(val[2,0,0,0]-0.0)<=self.RES_TOL,"wrong value for (2,0,0,0)")
         self.assertTrue(Lsup(val[2,0,0,1]-0.0)<=self.RES_TOL,"wrong value for (2,0,0,1)")
         self.assertTrue(Lsup(val[2,0,0,2]-0.0)<=self.RES_TOL,"wrong value for (2,0,0,2)")
         self.assertTrue(Lsup(val[2,0,1,0]-0.0)<=self.RES_TOL,"wrong value for (2,0,1,0)")
         self.assertTrue(Lsup(val[2,0,1,1]-0.0)<=self.RES_TOL,"wrong value for (2,0,1,1)")
         self.assertTrue(Lsup(val[2,0,1,2]-0.0)<=self.RES_TOL,"wrong value for (2,0,1,2)")
         self.assertTrue(Lsup(val[2,0,2,0]-1.0)<=self.RES_TOL,"wrong value for (2,0,2,0)")
         self.assertTrue(Lsup(val[2,0,2,1]-0.0)<=self.RES_TOL,"wrong value for (2,0,2,1)")
         self.assertTrue(Lsup(val[2,0,2,2]-0.0)<=self.RES_TOL,"wrong value for (2,0,2,2)")
         self.assertTrue(Lsup(val[2,1,0,0]-0.0)<=self.RES_TOL,"wrong value for (2,1,0,0)")
         self.assertTrue(Lsup(val[2,1,0,1]-0.0)<=self.RES_TOL,"wrong value for (2,1,0,1)")
         self.assertTrue(Lsup(val[2,1,0,2]-0.0)<=self.RES_TOL,"wrong value for (2,1,0,2)")
         self.assertTrue(Lsup(val[2,1,1,0]-0.0)<=self.RES_TOL,"wrong value for (2,1,1,0)")
         self.assertTrue(Lsup(val[2,1,1,1]-0.0)<=self.RES_TOL,"wrong value for (2,1,1,1)")
         self.assertTrue(Lsup(val[2,1,1,2]-0.0)<=self.RES_TOL,"wrong value for (2,1,1,2)")
         self.assertTrue(Lsup(val[2,1,2,0]-0.0)<=self.RES_TOL,"wrong value for (2,1,2,0)")
         self.assertTrue(Lsup(val[2,1,2,1]-1.0)<=self.RES_TOL,"wrong value for (2,1,2,1)")
         self.assertTrue(Lsup(val[2,1,2,2]-0.0)<=self.RES_TOL,"wrong value for (2,1,2,2)")
         self.assertTrue(Lsup(val[2,2,0,0]-0.0)<=self.RES_TOL,"wrong value for (2,2,0,0)")
         self.assertTrue(Lsup(val[2,2,0,1]-0.0)<=self.RES_TOL,"wrong value for (2,2,0,1)")
         self.assertTrue(Lsup(val[2,2,0,2]-0.0)<=self.RES_TOL,"wrong value for (2,2,0,2)")
         self.assertTrue(Lsup(val[2,2,1,0]-0.0)<=self.RES_TOL,"wrong value for (2,2,1,0)")
         self.assertTrue(Lsup(val[2,2,1,1]-0.0)<=self.RES_TOL,"wrong value for (2,2,1,1)")
         self.assertTrue(Lsup(val[2,2,1,2]-0.0)<=self.RES_TOL,"wrong value for (2,2,1,2)")
         self.assertTrue(Lsup(val[2,2,2,0]-0.0)<=self.RES_TOL,"wrong value for (2,2,2,0)")
         self.assertTrue(Lsup(val[2,2,2,1]-0.0)<=self.RES_TOL,"wrong value for (2,2,2,1)")
         self.assertTrue(Lsup(val[2,2,2,2]-1.0)<=self.RES_TOL,"wrong value for (2,2,2,2)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_unitVector_1(self):
      val=unitVector(i=0,d=1)
      self.assertEqual(val[0],1.0,"wrong value for 0 in the 0 vector")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_unitVector_2(self):
      val=unitVector(i=0,d=2)
      self.assertEqual(val[0],1.0,"wrong value for 0 in the 0 vector")
      self.assertEqual(val[1],0.0,"wrong value for 1 in the 0 vector")
      val=unitVector(i=1,d=2)
      self.assertEqual(val[0],0.0,"wrong value for 0 in the 1 vector")
      self.assertEqual(val[1],1.0,"wrong value for 1 in the 1 vector")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_unitVector_3(self):
      val=unitVector(i=0,d=3)
      self.assertEqual(val[0],1.0,"wrong value for 0 in the 0 vector")
      self.assertEqual(val[1],0.0,"wrong value for 1 in the 0 vector")
      self.assertEqual(val[2],0.0,"wrong value for 2 in the 0 vector")
      val=unitVector(i=1,d=3)
      self.assertEqual(val[0],0.0,"wrong value for 0 in the 1 vector")
      self.assertEqual(val[1],1.0,"wrong value for 1 in the 1 vector")
      self.assertEqual(val[2],0.0,"wrong value for 2 in the 1 vector")
      val=unitVector(i=2,d=3)
      self.assertEqual(val[0],0.0,"wrong value for 0 in the 2 vector")
      self.assertEqual(val[1],0.0,"wrong value for 1 in the 2 vector")
      self.assertEqual(val[2],1.0,"wrong value for 2 in the 2 vector")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_unitVector_domain(self):
      val=unitVector(i=0,d=self.domain)
      self.assertTrue(isinstance(val,numpy.ndarray),"wrong type of result.")
      self.assertEqual(val.shape,(self.domain.getDim(),),"wrong shape.")
      if self.domain.getDim()==2:
         self.assertEqual(val[0],1.0,"wrong value for 0 in the 0 vector")
         self.assertEqual(val[1],0.0,"wrong value for 1 in the 0 vector")
      if self.domain.getDim()==3:
         self.assertEqual(val[0],1.0,"wrong value for 0 in the 0 vector")
         self.assertEqual(val[1],0.0,"wrong value for 1 in the 0 vector")
         self.assertEqual(val[2],0.0,"wrong value for 2 in the 0 vector")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_unitVector_functionspace(self):
      val=unitVector(i=0,d=self.functionspace)
      self.assertTrue(isinstance(val,escript.Data),"wrong type of result.")
      self.assertEqual(val.getShape(),(self.functionspace.getDim(),),"wrong shape.")
      if self.domain.getDim()==2:
         self.assertTrue(Lsup(val[0]-1.0)<=self.RES_TOL,"wrong value for 0 in the 0 vector")
         self.assertTrue(Lsup(val[1]-0.0)<=self.RES_TOL,"wrong value for 1 in the 0 vector")
      if self.domain.getDim()==3:
         self.assertTrue(Lsup(val[0]-1.0)<=self.RES_TOL,"wrong value for 0 in the 0 vector")
         self.assertTrue(Lsup(val[1]-0.0)<=self.RES_TOL,"wrong value for 1 in the 0 vector")
         self.assertTrue(Lsup(val[2]-0.0)<=self.RES_TOL,"wrong value for 2 in the 0 vector")

