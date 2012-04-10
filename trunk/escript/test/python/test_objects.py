
########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Test suite for data objects. at the moment for dump and load only.

The tests must be linked with some function space class object in the setUp method:
to run the use:

   from esys.dudley import Brick
   class Test_DumpOnDudley(Test_Dump):
       def setUp(self):
          self.domain =Rectangle(NE,NE+1,2)
          self.domain_with_different_number_of_samples =Rectangle(2*NE,NE+1,2)
          self.domain_with_different_number_of_data_points_per_sample =Rectangle(2*NE,NE+1,2,integrationOrder=2)
          self.domain_with_different_sample_ordering =Rectangle(1,(NE+1)*NE,2)
          self.filename_base="."

   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_DumpOnDudley))
   unittest.TextTestRunner(verbosity=2).run(suite)

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

import unittest
import os
import numpy
from esys.escript import *

try:
     ESCRIPT_WORKDIR=os.environ['ESCRIPT_WORKDIR']
except KeyError:
     ESCRIPT_WORKDIR='.'

class Test_TableInterpolation(unittest.TestCase):
    RES_TOL=1.e-7 # RES_TOLerance to compare results    
        
        
    def test_NullFunctionSpace(self):
        arL=[[0, -1, -2, -3, -4], [1, 1, -2, -3, -4], [2, 2, 2, -3, -4], [3, 3, 3, 3, -4], [4, 4, 4, 4, 4]]
        arn=numpy.array(arL)
        ars=[arL,arn]
        d0=Data(0)
        d1=Data(1)
        d2=Data(2)
        d35=Data(3.5)
        d4=Data(4)
        dm05=Data(-0.5)
        d175=Data(1.75)
        d225=Data(2.25)
        for arr in ars:
            self.assertTrue(Lsup(d2.interpolateTable(arL,0, 1, d1, 0, 1, 100)+2)<self.RES_TOL)
            self.assertTrue(Lsup(d35.interpolateTable(arL,0, 1, d1, 0, 1, 100)+3.5)<self.RES_TOL)
            self.assertTrue(Lsup(d2.interpolateTable(arL,0,1, d35, 0, 1, 100)-3.5)<self.RES_TOL)
            self.assertTrue(Lsup(d225.interpolateTable(arL,0,1,d175,0,1, 100)-0)<self.RES_TOL)
               # Point out of bounds
            self.assertRaises(RuntimeError, d1.interpolateTable,arL,0, 1, d4, 0, 0.5, 100, check_boundaries=True )
            self.assertRaises(RuntimeError, d4.interpolateTable, arL,0, 0.5, d1, 0, 1, 100, check_boundaries=True  )
            self.assertRaises(RuntimeError, dm05.interpolateTable, arL,0,1, d1 , 0,1, 100, check_boundaries=True  )
            self.assertRaises(RuntimeError, d1.interpolateTable, arL,0,1, dm05 , 0,1, 100,check_boundaries=True  )
               #Test to ensure not check_boundaries does not throw in the above cases
            d1.interpolateTable(arL,0, 1, d4, 0, 0.5, 100, check_boundaries=False)
            d4.interpolateTable( arL,0, 0.5, d1, 0, 1, 100, check_boundaries=False  )
            dm05.interpolateTable( arL,0,1, d1 , 0,1, 100, check_boundaries=False  )
            d1.interpolateTable( arL,0,1, dm05 , 0,1, 100,check_boundaries=False  )
               
               # interpolated value too large
            self.assertRaises(RuntimeError, d2.interpolateTable, arL, 0, 1, d2, 0, 1, 1 )


    def test_FunctionSpace3D(self):
        vs=[(1,3,5,7,11,13,17,19), (-1,1,-1,1,-1,1,-1,1), (0.5, 17, 0.25, 42, 0.125, 35, 0.625, 49)]   #There is no particular significance to these numbers
        for fs in self.functionspaces:
            points=fs.getX()
            for t in vs:
                v0, v1, v2, v3, v4, v5, v6, v7 =t
                x=points[0]
                y=points[1]
                z=points[2]
                xmax=sup(x)
                xmin=inf(x)
                ymax=sup(y)
                ymin=inf(y)
                zmax=sup(z)
                zmin=inf(z)
                xwidth=(xmax-xmin)/(self.xn-1)
                ywidth=(ymax-ymin)/(self.yn-1)
                zwidth=(zmax-zmin)/(self.zn-1)          
                table=[]
                for k in range(self.zn):
                    face=[]
                    for j in range(self.yn):
                        row=[]
                        for i in range(self.xn):
                                row.append(v0+v1*xwidth*i+v2*ywidth*j+v3*i*j*xwidth*ywidth)
                        face.append(row)
                    table.append(face)
                ref=v0+v1*(x-xmin)+v2*(y-ymin)+v3*(x-xmin)*(y-ymin)
                lsupref=Lsup(ref)
                if lsupref!=0 and xwidth>0:             #This will happen in cases where there are no samples on this rank
                    res=interpolateTable(table, points, (xmin, ymin, zmin), (xwidth, ywidth, zwidth),900)
                    self.assertTrue(Lsup(res-ref)/Lsup(ref)<self.RES_TOL,"Failed for %s"%str(fs))
                    # Test undef
                    self.assertRaises(RuntimeError, interpolateTable, table, points, (xmin, ymin, zmin), (xwidth, ywidth, zwidth), -1)
                    # Test bounds checking
                    self.assertRaises(RuntimeError, interpolateTable, table, points, (xmin, ymin, zmin), (xwidth/3, ywidth, zwidth), 900,True)
                    self.assertRaises(RuntimeError, interpolateTable, table, points, (xmin, ymin, zmin), (xwidth, ywidth/3, zwidth), 900, True)
                    self.assertRaises(RuntimeError, interpolateTable, table, points, (xmin, ymin, zmin), (xwidth, ywidth, zwidth/3), 900, True)

    def test_FunctionSpace2D(self):
        vs=[(1,3,5,7), (-1,1,-1,1), (0.5, 17, 0.25, 42)]   #There is no particular significance to these numbers
        for fs in self.functionspaces:
            points=fs.getX()
            for t in vs:
                v0, v1, v2, v3 =t
                x=points[0]
                y=points[1]
                xmax=sup(x)
                xmin=inf(x)
                ymax=sup(y)
                ymin=inf(y)
                xwidth=(xmax-xmin)/(self.xn-1)
                ywidth=(ymax-ymin)/(self.yn-1)
                table=[]
                for j in range(self.yn):
                      row=[]
                      for i in range(self.xn):
                        row.append(v0+v1*xwidth*i+v2*ywidth*j+v3*i*j*xwidth*ywidth)
                      table.append(row)
                ref=v0+v1*(x-xmin)+v2*(y-ymin)+v3*(x-xmin)*(y-ymin) 
                lsupref=Lsup(ref)
                if lsupref!=0 and xwidth>0:             #This will happen in cases where there are no samples on this rank
                    res=x.interpolateTable(table,xmin,xwidth,y, ymin, ywidth,500)
                    self.assertTrue(Lsup(res-ref)/Lsup(ref)<self.RES_TOL,"Failed for %s"%str(fs))
                    #Now we try for the new interface
                    res=interpolateTable(table,points[0:2], (xmin, ymin), (xwidth, ywidth),500)
                    self.assertTrue(Lsup(res-ref)/Lsup(ref)<self.RES_TOL,"Failed for %s under unified call."%str(fs))
                    
            

    def test_FunctionSpace1D(self):
        vs=[(1,3), (-1,1), (0.5, 17)]     #There is no particular significance to these numbers
        for fs in self.functionspaces:
            points=fs.getX()
            for t in vs:
                v0, v1 =t
                x=points[0]
                xmax=sup(x)
                xmin=inf(x)
                xwidth=(xmax-xmin)/(self.xn-1)
                table=[]
                for i in range(self.xn):
                   table.append(v0+v1*xwidth*i)
                ref=v0+v1*(x-xmin)
                lsupref=Lsup(ref)
                if lsupref!=0 and xwidth>0:             #This will happen in cases where there are no samples on this rank
                   res=x.interpolateTable(table, xmin, xwidth,500)
                   self.assertTrue(Lsup(res-ref)/lsupref<self.RES_TOL,"Failed for %s"%str(fs))
                   #Now we try for the new interface
                   res=interpolateTable(table,points[0:1], (xmin,), (xwidth, ),500)
                   self.assertTrue(Lsup(res-ref)/Lsup(ref)<self.RES_TOL,"Failed for %s under unified call."%str(fs))
                   res=interpolateTable(table,points[0:1], xmin, xwidth,500)
                   self.assertTrue(Lsup(res-ref)/Lsup(ref)<self.RES_TOL,"Failed for %s under unified call (no tuple)."%str(fs))

class Test_saveCSV(unittest.TestCase):

   def test_save1(self):
        X=self.domain.getX()
        X0=X[0]
        fname=os.path.join(ESCRIPT_WORKDIR, "test_save1.csv")
        saveDataCSV(fname,C=X, D=X0)
        saveDataCSV(fname,append=True, J=X0, H=X)
        self.assertTrue(os.path.exists(fname), "test file not created")
        f=open(fname,'r')
        line=f.readline()
        self.assertEquals(line,"C_0, C_1, D\n")    #This tests both separator strings
        #Now we find out how many rows it has
        linecount=0
        while line != '':
                linecount+=1
                line=f.readline()
        self.assertEquals(linecount,self.linecount1*2)
        f.close()
        #Now to other output
        T2=Tensor(7,X.getFunctionSpace())
        T3=Tensor3(8,X.getFunctionSpace())
        T4=Tensor4(9,X.getFunctionSpace())
        saveDataCSV(fname,A=T2,B=T3,C=T4)
        f=open(fname,'r')
        line=f.readline()
        self.assertEquals(line, 'A_0_0, A_1_0, A_0_1, A_1_1, B_0_0_0, B_0_0_1, B_1_0_0, B_1_0_1, B_0_1_0, B_0_1_1, B_1_1_0, B_1_1_1, C_0_0_0_0, C_0_0_0_1, C_0_0_1_0, C_0_0_1_1, C_1_0_0_0, C_1_0_0_1, C_1_0_1_0, C_1_0_1_1, C_0_1_0_0, C_0_1_0_1, C_0_1_1_0, C_0_1_1_1, C_1_1_0_0, C_1_1_0_1, C_1_1_1_0, C_1_1_1_1\n')
        line=f.readline()
        line_expected=[7.]*4+[8.]*8+[9.]*16
        line_got=[float(elt) for elt in line.split(',')]
        self.assertEquals(line_got,line_expected)
        linecount=1
        while line != '':
                linecount+=1
                line=f.readline()
        self.assertEquals(linecount, self.linecount1)
        f.close()       
        #Now to test separators and mask
        saveDataCSV(fname, sep="|",csep="/", U=X, V=X0, mask=X0)
        f=open(fname,'r')
        line=f.readline()
        self.assertEquals(line, 'U/0|U/1|V\n')
        line=f.readline()
        line_got=[float(elt) for elt in line.split('|')]
        self.assertEquals(len(self.line_expected),len(line_got))
        for i in range(len(self.line_expected)):
                self.assertAlmostEquals(self.line_expected[i],line_got[i])
        linecount=1
        while line!='':
                linecount+=1
                line=f.readline()
        self.assertEquals(linecount,self.linecount2)
        
        
class Test_Domain(unittest.TestCase):

   def test_getListOfTags(self): # requires self.boundary_tag_list
       tags=FunctionOnBoundary(self.domain).getListOfTags()
       self.assertTrue(len(self.boundary_tag_list) == len(tags), "tag list length does not match")
       for i in self.boundary_tag_list:
           self.assertTrue(i in tags, "tag %s is missing."%i)

   def test_RandomData(self):
        fs=Function(self.domain)        # The choice of functionspace is arbitrary
        dat=RandomData((2,2,2,2),fs,8)  # Choice of seed is arbitrary
        self.assertTrue(Lsup(dat-1)<1.0001)

   def test_Factories(self):
        fs=Function(self.domain)        # The choice of functionspace is arbitrary
        dime=self.domain.getDim()
        if dime>0:
           z=[]
           bad=[]
           for i in range(dime):
                z+=[i]
                bad+=[i]
           bad+=[i]
           d=Vector(z,fs)
           self.assertTrue(d.getShape()==(dime,))
           self.assertRaises(RuntimeError, Vector, bad,fs)      #test wrong shape
           y=[]
           bad=[]
           for i in range(dime):
                y+=[z]
                bad+=[z]
           bad+=[z]
           z=y
           d=Tensor(z,fs)
           self.assertTrue(d.getShape()==(dime,dime))
           try:
                Tensor(bad,fs)
           except RuntimeError:
                pass
           else:
                self.fail("Tensor should have rejected bad shape")
           y=[]
           bad=[]
           for i in range(dime):
                y+=[z]
                bad+=[z]
           bad+=[z]
           z=y
           d=Tensor3(z,fs)
           self.assertTrue(d.getShape()==(dime,dime,dime))
           try:
                Tensor3(bad,fs)
           except RuntimeError:
                pass
           else:
                self.fail("Tensor3 should have rejected bad shape")
           y=[]
           bad=[]
           for i in range(dime):
                y+=[z]
                bad+=[z]
           bad+=[z]
           z=y
           d=Tensor4(z,fs)
           self.assertTrue(d.getShape()==(dime,dime,dime,dime))
           try:
                Tensor4(bad,fs)
           except RuntimeError:
                pass
           else:
                self.fail("Tensor4 should have rejected bad shape")


   def test_addTags(self):
        tag1="A"
        tag2="B"
        tag3="C"
        self.domain.setTagMap(tag1,1)
        self.assertTrue(self.domain.isValidTagName(tag1))
        self.assertTrue(not self.domain.isValidTagName(tag2))
        self.domain.setTagMap(tag2,2)
        self.assertTrue(self.domain.isValidTagName(tag1))
        self.assertTrue(self.domain.isValidTagName(tag2))
        self.assertTrue(not self.domain.isValidTagName(tag3))
        self.assertTrue(self.domain.getTag(tag1)==1)
        self.assertTrue(self.domain.getTag(tag2)==2)
        self.assertRaises(RuntimeError,self.domain.getTag,tag3)

        # set tag:
        s=Scalar(0,Function(self.domain))
        r=Scalar(0,Function(self.domain))
        s.setTaggedValue(tag1,1.)
        r.setTaggedValue(1,1.)
        s.setTaggedValue(tag2,2.)
        r.setTaggedValue(2,2.)
        self.assertRaises(RuntimeError,s.setTaggedValue,tag3,3.)        #tag3 does not exist
        self.assertTrue(Lsup(s-r)<=0.)
        # get tag:
        names=getTagNames(self.domain)
        self.assertTrue(len(names) == 6)
        self.assertTrue( tag1 in names )
        self.assertTrue( tag2 in names )
        self.assertTrue(self.domain.isValidTagName(tag1))
        self.assertTrue(self.domain.isValidTagName(tag2))
        # insert tag shortcut:
        s2=insertTaggedValues(Scalar(0,Function(self.domain)),**{ tag1 : 1., tag2 : 2.})
        self.assertTrue(Lsup(s2-r)<=0.)
   def test_functionspace_ContinuousFunction(self):
        fs=ContinuousFunction(self.domain)
        self.assertTrue(fs.getDomain()==self.domain)
        self.assertTrue(self.domain.getDim() == fs.getDim())
        x=fs.getX()
        self.assertTrue(x.getFunctionSpace() == fs)
        self.assertTrue(x.getShape() == (fs.getDim(),))
        self.assertTrue(inf(x[0])>=0.)
        if self.domain.getDim()>1: self.assertTrue(inf(x[1])>=0.)
        if self.domain.getDim()>2: self.assertTrue(inf(x[2])>=0.)
        self.assertTrue(sup(x[0])<=1.)
        if self.domain.getDim()>1: self.assertTrue(sup(x[1])<=1.)
        if self.domain.getDim()>2: self.assertTrue(sup(x[2])<=1.)

   def test_functionspace_Solution(self):
        fs=Solution(self.domain)
        self.assertTrue(fs.getDomain()==self.domain)
        self.assertTrue(self.domain.getDim() == fs.getDim())
        x=fs.getX()
        self.assertTrue(x.getFunctionSpace() == fs)
        self.assertTrue(x.getShape() == (fs.getDim(),))
        self.assertTrue(inf(x[0])>=0.)
        if self.domain.getDim()>1: self.assertTrue(inf(x[1])>=0.)
        if self.domain.getDim()>2: self.assertTrue(inf(x[2])>=0.)
        self.assertTrue(sup(x[0])<=1.)
        if self.domain.getDim()>1: self.assertTrue(sup(x[1])<=1.)
        if self.domain.getDim()>2: self.assertTrue(sup(x[2])<=1.)

   def test_functionspace_ReducedSolution(self):
        fs=ReducedSolution(self.domain)
        self.assertTrue(fs.getDomain()==self.domain)
        self.assertTrue(self.domain.getDim() == fs.getDim())
        x=fs.getX()
        self.assertTrue(x.getFunctionSpace() == fs)
        self.assertTrue(x.getShape() == (fs.getDim(),))
        self.assertTrue(inf(x[0])>=0.)
        if self.domain.getDim()>1: self.assertTrue(inf(x[1])>=0.)
        if self.domain.getDim()>2: self.assertTrue(inf(x[2])>=0.)
        self.assertTrue(sup(x[0])<=1.)
        if self.domain.getDim()>1: self.assertTrue(sup(x[1])<=1.)
        if self.domain.getDim()>2: self.assertTrue(sup(x[2])<=1.)

   def test_functionspace_Function(self):
        fs=Function(self.domain)
        self.assertTrue(fs.getDomain()==self.domain)
        self.assertTrue(self.domain.getDim() == fs.getDim())
        x=fs.getX()
        self.assertTrue(x.getFunctionSpace() == fs)
        self.assertTrue(x.getShape() == (fs.getDim(),))
        self.assertTrue(inf(x[0])>=0.)
        if self.domain.getDim()>1: self.assertTrue(inf(x[1])>=0.)
        if self.domain.getDim()>2: self.assertTrue(inf(x[2])>=0.)
        self.assertTrue(sup(x[0])<=1.)
        if self.domain.getDim()>1: self.assertTrue(sup(x[1])<=1.)
        if self.domain.getDim()>2: self.assertTrue(sup(x[2])<=1.)

   def test_functionspace_ReducedFunction(self):
        fs=ReducedFunction(self.domain)
        self.assertTrue(fs.getDomain()==self.domain)
        self.assertTrue(self.domain.getDim() == fs.getDim())
        x=fs.getX()
        self.assertTrue(x.getFunctionSpace() == fs)
        self.assertTrue(x.getShape() == (fs.getDim(),))
        self.assertTrue(inf(x[0])>=0.)
        if self.domain.getDim()>1: self.assertTrue(inf(x[1])>=0.)
        if self.domain.getDim()>2: self.assertTrue(inf(x[2])>=0.)
        self.assertTrue(sup(x[0])<=1.)
        if self.domain.getDim()>1: self.assertTrue(sup(x[1])<=1.)
        if self.domain.getDim()>2: self.assertTrue(sup(x[2])<=1.)
   def test_functionspace_FunctionOnBoundary(self):
        fs=FunctionOnBoundary(self.domain)
        self.assertTrue(fs.getDomain()==self.domain)
        self.assertTrue(self.domain.getDim() == fs.getDim())
        x=fs.getX()
        self.assertTrue(x.getFunctionSpace() == fs)
        self.assertTrue(x.getShape() == (fs.getDim(),))
        self.assertTrue(inf(x[0])>=0.)
        if self.domain.getDim()>1: self.assertTrue(inf(x[1])>=0.)
        if self.domain.getDim()>2: self.assertTrue(inf(x[2])>=0.)
        self.assertTrue(sup(x[0])<=1.)
        if self.domain.getDim()>1: self.assertTrue(sup(x[1])<=1.)
        if self.domain.getDim()>2: self.assertTrue(sup(x[2])<=1.)

   def test_functionspace_ReducedFunctionOnBoundary(self):
        fs=ReducedFunctionOnBoundary(self.domain)
        self.assertTrue(fs.getDomain()==self.domain)
        self.assertTrue(self.domain.getDim() == fs.getDim())
        x=fs.getX()
        self.assertTrue(x.getFunctionSpace() == fs)
        self.assertTrue(x.getShape() == (fs.getDim(),))
        self.assertTrue(inf(x[0])>=0.)
        if self.domain.getDim()>1: self.assertTrue(inf(x[1])>=0.)
        if self.domain.getDim()>2: self.assertTrue(inf(x[2])>=0.)
        self.assertTrue(sup(x[0])<=1.)
        if self.domain.getDim()>1: self.assertTrue(sup(x[1])<=1.)
        if self.domain.getDim()>2: self.assertTrue(sup(x[2])<=1.)
   #===========================================================================

class Test_GlobalMinMax(unittest.TestCase):
   def test_GlobalMinMax(self):
        myrank=getMPIRankWorld()
        d=Data(myrank,Function(self.domain))
        minproc=inf(d)
        maxproc=sup(d)          #This tells us where to expect values to be
        if d.getNumberOfDataPoints()>0:
                d.setValueOfDataPoint(0,myrank-0.001);
        p,n=d.minGlobalDataPoint()
        self.assertTrue(p==minproc,"Incorrect process indentified as holding min")
        self.assertTrue(n==0,"Incorrect position for min")
        if d.getNumberOfDataPoints()>0:
                d.setValueOfDataPoint(0,myrank+0.001)
        p,n=d.maxGlobalDataPoint()      
        self.assertTrue(p==maxproc,"Incorrect process indentified as holding min")
        self.assertTrue(n==0,"Incorrect position for min")

        

class Test_SetDataPointValue(unittest.TestCase):
   arg0=9.81
   arg1=numpy.array([3.098, -3.111])
   arg2=numpy.array([[3.82, -3.81, -0.957, 0.892, -1.367], [-4.589, -1.835, -2.679, -1.517, -4.2515], [-4.909, 1.634, -2.883, -2.135, 1.187], [0.6431, 4.638, -4.616, -0.196, -4.370]])
   arg3=numpy.array([[[-2.3667, -0.040], [-4.7398, -3.2412]], [[-2.125, -2.240], [2.237, -4.279]], [[0.68720, 2.4059], [-2.4964, 3.17453]], [[-4.907, -4.9431], [-0.3604, 0.4269]], [[1.4179, 3.326], [1.356, -0.4610]], [[3.378, 2.0902], [-2.6857, 1.3585]]])
   arg4=numpy.array([[[[-3.810, -1.3597, -1.5307, 1.099], [-1.828, 0.2526, -1.4429, 2.326], [4.9732, -2.063, 1.3153, -3.809]], [[-4.8902, -4.714, 1.520, -1.931], [-3.8847, 4.3867, 1.894030, 2.432], [-1.2082, -0.8304, 2.2612, 4.6399]]], [[[-4.5922, -3.309, -0.8171, -0.7210], [2.8051, -4.93047, 0.08450, 4.3824], [0.43204, 2.1908, 4.512633, -1.8218]], [[2.2493, -4.190, -2.3893, -4.147], [-2.104, -4.635, -4.2767, -3.53151], [-2.351, -1.6614, 2.9385, 4.099]]], [[[1.710, 0.2235, -3.4917, 0.8713], [-0.2881, 4.6278, 3.603, -2.1211], [-0.565, 4.294, -2.210827, -0.37651]], [[0.6578, -2.869, -2.490, -4.789], [3.232, 2.483, 0.9531, 2.260], [-1.785, 0.42156, -1.8379, 4.212]]]])
   def test_SetDataPointValue_Function_Rank0(self):
          d=Data(self.arg0,Function(self.domain))
          d.setValueOfDataPoint(0,self.arg0*2)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg0)
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg0*2)<=Lsup(self.arg0*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg0)<=Lsup(self.arg0), "wrong setting")
   def test_SetDataPointValue_Function_Rank1(self):
          d=Data(self.arg1,Function(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg2)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg1)
          d.setValueOfDataPoint(0,self.arg1*2)
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg1*2)<=Lsup(self.arg1*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg1)<=Lsup(self.arg1), "wrong setting")
   def test_SetDataPointValue_Function_Rank1_list(self):
          d=Data(self.arg1,Function(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg2)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg1)
          d.setValueOfDataPoint(0,(self.arg1*2).tolist())
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg1*2)<=Lsup(self.arg1*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg1)<=Lsup(self.arg1), "wrong setting")
   def test_SetDataPointValue_Function_Rank2(self):
          d=Data(self.arg2,Function(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg2)
          d.setValueOfDataPoint(0,self.arg2*2)
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg2*2)<=Lsup(self.arg2*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg2)<=Lsup(self.arg2), "wrong setting")
   def test_SetDataPointValue_Function_Rank2_list(self):
          d=Data(self.arg2,Function(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg2)
          d.setValueOfDataPoint(0,(self.arg2*2).tolist())
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg2*2)<=Lsup(self.arg2*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg2)<=Lsup(self.arg2), "wrong setting")
   def test_SetDataPointValue_Function_Rank3(self):
          d=Data(self.arg3,Function(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg3)
          d.setValueOfDataPoint(0,self.arg3*2)
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg3*2)<=Lsup(self.arg3*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg3)<=Lsup(self.arg3), "wrong setting")
   def test_SetDataPointValue_Function_Rank3_list(self):
          d=Data(self.arg3,Function(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg3)
          d.setValueOfDataPoint(0,(self.arg3*2).tolist())
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg3*2)<=Lsup(self.arg3*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg3)<=Lsup(self.arg3), "wrong setting")
   def test_SetDataPointValue_Function_Rank4(self):
          d=Data(self.arg4,Function(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg4)
          d.setValueOfDataPoint(0,self.arg4*2)
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg4*2)<=Lsup(self.arg4*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg4)<=Lsup(self.arg4), "wrong setting")
   def test_SetDataPointValue_Function_Rank4_list(self):
          d=Data(self.arg4,Function(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg4)
          d.setValueOfDataPoint(0,(self.arg4*2).tolist())
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg4*2)<=Lsup(self.arg4*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg4)<=Lsup(self.arg4), "wrong setting")
   #===========================================================================
   def test_SetDataPointValue_ReducedFunction_Rank0(self):
          d=Data(self.arg0,ReducedFunction(self.domain))
          d.setValueOfDataPoint(0,self.arg0*2)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg0)
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg0*2)<=Lsup(self.arg0*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg0)<=Lsup(self.arg0), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank1(self):
          d=Data(self.arg1,ReducedFunction(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg2)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg1)
          d.setValueOfDataPoint(0,self.arg1*2)
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg1*2)<=Lsup(self.arg1*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg1)<=Lsup(self.arg1), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank1_list(self):
          d=Data(self.arg1,ReducedFunction(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg2)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg1)
          d.setValueOfDataPoint(0,(self.arg1*2).tolist())
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg1*2)<=Lsup(self.arg1*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg1)<=Lsup(self.arg1), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank2(self):
          d=Data(self.arg2,ReducedFunction(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg2)
          d.setValueOfDataPoint(0,self.arg2*2)
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg2*2)<=Lsup(self.arg2*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg2)<=Lsup(self.arg2), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank2_list(self):
          d=Data(self.arg2,ReducedFunction(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg2)
          d.setValueOfDataPoint(0,(self.arg2*2).tolist())
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg2*2)<=Lsup(self.arg2*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg2)<=Lsup(self.arg2), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank3(self):
          d=Data(self.arg3,ReducedFunction(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg3)
          d.setValueOfDataPoint(0,self.arg3*2)
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg3*2)<=Lsup(self.arg3*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg3)<=Lsup(self.arg3), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank3_list(self):
          d=Data(self.arg3,ReducedFunction(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg3)
          d.setValueOfDataPoint(0,(self.arg3*2).tolist())
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg3*2)<=Lsup(self.arg3*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg3)<=Lsup(self.arg3), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank4(self):
          d=Data(self.arg4,ReducedFunction(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg4)
          d.setValueOfDataPoint(0,self.arg4*2)
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg4*2)<=Lsup(self.arg4*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg4)<=Lsup(self.arg4), "wrong setting")
   def test_SetDataPointValue_ReducedFunction_Rank4_list(self):
          d=Data(self.arg4,ReducedFunction(self.domain))
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, 0, self.arg1)
          self.assertRaises(RuntimeError, d.setValueOfDataPoint, -1, self.arg4)
          d.setValueOfDataPoint(0,(self.arg4*2).tolist())
          d_0=numpy.array(d.getTupleForDataPoint(0))
          d_1=numpy.array(d.getTupleForDataPoint(1))
          self.assertTrue(Lsup(d_0-self.arg4*2)<=Lsup(self.arg4*2), "wrong setting")
          self.assertTrue(Lsup(d_1-self.arg4)<=Lsup(self.arg4), "wrong setting")

class Test_Dump(unittest.TestCase):
   arg0=9.81
   arg1=numpy.array([3.098, -3.111])
   arg2=numpy.array([[3.82, -3.81, -0.957, 0.892, -1.367], [-4.589, -1.835, -2.679, -1.517, -4.2515], [-4.909, 1.634, -2.883, -2.135, 1.187], [0.6431, 4.638, -4.616, -0.196, -4.370]])
   arg3=numpy.array([[[-2.3667, -0.040], [-4.7398, -3.2412]], [[-2.125, -2.240], [2.237, -4.279]], [[0.68720, 2.4059], [-2.4964, 3.17453]], [[-4.907, -4.9431], [-0.3604, 0.4269]], [[1.4179, 3.326], [1.356, -0.4610]], [[3.378, 2.0902], [-2.6857, 1.3585]]])
   arg4=numpy.array([[[[-3.810, -1.3597, -1.5307, 1.099], [-1.828, 0.2526, -1.4429, 2.326], [4.9732, -2.063, 1.3153, -3.809]], [[-4.8902, -4.714, 1.520, -1.931], [-3.8847, 4.3867, 1.894030, 2.432], [-1.2082, -0.8304, 2.2612, 4.6399]]], [[[-4.5922, -3.309, -0.8171, -0.7210], [2.8051, -4.93047, 0.08450, 4.3824], [0.43204, 2.1908, 4.512633, -1.8218]], [[2.2493, -4.190, -2.3893, -4.147], [-2.104, -4.635, -4.2767, -3.53151], [-2.351, -1.6614, 2.9385, 4.099]]], [[[1.710, 0.2235, -3.4917, 0.8713], [-0.2881, 4.6278, 3.603, -2.1211], [-0.565, 4.294, -2.210827, -0.37651]], [[0.6578, -2.869, -2.490, -4.789], [3.232, 2.483, 0.9531, 2.260], [-1.785, 0.42156, -1.8379, 4.212]]]])

   def _diffDataObjects(self,d_ref,filemame, use_old_file=False):
       if not use_old_file: d_ref.dump(filemame)
       d=load(filemame, d_ref.getDomain())
       self.assertTrue(not d.isEmpty(),"data in %s are empty."%filemame)
       self.assertTrue(d_ref.getRank() == d.getRank(), "different rank in %s. "%filemame)
       self.assertTrue(d_ref.getShape() == d.getShape(), "different shape %s. "%filemame)
       self.assertTrue(d_ref.getFunctionSpace() == d.getFunctionSpace(), "wrong function space in %s."%filemame)
       self.assertTrue(Lsup(d_ref-d)<=0., "different entries %s."%filemame)

   #===========================================================================
   def test_DumpAndLoad_Constant_Solution_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_solution_rank0.nc")
          d=Data(self.arg0,Solution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_Solution_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_solution_rank1.nc")
          d=Data(self.arg1,Solution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_Solution_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_solution_rank2.nc")
          d=Data(self.arg2,Solution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_Solution_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_solution_rank3.nc")
          d=Data(self.arg3,Solution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_Solution_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_solution_rank4.nc")
          d=Data(self.arg4,Solution(self.domain))
          self._diffDataObjects(d,filemame)
   #===========================================================================
   def test_DumpAndLoad_Constant_ReducedSolution_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_solution_rank0.nc")
          d=Data(self.arg0,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_solution_rank1.nc")
          d=Data(self.arg1,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_solution_rank2.nc")
          d=Data(self.arg2,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_solution_rank3.nc")
          d=Data(self.arg3,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedSolution_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_solution_rank4.nc")
          d=Data(self.arg4,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)
   #===========================================================================
   def test_DumpAndLoad_Constant_ContinuousFunction_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_continuous_function_rank0.nc")
          d=Data(self.arg0,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_continuous_function_rank1.nc")
          d=Data(self.arg1,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_continuous_function_rank2.nc")
          d=Data(self.arg2,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_continuous_function_rank3.nc")
          d=Data(self.arg3,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ContinuousFunction_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_continuous_function_rank4.nc")
          d=Data(self.arg4,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Constant_Function_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_rank0.nc")
          d=Data(self.arg0,Function(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_Function_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_rank1.nc")
          d=Data(self.arg1,Function(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_Function_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_rank2.nc")
          d=Data(self.arg2,Function(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_Function_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_rank3.nc")
          d=Data(self.arg3,Function(self.domain))
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Constant_ReducedFunction_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_rank0.nc")
          d=Data(self.arg0,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedFunction_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_rank1.nc")
          d=Data(self.arg1,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedFunction_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_rank2.nc")
          d=Data(self.arg2,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedFunction_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_rank3.nc")
          d=Data(self.arg3,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)
   def test_DumpAndLoad_Constant_ReducedFunction_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_rank4.nc")
          d=Data(self.arg4,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_on_boundary_rank0.nc")
          d=Data(self.arg0,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_on_boundary_rank1.nc")
          d=Data(self.arg1,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_on_boundary_rank2.nc")
          d=Data(self.arg2,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_on_boundary_rank3.nc")
          d=Data(self.arg3,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_FunctionOnBoundary_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_function_on_boundary_rank4.nc")
          d=Data(self.arg4,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Constant_ReducedFunctionOnBoundary_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_on_boundary_rank0.nc")
          d=Data(self.arg0,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedFunctionOnBoundary_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_on_boundary_rank1.nc")
          d=Data(self.arg1,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedFunctionOnBoundary_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_on_boundary_rank2.nc")
          d=Data(self.arg2,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedFunctionOnBoundary_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_on_boundary_rank3.nc")
          d=Data(self.arg3,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Constant_ReducedFunctionOnBoundary_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"constant_reduced_function_on_boundary_rank4.nc")
          d=Data(self.arg4,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Expanded_Solution_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_solution_rank0.nc")
          d=Data(length(Solution(self.domain).getX())*self.arg0,Solution(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Solution(self.domain_with_different_sample_ordering).getX())*self.arg0,Solution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Solution_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_solution_rank1.nc")
          d=Data(length(Solution(self.domain).getX())*self.arg1,Solution(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Solution(self.domain_with_different_sample_ordering).getX())*self.arg1,Solution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Solution_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_solution_rank2.nc")
          d=Data(length(Solution(self.domain).getX())*self.arg2,Solution(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Solution(self.domain_with_different_sample_ordering).getX())*self.arg2,Solution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Solution_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_solution_rank3.nc")
          d=Data(length(Solution(self.domain).getX())*self.arg3,Solution(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Solution(self.domain_with_different_sample_ordering).getX())*self.arg3,Solution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Solution_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_solution_rank4.nc")
          d=Data(length(Solution(self.domain).getX())*self.arg4,Solution(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Solution(self.domain_with_different_sample_ordering).getX())*self.arg4,Solution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)
   #===========================================================================
   def test_DumpAndLoad_Expanded_ReducedSolution_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_solution_rank0.nc")
          d=Data(length(ReducedSolution(self.domain).getX())*self.arg0,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedSolution(self.domain_with_different_sample_ordering).getX())*self.arg0,ReducedSolution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_solution_rank1.nc")
          d=Data(length(ReducedSolution(self.domain).getX())*self.arg1,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedSolution(self.domain_with_different_sample_ordering).getX())*self.arg1,ReducedSolution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_solution_rank2.nc")
          d=Data(length(ReducedSolution(self.domain).getX())*self.arg2,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedSolution(self.domain_with_different_sample_ordering).getX())*self.arg2,ReducedSolution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_solution_rank3.nc")
          d=Data(length(ReducedSolution(self.domain).getX())*self.arg3,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedSolution(self.domain_with_different_sample_ordering).getX())*self.arg3,ReducedSolution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedSolution_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_solution_rank4.nc")
          d=Data(length(ReducedSolution(self.domain).getX())*self.arg4,ReducedSolution(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedSolution(self.domain_with_different_sample_ordering).getX())*self.arg4,ReducedSolution(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)
   #===========================================================================
   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_continuous_function_rank0.nc")
          d=Data(length(ContinuousFunction(self.domain).getX())*self.arg0,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ContinuousFunction(self.domain_with_different_sample_ordering).getX())*self.arg0,ContinuousFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_continuous_function_rank1.nc")
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          d=Data(length(ContinuousFunction(self.domain).getX())*self.arg1,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ContinuousFunction(self.domain_with_different_sample_ordering).getX())*self.arg1,ContinuousFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_continuous_function_rank2.nc")
          d=Data(length(ContinuousFunction(self.domain).getX())*self.arg2,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ContinuousFunction(self.domain_with_different_sample_ordering).getX())*self.arg2,ContinuousFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_continuous_function_rank3.nc")
          d=Data(length(ContinuousFunction(self.domain).getX())*self.arg3,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ContinuousFunction(self.domain_with_different_sample_ordering).getX())*self.arg3,ContinuousFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ContinuousFunction_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_continuous_function_rank4.nc")
          d=Data(length(ContinuousFunction(self.domain).getX())*self.arg4,ContinuousFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ContinuousFunction(self.domain_with_different_sample_ordering).getX())*self.arg4,ContinuousFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   #===========================================================================
   def test_DumpAndLoad_Expanded_Function_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_rank0.nc")
          d=Data(length(Function(self.domain).getX())*self.arg0,Function(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Function(self.domain_with_different_sample_ordering).getX())*self.arg0,Function(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Function_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_rank1.nc")
          d=Data(length(Function(self.domain).getX())*self.arg1,Function(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Function(self.domain_with_different_sample_ordering).getX())*self.arg1,Function(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Function_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_rank2.nc")
          d=Data(length(Function(self.domain).getX())*self.arg2,Function(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Function(self.domain_with_different_sample_ordering).getX())*self.arg2,Function(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Function_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_rank3.nc")
          d=Data(length(Function(self.domain).getX())*self.arg3,Function(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Function(self.domain_with_different_sample_ordering).getX())*self.arg3,Function(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_Function_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_rank4.nc")
          d=Data(length(Function(self.domain).getX())*self.arg4,Function(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(Function(self.domain_with_different_sample_ordering).getX())*self.arg4,Function(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   #===========================================================================
   def test_DumpAndLoad_Expanded_ReducedFunction_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_rank0.nc")
          d=Data(length(ReducedFunction(self.domain).getX())*self.arg0,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunction(self.domain_with_different_sample_ordering).getX())*self.arg0,ReducedFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunction_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_rank1.nc")
          d=Data(length(ReducedFunction(self.domain).getX())*self.arg1,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunction(self.domain_with_different_sample_ordering).getX())*self.arg1,ReducedFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunction_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_rank2.nc")
          d=Data(length(ReducedFunction(self.domain).getX())*self.arg2,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunction(self.domain_with_different_sample_ordering).getX())*self.arg2,ReducedFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunction_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_rank3.nc")
          d=Data(length(ReducedFunction(self.domain).getX())*self.arg3,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunction(self.domain_with_different_sample_ordering).getX())*self.arg3,ReducedFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunction_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_rank4.nc")
          d=Data(length(ReducedFunction(self.domain).getX())*self.arg4,ReducedFunction(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunction(self.domain_with_different_sample_ordering).getX())*self.arg4,ReducedFunction(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   #===========================================================================
   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_on_boundary_rank0.nc")
          d=Data(length(FunctionOnBoundary(self.domain).getX())*self.arg0,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(FunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg0,FunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_on_boundary_rank1.nc")
          d=Data(length(FunctionOnBoundary(self.domain).getX())*self.arg1,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(FunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg1,FunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_on_boundary_rank2.nc")
          d=Data(length(FunctionOnBoundary(self.domain).getX())*self.arg2,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(FunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg2,FunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_on_boundary_rank3.nc")
          d=Data(length(FunctionOnBoundary(self.domain).getX())*self.arg3,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(FunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg3,FunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_FunctionOnBoundary_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_function_on_boundary_rank4.nc")
          d=Data(length(FunctionOnBoundary(self.domain).getX())*self.arg4,FunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(FunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg4,FunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   #===========================================================================
   def test_DumpAndLoad_Expanded_ReducedFunctionOnBoundary_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_on_boundary_rank0.nc")
          d=Data(length(ReducedFunctionOnBoundary(self.domain).getX())*self.arg0,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg0,ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunctionOnBoundary_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_on_boundary_rank1.nc")
          d=Data(length(ReducedFunctionOnBoundary(self.domain).getX())*self.arg1,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg1,ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunctionOnBoundary_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_on_boundary_rank2.nc")
          d=Data(length(ReducedFunctionOnBoundary(self.domain).getX())*self.arg2,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg2,ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunctionOnBoundary_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_on_boundary_rank3.nc")
          d=Data(length(ReducedFunctionOnBoundary(self.domain).getX())*self.arg3,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg3,ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   def test_DumpAndLoad_Expanded_ReducedFunctionOnBoundary_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"expanded_reduced_function_on_boundary_rank4.nc")
          d=Data(length(ReducedFunctionOnBoundary(self.domain).getX())*self.arg4,ReducedFunctionOnBoundary(self.domain))
          self._diffDataObjects(d,filemame)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_samples)
          self.assertRaises(RuntimeError, load, filemame, self.domain_with_different_number_of_data_points_per_sample)
          if getMPISizeWorld() ==1:
             d=Data(length(ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering).getX())*self.arg4,ReducedFunctionOnBoundary(self.domain_with_different_sample_ordering))
             self._diffDataObjects(d,filemame, use_old_file=True)

   #===========================================================================
   ## This functionspace does not currently support tags.
   ## Instead, we test that the canTag() function throws in test_canTag_Failures.
   
   #def test_DumpAndLoad_Tagged_Solution_Rank0(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_solution_rank0.nc")
          #d=Data(self.arg0,Solution(self.domain))
          #d.setTaggedValue(1,self.arg0*2)
          #d.setTaggedValue(10,self.arg0*3)
          #d.setTaggedValue(100,self.arg0*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_Solution_Rank1(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_solution_rank1.nc")
          #d=Data(self.arg1,Solution(self.domain))
          #d.setTaggedValue(1,self.arg1*2)
          #d.setTaggedValue(10,self.arg1*3)
          #d.setTaggedValue(100,self.arg1*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_Solution_Rank2(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_solution_rank2.nc")
          #d=Data(self.arg2,Solution(self.domain))
          #d.setTaggedValue(1,self.arg2*2)
          #d.setTaggedValue(10,self.arg2*3)
          #d.setTaggedValue(100,self.arg2*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_Solution_Rank3(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_solution_rank3.nc")
          #d=Data(self.arg3,Solution(self.domain))
          #d.setTaggedValue(1,self.arg3*2)
          #d.setTaggedValue(10,self.arg3*3)
          #d.setTaggedValue(100,self.arg3*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_Solution_Rank4(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_solution_rank4.nc")
          #d=Data(self.arg4,Solution(self.domain))
          #d.setTaggedValue(1,self.arg4*2)
          #d.setTaggedValue(10,self.arg4*3)
          #d.setTaggedValue(100,self.arg4*4)
          #self._diffDataObjects(d,filemame)
   ##===========================================================================
   ## This functionspace does not currently support tags.
   ## Instead, we test that the canTag() function throws in test_canTag_Failures.
   
   #def test_DumpAndLoad_Tagged_ReducedSolution_Rank0(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_reduced_solution_rank0.nc")
          #d=Data(self.arg0,ReducedSolution(self.domain))
          #d.setTaggedValue(1,self.arg0*2)
          #d.setTaggedValue(10,self.arg0*3)
          #d.setTaggedValue(100,self.arg0*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_ReducedSolution_Rank1(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_reduced_solution_rank1.nc")
          #d=Data(self.arg1,ReducedSolution(self.domain))
          #d.setTaggedValue(1,self.arg1*2)
          #d.setTaggedValue(10,self.arg1*3)
          #d.setTaggedValue(100,self.arg1*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_ReducedSolution_Rank2(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_reduced_solution_rank2.nc")
          #d=Data(self.arg2,ReducedSolution(self.domain))
          #d.setTaggedValue(1,self.arg2*2)
          #d.setTaggedValue(10,self.arg2*3)
          #d.setTaggedValue(100,self.arg2*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_ReducedSolution_Rank3(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_reduced_solution_rank3.nc")
          #d=Data(self.arg3,ReducedSolution(self.domain))
          #d.setTaggedValue(1,self.arg3*2)
          #d.setTaggedValue(10,self.arg3*3)
          #d.setTaggedValue(100,self.arg3*4)
          #self._diffDataObjects(d,filemame)

   #def test_DumpAndLoad_Tagged_ReducedSolution_Rank4(self):
       #if loadIsConfigured():
          #filemame=os.path.join(self.filename_base,"tagged_reduced_solution_rank4.nc")
          #d=Data(self.arg4,ReducedSolution(self.domain))
          #d.setTaggedValue(1,self.arg4*2)
          #d.setTaggedValue(10,self.arg4*3)
          #d.setTaggedValue(100,self.arg4*4)
          #self._diffDataObjects(d,filemame)
   ##===========================================================================
   def test_DumpAndLoad_Tagged_ContinuousFunction_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_continuous_function_rank0.nc")
          d=Data(self.arg0,ContinuousFunction(self.domain))
          d.setTaggedValue(1,self.arg0*2)
          d.setTaggedValue(10,self.arg0*3)
          d.setTaggedValue(100,self.arg0*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ContinuousFunction_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_continuous_function_rank1.nc")
          d=Data(self.arg1,ContinuousFunction(self.domain))
          d.setTaggedValue(1,self.arg1*2)
          d.setTaggedValue(10,self.arg1*3)
          d.setTaggedValue(100,self.arg1*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ContinuousFunction_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_continuous_function_rank2.nc")
          d=Data(self.arg2,ContinuousFunction(self.domain))
          d.setTaggedValue(1,self.arg2*2)
          d.setTaggedValue(10,self.arg2*3)
          d.setTaggedValue(100,self.arg2*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ContinuousFunction_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_continuous_function_rank3.nc")
          d=Data(self.arg3,ContinuousFunction(self.domain))
          d.setTaggedValue(1,self.arg3*2)
          d.setTaggedValue(10,self.arg3*3)
          d.setTaggedValue(100,self.arg3*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ContinuousFunction_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_continuous_function_rank4.nc")
          d=Data(self.arg4,ContinuousFunction(self.domain))
          d.setTaggedValue(1,self.arg4*2)
          d.setTaggedValue(10,self.arg4*3)
          d.setTaggedValue(100,self.arg4*4)
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Tagged_Function_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_rank0.nc")
          d=Data(self.arg0,Function(self.domain))
          d.setTaggedValue(1,self.arg0*2)
          d.setTaggedValue(10,self.arg0*3)
          d.setTaggedValue(100,self.arg0*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_Function_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_rank1.nc")
          d=Data(self.arg1,Function(self.domain))
          d.setTaggedValue(1,self.arg1*2)
          d.setTaggedValue(10,self.arg1*3)
          d.setTaggedValue(100,self.arg1*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_Function_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_rank2.nc")
          d=Data(self.arg2,Function(self.domain))
          d.setTaggedValue(1,self.arg2*2)
          d.setTaggedValue(10,self.arg2*3)
          d.setTaggedValue(100,self.arg2*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_Function_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_rank3.nc")
          d=Data(self.arg3,Function(self.domain))
          d.setTaggedValue(1,self.arg3*2)
          d.setTaggedValue(10,self.arg3*3)
          d.setTaggedValue(100,self.arg3*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_Function_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_rank4.nc")
          d=Data(self.arg4,Function(self.domain))
          d.setTaggedValue(1,self.arg4*2)
          d.setTaggedValue(10,self.arg4*3)
          d.setTaggedValue(100,self.arg4*4)
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Tagged_FunctionOnBoundary_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_on_boundary_rank0.nc")
          d=Data(self.arg0,FunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg0*2)
          d.setTaggedValue(10,self.arg0*3)
          d.setTaggedValue(100,self.arg0*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_FunctionOnBoundary_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_on_boundary_rank1.nc")
          d=Data(self.arg1,FunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg1*2)
          d.setTaggedValue(10,self.arg1*3)
          d.setTaggedValue(100,self.arg1*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_FunctionOnBoundary_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_on_boundary_rank2.nc")
          d=Data(self.arg2,FunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg2*2)
          d.setTaggedValue(10,self.arg2*3)
          d.setTaggedValue(100,self.arg2*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_FunctionOnBoundary_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_on_boundary_rank3.nc")
          d=Data(self.arg3,FunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg3*2)
          d.setTaggedValue(10,self.arg3*3)
          d.setTaggedValue(100,self.arg3*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_FunctionOnBoundary_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_function_on_boundary_rank4.nc")
          d=Data(self.arg4,FunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg4*2)
          d.setTaggedValue(10,self.arg4*3)
          d.setTaggedValue(100,self.arg4*4)
          self._diffDataObjects(d,filemame)
   #===========================================================================
   def test_DumpAndLoad_Tagged_ReducedFunction_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_rank0.nc")
          d=Data(self.arg0,ReducedFunction(self.domain))
          d.setTaggedValue(1,self.arg0*2)
          d.setTaggedValue(10,self.arg0*3)
          d.setTaggedValue(100,self.arg0*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunction_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_rank1.nc")
          d=Data(self.arg1,ReducedFunction(self.domain))
          d.setTaggedValue(1,self.arg1*2)
          d.setTaggedValue(10,self.arg1*3)
          d.setTaggedValue(100,self.arg1*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunction_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_rank2.nc")
          d=Data(self.arg2,ReducedFunction(self.domain))
          d.setTaggedValue(1,self.arg2*2)
          d.setTaggedValue(10,self.arg2*3)
          d.setTaggedValue(100,self.arg2*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunction_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_rank3.nc")
          d=Data(self.arg3,ReducedFunction(self.domain))
          d.setTaggedValue(1,self.arg3*2)
          d.setTaggedValue(10,self.arg3*3)
          d.setTaggedValue(100,self.arg3*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunction_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_rank4.nc")
          d=Data(self.arg4,ReducedFunction(self.domain))
          d.setTaggedValue(1,self.arg4*2)
          d.setTaggedValue(10,self.arg4*3)
          d.setTaggedValue(100,self.arg4*4)
          self._diffDataObjects(d,filemame)

   #===========================================================================
   def test_DumpAndLoad_Tagged_ReducedFunctionOnBoundary_Rank0(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_on_boundary_rank0.nc")
          d=Data(self.arg0,ReducedFunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg0*2)
          d.setTaggedValue(10,self.arg0*3)
          d.setTaggedValue(100,self.arg0*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunctionOnBoundary_Rank1(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_on_boundary_rank1.nc")
          d=Data(self.arg1,ReducedFunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg1*2)
          d.setTaggedValue(10,self.arg1*3)
          d.setTaggedValue(100,self.arg1*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunctionOnBoundary_Rank2(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_on_boundary_rank2.nc")
          d=Data(self.arg2,ReducedFunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg2*2)
          d.setTaggedValue(10,self.arg2*3)
          d.setTaggedValue(100,self.arg2*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunctionOnBoundary_Rank3(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_on_boundary_rank3.nc")
          d=Data(self.arg3,ReducedFunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg3*2)
          d.setTaggedValue(10,self.arg3*3)
          d.setTaggedValue(100,self.arg3*4)
          self._diffDataObjects(d,filemame)

   def test_DumpAndLoad_Tagged_ReducedFunctionOnBoundary_Rank4(self):
       if loadIsConfigured():
          filemame=os.path.join(self.filename_base,"tagged_reduced_function_on_boundary_rank4.nc")
          d=Data(self.arg4,ReducedFunctionOnBoundary(self.domain))
          d.setTaggedValue(1,self.arg4*2)
          d.setTaggedValue(10,self.arg4*3)
          d.setTaggedValue(100,self.arg4*4)
          self._diffDataObjects(d,filemame)

class Test_Lazy(unittest.TestCase):
  def makeLazyObj(self):
        d=delay(Data(1,self.mainfs,True))
        e=delay(Data(2,self.mainfs,True))
        p=(d+e*d)/e
        q=p/(3*d)
        r1=q*q
        r2=q+q
        r3=q/4
        f=delay(Data(4,self.otherfs,True))
        t=Data(4,self.mainfs)
        t.tag()
        t=delay(t)
        t=t*2
        return r1,r2,r3,f,t
  
  def test_GroupRes(self):
        rr1,rr2,rr3,rf,rt=self.makeLazyObj()
        rr1=resolve(rr1)
        rr2=resolve(rr2)
        rr3=resolve(rr3)
        rf=resolve(rf)
        rt=resolve(rt)
        r1,r2,r3,f,t=self.makeLazyObj()
        resolveGroup((r1,r2,r3))
        err=Lsup(rr1-r1)+Lsup(rr2-r2)+Lsup(rr3-r3)
        self.assertTrue(err<0.001, "Same functionspace group resolve")
        r1,r2,r3,f,t=self.makeLazyObj()
        resolveGroup((r1,r2,r3,rt))
        err=Lsup(rr1-r1)+Lsup(rr2-r2)+Lsup(rr3-r3)+Lsup(rt-t)
        self.assertTrue(err<0.001, "Same functionspace group resolve with early collapse")
        r1,r2,r3,f,t=self.makeLazyObj()
        err=Lsup(rr1-r1)+Lsup(rr2-r2)+Lsup(rr3-r3)+Lsup(rt-t)+Lsup(rf-f)
        self.assertTrue(err<0.001, "Same functionspace group resolve with mixed functionspaces")
        
        
