
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
Test suite for data objects.

The tests must be linked with some function space class object in the setUp
method to run:

   from esys.finley import Brick
   class Test_DumpOnFinley(Test_Dump):
       def setUp(self):
          self.domain=Rectangle(NE,NE+1,2)
          self.domain_with_different_number_of_samples=Rectangle(2*NE,NE+1,2)
          self.domain_with_different_number_of_data_points_per_sample=Rectangle(2*NE,NE+1,2,integrationOrder=2)
          self.domain_with_different_sample_ordering=Rectangle(1,(NE+1)*NE,2)
          self.filename_base="."

   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_DumpOnFinley))
   unittest.TextTestRunner(verbosity=2).run(suite)

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

from esys.escript import *
import esys.escriptcore.utestselect as unittest
import os


try:
    import numpy
    HAVE_NUMPY = True
except:
    HAVE_NUMPY = False

try:
     ESCRIPT_WORKDIR=os.environ['ESCRIPT_WORKDIR']
except KeyError:
     ESCRIPT_WORKDIR='.'

class Test_tagMap(unittest.TestCase):

    def test_makeTagMap(self):
        for fs in self.functionspaces:
            d=makeTagMap(fs)

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
                    self.assertRaises(RuntimeError, interpolateTable, table, points, (xmin, ymin, zmin), (xwidth, ywidth,
zwidth), -1)
                    # Test bounds checking
                    self.assertRaises(RuntimeError, interpolateTable, table, points, (xmin, ymin, zmin), (xwidth/3, ywidth,
zwidth), 900,True)
                    self.assertRaises(RuntimeError, interpolateTable, table, points, (xmin, ymin, zmin), (xwidth, ywidth/3,
zwidth), 900, True)
                    self.assertRaises(RuntimeError, interpolateTable, table, points, (xmin, ymin, zmin), (xwidth, ywidth,
zwidth/3), 900, True)

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
   def setUp(self):
        self.workdir=ESCRIPT_WORKDIR

   def test_csv_header_separator_and_append(self):
        X=self.domain.getX()
        X0=X[0]
        fname=os.path.join(self.workdir, "test_save1.csv")
        saveDataCSV(fname, C=X, D=X0)
        self.assertTrue(os.path.exists(fname), "test file not created")
        saveDataCSV(fname,append=True, J=X0, H=X)
        self.assertTrue(os.path.exists(fname), "test file not created")
        line=open(fname, 'r').readline()
        # Test both separator strings for vector, scalar
        self.assertEqual(line, "C_0, C_1, D\n")

        # Test Tensor header
        T2=Tensor(7,X.getFunctionSpace())
        T3=Tensor3(8,X.getFunctionSpace())
        T4=Tensor4(9,X.getFunctionSpace())
        saveDataCSV(fname,A=T2,B=T3,C=T4)
        line=open(fname,'r').readline()
        self.assertEqual(line, 'A_0_0, A_1_0, A_0_1, A_1_1, B_0_0_0, B_0_0_1, B_1_0_0, B_1_0_1, B_0_1_0, B_0_1_1, B_1_1_0, B_1_1_1, C_0_0_0_0, C_0_0_0_1, C_0_0_1_0, C_0_0_1_1, C_1_0_0_0, C_1_0_0_1, C_1_0_1_0, C_1_0_1_1, C_0_1_0_0, C_0_1_0_1, C_0_1_1_0, C_0_1_1_1, C_1_1_0_0, C_1_1_0_1, C_1_1_1_0, C_1_1_1_1\n')

        # test different separator
        saveDataCSV(fname, sep="|",csep="/", U=X, V=X0)
        line=open(fname,'r').readline()
        self.assertEqual(line, 'U/0|U/1|V\n')
        MPIBarrierWorld()
        if getMPIRankWorld()==0: os.unlink(fname)

   def test_saveCSV_functionspaces(self):
        for i in range(len(self.functionspaces)):
            FS=self.functionspaces[i]
            X=FS(self.domain).getX()
            X0=X[0]
            fname=os.path.join(self.workdir, "test_save2.csv")
            saveDataCSV(fname,C=X, D=X0)
            f=open(fname,'r')
            # test number of rows written
            linecount=0
            line=f.readline() # skip header
            while line != '':
                linecount+=1
                line=f.readline()
            self.assertEqual(linecount, self.linecounts[i])
            f.close()

            # Now check data content
            T2=Tensor(7, X.getFunctionSpace())
            T3=Tensor3(8, X.getFunctionSpace())
            T4=Tensor4(9, X.getFunctionSpace())
            expected=[7.]*4+[8.]*8+[9.]*16
            saveDataCSV(fname, A=T2, B=T3, C=T4)
            f=open(fname, 'r')
            f.readline() # skip header
            line=f.readline()
            linecount=1
            while line != '':
                line_got=[float(elt) for elt in line.split(',')]
                self.assertEqual(line_got, expected)
                linecount+=1
                line=f.readline()
            self.assertEqual(linecount, self.linecounts[i])
            f.close()

            # As above but with mask and variable data
            saveDataCSV(fname, U=X, V=X0, mask=X0)
            f=open(fname, 'r')
            f.readline() # skip header
            line=f.readline()
            line_got=[float(elt) for elt in line.split(',')]
            self.assertEqual(len(self.firstline[i]),len(line_got))
            for j in range(len(self.firstline[i])):
                if self.firstline[i][j] is not None:
                    self.assertAlmostEqual(self.firstline[i][j],line_got[j])
            linecount=1
            while line!='':
                linecount+=1
                line=f.readline()
            self.assertEqual(linecount, self.linecounts_masked[i])
            f.close()
            MPIBarrierWorld()
            if getMPIRankWorld()==0: os.unlink(fname)


class Test_Domain(unittest.TestCase):

   def test_getListOfTags(self): # requires self.boundary_tag_list
       tags=FunctionOnBoundary(self.domain).getListOfTags()
       self.assertTrue(len(self.boundary_tag_list) == len(tags), "tag list length does not match")
       for i in self.boundary_tag_list:
           self.assertTrue(i in tags, "tag %s is missing."%i)

   def test_RandomData(self):
        fs=Function(self.rdomain)        # The choice of functionspace is arbitrary
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
        self.assertRaises(ValueError,self.domain.getTag,tag3)

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

   @unittest.skipIf(HAVE_NUMPY is False, "Numpy is not installed")
   @unittest.skipIf(getMPISizeWorld()>1, "more than 1 MPI rank")
   def test_getNumpyX(self):
      if hasFeature("boostnumpy"):
         tups=self.domain.getX().toListOfTuples()
         numps=self.domain.getNumpyX()
         for i in range(0,tups.__len__()):
            for x in range(0, self.domain.getDim()):
               self.assertEqual(float(tups[i][x]),float(numps[x][i]))
   #===========================================================================

class Test_SetDataPointValue(unittest.TestCase):
    args=[9.81,
        numpy.array([3.098, -3.111]),
        numpy.array([[3.82, -3.81, -0.957, 0.892, -1.367], [-4.589, -1.835, -2.679, -1.517, -4.2515], [-4.909, 1.634, -2.883,
-2.135, 1.187], [0.6431, 4.638, -4.616, -0.196, -4.370]]),
        numpy.array([[[-2.3667, -0.040], [-4.7398, -3.2412]], [[-2.125, -2.240], [2.237, -4.279]], [[0.68720, 2.4059],
[-2.4964, 3.17453]], [[-4.907, -4.9431], [-0.3604, 0.4269]], [[1.4179, 3.326], [1.356, -0.4610]], [[3.378, 2.0902], [-2.6857,
1.3585]]]),
        numpy.array([[[[-3.810, -1.3597, -1.5307, 1.099], [-1.828, 0.2526, -1.4429, 2.326], [4.9732, -2.063, 1.3153, -3.809]],
[[-4.8902, -4.714, 1.520, -1.931], [-3.8847, 4.3867, 1.894030, 2.432], [-1.2082, -0.8304, 2.2612, 4.6399]]], [[[-4.5922,
-3.309, -0.8171, -0.7210], [2.8051, -4.93047, 0.08450, 4.3824], [0.43204, 2.1908, 4.512633, -1.8218]], [[2.2493, -4.190,
-2.3893, -4.147], [-2.104, -4.635, -4.2767, -3.53151], [-2.351, -1.6614, 2.9385, 4.099]]], [[[1.710, 0.2235, -3.4917, 0.8713],
[-0.2881, 4.6278, 3.603, -2.1211], [-0.565, 4.294, -2.210827, -0.37651]], [[0.6578, -2.869, -2.490, -4.789], [3.232, 2.483,
0.9531, 2.260], [-1.785, 0.42156, -1.8379, 4.212]]]])
        ]
    def test_SetDataPointValue_Function(self):
        for use_list in [False, True]:
            for rank in range(5):
                d=Data(self.args[rank], Function(self.domain))
                d.setValueOfDataPoint(0, self.args[rank]*2)
                param = self.args[1]
                if rank == 1:
                    param = self.args[2]
                for value in [-1, 0]:
                    try:
                        d.setValueOfDataPoint(0, param)
                        self.fail("setting value to %d should have thrown an exception for rank %d"%(value, rank))
                    except RuntimeError:
                        pass
                if rank > 0 and use_list:
                    d.setValueOfDataPoint(0,(self.args[rank]*2).tolist())
                d_0=numpy.array(d.getTupleForDataPoint(0))
                d_1=numpy.array(d.getTupleForDataPoint(1))
                errorstring = "wrong setting for data of rank {0}{1}".format(rank,
                        ", using a list" if rank > 0 and use_list else "")
                self.assertLessEqual(Lsup(d_0-self.args[rank]*2),
                        Lsup(self.args[rank]*2), errorstring)
                self.assertLessEqual(Lsup(d_1-self.args[rank]),
                        Lsup(self.args[rank]), errorstring)


    def test_SetDataPointValue_ReducedFunction(self):
        for use_list in [False, True]:
            for rank in range(5):
                d=Data(self.args[rank], ReducedFunction(self.domain))
                d.setValueOfDataPoint(0,self.args[rank]*2)
                param = self.args[1]
                if rank == 1:
                    param = self.args[2]
                for value in [-1, 0]:
                    try:
                        d.setValueOfDataPoint(0, param)
                        self.fail("setting value to %d should have thrown an exception for rank %d"%(value, rank))
                    except RuntimeError:
                        pass
                if rank > 0 and use_list:
                    d.setValueOfDataPoint(0,(self.args[rank]*2).tolist())
                d_0=numpy.array(d.getTupleForDataPoint(0))
                d_1=numpy.array(d.getTupleForDataPoint(1))
                errorstring = "wrong setting for data of rank {0}{1}".format(rank,
                        ", using a list" if rank > 0 and use_list else "")
                self.assertLessEqual(Lsup(d_0-self.args[rank]*2),
                        Lsup(self.args[rank]*2), errorstring)
                self.assertLessEqual(Lsup(d_1-self.args[rank]),
                        Lsup(self.args[rank]), errorstring)

@unittest.skipIf(not loadIsConfigured(), "load not configured")
class Test_Dump(unittest.TestCase):
   args=[9.81,
        numpy.array([3.098, -3.111]),
        numpy.array([[3.82, -3.81, -0.957, 0.892, -1.367], [-4.589, -1.835, -2.679, -1.517, -4.2515], [-4.909, 1.634, -2.883,
-2.135, 1.187], [0.6431, 4.638, -4.616, -0.196, -4.370]]),
        numpy.array([[[-2.3667, -0.040], [-4.7398, -3.2412]], [[-2.125, -2.240], [2.237, -4.279]], [[0.68720, 2.4059],
[-2.4964, 3.17453]], [[-4.907, -4.9431], [-0.3604, 0.4269]], [[1.4179, 3.326], [1.356, -0.4610]], [[3.378, 2.0902], [-2.6857,
1.3585]]]),
        numpy.array([[[[-3.810, -1.3597, -1.5307, 1.099], [-1.828, 0.2526, -1.4429, 2.326], [4.9732, -2.063, 1.3153, -3.809]],
[[-4.8902, -4.714, 1.520, -1.931], [-3.8847, 4.3867, 1.894030, 2.432], [-1.2082, -0.8304, 2.2612, 4.6399]]], [[[-4.5922,
-3.309, -0.8171, -0.7210], [2.8051, -4.93047, 0.08450, 4.3824], [0.43204, 2.1908, 4.512633, -1.8218]], [[2.2493, -4.190,
-2.3893, -4.147], [-2.104, -4.635, -4.2767, -3.53151], [-2.351, -1.6614, 2.9385, 4.099]]], [[[1.710, 0.2235, -3.4917, 0.8713],
[-0.2881, 4.6278, 3.603, -2.1211], [-0.565, 4.294, -2.210827, -0.37651]], [[0.6578, -2.869, -2.490, -4.789], [3.232, 2.483,
0.9531, 2.260], [-1.785, 0.42156, -1.8379, 4.212]]]])
        ]

   def _diffDataObjects(self,d_ref,filename, use_old_file=False):
       if not use_old_file:
            d_ref.dump(filename)
       d=load(filename, d_ref.getDomain())
       self.assertTrue(not d.isEmpty(),"data in %s are empty."%filename)
       self.assertTrue(d_ref.getRank() == d.getRank(), "different rank in %s. "%filename)
       self.assertTrue(d_ref.getShape() == d.getShape(), "different shape %s. "%filename)
       self.assertTrue(d_ref.getFunctionSpace() == d.getFunctionSpace(), "wrong function space in %s."%filename)
       self.assertTrue(Lsup(d_ref-d)<=0., "different entries %s."%filename)

   #===========================================================================
   def test_DumpAndLoad_Constant(self):
        for functionspace, spacename in [
                (Solution, "solution"),
                (ReducedSolution, "reduced_solution"),
                (ContinuousFunction, "continuous_function"),
                (Function, "function"),
                (ReducedFunction, "reduced_function"),
                (FunctionOnBoundary, "function_on_boundary"),
                (ReducedFunctionOnBoundary, "reduced_function_on_boundary")
            ]:

            for rank in range(5):
                print(spacename, rank)
                filename=os.path.join(self.filename_base,
                        "constant_{0}_rank{1}.h5".format(spacename, rank))
                d=Data(self.args[rank], functionspace(self.domain))
                self._diffDataObjects(d,filename)

   def test_DumpAndLoad_Expanded(self):
        for functionspace, spacename in [
                (Solution, "solution"),
                (ReducedSolution, "reduced_solution"),
                (ContinuousFunction, "continuous_function"),
                (Function, "function"),
                (ReducedFunction, "reduced_function"),
                (FunctionOnBoundary, "function_on_boundary"),
                (ReducedFunctionOnBoundary, "reduced_function_on_boundary")
            ]:

            for rank in range(5):
                filename=os.path.join(self.filename_base,
                        "expanded_{0}_rank{1}.h5".format(spacename, rank))
                d=Data(length(functionspace(self.domain).getX()) * self.args[0],
                        functionspace(self.domain))
                self._diffDataObjects(d,filename)
                self.assertRaises(RuntimeError, load, filename,
                        self.domain_with_different_number_of_samples)
                self.assertRaises(RuntimeError, load, filename,
                        self.domain_with_different_number_of_data_points_per_sample)
                if getMPISizeWorld() ==1:
                    d=Data(length(functionspace(self.domain_with_different_sample_ordering).getX()) *
                        self.args[0], functionspace(self.domain_with_different_sample_ordering))
                    self._diffDataObjects(d, filename, use_old_file=True)

   def test_DumpAndLoad_Tagged(self):
        for functionspace, spacename in [
                #(Solution, "solution"),                #commented in original
                #(ReducedSolution, "reduced_solution"), #commented in original
                (ContinuousFunction, "continuous_function"),
                (Function, "function"),
                (ReducedFunction, "reduced_function"),
                (FunctionOnBoundary, "function_on_boundary"),
                (ReducedFunctionOnBoundary, "reduced_function_on_boundary")
            ]:

            for rank in range(5):
                filename=os.path.join(self.filename_base,
                        "tagged_{0}_rank{1}.h5".format(spacename, rank))
                d=Data(self.args[rank],ContinuousFunction(self.domain))
                d.setTaggedValue(1,self.args[rank]*2)
                d.setTaggedValue(10,self.args[rank]*3)
                d.setTaggedValue(100,self.args[rank]*4)
                self._diffDataObjects(d,filename)

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

  def test_data_getX_Scalar(self): # This tests the Data getXFromFunctionSpace function
        s = Scalar(0, ContinuousFunction(self.domain))
        x1=s.getX()
        x2=s.getFunctionSpace().getX()
        y=x1-x2
        self.assertTrue(inf(y)==0,'Fail test_data_getX_Scalar')
        self.assertTrue(sup(y)==0,'Fail test_data_getX_Scalar')

  def test_complex_data(self):
        s1=Scalar(0,ContinuousFunction(self.domain))
        s2=ComplexScalar(0,ContinuousFunction(self.domain))
        self.assertTrue(s1.isComplex() == False, 'Fail test_complex_data')
        self.assertTrue(s2.isComplex() == True, 'Fail test_complex_data')
        v1=Vector(0,ContinuousFunction(self.domain))
        v2=ComplexVector(0,ContinuousFunction(self.domain))
        self.assertTrue(v1.isComplex() == False, 'Fail test_complex_data')
        self.assertTrue(v2.isComplex() == True, 'Fail test_complex_data')
        t1=Tensor(0,ContinuousFunction(self.domain))
        t2=ComplexTensor(0,ContinuousFunction(self.domain))
        self.assertTrue(t1.isComplex() == False, 'Fail test_complex_data')
        self.assertTrue(t2.isComplex() == True, 'Fail test_complex_data')
        t3=Tensor3(0,ContinuousFunction(self.domain))
        t4=ComplexTensor3(0,ContinuousFunction(self.domain))
        self.assertTrue(t3.isComplex() == False, 'Fail test_complex_data')
        self.assertTrue(t4.isComplex() == True, 'Fail test_complex_data')
        t5=Tensor4(0,ContinuousFunction(self.domain))
        t6=ComplexTensor4(0,ContinuousFunction(self.domain))
        self.assertTrue(t5.isComplex() == False, 'Fail test_complex_data')
        self.assertTrue(t6.isComplex() == True, 'Fail test_complex_data')
        d1=Data(0,ContinuousFunction(self.domain))
        d2=ComplexData(0,ContinuousFunction(self.domain))
        d3=ComplexData(1j,ContinuousFunction(self.domain))
        self.assertTrue(d1.isComplex() == False, 'Fail test_complex_data')
        self.assertTrue(d2.isComplex() == True, 'Fail test_complex_data')
        self.assertTrue(d3.isComplex() == True, 'Fail test_complex_data')
