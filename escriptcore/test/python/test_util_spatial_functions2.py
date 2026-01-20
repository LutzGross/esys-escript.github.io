
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
basic tests for functions in util.py effecting the spatial distribution 

it is assumed that the domain is the usint square/cube 

not all these test will run for all domains. check the doc string for the assumptions of a particular test

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""
__author__="Lutz Gross, l.gross@uq.edu.au"

import esys.escriptcore.utestselect as unittest
from esys.escript import *
from numpy import array
import numpy

from test_util_spatial_functions1 import Test_Util_SpatialFunctions_noGradOnBoundary_noContact
from test_util_integrals import Test_Util_Integration
from test_util_interpolation import Test_Util_Interpolation

#begin

class Test_Util_SpatialFunctions_noGradOnBoundary(Test_Util_SpatialFunctions_noGradOnBoundary_noContact, Test_Util_Integration, Test_Util_Interpolation):
   RES_TOL=1.e-8

   def test_normal_onFunctionOnContactZero(self):
     """
     test getNormal() on contact side 0

     assumptions: FunctionOnContactZero(self.domain) exists
     """
     dim=self.domain.getDim()
     f=FunctionOnContactZero(self.domain)
     x=f.getX()
     ref=Vector(0.,what=f)
     if dim==3:
         ref+=whereZero(x[0]-0.5,tol=self.RES_TOL)*[1,0,0]
     else:
         ref+=whereZero(x[0]-0.5,tol=self.RES_TOL)*[1,0]

     res=f.getNormal()
     self.assertEqual(res.getShape(),(dim,),"wrong shape of result.")
     self.assertEqual(res.getFunctionSpace(),f,"wrong functionspace of result.")
     self.assertTrue(Lsup(length(res)-1)<=self.RES_TOL,"wrong length")
     self.assertTrue(Lsup(abs(inner(ref,res))-1)<=self.RES_TOL,"wrong direction")
   def test_normal_onReducedFunctionOnContactZero(self):
     """
     test getNormal() on contact side 0

     assumptions: FunctionOnContactZero(self.domain) exists
     """
     dim=self.domain.getDim()
     f=ReducedFunctionOnContactZero(self.domain)
     x=f.getX()
     ref=Vector(0.,what=f)
     if dim==3:
         ref+=whereZero(x[0]-0.5,tol=self.RES_TOL)*[1,0,0]
     else:
         ref+=whereZero(x[0]-0.5,tol=self.RES_TOL)*[1,0]

     res=f.getNormal()
     self.assertEqual(res.getShape(),(dim,),"wrong shape of result.")
     self.assertEqual(res.getFunctionSpace(),f,"wrong functionspace of result.")
     self.assertTrue(Lsup(length(res)-1)<=self.RES_TOL,"wrong length")
     self.assertTrue(Lsup(abs(inner(ref,res))-1)<=self.RES_TOL,"wrong direction")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_normal_onFunctionOnContactOne(self):
     """
     test getNormal() on contact side 1

     assumptions: FunctionOnContactOne(self.domain) exists
     """
     dim=self.domain.getDim()
     f=FunctionOnContactOne(self.domain)
     x=f.getX()
     ref=Vector(0.,what=f)
     if dim==3:
         ref+=whereZero(x[0]-0.5,tol=self.RES_TOL)*[-1,0,0]
     else:
         ref+=whereZero(x[0]-0.5,tol=self.RES_TOL)*[-1,0]

     res=f.getNormal()
     self.assertEqual(res.getShape(),(dim,),"wrong shape of result.")
     self.assertEqual(res.getFunctionSpace(),f,"wrong functionspace of result.")
     self.assertTrue(Lsup(length(res)-1)<=self.RES_TOL,"wrong length")
     self.assertTrue(Lsup(abs(inner(ref,res))-1)<=self.RES_TOL,"wrong direction")
   def test_normal_onReducedFunctionOnContactOne(self):
     """
     test getNormal() on contact side 1

     assumptions: FunctionOnContactOne(self.domain) exists
     """
     dim=self.domain.getDim()
     f=ReducedFunctionOnContactOne(self.domain)
     x=f.getX()
     ref=Vector(0.,what=f)
     if dim==3:
         ref+=whereZero(x[0]-0.5,tol=self.RES_TOL)*[-1,0,0]
     else:
         ref+=whereZero(x[0]-0.5,tol=self.RES_TOL)*[-1,0]

     res=f.getNormal()
     self.assertEqual(res.getShape(),(dim,),"wrong shape of result.")
     self.assertEqual(res.getFunctionSpace(),f,"wrong functionspace of result.")
     self.assertTrue(Lsup(length(res)-1)<=self.RES_TOL,"wrong length")
     self.assertTrue(Lsup(abs(inner(ref,res))-1)<=self.RES_TOL,"wrong direction")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunctionOnContactZero_fromData_rank0(self):
      """
      tests L2-norm of Data on the FunctionOnContactZero

      assumptions: self.domain supports integration on FunctionOnContactZero
      """
      dim=self.domain.getDim()
      w=FunctionOnContactZero(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(),w)
        arg=(-0.908822474879)*x[0]
        ref=sqrt((0.0)+(0.206489572711))

      else:
        arg=Data(0,(),w)
        arg=(-0.893168521371)*x[0]
        ref=sqrt((0.0)+(0.199437501892))

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunctionOnContactZero_fromData_rank1(self):
      """
      tests L2-norm of Data on the FunctionOnContactZero

      assumptions: self.domain supports integration on FunctionOnContactZero
      """
      dim=self.domain.getDim()
      w=FunctionOnContactZero(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(2,),w)
        arg[0]=(-0.30839170903)*x[0]
        arg[1]=(0.108352504587)*x[1]
        ref=sqrt((0.00391342175007)+(0.0237763615497))

      else:
        arg=Data(0,(3,),w)
        arg[0]=(-0.260572084359)*x[0]
        arg[1]=(0.250696802346)*x[1]
        arg[2]=(-0.498857619202)*x[2]
        ref=sqrt((0.103902603647)+(0.0169744527868))

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunctionOnContactZero_fromData_rank2(self):
      """
      tests L2-norm of Data on the FunctionOnContactZero

      assumptions: self.domain supports integration on FunctionOnContactZero
      """
      dim=self.domain.getDim()
      w=FunctionOnContactZero(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 2),w)
        arg[0,0]=(0.0107777943992)*x[0]
        arg[0,1]=(0.0188929341636)*x[1]
        arg[1,0]=(-0.00268267834861)*x[0]
        arg[1,1]=(-0.656039971477)*x[1]
        arg[2,0]=(-0.997645982047)*x[0]
        arg[2,1]=(0.993049748017)*x[1]
        arg[3,0]=(0.938636169945)*x[0]
        arg[3,1]=(-0.609892243528)*x[1]
        ref=sqrt((0.596287245963)+(0.46911468066))

      else:
        arg=Data(0,(4, 3),w)
        arg[0,0]=(0.128987155564)*x[0]
        arg[0,1]=(-0.582660210551)*x[1]
        arg[0,2]=(-0.605040149822)*x[2]
        arg[1,0]=(0.494540007291)*x[0]
        arg[1,1]=(0.563900397935)*x[1]
        arg[1,2]=(-0.806214220737)*x[2]
        arg[2,0]=(0.878027993704)*x[0]
        arg[2,1]=(0.597165310427)*x[1]
        arg[2,2]=(-0.444798130742)*x[2]
        arg[3,0]=(0.982067213342)*x[0]
        arg[3,1]=(-0.160013302763)*x[1]
        arg[3,2]=(-0.490830124654)*x[2]
        ref=sqrt((0.831500595262)+(0.49914916859))

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunctionOnContactZero_fromData_rank3(self):
      """
      tests L2-norm of Data on the FunctionOnContactZero

      assumptions: self.domain supports integration on FunctionOnContactZero
      """
      dim=self.domain.getDim()
      w=FunctionOnContactZero(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(6, 2, 2),w)
        arg[0,0,0]=(0.83264345248)*x[0]
        arg[0,0,1]=(-0.29697614845)*x[1]
        arg[0,1,0]=(-0.903150797089)*x[0]
        arg[0,1,1]=(-0.569048463149)*x[1]
        arg[1,0,0]=(0.317087339385)*x[0]
        arg[1,0,1]=(0.0120740865344)*x[1]
        arg[1,1,0]=(-0.221145803457)*x[0]
        arg[1,1,1]=(-0.15362409313)*x[1]
        arg[2,0,0]=(0.459428065932)*x[0]
        arg[2,0,1]=(-0.495488213623)*x[1]
        arg[2,1,0]=(0.919114285084)*x[0]
        arg[2,1,1]=(0.54511942052)*x[1]
        arg[3,0,0]=(0.442207797529)*x[0]
        arg[3,0,1]=(0.958341188449)*x[1]
        arg[3,1,0]=(-0.0915205976682)*x[0]
        arg[3,1,1]=(-0.410964460336)*x[1]
        arg[4,0,0]=(-0.243030095673)*x[0]
        arg[4,0,1]=(-0.413661756935)*x[1]
        arg[4,1,0]=(-0.0973014277464)*x[0]
        arg[4,1,1]=(0.204130676874)*x[1]
        arg[5,0,0]=(0.38083750131)*x[0]
        arg[5,0,1]=(0.286331574579)*x[1]
        arg[5,1,0]=(0.583320730408)*x[0]
        arg[5,1,1]=(-0.0899419806181)*x[1]
        ref=sqrt((0.789530406063)+(0.868006693351))

      else:
        arg=Data(0,(6, 2, 3),w)
        arg[0,0,0]=(-0.98411764268)*x[0]
        arg[0,0,1]=(-0.180549567137)*x[1]
        arg[0,0,2]=(0.357923974709)*x[2]
        arg[0,1,0]=(0.67691778218)*x[0]
        arg[0,1,1]=(-0.643267669153)*x[1]
        arg[0,1,2]=(0.452374328301)*x[2]
        arg[1,0,0]=(0.494303694439)*x[0]
        arg[1,0,1]=(0.00767629822214)*x[1]
        arg[1,0,2]=(0.232436713422)*x[2]
        arg[1,1,0]=(-0.698568975983)*x[0]
        arg[1,1,1]=(-0.08483720353)*x[1]
        arg[1,1,2]=(-0.723935799719)*x[2]
        arg[2,0,0]=(0.931313239646)*x[0]
        arg[2,0,1]=(-0.0988070442341)*x[1]
        arg[2,0,2]=(-0.166201798278)*x[2]
        arg[2,1,0]=(-0.733172171406)*x[0]
        arg[2,1,1]=(0.24542706584)*x[1]
        arg[2,1,2]=(0.49170717592)*x[2]
        arg[3,0,0]=(-0.718346096869)*x[0]
        arg[3,0,1]=(0.297442410539)*x[1]
        arg[3,0,2]=(0.0277708718174)*x[2]
        arg[3,1,0]=(0.735100702194)*x[0]
        arg[3,1,1]=(-0.808007526296)*x[1]
        arg[3,1,2]=(-0.423023238313)*x[2]
        arg[4,0,0]=(0.383188529876)*x[0]
        arg[4,0,1]=(-0.190101601814)*x[1]
        arg[4,0,2]=(-0.404380990215)*x[2]
        arg[4,1,0]=(-0.31098204217)*x[0]
        arg[4,1,1]=(-0.715247945432)*x[1]
        arg[4,1,2]=(0.435774772982)*x[2]
        arg[5,0,0]=(-0.0662056096206)*x[0]
        arg[5,0,1]=(0.746710247367)*x[1]
        arg[5,0,2]=(0.753267249301)*x[2]
        arg[5,1,0]=(0.258102962518)*x[0]
        arg[5,1,1]=(-0.0270028483187)*x[1]
        arg[5,1,2]=(-0.971060107954)*x[2]
        ref=sqrt((1.86493016409)+(1.23371587938))

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunctionOnContactZero_fromData_rank4(self):
      """
      tests L2-norm of Data on the FunctionOnContactZero

      assumptions: self.domain supports integration on FunctionOnContactZero
      """
      dim=self.domain.getDim()
      w=FunctionOnContactZero(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 5, 3, 2),w)
        arg[0,0,0,0]=(-0.144741217139)*x[0]
        arg[0,0,0,1]=(-0.188187241012)*x[1]
        arg[0,0,1,0]=(0.056339373733)*x[0]
        arg[0,0,1,1]=(0.222177130759)*x[1]
        arg[0,0,2,0]=(-0.437240978497)*x[0]
        arg[0,0,2,1]=(-0.479213082014)*x[1]
        arg[0,1,0,0]=(0.343204141181)*x[0]
        arg[0,1,0,1]=(0.099683148208)*x[1]
        arg[0,1,1,0]=(0.167657906331)*x[0]
        arg[0,1,1,1]=(0.371064917656)*x[1]
        arg[0,1,2,0]=(-0.175796118619)*x[0]
        arg[0,1,2,1]=(-0.523477567832)*x[1]
        arg[0,2,0,0]=(0.188832282277)*x[0]
        arg[0,2,0,1]=(0.948618356423)*x[1]
        arg[0,2,1,0]=(-0.3835066594)*x[0]
        arg[0,2,1,1]=(0.872875257641)*x[1]
        arg[0,2,2,0]=(-0.98414676194)*x[0]
        arg[0,2,2,1]=(0.59170993756)*x[1]
        arg[0,3,0,0]=(0.875820385471)*x[0]
        arg[0,3,0,1]=(-0.209139794169)*x[1]
        arg[0,3,1,0]=(0.783312979835)*x[0]
        arg[0,3,1,1]=(-0.440048708153)*x[1]
        arg[0,3,2,0]=(-0.906741480158)*x[0]
        arg[0,3,2,1]=(0.63834196559)*x[1]
        arg[0,4,0,0]=(0.0142616166715)*x[0]
        arg[0,4,0,1]=(-0.0745581409622)*x[1]
        arg[0,4,1,0]=(-0.779966728505)*x[0]
        arg[0,4,1,1]=(0.734744791584)*x[1]
        arg[0,4,2,0]=(0.0407627699446)*x[0]
        arg[0,4,2,1]=(-0.629666423098)*x[1]
        arg[1,0,0,0]=(-0.782291768803)*x[0]
        arg[1,0,0,1]=(-0.197799867335)*x[1]
        arg[1,0,1,0]=(-0.487323212409)*x[0]
        arg[1,0,1,1]=(-0.789940409633)*x[1]
        arg[1,0,2,0]=(-0.523871448728)*x[0]
        arg[1,0,2,1]=(-0.622272488658)*x[1]
        arg[1,1,0,0]=(0.370095832941)*x[0]
        arg[1,1,0,1]=(-0.0683624906905)*x[1]
        arg[1,1,1,0]=(-0.180583695959)*x[0]
        arg[1,1,1,1]=(-0.0529097308814)*x[1]
        arg[1,1,2,0]=(0.440088045764)*x[0]
        arg[1,1,2,1]=(-0.207665675854)*x[1]
        arg[1,2,0,0]=(-0.360120735448)*x[0]
        arg[1,2,0,1]=(0.877797020694)*x[1]
        arg[1,2,1,0]=(0.655192721419)*x[0]
        arg[1,2,1,1]=(0.978781227463)*x[1]
        arg[1,2,2,0]=(-0.205886258554)*x[0]
        arg[1,2,2,1]=(-0.0050396100363)*x[1]
        arg[1,3,0,0]=(-0.484332943177)*x[0]
        arg[1,3,0,1]=(0.964697361018)*x[1]
        arg[1,3,1,0]=(0.0815652906824)*x[0]
        arg[1,3,1,1]=(-0.531163025802)*x[1]
        arg[1,3,2,0]=(-0.229615479833)*x[0]
        arg[1,3,2,1]=(0.139971291609)*x[1]
        arg[1,4,0,0]=(0.960839193726)*x[0]
        arg[1,4,0,1]=(-0.9744930023)*x[1]
        arg[1,4,1,0]=(-0.0919000363634)*x[0]
        arg[1,4,1,1]=(0.483909186003)*x[1]
        arg[1,4,2,0]=(0.248597552464)*x[0]
        arg[1,4,2,1]=(0.965872402486)*x[1]
        arg[2,0,0,0]=(-0.935049188033)*x[0]
        arg[2,0,0,1]=(-0.823239964421)*x[1]
        arg[2,0,1,0]=(0.30901193713)*x[0]
        arg[2,0,1,1]=(-0.0538396976182)*x[1]
        arg[2,0,2,0]=(-0.449416473099)*x[0]
        arg[2,0,2,1]=(-0.857785267819)*x[1]
        arg[2,1,0,0]=(-0.167802578334)*x[0]
        arg[2,1,0,1]=(-0.7016101186)*x[1]
        arg[2,1,1,0]=(0.995848317929)*x[0]
        arg[2,1,1,1]=(-0.824595018011)*x[1]
        arg[2,1,2,0]=(-0.826030527686)*x[0]
        arg[2,1,2,1]=(-0.856151725651)*x[1]
        arg[2,2,0,0]=(0.342047605028)*x[0]
        arg[2,2,0,1]=(-0.342370307931)*x[1]
        arg[2,2,1,0]=(-0.540693890831)*x[0]
        arg[2,2,1,1]=(-0.720754897785)*x[1]
        arg[2,2,2,0]=(-0.0941215154824)*x[0]
        arg[2,2,2,1]=(0.115887914141)*x[1]
        arg[2,3,0,0]=(-0.0920059639279)*x[0]
        arg[2,3,0,1]=(0.910948619784)*x[1]
        arg[2,3,1,0]=(-0.40846609126)*x[0]
        arg[2,3,1,1]=(0.542728733954)*x[1]
        arg[2,3,2,0]=(0.451361830368)*x[0]
        arg[2,3,2,1]=(0.0767538617936)*x[1]
        arg[2,4,0,0]=(-0.431757008766)*x[0]
        arg[2,4,0,1]=(0.258988103513)*x[1]
        arg[2,4,1,0]=(-0.584559859464)*x[0]
        arg[2,4,1,1]=(-0.776827606804)*x[1]
        arg[2,4,2,0]=(-0.867078960637)*x[0]
        arg[2,4,2,1]=(-0.0972960135129)*x[1]
        arg[3,0,0,0]=(0.999511727022)*x[0]
        arg[3,0,0,1]=(-0.753218704622)*x[1]
        arg[3,0,1,0]=(0.506439354885)*x[0]
        arg[3,0,1,1]=(-0.0845853282297)*x[1]
        arg[3,0,2,0]=(-0.871475028801)*x[0]
        arg[3,0,2,1]=(-0.69261642657)*x[1]
        arg[3,1,0,0]=(0.549778399933)*x[0]
        arg[3,1,0,1]=(0.246187536149)*x[1]
        arg[3,1,1,0]=(-0.620454200676)*x[0]
        arg[3,1,1,1]=(0.372738564315)*x[1]
        arg[3,1,2,0]=(0.0544097075138)*x[0]
        arg[3,1,2,1]=(-0.0883198944676)*x[1]
        arg[3,2,0,0]=(-0.671057180545)*x[0]
        arg[3,2,0,1]=(-0.118991797907)*x[1]
        arg[3,2,1,0]=(-0.27196730081)*x[0]
        arg[3,2,1,1]=(-0.458807068968)*x[1]
        arg[3,2,2,0]=(-0.31485399083)*x[0]
        arg[3,2,2,1]=(0.487291960328)*x[1]
        arg[3,3,0,0]=(-0.750302723531)*x[0]
        arg[3,3,0,1]=(-0.486428789771)*x[1]
        arg[3,3,1,0]=(-0.295909324594)*x[0]
        arg[3,3,1,1]=(0.325702372596)*x[1]
        arg[3,3,2,0]=(-0.512109540631)*x[0]
        arg[3,3,2,1]=(0.563284721908)*x[1]
        arg[3,4,0,0]=(0.53579406798)*x[0]
        arg[3,4,0,1]=(-0.468431927975)*x[1]
        arg[3,4,1,0]=(0.570851401329)*x[0]
        arg[3,4,1,1]=(0.107791149865)*x[1]
        arg[3,4,2,0]=(-0.211543670157)*x[0]
        arg[3,4,2,1]=(0.848189550468)*x[1]
        ref=sqrt((6.54316712884)+(4.39564064326))

      else:
        arg=Data(0,(4, 5, 3, 3),w)
        arg[0,0,0,0]=(-0.392256494872)*x[0]
        arg[0,0,0,1]=(0.271488478841)*x[1]
        arg[0,0,0,2]=(-0.878268531987)*x[2]
        arg[0,0,1,0]=(0.106903787643)*x[0]
        arg[0,0,1,1]=(-0.117986434516)*x[1]
        arg[0,0,1,2]=(0.912673598522)*x[2]
        arg[0,0,2,0]=(-0.652598945824)*x[0]
        arg[0,0,2,1]=(-0.984195895647)*x[1]
        arg[0,0,2,2]=(-0.137845459589)*x[2]
        arg[0,1,0,0]=(0.696713200721)*x[0]
        arg[0,1,0,1]=(0.480337929223)*x[1]
        arg[0,1,0,2]=(-0.628726114325)*x[2]
        arg[0,1,1,0]=(-0.356503461347)*x[0]
        arg[0,1,1,1]=(-0.64670584045)*x[1]
        arg[0,1,1,2]=(-0.737787618908)*x[2]
        arg[0,1,2,0]=(-0.14207595302)*x[0]
        arg[0,1,2,1]=(-0.573485525334)*x[1]
        arg[0,1,2,2]=(-0.955542178481)*x[2]
        arg[0,2,0,0]=(0.371952883975)*x[0]
        arg[0,2,0,1]=(0.114016178762)*x[1]
        arg[0,2,0,2]=(0.552721218169)*x[2]
        arg[0,2,1,0]=(0.318426742113)*x[0]
        arg[0,2,1,1]=(0.862220239384)*x[1]
        arg[0,2,1,2]=(0.887376889978)*x[2]
        arg[0,2,2,0]=(0.203656742981)*x[0]
        arg[0,2,2,1]=(0.350541335144)*x[1]
        arg[0,2,2,2]=(-0.448401957447)*x[2]
        arg[0,3,0,0]=(-0.349200084042)*x[0]
        arg[0,3,0,1]=(-0.546012602125)*x[1]
        arg[0,3,0,2]=(-0.931488270979)*x[2]
        arg[0,3,1,0]=(0.207457831058)*x[0]
        arg[0,3,1,1]=(0.557715840549)*x[1]
        arg[0,3,1,2]=(-0.978130146744)*x[2]
        arg[0,3,2,0]=(-0.55088967957)*x[0]
        arg[0,3,2,1]=(0.0490328838057)*x[1]
        arg[0,3,2,2]=(0.150209929122)*x[2]
        arg[0,4,0,0]=(-0.484145614698)*x[0]
        arg[0,4,0,1]=(0.393078411279)*x[1]
        arg[0,4,0,2]=(0.0678631863917)*x[2]
        arg[0,4,1,0]=(-0.350980464628)*x[0]
        arg[0,4,1,1]=(-0.784203839564)*x[1]
        arg[0,4,1,2]=(0.636960296147)*x[2]
        arg[0,4,2,0]=(0.592799503581)*x[0]
        arg[0,4,2,1]=(-0.672104833683)*x[1]
        arg[0,4,2,2]=(0.0366914082467)*x[2]
        arg[1,0,0,0]=(-0.147168019774)*x[0]
        arg[1,0,0,1]=(-0.0823637938956)*x[1]
        arg[1,0,0,2]=(-0.852729690176)*x[2]
        arg[1,0,1,0]=(-0.612338430408)*x[0]
        arg[1,0,1,1]=(-0.85820035747)*x[1]
        arg[1,0,1,2]=(-0.463664966162)*x[2]
        arg[1,0,2,0]=(0.274600720491)*x[0]
        arg[1,0,2,1]=(-0.488508234093)*x[1]
        arg[1,0,2,2]=(-0.28251466519)*x[2]
        arg[1,1,0,0]=(-0.0196532462794)*x[0]
        arg[1,1,0,1]=(0.239352528871)*x[1]
        arg[1,1,0,2]=(-0.17054773873)*x[2]
        arg[1,1,1,0]=(-0.2396627789)*x[0]
        arg[1,1,1,1]=(0.868970323003)*x[1]
        arg[1,1,1,2]=(0.401487430312)*x[2]
        arg[1,1,2,0]=(-0.624411449783)*x[0]
        arg[1,1,2,1]=(0.0036634266684)*x[1]
        arg[1,1,2,2]=(0.736129120967)*x[2]
        arg[1,2,0,0]=(0.183756511707)*x[0]
        arg[1,2,0,1]=(-0.288651848639)*x[1]
        arg[1,2,0,2]=(-0.0672121537447)*x[2]
        arg[1,2,1,0]=(-0.323274725936)*x[0]
        arg[1,2,1,1]=(0.298001016025)*x[1]
        arg[1,2,1,2]=(-0.976052460675)*x[2]
        arg[1,2,2,0]=(0.596504441096)*x[0]
        arg[1,2,2,1]=(0.873776068983)*x[1]
        arg[1,2,2,2]=(-0.994068273196)*x[2]
        arg[1,3,0,0]=(-0.495387299681)*x[0]
        arg[1,3,0,1]=(-0.123674756551)*x[1]
        arg[1,3,0,2]=(0.581213818577)*x[2]
        arg[1,3,1,0]=(0.146405749701)*x[0]
        arg[1,3,1,1]=(-0.594994686675)*x[1]
        arg[1,3,1,2]=(-0.059093568436)*x[2]
        arg[1,3,2,0]=(0.651004255104)*x[0]
        arg[1,3,2,1]=(-0.977880706193)*x[1]
        arg[1,3,2,2]=(0.370344651319)*x[2]
        arg[1,4,0,0]=(-0.503657215247)*x[0]
        arg[1,4,0,1]=(-0.170885297253)*x[1]
        arg[1,4,0,2]=(0.533424480956)*x[2]
        arg[1,4,1,0]=(-0.220533193308)*x[0]
        arg[1,4,1,1]=(0.344537611882)*x[1]
        arg[1,4,1,2]=(0.861473877282)*x[2]
        arg[1,4,2,0]=(-0.0923010438884)*x[0]
        arg[1,4,2,1]=(-0.338256780498)*x[1]
        arg[1,4,2,2]=(0.528567959345)*x[2]
        arg[2,0,0,0]=(-0.0423053381485)*x[0]
        arg[2,0,0,1]=(0.856798579151)*x[1]
        arg[2,0,0,2]=(0.383258153853)*x[2]
        arg[2,0,1,0]=(0.350994872736)*x[0]
        arg[2,0,1,1]=(-0.78055158106)*x[1]
        arg[2,0,1,2]=(-0.770876699915)*x[2]
        arg[2,0,2,0]=(-0.935133287106)*x[0]
        arg[2,0,2,1]=(0.618238076989)*x[1]
        arg[2,0,2,2]=(-0.846783087949)*x[2]
        arg[2,1,0,0]=(-0.657394511405)*x[0]
        arg[2,1,0,1]=(0.576218821654)*x[1]
        arg[2,1,0,2]=(0.0269446356493)*x[2]
        arg[2,1,1,0]=(-0.310710230949)*x[0]
        arg[2,1,1,1]=(0.425412515598)*x[1]
        arg[2,1,1,2]=(-0.225225290862)*x[2]
        arg[2,1,2,0]=(-0.539928589495)*x[0]
        arg[2,1,2,1]=(-0.348082121765)*x[1]
        arg[2,1,2,2]=(-0.0287274646233)*x[2]
        arg[2,2,0,0]=(0.746132865923)*x[0]
        arg[2,2,0,1]=(-0.0234203693548)*x[1]
        arg[2,2,0,2]=(0.517411821941)*x[2]
        arg[2,2,1,0]=(-0.183204217349)*x[0]
        arg[2,2,1,1]=(0.714988861836)*x[1]
        arg[2,2,1,2]=(0.829083318937)*x[2]
        arg[2,2,2,0]=(0.458067555841)*x[0]
        arg[2,2,2,1]=(0.639317125869)*x[1]
        arg[2,2,2,2]=(0.104611520408)*x[2]
        arg[2,3,0,0]=(-0.420513461135)*x[0]
        arg[2,3,0,1]=(0.888686162754)*x[1]
        arg[2,3,0,2]=(0.939305777879)*x[2]
        arg[2,3,1,0]=(0.856795132015)*x[0]
        arg[2,3,1,1]=(0.817593141895)*x[1]
        arg[2,3,1,2]=(0.962503342535)*x[2]
        arg[2,3,2,0]=(-0.334586700245)*x[0]
        arg[2,3,2,1]=(0.182696129528)*x[1]
        arg[2,3,2,2]=(-0.707271571206)*x[2]
        arg[2,4,0,0]=(-0.594795981069)*x[0]
        arg[2,4,0,1]=(0.451239168073)*x[1]
        arg[2,4,0,2]=(0.191212211556)*x[2]
        arg[2,4,1,0]=(-0.503465984944)*x[0]
        arg[2,4,1,1]=(0.725377884208)*x[1]
        arg[2,4,1,2]=(-0.40719255752)*x[2]
        arg[2,4,2,0]=(0.268071476451)*x[0]
        arg[2,4,2,1]=(0.85066639942)*x[1]
        arg[2,4,2,2]=(-0.906021406945)*x[2]
        arg[3,0,0,0]=(-0.922179152122)*x[0]
        arg[3,0,0,1]=(-0.0903841240007)*x[1]
        arg[3,0,0,2]=(-0.751482803516)*x[2]
        arg[3,0,1,0]=(0.960697809119)*x[0]
        arg[3,0,1,1]=(0.638878873158)*x[1]
        arg[3,0,1,2]=(0.390932234724)*x[2]
        arg[3,0,2,0]=(-0.925078301694)*x[0]
        arg[3,0,2,1]=(0.793590580665)*x[1]
        arg[3,0,2,2]=(0.535478366911)*x[2]
        arg[3,1,0,0]=(-0.431951993217)*x[0]
        arg[3,1,0,1]=(0.211750261417)*x[1]
        arg[3,1,0,2]=(-0.930706580442)*x[2]
        arg[3,1,1,0]=(0.330979313323)*x[0]
        arg[3,1,1,1]=(-0.838919076081)*x[1]
        arg[3,1,1,2]=(0.134250050168)*x[2]
        arg[3,1,2,0]=(0.414922811301)*x[0]
        arg[3,1,2,1]=(-0.663692878121)*x[1]
        arg[3,1,2,2]=(0.88499278543)*x[2]
        arg[3,2,0,0]=(0.0884662742233)*x[0]
        arg[3,2,0,1]=(-0.412630722821)*x[1]
        arg[3,2,0,2]=(-0.730850884928)*x[2]
        arg[3,2,1,0]=(0.722207547366)*x[0]
        arg[3,2,1,1]=(-0.260067950749)*x[1]
        arg[3,2,1,2]=(0.426259201494)*x[2]
        arg[3,2,2,0]=(0.0516111795322)*x[0]
        arg[3,2,2,1]=(0.922853710048)*x[1]
        arg[3,2,2,2]=(-0.991912758116)*x[2]
        arg[3,3,0,0]=(0.263933965905)*x[0]
        arg[3,3,0,1]=(0.840541758799)*x[1]
        arg[3,3,0,2]=(0.417658511125)*x[2]
        arg[3,3,1,0]=(-0.901745614723)*x[0]
        arg[3,3,1,1]=(-0.623608908699)*x[1]
        arg[3,3,1,2]=(0.0522167208784)*x[2]
        arg[3,3,2,0]=(-0.549431264931)*x[0]
        arg[3,3,2,1]=(0.919971855457)*x[1]
        arg[3,3,2,2]=(0.142757773397)*x[2]
        arg[3,4,0,0]=(-0.258875259824)*x[0]
        arg[3,4,0,1]=(-0.0373872187041)*x[1]
        arg[3,4,0,2]=(0.445989164864)*x[2]
        arg[3,4,1,0]=(0.0273711397038)*x[0]
        arg[3,4,1,1]=(-0.85522629)*x[1]
        arg[3,4,1,2]=(-0.835392581226)*x[2]
        arg[3,4,2,0]=(0.402929450189)*x[0]
        arg[3,4,2,1]=(0.549623033221)*x[1]
        arg[3,4,2,2]=(0.237803212109)*x[2]
        ref=sqrt((15.1448100099)+(3.61598249819))

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunctionOnContactZero_fromData_rank0(self):
      """
      tests L2-norm of Data on the ReducedFunctionOnContactZero

      assumptions: self.domain supports integration on ReducedFunctionOnContactZero
      """
      dim=self.domain.getDim()
      w=ReducedFunctionOnContactZero(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(),w)
        arg=(-0.908822474879)*sqrt(x[0])
        ref=sqrt((0.0)*3./2.+(0.206489572711)/0.5)

      else:
        arg=Data(0,(),w)
        arg=(-0.893168521371)*sqrt(x[0])
        ref=sqrt((0.0)*3./2.+(0.199437501892)/0.5)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunctionOnContactZero_fromData_rank1(self):
      """
      tests L2-norm of Data on the ReducedFunctionOnContactZero

      assumptions: self.domain supports integration on ReducedFunctionOnContactZero
      """
      dim=self.domain.getDim()
      w=ReducedFunctionOnContactZero(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(2,),w)
        arg[0]=(-0.30839170903)*sqrt(x[0])
        arg[1]=(0.108352504587)*sqrt(x[1])
        ref=sqrt((0.00391342175007)*3./2.+(0.0237763615497)/0.5)

      else:
        arg=Data(0,(3,),w)
        arg[0]=(-0.260572084359)*sqrt(x[0])
        arg[1]=(0.250696802346)*sqrt(x[1])
        arg[2]=(-0.498857619202)*sqrt(x[2])
        ref=sqrt((0.103902603647)*3./2.+(0.0169744527868)/0.5)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunctionOnContactZero_fromData_rank2(self):
      """
      tests L2-norm of Data on the ReducedFunctionOnContactZero

      assumptions: self.domain supports integration on ReducedFunctionOnContactZero
      """
      dim=self.domain.getDim()
      w=ReducedFunctionOnContactZero(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 2),w)
        arg[0,0]=(0.0107777943992)*sqrt(x[0])
        arg[0,1]=(0.0188929341636)*sqrt(x[1])
        arg[1,0]=(-0.00268267834861)*sqrt(x[0])
        arg[1,1]=(-0.656039971477)*sqrt(x[1])
        arg[2,0]=(-0.997645982047)*sqrt(x[0])
        arg[2,1]=(0.993049748017)*sqrt(x[1])
        arg[3,0]=(0.938636169945)*sqrt(x[0])
        arg[3,1]=(-0.609892243528)*sqrt(x[1])
        ref=sqrt((0.596287245963)*3./2.+(0.46911468066)/0.5)

      else:
        arg=Data(0,(4, 3),w)
        arg[0,0]=(0.128987155564)*sqrt(x[0])
        arg[0,1]=(-0.582660210551)*sqrt(x[1])
        arg[0,2]=(-0.605040149822)*sqrt(x[2])
        arg[1,0]=(0.494540007291)*sqrt(x[0])
        arg[1,1]=(0.563900397935)*sqrt(x[1])
        arg[1,2]=(-0.806214220737)*sqrt(x[2])
        arg[2,0]=(0.878027993704)*sqrt(x[0])
        arg[2,1]=(0.597165310427)*sqrt(x[1])
        arg[2,2]=(-0.444798130742)*sqrt(x[2])
        arg[3,0]=(0.982067213342)*sqrt(x[0])
        arg[3,1]=(-0.160013302763)*sqrt(x[1])
        arg[3,2]=(-0.490830124654)*sqrt(x[2])
        ref=sqrt((0.831500595262)*3./2.+(0.49914916859)/0.5)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunctionOnContactZero_fromData_rank3(self):
      """
      tests L2-norm of Data on the ReducedFunctionOnContactZero

      assumptions: self.domain supports integration on ReducedFunctionOnContactZero
      """
      dim=self.domain.getDim()
      w=ReducedFunctionOnContactZero(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(6, 2, 2),w)
        arg[0,0,0]=(0.83264345248)*sqrt(x[0])
        arg[0,0,1]=(-0.29697614845)*sqrt(x[1])
        arg[0,1,0]=(-0.903150797089)*sqrt(x[0])
        arg[0,1,1]=(-0.569048463149)*sqrt(x[1])
        arg[1,0,0]=(0.317087339385)*sqrt(x[0])
        arg[1,0,1]=(0.0120740865344)*sqrt(x[1])
        arg[1,1,0]=(-0.221145803457)*sqrt(x[0])
        arg[1,1,1]=(-0.15362409313)*sqrt(x[1])
        arg[2,0,0]=(0.459428065932)*sqrt(x[0])
        arg[2,0,1]=(-0.495488213623)*sqrt(x[1])
        arg[2,1,0]=(0.919114285084)*sqrt(x[0])
        arg[2,1,1]=(0.54511942052)*sqrt(x[1])
        arg[3,0,0]=(0.442207797529)*sqrt(x[0])
        arg[3,0,1]=(0.958341188449)*sqrt(x[1])
        arg[3,1,0]=(-0.0915205976682)*sqrt(x[0])
        arg[3,1,1]=(-0.410964460336)*sqrt(x[1])
        arg[4,0,0]=(-0.243030095673)*sqrt(x[0])
        arg[4,0,1]=(-0.413661756935)*sqrt(x[1])
        arg[4,1,0]=(-0.0973014277464)*sqrt(x[0])
        arg[4,1,1]=(0.204130676874)*sqrt(x[1])
        arg[5,0,0]=(0.38083750131)*sqrt(x[0])
        arg[5,0,1]=(0.286331574579)*sqrt(x[1])
        arg[5,1,0]=(0.583320730408)*sqrt(x[0])
        arg[5,1,1]=(-0.0899419806181)*sqrt(x[1])
        ref=sqrt((0.789530406063)*3./2.+(0.868006693351)/0.5)

      else:
        arg=Data(0,(6, 2, 3),w)
        arg[0,0,0]=(-0.98411764268)*sqrt(x[0])
        arg[0,0,1]=(-0.180549567137)*sqrt(x[1])
        arg[0,0,2]=(0.357923974709)*sqrt(x[2])
        arg[0,1,0]=(0.67691778218)*sqrt(x[0])
        arg[0,1,1]=(-0.643267669153)*sqrt(x[1])
        arg[0,1,2]=(0.452374328301)*sqrt(x[2])
        arg[1,0,0]=(0.494303694439)*sqrt(x[0])
        arg[1,0,1]=(0.00767629822214)*sqrt(x[1])
        arg[1,0,2]=(0.232436713422)*sqrt(x[2])
        arg[1,1,0]=(-0.698568975983)*sqrt(x[0])
        arg[1,1,1]=(-0.08483720353)*sqrt(x[1])
        arg[1,1,2]=(-0.723935799719)*sqrt(x[2])
        arg[2,0,0]=(0.931313239646)*sqrt(x[0])
        arg[2,0,1]=(-0.0988070442341)*sqrt(x[1])
        arg[2,0,2]=(-0.166201798278)*sqrt(x[2])
        arg[2,1,0]=(-0.733172171406)*sqrt(x[0])
        arg[2,1,1]=(0.24542706584)*sqrt(x[1])
        arg[2,1,2]=(0.49170717592)*sqrt(x[2])
        arg[3,0,0]=(-0.718346096869)*sqrt(x[0])
        arg[3,0,1]=(0.297442410539)*sqrt(x[1])
        arg[3,0,2]=(0.0277708718174)*sqrt(x[2])
        arg[3,1,0]=(0.735100702194)*sqrt(x[0])
        arg[3,1,1]=(-0.808007526296)*sqrt(x[1])
        arg[3,1,2]=(-0.423023238313)*sqrt(x[2])
        arg[4,0,0]=(0.383188529876)*sqrt(x[0])
        arg[4,0,1]=(-0.190101601814)*sqrt(x[1])
        arg[4,0,2]=(-0.404380990215)*sqrt(x[2])
        arg[4,1,0]=(-0.31098204217)*sqrt(x[0])
        arg[4,1,1]=(-0.715247945432)*sqrt(x[1])
        arg[4,1,2]=(0.435774772982)*sqrt(x[2])
        arg[5,0,0]=(-0.0662056096206)*sqrt(x[0])
        arg[5,0,1]=(0.746710247367)*sqrt(x[1])
        arg[5,0,2]=(0.753267249301)*sqrt(x[2])
        arg[5,1,0]=(0.258102962518)*sqrt(x[0])
        arg[5,1,1]=(-0.0270028483187)*sqrt(x[1])
        arg[5,1,2]=(-0.971060107954)*sqrt(x[2])
        ref=sqrt((1.86493016409)*3./2.+(1.23371587938)/0.5)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunctionOnContactZero_fromData_rank4(self):
      """
      tests L2-norm of Data on the ReducedFunctionOnContactZero

      assumptions: self.domain supports integration on ReducedFunctionOnContactZero
      """
      dim=self.domain.getDim()
      w=ReducedFunctionOnContactZero(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 5, 3, 2),w)
        arg[0,0,0,0]=(-0.144741217139)*sqrt(x[0])
        arg[0,0,0,1]=(-0.188187241012)*sqrt(x[1])
        arg[0,0,1,0]=(0.056339373733)*sqrt(x[0])
        arg[0,0,1,1]=(0.222177130759)*sqrt(x[1])
        arg[0,0,2,0]=(-0.437240978497)*sqrt(x[0])
        arg[0,0,2,1]=(-0.479213082014)*sqrt(x[1])
        arg[0,1,0,0]=(0.343204141181)*sqrt(x[0])
        arg[0,1,0,1]=(0.099683148208)*sqrt(x[1])
        arg[0,1,1,0]=(0.167657906331)*sqrt(x[0])
        arg[0,1,1,1]=(0.371064917656)*sqrt(x[1])
        arg[0,1,2,0]=(-0.175796118619)*sqrt(x[0])
        arg[0,1,2,1]=(-0.523477567832)*sqrt(x[1])
        arg[0,2,0,0]=(0.188832282277)*sqrt(x[0])
        arg[0,2,0,1]=(0.948618356423)*sqrt(x[1])
        arg[0,2,1,0]=(-0.3835066594)*sqrt(x[0])
        arg[0,2,1,1]=(0.872875257641)*sqrt(x[1])
        arg[0,2,2,0]=(-0.98414676194)*sqrt(x[0])
        arg[0,2,2,1]=(0.59170993756)*sqrt(x[1])
        arg[0,3,0,0]=(0.875820385471)*sqrt(x[0])
        arg[0,3,0,1]=(-0.209139794169)*sqrt(x[1])
        arg[0,3,1,0]=(0.783312979835)*sqrt(x[0])
        arg[0,3,1,1]=(-0.440048708153)*sqrt(x[1])
        arg[0,3,2,0]=(-0.906741480158)*sqrt(x[0])
        arg[0,3,2,1]=(0.63834196559)*sqrt(x[1])
        arg[0,4,0,0]=(0.0142616166715)*sqrt(x[0])
        arg[0,4,0,1]=(-0.0745581409622)*sqrt(x[1])
        arg[0,4,1,0]=(-0.779966728505)*sqrt(x[0])
        arg[0,4,1,1]=(0.734744791584)*sqrt(x[1])
        arg[0,4,2,0]=(0.0407627699446)*sqrt(x[0])
        arg[0,4,2,1]=(-0.629666423098)*sqrt(x[1])
        arg[1,0,0,0]=(-0.782291768803)*sqrt(x[0])
        arg[1,0,0,1]=(-0.197799867335)*sqrt(x[1])
        arg[1,0,1,0]=(-0.487323212409)*sqrt(x[0])
        arg[1,0,1,1]=(-0.789940409633)*sqrt(x[1])
        arg[1,0,2,0]=(-0.523871448728)*sqrt(x[0])
        arg[1,0,2,1]=(-0.622272488658)*sqrt(x[1])
        arg[1,1,0,0]=(0.370095832941)*sqrt(x[0])
        arg[1,1,0,1]=(-0.0683624906905)*sqrt(x[1])
        arg[1,1,1,0]=(-0.180583695959)*sqrt(x[0])
        arg[1,1,1,1]=(-0.0529097308814)*sqrt(x[1])
        arg[1,1,2,0]=(0.440088045764)*sqrt(x[0])
        arg[1,1,2,1]=(-0.207665675854)*sqrt(x[1])
        arg[1,2,0,0]=(-0.360120735448)*sqrt(x[0])
        arg[1,2,0,1]=(0.877797020694)*sqrt(x[1])
        arg[1,2,1,0]=(0.655192721419)*sqrt(x[0])
        arg[1,2,1,1]=(0.978781227463)*sqrt(x[1])
        arg[1,2,2,0]=(-0.205886258554)*sqrt(x[0])
        arg[1,2,2,1]=(-0.0050396100363)*sqrt(x[1])
        arg[1,3,0,0]=(-0.484332943177)*sqrt(x[0])
        arg[1,3,0,1]=(0.964697361018)*sqrt(x[1])
        arg[1,3,1,0]=(0.0815652906824)*sqrt(x[0])
        arg[1,3,1,1]=(-0.531163025802)*sqrt(x[1])
        arg[1,3,2,0]=(-0.229615479833)*sqrt(x[0])
        arg[1,3,2,1]=(0.139971291609)*sqrt(x[1])
        arg[1,4,0,0]=(0.960839193726)*sqrt(x[0])
        arg[1,4,0,1]=(-0.9744930023)*sqrt(x[1])
        arg[1,4,1,0]=(-0.0919000363634)*sqrt(x[0])
        arg[1,4,1,1]=(0.483909186003)*sqrt(x[1])
        arg[1,4,2,0]=(0.248597552464)*sqrt(x[0])
        arg[1,4,2,1]=(0.965872402486)*sqrt(x[1])
        arg[2,0,0,0]=(-0.935049188033)*sqrt(x[0])
        arg[2,0,0,1]=(-0.823239964421)*sqrt(x[1])
        arg[2,0,1,0]=(0.30901193713)*sqrt(x[0])
        arg[2,0,1,1]=(-0.0538396976182)*sqrt(x[1])
        arg[2,0,2,0]=(-0.449416473099)*sqrt(x[0])
        arg[2,0,2,1]=(-0.857785267819)*sqrt(x[1])
        arg[2,1,0,0]=(-0.167802578334)*sqrt(x[0])
        arg[2,1,0,1]=(-0.7016101186)*sqrt(x[1])
        arg[2,1,1,0]=(0.995848317929)*sqrt(x[0])
        arg[2,1,1,1]=(-0.824595018011)*sqrt(x[1])
        arg[2,1,2,0]=(-0.826030527686)*sqrt(x[0])
        arg[2,1,2,1]=(-0.856151725651)*sqrt(x[1])
        arg[2,2,0,0]=(0.342047605028)*sqrt(x[0])
        arg[2,2,0,1]=(-0.342370307931)*sqrt(x[1])
        arg[2,2,1,0]=(-0.540693890831)*sqrt(x[0])
        arg[2,2,1,1]=(-0.720754897785)*sqrt(x[1])
        arg[2,2,2,0]=(-0.0941215154824)*sqrt(x[0])
        arg[2,2,2,1]=(0.115887914141)*sqrt(x[1])
        arg[2,3,0,0]=(-0.0920059639279)*sqrt(x[0])
        arg[2,3,0,1]=(0.910948619784)*sqrt(x[1])
        arg[2,3,1,0]=(-0.40846609126)*sqrt(x[0])
        arg[2,3,1,1]=(0.542728733954)*sqrt(x[1])
        arg[2,3,2,0]=(0.451361830368)*sqrt(x[0])
        arg[2,3,2,1]=(0.0767538617936)*sqrt(x[1])
        arg[2,4,0,0]=(-0.431757008766)*sqrt(x[0])
        arg[2,4,0,1]=(0.258988103513)*sqrt(x[1])
        arg[2,4,1,0]=(-0.584559859464)*sqrt(x[0])
        arg[2,4,1,1]=(-0.776827606804)*sqrt(x[1])
        arg[2,4,2,0]=(-0.867078960637)*sqrt(x[0])
        arg[2,4,2,1]=(-0.0972960135129)*sqrt(x[1])
        arg[3,0,0,0]=(0.999511727022)*sqrt(x[0])
        arg[3,0,0,1]=(-0.753218704622)*sqrt(x[1])
        arg[3,0,1,0]=(0.506439354885)*sqrt(x[0])
        arg[3,0,1,1]=(-0.0845853282297)*sqrt(x[1])
        arg[3,0,2,0]=(-0.871475028801)*sqrt(x[0])
        arg[3,0,2,1]=(-0.69261642657)*sqrt(x[1])
        arg[3,1,0,0]=(0.549778399933)*sqrt(x[0])
        arg[3,1,0,1]=(0.246187536149)*sqrt(x[1])
        arg[3,1,1,0]=(-0.620454200676)*sqrt(x[0])
        arg[3,1,1,1]=(0.372738564315)*sqrt(x[1])
        arg[3,1,2,0]=(0.0544097075138)*sqrt(x[0])
        arg[3,1,2,1]=(-0.0883198944676)*sqrt(x[1])
        arg[3,2,0,0]=(-0.671057180545)*sqrt(x[0])
        arg[3,2,0,1]=(-0.118991797907)*sqrt(x[1])
        arg[3,2,1,0]=(-0.27196730081)*sqrt(x[0])
        arg[3,2,1,1]=(-0.458807068968)*sqrt(x[1])
        arg[3,2,2,0]=(-0.31485399083)*sqrt(x[0])
        arg[3,2,2,1]=(0.487291960328)*sqrt(x[1])
        arg[3,3,0,0]=(-0.750302723531)*sqrt(x[0])
        arg[3,3,0,1]=(-0.486428789771)*sqrt(x[1])
        arg[3,3,1,0]=(-0.295909324594)*sqrt(x[0])
        arg[3,3,1,1]=(0.325702372596)*sqrt(x[1])
        arg[3,3,2,0]=(-0.512109540631)*sqrt(x[0])
        arg[3,3,2,1]=(0.563284721908)*sqrt(x[1])
        arg[3,4,0,0]=(0.53579406798)*sqrt(x[0])
        arg[3,4,0,1]=(-0.468431927975)*sqrt(x[1])
        arg[3,4,1,0]=(0.570851401329)*sqrt(x[0])
        arg[3,4,1,1]=(0.107791149865)*sqrt(x[1])
        arg[3,4,2,0]=(-0.211543670157)*sqrt(x[0])
        arg[3,4,2,1]=(0.848189550468)*sqrt(x[1])
        ref=sqrt((6.54316712884)*3./2.+(4.39564064326)/0.5)

      else:
        arg=Data(0,(4, 5, 3, 3),w)
        arg[0,0,0,0]=(-0.392256494872)*sqrt(x[0])
        arg[0,0,0,1]=(0.271488478841)*sqrt(x[1])
        arg[0,0,0,2]=(-0.878268531987)*sqrt(x[2])
        arg[0,0,1,0]=(0.106903787643)*sqrt(x[0])
        arg[0,0,1,1]=(-0.117986434516)*sqrt(x[1])
        arg[0,0,1,2]=(0.912673598522)*sqrt(x[2])
        arg[0,0,2,0]=(-0.652598945824)*sqrt(x[0])
        arg[0,0,2,1]=(-0.984195895647)*sqrt(x[1])
        arg[0,0,2,2]=(-0.137845459589)*sqrt(x[2])
        arg[0,1,0,0]=(0.696713200721)*sqrt(x[0])
        arg[0,1,0,1]=(0.480337929223)*sqrt(x[1])
        arg[0,1,0,2]=(-0.628726114325)*sqrt(x[2])
        arg[0,1,1,0]=(-0.356503461347)*sqrt(x[0])
        arg[0,1,1,1]=(-0.64670584045)*sqrt(x[1])
        arg[0,1,1,2]=(-0.737787618908)*sqrt(x[2])
        arg[0,1,2,0]=(-0.14207595302)*sqrt(x[0])
        arg[0,1,2,1]=(-0.573485525334)*sqrt(x[1])
        arg[0,1,2,2]=(-0.955542178481)*sqrt(x[2])
        arg[0,2,0,0]=(0.371952883975)*sqrt(x[0])
        arg[0,2,0,1]=(0.114016178762)*sqrt(x[1])
        arg[0,2,0,2]=(0.552721218169)*sqrt(x[2])
        arg[0,2,1,0]=(0.318426742113)*sqrt(x[0])
        arg[0,2,1,1]=(0.862220239384)*sqrt(x[1])
        arg[0,2,1,2]=(0.887376889978)*sqrt(x[2])
        arg[0,2,2,0]=(0.203656742981)*sqrt(x[0])
        arg[0,2,2,1]=(0.350541335144)*sqrt(x[1])
        arg[0,2,2,2]=(-0.448401957447)*sqrt(x[2])
        arg[0,3,0,0]=(-0.349200084042)*sqrt(x[0])
        arg[0,3,0,1]=(-0.546012602125)*sqrt(x[1])
        arg[0,3,0,2]=(-0.931488270979)*sqrt(x[2])
        arg[0,3,1,0]=(0.207457831058)*sqrt(x[0])
        arg[0,3,1,1]=(0.557715840549)*sqrt(x[1])
        arg[0,3,1,2]=(-0.978130146744)*sqrt(x[2])
        arg[0,3,2,0]=(-0.55088967957)*sqrt(x[0])
        arg[0,3,2,1]=(0.0490328838057)*sqrt(x[1])
        arg[0,3,2,2]=(0.150209929122)*sqrt(x[2])
        arg[0,4,0,0]=(-0.484145614698)*sqrt(x[0])
        arg[0,4,0,1]=(0.393078411279)*sqrt(x[1])
        arg[0,4,0,2]=(0.0678631863917)*sqrt(x[2])
        arg[0,4,1,0]=(-0.350980464628)*sqrt(x[0])
        arg[0,4,1,1]=(-0.784203839564)*sqrt(x[1])
        arg[0,4,1,2]=(0.636960296147)*sqrt(x[2])
        arg[0,4,2,0]=(0.592799503581)*sqrt(x[0])
        arg[0,4,2,1]=(-0.672104833683)*sqrt(x[1])
        arg[0,4,2,2]=(0.0366914082467)*sqrt(x[2])
        arg[1,0,0,0]=(-0.147168019774)*sqrt(x[0])
        arg[1,0,0,1]=(-0.0823637938956)*sqrt(x[1])
        arg[1,0,0,2]=(-0.852729690176)*sqrt(x[2])
        arg[1,0,1,0]=(-0.612338430408)*sqrt(x[0])
        arg[1,0,1,1]=(-0.85820035747)*sqrt(x[1])
        arg[1,0,1,2]=(-0.463664966162)*sqrt(x[2])
        arg[1,0,2,0]=(0.274600720491)*sqrt(x[0])
        arg[1,0,2,1]=(-0.488508234093)*sqrt(x[1])
        arg[1,0,2,2]=(-0.28251466519)*sqrt(x[2])
        arg[1,1,0,0]=(-0.0196532462794)*sqrt(x[0])
        arg[1,1,0,1]=(0.239352528871)*sqrt(x[1])
        arg[1,1,0,2]=(-0.17054773873)*sqrt(x[2])
        arg[1,1,1,0]=(-0.2396627789)*sqrt(x[0])
        arg[1,1,1,1]=(0.868970323003)*sqrt(x[1])
        arg[1,1,1,2]=(0.401487430312)*sqrt(x[2])
        arg[1,1,2,0]=(-0.624411449783)*sqrt(x[0])
        arg[1,1,2,1]=(0.0036634266684)*sqrt(x[1])
        arg[1,1,2,2]=(0.736129120967)*sqrt(x[2])
        arg[1,2,0,0]=(0.183756511707)*sqrt(x[0])
        arg[1,2,0,1]=(-0.288651848639)*sqrt(x[1])
        arg[1,2,0,2]=(-0.0672121537447)*sqrt(x[2])
        arg[1,2,1,0]=(-0.323274725936)*sqrt(x[0])
        arg[1,2,1,1]=(0.298001016025)*sqrt(x[1])
        arg[1,2,1,2]=(-0.976052460675)*sqrt(x[2])
        arg[1,2,2,0]=(0.596504441096)*sqrt(x[0])
        arg[1,2,2,1]=(0.873776068983)*sqrt(x[1])
        arg[1,2,2,2]=(-0.994068273196)*sqrt(x[2])
        arg[1,3,0,0]=(-0.495387299681)*sqrt(x[0])
        arg[1,3,0,1]=(-0.123674756551)*sqrt(x[1])
        arg[1,3,0,2]=(0.581213818577)*sqrt(x[2])
        arg[1,3,1,0]=(0.146405749701)*sqrt(x[0])
        arg[1,3,1,1]=(-0.594994686675)*sqrt(x[1])
        arg[1,3,1,2]=(-0.059093568436)*sqrt(x[2])
        arg[1,3,2,0]=(0.651004255104)*sqrt(x[0])
        arg[1,3,2,1]=(-0.977880706193)*sqrt(x[1])
        arg[1,3,2,2]=(0.370344651319)*sqrt(x[2])
        arg[1,4,0,0]=(-0.503657215247)*sqrt(x[0])
        arg[1,4,0,1]=(-0.170885297253)*sqrt(x[1])
        arg[1,4,0,2]=(0.533424480956)*sqrt(x[2])
        arg[1,4,1,0]=(-0.220533193308)*sqrt(x[0])
        arg[1,4,1,1]=(0.344537611882)*sqrt(x[1])
        arg[1,4,1,2]=(0.861473877282)*sqrt(x[2])
        arg[1,4,2,0]=(-0.0923010438884)*sqrt(x[0])
        arg[1,4,2,1]=(-0.338256780498)*sqrt(x[1])
        arg[1,4,2,2]=(0.528567959345)*sqrt(x[2])
        arg[2,0,0,0]=(-0.0423053381485)*sqrt(x[0])
        arg[2,0,0,1]=(0.856798579151)*sqrt(x[1])
        arg[2,0,0,2]=(0.383258153853)*sqrt(x[2])
        arg[2,0,1,0]=(0.350994872736)*sqrt(x[0])
        arg[2,0,1,1]=(-0.78055158106)*sqrt(x[1])
        arg[2,0,1,2]=(-0.770876699915)*sqrt(x[2])
        arg[2,0,2,0]=(-0.935133287106)*sqrt(x[0])
        arg[2,0,2,1]=(0.618238076989)*sqrt(x[1])
        arg[2,0,2,2]=(-0.846783087949)*sqrt(x[2])
        arg[2,1,0,0]=(-0.657394511405)*sqrt(x[0])
        arg[2,1,0,1]=(0.576218821654)*sqrt(x[1])
        arg[2,1,0,2]=(0.0269446356493)*sqrt(x[2])
        arg[2,1,1,0]=(-0.310710230949)*sqrt(x[0])
        arg[2,1,1,1]=(0.425412515598)*sqrt(x[1])
        arg[2,1,1,2]=(-0.225225290862)*sqrt(x[2])
        arg[2,1,2,0]=(-0.539928589495)*sqrt(x[0])
        arg[2,1,2,1]=(-0.348082121765)*sqrt(x[1])
        arg[2,1,2,2]=(-0.0287274646233)*sqrt(x[2])
        arg[2,2,0,0]=(0.746132865923)*sqrt(x[0])
        arg[2,2,0,1]=(-0.0234203693548)*sqrt(x[1])
        arg[2,2,0,2]=(0.517411821941)*sqrt(x[2])
        arg[2,2,1,0]=(-0.183204217349)*sqrt(x[0])
        arg[2,2,1,1]=(0.714988861836)*sqrt(x[1])
        arg[2,2,1,2]=(0.829083318937)*sqrt(x[2])
        arg[2,2,2,0]=(0.458067555841)*sqrt(x[0])
        arg[2,2,2,1]=(0.639317125869)*sqrt(x[1])
        arg[2,2,2,2]=(0.104611520408)*sqrt(x[2])
        arg[2,3,0,0]=(-0.420513461135)*sqrt(x[0])
        arg[2,3,0,1]=(0.888686162754)*sqrt(x[1])
        arg[2,3,0,2]=(0.939305777879)*sqrt(x[2])
        arg[2,3,1,0]=(0.856795132015)*sqrt(x[0])
        arg[2,3,1,1]=(0.817593141895)*sqrt(x[1])
        arg[2,3,1,2]=(0.962503342535)*sqrt(x[2])
        arg[2,3,2,0]=(-0.334586700245)*sqrt(x[0])
        arg[2,3,2,1]=(0.182696129528)*sqrt(x[1])
        arg[2,3,2,2]=(-0.707271571206)*sqrt(x[2])
        arg[2,4,0,0]=(-0.594795981069)*sqrt(x[0])
        arg[2,4,0,1]=(0.451239168073)*sqrt(x[1])
        arg[2,4,0,2]=(0.191212211556)*sqrt(x[2])
        arg[2,4,1,0]=(-0.503465984944)*sqrt(x[0])
        arg[2,4,1,1]=(0.725377884208)*sqrt(x[1])
        arg[2,4,1,2]=(-0.40719255752)*sqrt(x[2])
        arg[2,4,2,0]=(0.268071476451)*sqrt(x[0])
        arg[2,4,2,1]=(0.85066639942)*sqrt(x[1])
        arg[2,4,2,2]=(-0.906021406945)*sqrt(x[2])
        arg[3,0,0,0]=(-0.922179152122)*sqrt(x[0])
        arg[3,0,0,1]=(-0.0903841240007)*sqrt(x[1])
        arg[3,0,0,2]=(-0.751482803516)*sqrt(x[2])
        arg[3,0,1,0]=(0.960697809119)*sqrt(x[0])
        arg[3,0,1,1]=(0.638878873158)*sqrt(x[1])
        arg[3,0,1,2]=(0.390932234724)*sqrt(x[2])
        arg[3,0,2,0]=(-0.925078301694)*sqrt(x[0])
        arg[3,0,2,1]=(0.793590580665)*sqrt(x[1])
        arg[3,0,2,2]=(0.535478366911)*sqrt(x[2])
        arg[3,1,0,0]=(-0.431951993217)*sqrt(x[0])
        arg[3,1,0,1]=(0.211750261417)*sqrt(x[1])
        arg[3,1,0,2]=(-0.930706580442)*sqrt(x[2])
        arg[3,1,1,0]=(0.330979313323)*sqrt(x[0])
        arg[3,1,1,1]=(-0.838919076081)*sqrt(x[1])
        arg[3,1,1,2]=(0.134250050168)*sqrt(x[2])
        arg[3,1,2,0]=(0.414922811301)*sqrt(x[0])
        arg[3,1,2,1]=(-0.663692878121)*sqrt(x[1])
        arg[3,1,2,2]=(0.88499278543)*sqrt(x[2])
        arg[3,2,0,0]=(0.0884662742233)*sqrt(x[0])
        arg[3,2,0,1]=(-0.412630722821)*sqrt(x[1])
        arg[3,2,0,2]=(-0.730850884928)*sqrt(x[2])
        arg[3,2,1,0]=(0.722207547366)*sqrt(x[0])
        arg[3,2,1,1]=(-0.260067950749)*sqrt(x[1])
        arg[3,2,1,2]=(0.426259201494)*sqrt(x[2])
        arg[3,2,2,0]=(0.0516111795322)*sqrt(x[0])
        arg[3,2,2,1]=(0.922853710048)*sqrt(x[1])
        arg[3,2,2,2]=(-0.991912758116)*sqrt(x[2])
        arg[3,3,0,0]=(0.263933965905)*sqrt(x[0])
        arg[3,3,0,1]=(0.840541758799)*sqrt(x[1])
        arg[3,3,0,2]=(0.417658511125)*sqrt(x[2])
        arg[3,3,1,0]=(-0.901745614723)*sqrt(x[0])
        arg[3,3,1,1]=(-0.623608908699)*sqrt(x[1])
        arg[3,3,1,2]=(0.0522167208784)*sqrt(x[2])
        arg[3,3,2,0]=(-0.549431264931)*sqrt(x[0])
        arg[3,3,2,1]=(0.919971855457)*sqrt(x[1])
        arg[3,3,2,2]=(0.142757773397)*sqrt(x[2])
        arg[3,4,0,0]=(-0.258875259824)*sqrt(x[0])
        arg[3,4,0,1]=(-0.0373872187041)*sqrt(x[1])
        arg[3,4,0,2]=(0.445989164864)*sqrt(x[2])
        arg[3,4,1,0]=(0.0273711397038)*sqrt(x[0])
        arg[3,4,1,1]=(-0.85522629)*sqrt(x[1])
        arg[3,4,1,2]=(-0.835392581226)*sqrt(x[2])
        arg[3,4,2,0]=(0.402929450189)*sqrt(x[0])
        arg[3,4,2,1]=(0.549623033221)*sqrt(x[1])
        arg[3,4,2,2]=(0.237803212109)*sqrt(x[2])
        ref=sqrt((15.1448100099)*3./2.+(3.61598249819)/0.5)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunctionOnContactOne_fromData_rank0(self):
      """
      tests L2-norm of Data on the FunctionOnContactOne

      assumptions: self.domain supports integration on FunctionOnContactOne
      """
      dim=self.domain.getDim()
      w=FunctionOnContactOne(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(),w)
        arg=(0.39039157828)*x[0]
        ref=sqrt((0.0)+(0.038101396098))

      else:
        arg=Data(0,(),w)
        arg=(0.893204658234)*x[0]
        ref=sqrt((0.0)+(0.199453640373))

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunctionOnContactOne_fromData_rank1(self):
      """
      tests L2-norm of Data on the FunctionOnContactOne

      assumptions: self.domain supports integration on FunctionOnContactOne
      """
      dim=self.domain.getDim()
      w=FunctionOnContactOne(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(2,),w)
        arg[0]=(0.697576620054)*x[0]
        arg[1]=(-0.302037273777)*x[1]
        ref=sqrt((0.0304088382502)+(0.121653285212))

      else:
        arg=Data(0,(3,),w)
        arg[0]=(0.649308272675)*x[0]
        arg[1]=(-0.408829699836)*x[1]
        arg[2]=(-0.984892977949)*x[2]
        ref=sqrt((0.37905196716)+(0.105400308241))

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunctionOnContactOne_fromData_rank2(self):
      """
      tests L2-norm of Data on the FunctionOnContactOne

      assumptions: self.domain supports integration on FunctionOnContactOne
      """
      dim=self.domain.getDim()
      w=FunctionOnContactOne(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 2),w)
        arg[0,0]=(0.504133645896)*x[0]
        arg[0,1]=(0.728228500775)*x[1]
        arg[1,0]=(0.968237164633)*x[0]
        arg[1,1]=(-0.205351628957)*x[1]
        arg[2,0]=(-0.154008228702)*x[0]
        arg[2,1]=(-0.972088404566)*x[1]
        arg[3,0]=(0.585153521593)*x[0]
        arg[3,1]=(0.371861820001)*x[1]
        ref=sqrt((0.551907706774)+(0.389439279561))

      else:
        arg=Data(0,(4, 3),w)
        arg[0,0]=(-0.4847970605)*x[0]
        arg[0,1]=(-0.664433021479)*x[1]
        arg[0,2]=(-0.93730914676)*x[2]
        arg[1,0]=(0.712624065023)*x[0]
        arg[1,1]=(0.105431976424)*x[1]
        arg[1,2]=(-0.452636806104)*x[2]
        arg[2,0]=(0.488406064176)*x[0]
        arg[2,1]=(0.0594879976584)*x[1]
        arg[2,2]=(0.596963883202)*x[2]
        arg[3,0]=(0.619366395328)*x[0]
        arg[3,1]=(0.719965901041)*x[1]
        arg[3,2]=(-0.568971277169)*x[2]
        ref=sqrt((0.912666523048)+(0.341254115776))

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunctionOnContactOne_fromData_rank3(self):
      """
      tests L2-norm of Data on the FunctionOnContactOne

      assumptions: self.domain supports integration on FunctionOnContactOne
      """
      dim=self.domain.getDim()
      w=FunctionOnContactOne(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(6, 2, 2),w)
        arg[0,0,0]=(-0.481233278183)*x[0]
        arg[0,0,1]=(-0.883447576048)*x[1]
        arg[0,1,0]=(0.525938652633)*x[0]
        arg[0,1,1]=(-0.694971850651)*x[1]
        arg[1,0,0]=(-0.935880633291)*x[0]
        arg[1,0,1]=(0.970445286468)*x[1]
        arg[1,1,0]=(-0.302571977945)*x[0]
        arg[1,1,1]=(-0.237750319821)*x[1]
        arg[2,0,0]=(0.922940652883)*x[0]
        arg[2,0,1]=(-0.624756614078)*x[1]
        arg[2,1,0]=(-0.146839388764)*x[0]
        arg[2,1,1]=(-0.995551211092)*x[1]
        arg[3,0,0]=(0.719184541649)*x[0]
        arg[3,0,1]=(-0.629461063375)*x[1]
        arg[3,1,0]=(-0.925581215276)*x[0]
        arg[3,1,1]=(-0.564418658495)*x[1]
        arg[4,0,0]=(-0.724035884026)*x[0]
        arg[4,0,1]=(0.78393984221)*x[1]
        arg[4,1,0]=(-0.722375249317)*x[0]
        arg[4,1,1]=(0.27571990944)*x[1]
        arg[5,0,0]=(-0.197579392553)*x[0]
        arg[5,0,1]=(0.341717398275)*x[1]
        arg[5,1,0]=(0.926327822859)*x[0]
        arg[5,1,1]=(0.44114488285)*x[1]
        ref=sqrt((1.78665006238)+(1.41652558894))

      else:
        arg=Data(0,(6, 2, 3),w)
        arg[0,0,0]=(-0.348448140437)*x[0]
        arg[0,0,1]=(-0.530793589782)*x[1]
        arg[0,0,2]=(0.386511493376)*x[2]
        arg[0,1,0]=(0.00776536299593)*x[0]
        arg[0,1,1]=(-0.0696129074064)*x[1]
        arg[0,1,2]=(0.114179504216)*x[2]
        arg[1,0,0]=(-0.681009917406)*x[0]
        arg[1,0,1]=(0.763595403719)*x[1]
        arg[1,0,2]=(-0.54459342544)*x[2]
        arg[1,1,0]=(0.0839255723157)*x[0]
        arg[1,1,1]=(0.269389515003)*x[1]
        arg[1,1,2]=(0.612446509453)*x[2]
        arg[2,0,0]=(0.600070579238)*x[0]
        arg[2,0,1]=(-0.49585867531)*x[1]
        arg[2,0,2]=(-0.766857395141)*x[2]
        arg[2,1,0]=(-0.84515816441)*x[0]
        arg[2,1,1]=(-0.594337789049)*x[1]
        arg[2,1,2]=(0.774596009811)*x[2]
        arg[3,0,0]=(-0.833982234688)*x[0]
        arg[3,0,1]=(0.747667281218)*x[1]
        arg[3,0,2]=(-0.77025764051)*x[2]
        arg[3,1,0]=(0.0367023286497)*x[0]
        arg[3,1,1]=(-0.503080137993)*x[1]
        arg[3,1,2]=(-0.129478323797)*x[2]
        arg[4,0,0]=(0.304123530695)*x[0]
        arg[4,0,1]=(-0.873023368289)*x[1]
        arg[4,0,2]=(-0.712740806672)*x[2]
        arg[4,1,0]=(-0.0629132979843)*x[0]
        arg[4,1,1]=(0.096843438525)*x[1]
        arg[4,1,2]=(-0.608716559621)*x[2]
        arg[5,0,0]=(-0.642902822517)*x[0]
        arg[5,0,1]=(0.478134785533)*x[1]
        arg[5,0,2]=(0.757891512092)*x[2]
        arg[5,1,0]=(-0.144492755065)*x[0]
        arg[5,1,1]=(0.286775195938)*x[1]
        arg[5,1,2]=(-0.063580488638)*x[2]
        ref=sqrt((2.50835195543)+(0.723549067206))

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunctionOnContactOne_fromData_rank4(self):
      """
      tests L2-norm of Data on the FunctionOnContactOne

      assumptions: self.domain supports integration on FunctionOnContactOne
      """
      dim=self.domain.getDim()
      w=FunctionOnContactOne(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 5, 3, 2),w)
        arg[0,0,0,0]=(0.523372259748)*x[0]
        arg[0,0,0,1]=(-0.732795389957)*x[1]
        arg[0,0,1,0]=(0.0571737128606)*x[0]
        arg[0,0,1,1]=(-0.627260438721)*x[1]
        arg[0,0,2,0]=(-0.945159375045)*x[0]
        arg[0,0,2,1]=(-0.145103328237)*x[1]
        arg[0,1,0,0]=(-0.0586985419455)*x[0]
        arg[0,1,0,1]=(-0.788101334389)*x[1]
        arg[0,1,1,0]=(0.615023802803)*x[0]
        arg[0,1,1,1]=(0.965381283357)*x[1]
        arg[0,1,2,0]=(0.173979432067)*x[0]
        arg[0,1,2,1]=(-0.980483156424)*x[1]
        arg[0,2,0,0]=(-0.345940909471)*x[0]
        arg[0,2,0,1]=(0.00162430279595)*x[1]
        arg[0,2,1,0]=(-0.0847413350082)*x[0]
        arg[0,2,1,1]=(-0.911248036753)*x[1]
        arg[0,2,2,0]=(-0.817295589942)*x[0]
        arg[0,2,2,1]=(-0.0514569401185)*x[1]
        arg[0,3,0,0]=(0.329704440159)*x[0]
        arg[0,3,0,1]=(0.580379057442)*x[1]
        arg[0,3,1,0]=(0.737014729605)*x[0]
        arg[0,3,1,1]=(0.535521740242)*x[1]
        arg[0,3,2,0]=(-0.943554721405)*x[0]
        arg[0,3,2,1]=(0.798294295959)*x[1]
        arg[0,4,0,0]=(-0.1814410847)*x[0]
        arg[0,4,0,1]=(0.968476629577)*x[1]
        arg[0,4,1,0]=(0.26595250332)*x[0]
        arg[0,4,1,1]=(-0.493516007909)*x[1]
        arg[0,4,2,0]=(0.374323739856)*x[0]
        arg[0,4,2,1]=(-0.914413601829)*x[1]
        arg[1,0,0,0]=(-0.701090812201)*x[0]
        arg[1,0,0,1]=(-0.403596866322)*x[1]
        arg[1,0,1,0]=(0.587183778075)*x[0]
        arg[1,0,1,1]=(-0.218129059558)*x[1]
        arg[1,0,2,0]=(-0.0205642762349)*x[0]
        arg[1,0,2,1]=(0.441382235525)*x[1]
        arg[1,1,0,0]=(0.909683859519)*x[0]
        arg[1,1,0,1]=(0.960088832754)*x[1]
        arg[1,1,1,0]=(0.0450355911985)*x[0]
        arg[1,1,1,1]=(0.366849576658)*x[1]
        arg[1,1,2,0]=(0.315490164497)*x[0]
        arg[1,1,2,1]=(0.533746946808)*x[1]
        arg[1,2,0,0]=(-0.552919155957)*x[0]
        arg[1,2,0,1]=(-0.647373009727)*x[1]
        arg[1,2,1,0]=(0.32617079311)*x[0]
        arg[1,2,1,1]=(-0.3424263412)*x[1]
        arg[1,2,2,0]=(0.848733118893)*x[0]
        arg[1,2,2,1]=(-0.65341712364)*x[1]
        arg[1,3,0,0]=(-0.698538235208)*x[0]
        arg[1,3,0,1]=(-0.176858090319)*x[1]
        arg[1,3,1,0]=(-0.956962782849)*x[0]
        arg[1,3,1,1]=(0.53041598553)*x[1]
        arg[1,3,2,0]=(-0.241715236701)*x[0]
        arg[1,3,2,1]=(0.0340608544103)*x[1]
        arg[1,4,0,0]=(-0.399279751219)*x[0]
        arg[1,4,0,1]=(0.383860128496)*x[1]
        arg[1,4,1,0]=(0.229189108307)*x[0]
        arg[1,4,1,1]=(0.963667605308)*x[1]
        arg[1,4,2,0]=(0.731912848998)*x[0]
        arg[1,4,2,1]=(0.0534970206048)*x[1]
        arg[2,0,0,0]=(0.422292915614)*x[0]
        arg[2,0,0,1]=(0.915607515426)*x[1]
        arg[2,0,1,0]=(0.451054723135)*x[0]
        arg[2,0,1,1]=(0.0400870534188)*x[1]
        arg[2,0,2,0]=(-0.851047453002)*x[0]
        arg[2,0,2,1]=(0.300957420969)*x[1]
        arg[2,1,0,0]=(0.636801515211)*x[0]
        arg[2,1,0,1]=(0.386798062533)*x[1]
        arg[2,1,1,0]=(-0.0501653606261)*x[0]
        arg[2,1,1,1]=(0.46515271421)*x[1]
        arg[2,1,2,0]=(-0.416436918527)*x[0]
        arg[2,1,2,1]=(0.214018855333)*x[1]
        arg[2,2,0,0]=(0.0516439357756)*x[0]
        arg[2,2,0,1]=(-0.237707451356)*x[1]
        arg[2,2,1,0]=(0.117559050121)*x[0]
        arg[2,2,1,1]=(0.164349348553)*x[1]
        arg[2,2,2,0]=(0.0524038362463)*x[0]
        arg[2,2,2,1]=(0.690751059081)*x[1]
        arg[2,3,0,0]=(-0.942079336547)*x[0]
        arg[2,3,0,1]=(0.340143899875)*x[1]
        arg[2,3,1,0]=(-0.630411753789)*x[0]
        arg[2,3,1,1]=(-0.722078100418)*x[1]
        arg[2,3,2,0]=(-0.310606406168)*x[0]
        arg[2,3,2,1]=(0.728069014178)*x[1]
        arg[2,4,0,0]=(0.981678768954)*x[0]
        arg[2,4,0,1]=(-0.273023706956)*x[1]
        arg[2,4,1,0]=(0.307187251211)*x[0]
        arg[2,4,1,1]=(-0.885309661176)*x[1]
        arg[2,4,2,0]=(-0.808022645995)*x[0]
        arg[2,4,2,1]=(-0.21488080753)*x[1]
        arg[3,0,0,0]=(-0.785416067695)*x[0]
        arg[3,0,0,1]=(0.344457091999)*x[1]
        arg[3,0,1,0]=(-0.845835855973)*x[0]
        arg[3,0,1,1]=(0.365474280353)*x[1]
        arg[3,0,2,0]=(-0.250535373488)*x[0]
        arg[3,0,2,1]=(-0.888004656815)*x[1]
        arg[3,1,0,0]=(-0.0652508248717)*x[0]
        arg[3,1,0,1]=(0.402069695186)*x[1]
        arg[3,1,1,0]=(-0.229714763188)*x[0]
        arg[3,1,1,1]=(-0.866938136156)*x[1]
        arg[3,1,2,0]=(-0.638603745536)*x[0]
        arg[3,1,2,1]=(-0.997108734664)*x[1]
        arg[3,2,0,0]=(-0.182784593825)*x[0]
        arg[3,2,0,1]=(0.426951207538)*x[1]
        arg[3,2,1,0]=(-0.958149736654)*x[0]
        arg[3,2,1,1]=(-0.58458370587)*x[1]
        arg[3,2,2,0]=(0.410243417149)*x[0]
        arg[3,2,2,1]=(0.764474507608)*x[1]
        arg[3,3,0,0]=(0.235623988334)*x[0]
        arg[3,3,0,1]=(-0.105137401903)*x[1]
        arg[3,3,1,0]=(-0.957743319072)*x[0]
        arg[3,3,1,1]=(0.285680328399)*x[1]
        arg[3,3,2,0]=(-0.525269717395)*x[0]
        arg[3,3,2,1]=(0.655974370871)*x[1]
        arg[3,4,0,0]=(0.926779313365)*x[0]
        arg[3,4,0,1]=(-0.111556488031)*x[1]
        arg[3,4,1,0]=(-0.672785588938)*x[0]
        arg[3,4,1,1]=(-0.687907843899)*x[1]
        arg[3,4,2,0]=(0.503984478085)*x[0]
        arg[3,4,2,1]=(0.213194091516)*x[1]
        ref=sqrt((6.92178561091)+(4.96610502195))

      else:
        arg=Data(0,(4, 5, 3, 3),w)
        arg[0,0,0,0]=(-0.641108959812)*x[0]
        arg[0,0,0,1]=(0.395837684631)*x[1]
        arg[0,0,0,2]=(0.36258240616)*x[2]
        arg[0,0,1,0]=(-0.0543293722691)*x[0]
        arg[0,0,1,1]=(0.850820297369)*x[1]
        arg[0,0,1,2]=(0.0910730209593)*x[2]
        arg[0,0,2,0]=(-0.831108730829)*x[0]
        arg[0,0,2,1]=(-0.653018722047)*x[1]
        arg[0,0,2,2]=(0.29052921548)*x[2]
        arg[0,1,0,0]=(-0.115846569271)*x[0]
        arg[0,1,0,1]=(-0.828080308748)*x[1]
        arg[0,1,0,2]=(0.715255255273)*x[2]
        arg[0,1,1,0]=(-0.0106077934044)*x[0]
        arg[0,1,1,1]=(0.358323295554)*x[1]
        arg[0,1,1,2]=(-0.916496266702)*x[2]
        arg[0,1,2,0]=(0.740969432327)*x[0]
        arg[0,1,2,1]=(0.0631507192754)*x[1]
        arg[0,1,2,2]=(0.440009320758)*x[2]
        arg[0,2,0,0]=(0.689964376385)*x[0]
        arg[0,2,0,1]=(0.58199412545)*x[1]
        arg[0,2,0,2]=(0.649058556936)*x[2]
        arg[0,2,1,0]=(0.811660544475)*x[0]
        arg[0,2,1,1]=(-0.121902248927)*x[1]
        arg[0,2,1,2]=(-0.0750232179247)*x[2]
        arg[0,2,2,0]=(0.779166450882)*x[0]
        arg[0,2,2,1]=(0.298560101475)*x[1]
        arg[0,2,2,2]=(-0.776198356206)*x[2]
        arg[0,3,0,0]=(0.79804128408)*x[0]
        arg[0,3,0,1]=(-0.378223308322)*x[1]
        arg[0,3,0,2]=(0.191534109359)*x[2]
        arg[0,3,1,0]=(-0.59973117885)*x[0]
        arg[0,3,1,1]=(0.149042766249)*x[1]
        arg[0,3,1,2]=(-0.357473233062)*x[2]
        arg[0,3,2,0]=(-0.343758339283)*x[0]
        arg[0,3,2,1]=(-0.140764743713)*x[1]
        arg[0,3,2,2]=(-0.506377878901)*x[2]
        arg[0,4,0,0]=(-0.967765607304)*x[0]
        arg[0,4,0,1]=(-0.121373931498)*x[1]
        arg[0,4,0,2]=(-0.300794608584)*x[2]
        arg[0,4,1,0]=(-0.604537150389)*x[0]
        arg[0,4,1,1]=(-0.290834102056)*x[1]
        arg[0,4,1,2]=(0.509331591828)*x[2]
        arg[0,4,2,0]=(0.643653985564)*x[0]
        arg[0,4,2,1]=(-0.599746831159)*x[1]
        arg[0,4,2,2]=(0.843193710828)*x[2]
        arg[1,0,0,0]=(-0.598249432719)*x[0]
        arg[1,0,0,1]=(0.840724468488)*x[1]
        arg[1,0,0,2]=(0.558893493847)*x[2]
        arg[1,0,1,0]=(0.843245438555)*x[0]
        arg[1,0,1,1]=(0.280919545563)*x[1]
        arg[1,0,1,2]=(-0.0549504830615)*x[2]
        arg[1,0,2,0]=(0.295812998973)*x[0]
        arg[1,0,2,1]=(0.513221249344)*x[1]
        arg[1,0,2,2]=(-0.432253600492)*x[2]
        arg[1,1,0,0]=(0.300696333449)*x[0]
        arg[1,1,0,1]=(0.964202717779)*x[1]
        arg[1,1,0,2]=(0.131730374204)*x[2]
        arg[1,1,1,0]=(0.209135189207)*x[0]
        arg[1,1,1,1]=(0.301389564979)*x[1]
        arg[1,1,1,2]=(0.378481289816)*x[2]
        arg[1,1,2,0]=(0.407413524533)*x[0]
        arg[1,1,2,1]=(0.184345610049)*x[1]
        arg[1,1,2,2]=(-0.884401323556)*x[2]
        arg[1,2,0,0]=(0.886883417187)*x[0]
        arg[1,2,0,1]=(-0.0524504232352)*x[1]
        arg[1,2,0,2]=(0.049116549681)*x[2]
        arg[1,2,1,0]=(-0.41240815391)*x[0]
        arg[1,2,1,1]=(-0.523033558623)*x[1]
        arg[1,2,1,2]=(-0.920238404199)*x[2]
        arg[1,2,2,0]=(-0.635169424461)*x[0]
        arg[1,2,2,1]=(0.773163888164)*x[1]
        arg[1,2,2,2]=(0.554400056145)*x[2]
        arg[1,3,0,0]=(0.573361031392)*x[0]
        arg[1,3,0,1]=(-0.107749750058)*x[1]
        arg[1,3,0,2]=(-0.971799077244)*x[2]
        arg[1,3,1,0]=(-0.103579714499)*x[0]
        arg[1,3,1,1]=(0.93329275867)*x[1]
        arg[1,3,1,2]=(0.994409555402)*x[2]
        arg[1,3,2,0]=(0.571961012279)*x[0]
        arg[1,3,2,1]=(-0.54961086377)*x[1]
        arg[1,3,2,2]=(0.563815465418)*x[2]
        arg[1,4,0,0]=(-0.820638948871)*x[0]
        arg[1,4,0,1]=(0.407750087109)*x[1]
        arg[1,4,0,2]=(0.0911636159786)*x[2]
        arg[1,4,1,0]=(0.371133349901)*x[0]
        arg[1,4,1,1]=(0.688167208194)*x[1]
        arg[1,4,1,2]=(0.188789217725)*x[2]
        arg[1,4,2,0]=(0.961439941063)*x[0]
        arg[1,4,2,1]=(0.524326385264)*x[1]
        arg[1,4,2,2]=(-0.0866064394717)*x[2]
        arg[2,0,0,0]=(0.392756517944)*x[0]
        arg[2,0,0,1]=(0.713723092455)*x[1]
        arg[2,0,0,2]=(0.68441988583)*x[2]
        arg[2,0,1,0]=(-0.471505369012)*x[0]
        arg[2,0,1,1]=(0.152712812211)*x[1]
        arg[2,0,1,2]=(-0.579106568625)*x[2]
        arg[2,0,2,0]=(0.637759469489)*x[0]
        arg[2,0,2,1]=(-0.986405270597)*x[1]
        arg[2,0,2,2]=(0.498109007475)*x[2]
        arg[2,1,0,0]=(0.870185247694)*x[0]
        arg[2,1,0,1]=(0.738095675991)*x[1]
        arg[2,1,0,2]=(-0.613836737976)*x[2]
        arg[2,1,1,0]=(-0.507009374688)*x[0]
        arg[2,1,1,1]=(0.00752595652646)*x[1]
        arg[2,1,1,2]=(-0.389118104751)*x[2]
        arg[2,1,2,0]=(-0.466733433552)*x[0]
        arg[2,1,2,1]=(-0.356829687828)*x[1]
        arg[2,1,2,2]=(0.691337110432)*x[2]
        arg[2,2,0,0]=(-0.684116912917)*x[0]
        arg[2,2,0,1]=(0.705774798227)*x[1]
        arg[2,2,0,2]=(0.522279733527)*x[2]
        arg[2,2,1,0]=(0.390315182318)*x[0]
        arg[2,2,1,1]=(-0.450165044241)*x[1]
        arg[2,2,1,2]=(0.436062806791)*x[2]
        arg[2,2,2,0]=(0.973600590724)*x[0]
        arg[2,2,2,1]=(-0.89116205465)*x[1]
        arg[2,2,2,2]=(-0.575262591791)*x[2]
        arg[2,3,0,0]=(0.842833949385)*x[0]
        arg[2,3,0,1]=(-0.69939817904)*x[1]
        arg[2,3,0,2]=(0.343679860971)*x[2]
        arg[2,3,1,0]=(-0.387030334358)*x[0]
        arg[2,3,1,1]=(0.901637622406)*x[1]
        arg[2,3,1,2]=(0.054184099786)*x[2]
        arg[2,3,2,0]=(0.652712575782)*x[0]
        arg[2,3,2,1]=(0.286474075274)*x[1]
        arg[2,3,2,2]=(-0.745166160145)*x[2]
        arg[2,4,0,0]=(-0.104591304669)*x[0]
        arg[2,4,0,1]=(0.417496371683)*x[1]
        arg[2,4,0,2]=(0.0579341684769)*x[2]
        arg[2,4,1,0]=(0.24486914915)*x[0]
        arg[2,4,1,1]=(-0.166211416717)*x[1]
        arg[2,4,1,2]=(0.651515760357)*x[2]
        arg[2,4,2,0]=(0.934693067557)*x[0]
        arg[2,4,2,1]=(0.628486602755)*x[1]
        arg[2,4,2,2]=(0.855313741494)*x[2]
        arg[3,0,0,0]=(0.2241978058)*x[0]
        arg[3,0,0,1]=(0.389284487657)*x[1]
        arg[3,0,0,2]=(0.711521372962)*x[2]
        arg[3,0,1,0]=(0.653184309521)*x[0]
        arg[3,0,1,1]=(0.270052840275)*x[1]
        arg[3,0,1,2]=(0.388881294549)*x[2]
        arg[3,0,2,0]=(0.973583452692)*x[0]
        arg[3,0,2,1]=(0.508175099792)*x[1]
        arg[3,0,2,2]=(-0.989999307019)*x[2]
        arg[3,1,0,0]=(0.935765996464)*x[0]
        arg[3,1,0,1]=(-0.186350920716)*x[1]
        arg[3,1,0,2]=(-0.809013452979)*x[2]
        arg[3,1,1,0]=(0.467966757364)*x[0]
        arg[3,1,1,1]=(-0.347420416508)*x[1]
        arg[3,1,1,2]=(0.686386725318)*x[2]
        arg[3,1,2,0]=(-0.8117641465)*x[0]
        arg[3,1,2,1]=(-0.357390874294)*x[1]
        arg[3,1,2,2]=(-0.0202262660381)*x[2]
        arg[3,2,0,0]=(-0.181594076044)*x[0]
        arg[3,2,0,1]=(0.0655758078072)*x[1]
        arg[3,2,0,2]=(0.926686460908)*x[2]
        arg[3,2,1,0]=(-0.700426183586)*x[0]
        arg[3,2,1,1]=(0.600692261192)*x[1]
        arg[3,2,1,2]=(0.757146807006)*x[2]
        arg[3,2,2,0]=(0.588097603826)*x[0]
        arg[3,2,2,1]=(-0.748732775295)*x[1]
        arg[3,2,2,2]=(0.628960378488)*x[2]
        arg[3,3,0,0]=(-0.0612899752684)*x[0]
        arg[3,3,0,1]=(0.567198903535)*x[1]
        arg[3,3,0,2]=(0.743996205698)*x[2]
        arg[3,3,1,0]=(0.333763157477)*x[0]
        arg[3,3,1,1]=(0.795773692731)*x[1]
        arg[3,3,1,2]=(-0.676301506816)*x[2]
        arg[3,3,2,0]=(-0.255504087739)*x[0]
        arg[3,3,2,1]=(-0.556263970561)*x[1]
        arg[3,3,2,2]=(0.721245074469)*x[2]
        arg[3,4,0,0]=(0.384551582302)*x[0]
        arg[3,4,0,1]=(0.849519144853)*x[1]
        arg[3,4,0,2]=(0.101881023262)*x[2]
        arg[3,4,1,0]=(-0.685891926653)*x[0]
        arg[3,4,1,1]=(0.390887672448)*x[1]
        arg[3,4,1,2]=(0.0911617386398)*x[2]
        arg[3,4,2,0]=(-0.889618312513)*x[0]
        arg[3,4,2,1]=(0.130310876981)*x[1]
        arg[3,4,2,2]=(0.958562278204)*x[2]
        ref=sqrt((12.9075635889)+(5.72789470359))

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunctionOnContactOne_fromData_rank0(self):
      """
      tests L2-norm of Data on the ReducedFunctionOnContactOne

      assumptions: self.domain supports integration on ReducedFunctionOnContactOne
      """
      dim=self.domain.getDim()
      w=ReducedFunctionOnContactOne(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(),w)
        arg=(0.39039157828)*sqrt(x[0])
        ref=sqrt((0.0)*3./2.+(0.038101396098)/0.5)

      else:
        arg=Data(0,(),w)
        arg=(0.893204658234)*sqrt(x[0])
        ref=sqrt((0.0)*3./2.+(0.199453640373)/0.5)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunctionOnContactOne_fromData_rank1(self):
      """
      tests L2-norm of Data on the ReducedFunctionOnContactOne

      assumptions: self.domain supports integration on ReducedFunctionOnContactOne
      """
      dim=self.domain.getDim()
      w=ReducedFunctionOnContactOne(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(2,),w)
        arg[0]=(0.697576620054)*sqrt(x[0])
        arg[1]=(-0.302037273777)*sqrt(x[1])
        ref=sqrt((0.0304088382502)*3./2.+(0.121653285212)/0.5)

      else:
        arg=Data(0,(3,),w)
        arg[0]=(0.649308272675)*sqrt(x[0])
        arg[1]=(-0.408829699836)*sqrt(x[1])
        arg[2]=(-0.984892977949)*sqrt(x[2])
        ref=sqrt((0.37905196716)*3./2.+(0.105400308241)/0.5)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunctionOnContactOne_fromData_rank2(self):
      """
      tests L2-norm of Data on the ReducedFunctionOnContactOne

      assumptions: self.domain supports integration on ReducedFunctionOnContactOne
      """
      dim=self.domain.getDim()
      w=ReducedFunctionOnContactOne(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 2),w)
        arg[0,0]=(0.504133645896)*sqrt(x[0])
        arg[0,1]=(0.728228500775)*sqrt(x[1])
        arg[1,0]=(0.968237164633)*sqrt(x[0])
        arg[1,1]=(-0.205351628957)*sqrt(x[1])
        arg[2,0]=(-0.154008228702)*sqrt(x[0])
        arg[2,1]=(-0.972088404566)*sqrt(x[1])
        arg[3,0]=(0.585153521593)*sqrt(x[0])
        arg[3,1]=(0.371861820001)*sqrt(x[1])
        ref=sqrt((0.551907706774)*3./2.+(0.389439279561)/0.5)

      else:
        arg=Data(0,(4, 3),w)
        arg[0,0]=(-0.4847970605)*sqrt(x[0])
        arg[0,1]=(-0.664433021479)*sqrt(x[1])
        arg[0,2]=(-0.93730914676)*sqrt(x[2])
        arg[1,0]=(0.712624065023)*sqrt(x[0])
        arg[1,1]=(0.105431976424)*sqrt(x[1])
        arg[1,2]=(-0.452636806104)*sqrt(x[2])
        arg[2,0]=(0.488406064176)*sqrt(x[0])
        arg[2,1]=(0.0594879976584)*sqrt(x[1])
        arg[2,2]=(0.596963883202)*sqrt(x[2])
        arg[3,0]=(0.619366395328)*sqrt(x[0])
        arg[3,1]=(0.719965901041)*sqrt(x[1])
        arg[3,2]=(-0.568971277169)*sqrt(x[2])
        ref=sqrt((0.912666523048)*3./2.+(0.341254115776)/0.5)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunctionOnContactOne_fromData_rank3(self):
      """
      tests L2-norm of Data on the ReducedFunctionOnContactOne

      assumptions: self.domain supports integration on ReducedFunctionOnContactOne
      """
      dim=self.domain.getDim()
      w=ReducedFunctionOnContactOne(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(6, 2, 2),w)
        arg[0,0,0]=(-0.481233278183)*sqrt(x[0])
        arg[0,0,1]=(-0.883447576048)*sqrt(x[1])
        arg[0,1,0]=(0.525938652633)*sqrt(x[0])
        arg[0,1,1]=(-0.694971850651)*sqrt(x[1])
        arg[1,0,0]=(-0.935880633291)*sqrt(x[0])
        arg[1,0,1]=(0.970445286468)*sqrt(x[1])
        arg[1,1,0]=(-0.302571977945)*sqrt(x[0])
        arg[1,1,1]=(-0.237750319821)*sqrt(x[1])
        arg[2,0,0]=(0.922940652883)*sqrt(x[0])
        arg[2,0,1]=(-0.624756614078)*sqrt(x[1])
        arg[2,1,0]=(-0.146839388764)*sqrt(x[0])
        arg[2,1,1]=(-0.995551211092)*sqrt(x[1])
        arg[3,0,0]=(0.719184541649)*sqrt(x[0])
        arg[3,0,1]=(-0.629461063375)*sqrt(x[1])
        arg[3,1,0]=(-0.925581215276)*sqrt(x[0])
        arg[3,1,1]=(-0.564418658495)*sqrt(x[1])
        arg[4,0,0]=(-0.724035884026)*sqrt(x[0])
        arg[4,0,1]=(0.78393984221)*sqrt(x[1])
        arg[4,1,0]=(-0.722375249317)*sqrt(x[0])
        arg[4,1,1]=(0.27571990944)*sqrt(x[1])
        arg[5,0,0]=(-0.197579392553)*sqrt(x[0])
        arg[5,0,1]=(0.341717398275)*sqrt(x[1])
        arg[5,1,0]=(0.926327822859)*sqrt(x[0])
        arg[5,1,1]=(0.44114488285)*sqrt(x[1])
        ref=sqrt((1.78665006238)*3./2.+(1.41652558894)/0.5)

      else:
        arg=Data(0,(6, 2, 3),w)
        arg[0,0,0]=(-0.348448140437)*sqrt(x[0])
        arg[0,0,1]=(-0.530793589782)*sqrt(x[1])
        arg[0,0,2]=(0.386511493376)*sqrt(x[2])
        arg[0,1,0]=(0.00776536299593)*sqrt(x[0])
        arg[0,1,1]=(-0.0696129074064)*sqrt(x[1])
        arg[0,1,2]=(0.114179504216)*sqrt(x[2])
        arg[1,0,0]=(-0.681009917406)*sqrt(x[0])
        arg[1,0,1]=(0.763595403719)*sqrt(x[1])
        arg[1,0,2]=(-0.54459342544)*sqrt(x[2])
        arg[1,1,0]=(0.0839255723157)*sqrt(x[0])
        arg[1,1,1]=(0.269389515003)*sqrt(x[1])
        arg[1,1,2]=(0.612446509453)*sqrt(x[2])
        arg[2,0,0]=(0.600070579238)*sqrt(x[0])
        arg[2,0,1]=(-0.49585867531)*sqrt(x[1])
        arg[2,0,2]=(-0.766857395141)*sqrt(x[2])
        arg[2,1,0]=(-0.84515816441)*sqrt(x[0])
        arg[2,1,1]=(-0.594337789049)*sqrt(x[1])
        arg[2,1,2]=(0.774596009811)*sqrt(x[2])
        arg[3,0,0]=(-0.833982234688)*sqrt(x[0])
        arg[3,0,1]=(0.747667281218)*sqrt(x[1])
        arg[3,0,2]=(-0.77025764051)*sqrt(x[2])
        arg[3,1,0]=(0.0367023286497)*sqrt(x[0])
        arg[3,1,1]=(-0.503080137993)*sqrt(x[1])
        arg[3,1,2]=(-0.129478323797)*sqrt(x[2])
        arg[4,0,0]=(0.304123530695)*sqrt(x[0])
        arg[4,0,1]=(-0.873023368289)*sqrt(x[1])
        arg[4,0,2]=(-0.712740806672)*sqrt(x[2])
        arg[4,1,0]=(-0.0629132979843)*sqrt(x[0])
        arg[4,1,1]=(0.096843438525)*sqrt(x[1])
        arg[4,1,2]=(-0.608716559621)*sqrt(x[2])
        arg[5,0,0]=(-0.642902822517)*sqrt(x[0])
        arg[5,0,1]=(0.478134785533)*sqrt(x[1])
        arg[5,0,2]=(0.757891512092)*sqrt(x[2])
        arg[5,1,0]=(-0.144492755065)*sqrt(x[0])
        arg[5,1,1]=(0.286775195938)*sqrt(x[1])
        arg[5,1,2]=(-0.063580488638)*sqrt(x[2])
        ref=sqrt((2.50835195543)*3./2.+(0.723549067206)/0.5)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunctionOnContactOne_fromData_rank4(self):
      """
      tests L2-norm of Data on the ReducedFunctionOnContactOne

      assumptions: self.domain supports integration on ReducedFunctionOnContactOne
      """
      dim=self.domain.getDim()
      w=ReducedFunctionOnContactOne(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 5, 3, 2),w)
        arg[0,0,0,0]=(0.523372259748)*sqrt(x[0])
        arg[0,0,0,1]=(-0.732795389957)*sqrt(x[1])
        arg[0,0,1,0]=(0.0571737128606)*sqrt(x[0])
        arg[0,0,1,1]=(-0.627260438721)*sqrt(x[1])
        arg[0,0,2,0]=(-0.945159375045)*sqrt(x[0])
        arg[0,0,2,1]=(-0.145103328237)*sqrt(x[1])
        arg[0,1,0,0]=(-0.0586985419455)*sqrt(x[0])
        arg[0,1,0,1]=(-0.788101334389)*sqrt(x[1])
        arg[0,1,1,0]=(0.615023802803)*sqrt(x[0])
        arg[0,1,1,1]=(0.965381283357)*sqrt(x[1])
        arg[0,1,2,0]=(0.173979432067)*sqrt(x[0])
        arg[0,1,2,1]=(-0.980483156424)*sqrt(x[1])
        arg[0,2,0,0]=(-0.345940909471)*sqrt(x[0])
        arg[0,2,0,1]=(0.00162430279595)*sqrt(x[1])
        arg[0,2,1,0]=(-0.0847413350082)*sqrt(x[0])
        arg[0,2,1,1]=(-0.911248036753)*sqrt(x[1])
        arg[0,2,2,0]=(-0.817295589942)*sqrt(x[0])
        arg[0,2,2,1]=(-0.0514569401185)*sqrt(x[1])
        arg[0,3,0,0]=(0.329704440159)*sqrt(x[0])
        arg[0,3,0,1]=(0.580379057442)*sqrt(x[1])
        arg[0,3,1,0]=(0.737014729605)*sqrt(x[0])
        arg[0,3,1,1]=(0.535521740242)*sqrt(x[1])
        arg[0,3,2,0]=(-0.943554721405)*sqrt(x[0])
        arg[0,3,2,1]=(0.798294295959)*sqrt(x[1])
        arg[0,4,0,0]=(-0.1814410847)*sqrt(x[0])
        arg[0,4,0,1]=(0.968476629577)*sqrt(x[1])
        arg[0,4,1,0]=(0.26595250332)*sqrt(x[0])
        arg[0,4,1,1]=(-0.493516007909)*sqrt(x[1])
        arg[0,4,2,0]=(0.374323739856)*sqrt(x[0])
        arg[0,4,2,1]=(-0.914413601829)*sqrt(x[1])
        arg[1,0,0,0]=(-0.701090812201)*sqrt(x[0])
        arg[1,0,0,1]=(-0.403596866322)*sqrt(x[1])
        arg[1,0,1,0]=(0.587183778075)*sqrt(x[0])
        arg[1,0,1,1]=(-0.218129059558)*sqrt(x[1])
        arg[1,0,2,0]=(-0.0205642762349)*sqrt(x[0])
        arg[1,0,2,1]=(0.441382235525)*sqrt(x[1])
        arg[1,1,0,0]=(0.909683859519)*sqrt(x[0])
        arg[1,1,0,1]=(0.960088832754)*sqrt(x[1])
        arg[1,1,1,0]=(0.0450355911985)*sqrt(x[0])
        arg[1,1,1,1]=(0.366849576658)*sqrt(x[1])
        arg[1,1,2,0]=(0.315490164497)*sqrt(x[0])
        arg[1,1,2,1]=(0.533746946808)*sqrt(x[1])
        arg[1,2,0,0]=(-0.552919155957)*sqrt(x[0])
        arg[1,2,0,1]=(-0.647373009727)*sqrt(x[1])
        arg[1,2,1,0]=(0.32617079311)*sqrt(x[0])
        arg[1,2,1,1]=(-0.3424263412)*sqrt(x[1])
        arg[1,2,2,0]=(0.848733118893)*sqrt(x[0])
        arg[1,2,2,1]=(-0.65341712364)*sqrt(x[1])
        arg[1,3,0,0]=(-0.698538235208)*sqrt(x[0])
        arg[1,3,0,1]=(-0.176858090319)*sqrt(x[1])
        arg[1,3,1,0]=(-0.956962782849)*sqrt(x[0])
        arg[1,3,1,1]=(0.53041598553)*sqrt(x[1])
        arg[1,3,2,0]=(-0.241715236701)*sqrt(x[0])
        arg[1,3,2,1]=(0.0340608544103)*sqrt(x[1])
        arg[1,4,0,0]=(-0.399279751219)*sqrt(x[0])
        arg[1,4,0,1]=(0.383860128496)*sqrt(x[1])
        arg[1,4,1,0]=(0.229189108307)*sqrt(x[0])
        arg[1,4,1,1]=(0.963667605308)*sqrt(x[1])
        arg[1,4,2,0]=(0.731912848998)*sqrt(x[0])
        arg[1,4,2,1]=(0.0534970206048)*sqrt(x[1])
        arg[2,0,0,0]=(0.422292915614)*sqrt(x[0])
        arg[2,0,0,1]=(0.915607515426)*sqrt(x[1])
        arg[2,0,1,0]=(0.451054723135)*sqrt(x[0])
        arg[2,0,1,1]=(0.0400870534188)*sqrt(x[1])
        arg[2,0,2,0]=(-0.851047453002)*sqrt(x[0])
        arg[2,0,2,1]=(0.300957420969)*sqrt(x[1])
        arg[2,1,0,0]=(0.636801515211)*sqrt(x[0])
        arg[2,1,0,1]=(0.386798062533)*sqrt(x[1])
        arg[2,1,1,0]=(-0.0501653606261)*sqrt(x[0])
        arg[2,1,1,1]=(0.46515271421)*sqrt(x[1])
        arg[2,1,2,0]=(-0.416436918527)*sqrt(x[0])
        arg[2,1,2,1]=(0.214018855333)*sqrt(x[1])
        arg[2,2,0,0]=(0.0516439357756)*sqrt(x[0])
        arg[2,2,0,1]=(-0.237707451356)*sqrt(x[1])
        arg[2,2,1,0]=(0.117559050121)*sqrt(x[0])
        arg[2,2,1,1]=(0.164349348553)*sqrt(x[1])
        arg[2,2,2,0]=(0.0524038362463)*sqrt(x[0])
        arg[2,2,2,1]=(0.690751059081)*sqrt(x[1])
        arg[2,3,0,0]=(-0.942079336547)*sqrt(x[0])
        arg[2,3,0,1]=(0.340143899875)*sqrt(x[1])
        arg[2,3,1,0]=(-0.630411753789)*sqrt(x[0])
        arg[2,3,1,1]=(-0.722078100418)*sqrt(x[1])
        arg[2,3,2,0]=(-0.310606406168)*sqrt(x[0])
        arg[2,3,2,1]=(0.728069014178)*sqrt(x[1])
        arg[2,4,0,0]=(0.981678768954)*sqrt(x[0])
        arg[2,4,0,1]=(-0.273023706956)*sqrt(x[1])
        arg[2,4,1,0]=(0.307187251211)*sqrt(x[0])
        arg[2,4,1,1]=(-0.885309661176)*sqrt(x[1])
        arg[2,4,2,0]=(-0.808022645995)*sqrt(x[0])
        arg[2,4,2,1]=(-0.21488080753)*sqrt(x[1])
        arg[3,0,0,0]=(-0.785416067695)*sqrt(x[0])
        arg[3,0,0,1]=(0.344457091999)*sqrt(x[1])
        arg[3,0,1,0]=(-0.845835855973)*sqrt(x[0])
        arg[3,0,1,1]=(0.365474280353)*sqrt(x[1])
        arg[3,0,2,0]=(-0.250535373488)*sqrt(x[0])
        arg[3,0,2,1]=(-0.888004656815)*sqrt(x[1])
        arg[3,1,0,0]=(-0.0652508248717)*sqrt(x[0])
        arg[3,1,0,1]=(0.402069695186)*sqrt(x[1])
        arg[3,1,1,0]=(-0.229714763188)*sqrt(x[0])
        arg[3,1,1,1]=(-0.866938136156)*sqrt(x[1])
        arg[3,1,2,0]=(-0.638603745536)*sqrt(x[0])
        arg[3,1,2,1]=(-0.997108734664)*sqrt(x[1])
        arg[3,2,0,0]=(-0.182784593825)*sqrt(x[0])
        arg[3,2,0,1]=(0.426951207538)*sqrt(x[1])
        arg[3,2,1,0]=(-0.958149736654)*sqrt(x[0])
        arg[3,2,1,1]=(-0.58458370587)*sqrt(x[1])
        arg[3,2,2,0]=(0.410243417149)*sqrt(x[0])
        arg[3,2,2,1]=(0.764474507608)*sqrt(x[1])
        arg[3,3,0,0]=(0.235623988334)*sqrt(x[0])
        arg[3,3,0,1]=(-0.105137401903)*sqrt(x[1])
        arg[3,3,1,0]=(-0.957743319072)*sqrt(x[0])
        arg[3,3,1,1]=(0.285680328399)*sqrt(x[1])
        arg[3,3,2,0]=(-0.525269717395)*sqrt(x[0])
        arg[3,3,2,1]=(0.655974370871)*sqrt(x[1])
        arg[3,4,0,0]=(0.926779313365)*sqrt(x[0])
        arg[3,4,0,1]=(-0.111556488031)*sqrt(x[1])
        arg[3,4,1,0]=(-0.672785588938)*sqrt(x[0])
        arg[3,4,1,1]=(-0.687907843899)*sqrt(x[1])
        arg[3,4,2,0]=(0.503984478085)*sqrt(x[0])
        arg[3,4,2,1]=(0.213194091516)*sqrt(x[1])
        ref=sqrt((6.92178561091)*3./2.+(4.96610502195)/0.5)

      else:
        arg=Data(0,(4, 5, 3, 3),w)
        arg[0,0,0,0]=(-0.641108959812)*sqrt(x[0])
        arg[0,0,0,1]=(0.395837684631)*sqrt(x[1])
        arg[0,0,0,2]=(0.36258240616)*sqrt(x[2])
        arg[0,0,1,0]=(-0.0543293722691)*sqrt(x[0])
        arg[0,0,1,1]=(0.850820297369)*sqrt(x[1])
        arg[0,0,1,2]=(0.0910730209593)*sqrt(x[2])
        arg[0,0,2,0]=(-0.831108730829)*sqrt(x[0])
        arg[0,0,2,1]=(-0.653018722047)*sqrt(x[1])
        arg[0,0,2,2]=(0.29052921548)*sqrt(x[2])
        arg[0,1,0,0]=(-0.115846569271)*sqrt(x[0])
        arg[0,1,0,1]=(-0.828080308748)*sqrt(x[1])
        arg[0,1,0,2]=(0.715255255273)*sqrt(x[2])
        arg[0,1,1,0]=(-0.0106077934044)*sqrt(x[0])
        arg[0,1,1,1]=(0.358323295554)*sqrt(x[1])
        arg[0,1,1,2]=(-0.916496266702)*sqrt(x[2])
        arg[0,1,2,0]=(0.740969432327)*sqrt(x[0])
        arg[0,1,2,1]=(0.0631507192754)*sqrt(x[1])
        arg[0,1,2,2]=(0.440009320758)*sqrt(x[2])
        arg[0,2,0,0]=(0.689964376385)*sqrt(x[0])
        arg[0,2,0,1]=(0.58199412545)*sqrt(x[1])
        arg[0,2,0,2]=(0.649058556936)*sqrt(x[2])
        arg[0,2,1,0]=(0.811660544475)*sqrt(x[0])
        arg[0,2,1,1]=(-0.121902248927)*sqrt(x[1])
        arg[0,2,1,2]=(-0.0750232179247)*sqrt(x[2])
        arg[0,2,2,0]=(0.779166450882)*sqrt(x[0])
        arg[0,2,2,1]=(0.298560101475)*sqrt(x[1])
        arg[0,2,2,2]=(-0.776198356206)*sqrt(x[2])
        arg[0,3,0,0]=(0.79804128408)*sqrt(x[0])
        arg[0,3,0,1]=(-0.378223308322)*sqrt(x[1])
        arg[0,3,0,2]=(0.191534109359)*sqrt(x[2])
        arg[0,3,1,0]=(-0.59973117885)*sqrt(x[0])
        arg[0,3,1,1]=(0.149042766249)*sqrt(x[1])
        arg[0,3,1,2]=(-0.357473233062)*sqrt(x[2])
        arg[0,3,2,0]=(-0.343758339283)*sqrt(x[0])
        arg[0,3,2,1]=(-0.140764743713)*sqrt(x[1])
        arg[0,3,2,2]=(-0.506377878901)*sqrt(x[2])
        arg[0,4,0,0]=(-0.967765607304)*sqrt(x[0])
        arg[0,4,0,1]=(-0.121373931498)*sqrt(x[1])
        arg[0,4,0,2]=(-0.300794608584)*sqrt(x[2])
        arg[0,4,1,0]=(-0.604537150389)*sqrt(x[0])
        arg[0,4,1,1]=(-0.290834102056)*sqrt(x[1])
        arg[0,4,1,2]=(0.509331591828)*sqrt(x[2])
        arg[0,4,2,0]=(0.643653985564)*sqrt(x[0])
        arg[0,4,2,1]=(-0.599746831159)*sqrt(x[1])
        arg[0,4,2,2]=(0.843193710828)*sqrt(x[2])
        arg[1,0,0,0]=(-0.598249432719)*sqrt(x[0])
        arg[1,0,0,1]=(0.840724468488)*sqrt(x[1])
        arg[1,0,0,2]=(0.558893493847)*sqrt(x[2])
        arg[1,0,1,0]=(0.843245438555)*sqrt(x[0])
        arg[1,0,1,1]=(0.280919545563)*sqrt(x[1])
        arg[1,0,1,2]=(-0.0549504830615)*sqrt(x[2])
        arg[1,0,2,0]=(0.295812998973)*sqrt(x[0])
        arg[1,0,2,1]=(0.513221249344)*sqrt(x[1])
        arg[1,0,2,2]=(-0.432253600492)*sqrt(x[2])
        arg[1,1,0,0]=(0.300696333449)*sqrt(x[0])
        arg[1,1,0,1]=(0.964202717779)*sqrt(x[1])
        arg[1,1,0,2]=(0.131730374204)*sqrt(x[2])
        arg[1,1,1,0]=(0.209135189207)*sqrt(x[0])
        arg[1,1,1,1]=(0.301389564979)*sqrt(x[1])
        arg[1,1,1,2]=(0.378481289816)*sqrt(x[2])
        arg[1,1,2,0]=(0.407413524533)*sqrt(x[0])
        arg[1,1,2,1]=(0.184345610049)*sqrt(x[1])
        arg[1,1,2,2]=(-0.884401323556)*sqrt(x[2])
        arg[1,2,0,0]=(0.886883417187)*sqrt(x[0])
        arg[1,2,0,1]=(-0.0524504232352)*sqrt(x[1])
        arg[1,2,0,2]=(0.049116549681)*sqrt(x[2])
        arg[1,2,1,0]=(-0.41240815391)*sqrt(x[0])
        arg[1,2,1,1]=(-0.523033558623)*sqrt(x[1])
        arg[1,2,1,2]=(-0.920238404199)*sqrt(x[2])
        arg[1,2,2,0]=(-0.635169424461)*sqrt(x[0])
        arg[1,2,2,1]=(0.773163888164)*sqrt(x[1])
        arg[1,2,2,2]=(0.554400056145)*sqrt(x[2])
        arg[1,3,0,0]=(0.573361031392)*sqrt(x[0])
        arg[1,3,0,1]=(-0.107749750058)*sqrt(x[1])
        arg[1,3,0,2]=(-0.971799077244)*sqrt(x[2])
        arg[1,3,1,0]=(-0.103579714499)*sqrt(x[0])
        arg[1,3,1,1]=(0.93329275867)*sqrt(x[1])
        arg[1,3,1,2]=(0.994409555402)*sqrt(x[2])
        arg[1,3,2,0]=(0.571961012279)*sqrt(x[0])
        arg[1,3,2,1]=(-0.54961086377)*sqrt(x[1])
        arg[1,3,2,2]=(0.563815465418)*sqrt(x[2])
        arg[1,4,0,0]=(-0.820638948871)*sqrt(x[0])
        arg[1,4,0,1]=(0.407750087109)*sqrt(x[1])
        arg[1,4,0,2]=(0.0911636159786)*sqrt(x[2])
        arg[1,4,1,0]=(0.371133349901)*sqrt(x[0])
        arg[1,4,1,1]=(0.688167208194)*sqrt(x[1])
        arg[1,4,1,2]=(0.188789217725)*sqrt(x[2])
        arg[1,4,2,0]=(0.961439941063)*sqrt(x[0])
        arg[1,4,2,1]=(0.524326385264)*sqrt(x[1])
        arg[1,4,2,2]=(-0.0866064394717)*sqrt(x[2])
        arg[2,0,0,0]=(0.392756517944)*sqrt(x[0])
        arg[2,0,0,1]=(0.713723092455)*sqrt(x[1])
        arg[2,0,0,2]=(0.68441988583)*sqrt(x[2])
        arg[2,0,1,0]=(-0.471505369012)*sqrt(x[0])
        arg[2,0,1,1]=(0.152712812211)*sqrt(x[1])
        arg[2,0,1,2]=(-0.579106568625)*sqrt(x[2])
        arg[2,0,2,0]=(0.637759469489)*sqrt(x[0])
        arg[2,0,2,1]=(-0.986405270597)*sqrt(x[1])
        arg[2,0,2,2]=(0.498109007475)*sqrt(x[2])
        arg[2,1,0,0]=(0.870185247694)*sqrt(x[0])
        arg[2,1,0,1]=(0.738095675991)*sqrt(x[1])
        arg[2,1,0,2]=(-0.613836737976)*sqrt(x[2])
        arg[2,1,1,0]=(-0.507009374688)*sqrt(x[0])
        arg[2,1,1,1]=(0.00752595652646)*sqrt(x[1])
        arg[2,1,1,2]=(-0.389118104751)*sqrt(x[2])
        arg[2,1,2,0]=(-0.466733433552)*sqrt(x[0])
        arg[2,1,2,1]=(-0.356829687828)*sqrt(x[1])
        arg[2,1,2,2]=(0.691337110432)*sqrt(x[2])
        arg[2,2,0,0]=(-0.684116912917)*sqrt(x[0])
        arg[2,2,0,1]=(0.705774798227)*sqrt(x[1])
        arg[2,2,0,2]=(0.522279733527)*sqrt(x[2])
        arg[2,2,1,0]=(0.390315182318)*sqrt(x[0])
        arg[2,2,1,1]=(-0.450165044241)*sqrt(x[1])
        arg[2,2,1,2]=(0.436062806791)*sqrt(x[2])
        arg[2,2,2,0]=(0.973600590724)*sqrt(x[0])
        arg[2,2,2,1]=(-0.89116205465)*sqrt(x[1])
        arg[2,2,2,2]=(-0.575262591791)*sqrt(x[2])
        arg[2,3,0,0]=(0.842833949385)*sqrt(x[0])
        arg[2,3,0,1]=(-0.69939817904)*sqrt(x[1])
        arg[2,3,0,2]=(0.343679860971)*sqrt(x[2])
        arg[2,3,1,0]=(-0.387030334358)*sqrt(x[0])
        arg[2,3,1,1]=(0.901637622406)*sqrt(x[1])
        arg[2,3,1,2]=(0.054184099786)*sqrt(x[2])
        arg[2,3,2,0]=(0.652712575782)*sqrt(x[0])
        arg[2,3,2,1]=(0.286474075274)*sqrt(x[1])
        arg[2,3,2,2]=(-0.745166160145)*sqrt(x[2])
        arg[2,4,0,0]=(-0.104591304669)*sqrt(x[0])
        arg[2,4,0,1]=(0.417496371683)*sqrt(x[1])
        arg[2,4,0,2]=(0.0579341684769)*sqrt(x[2])
        arg[2,4,1,0]=(0.24486914915)*sqrt(x[0])
        arg[2,4,1,1]=(-0.166211416717)*sqrt(x[1])
        arg[2,4,1,2]=(0.651515760357)*sqrt(x[2])
        arg[2,4,2,0]=(0.934693067557)*sqrt(x[0])
        arg[2,4,2,1]=(0.628486602755)*sqrt(x[1])
        arg[2,4,2,2]=(0.855313741494)*sqrt(x[2])
        arg[3,0,0,0]=(0.2241978058)*sqrt(x[0])
        arg[3,0,0,1]=(0.389284487657)*sqrt(x[1])
        arg[3,0,0,2]=(0.711521372962)*sqrt(x[2])
        arg[3,0,1,0]=(0.653184309521)*sqrt(x[0])
        arg[3,0,1,1]=(0.270052840275)*sqrt(x[1])
        arg[3,0,1,2]=(0.388881294549)*sqrt(x[2])
        arg[3,0,2,0]=(0.973583452692)*sqrt(x[0])
        arg[3,0,2,1]=(0.508175099792)*sqrt(x[1])
        arg[3,0,2,2]=(-0.989999307019)*sqrt(x[2])
        arg[3,1,0,0]=(0.935765996464)*sqrt(x[0])
        arg[3,1,0,1]=(-0.186350920716)*sqrt(x[1])
        arg[3,1,0,2]=(-0.809013452979)*sqrt(x[2])
        arg[3,1,1,0]=(0.467966757364)*sqrt(x[0])
        arg[3,1,1,1]=(-0.347420416508)*sqrt(x[1])
        arg[3,1,1,2]=(0.686386725318)*sqrt(x[2])
        arg[3,1,2,0]=(-0.8117641465)*sqrt(x[0])
        arg[3,1,2,1]=(-0.357390874294)*sqrt(x[1])
        arg[3,1,2,2]=(-0.0202262660381)*sqrt(x[2])
        arg[3,2,0,0]=(-0.181594076044)*sqrt(x[0])
        arg[3,2,0,1]=(0.0655758078072)*sqrt(x[1])
        arg[3,2,0,2]=(0.926686460908)*sqrt(x[2])
        arg[3,2,1,0]=(-0.700426183586)*sqrt(x[0])
        arg[3,2,1,1]=(0.600692261192)*sqrt(x[1])
        arg[3,2,1,2]=(0.757146807006)*sqrt(x[2])
        arg[3,2,2,0]=(0.588097603826)*sqrt(x[0])
        arg[3,2,2,1]=(-0.748732775295)*sqrt(x[1])
        arg[3,2,2,2]=(0.628960378488)*sqrt(x[2])
        arg[3,3,0,0]=(-0.0612899752684)*sqrt(x[0])
        arg[3,3,0,1]=(0.567198903535)*sqrt(x[1])
        arg[3,3,0,2]=(0.743996205698)*sqrt(x[2])
        arg[3,3,1,0]=(0.333763157477)*sqrt(x[0])
        arg[3,3,1,1]=(0.795773692731)*sqrt(x[1])
        arg[3,3,1,2]=(-0.676301506816)*sqrt(x[2])
        arg[3,3,2,0]=(-0.255504087739)*sqrt(x[0])
        arg[3,3,2,1]=(-0.556263970561)*sqrt(x[1])
        arg[3,3,2,2]=(0.721245074469)*sqrt(x[2])
        arg[3,4,0,0]=(0.384551582302)*sqrt(x[0])
        arg[3,4,0,1]=(0.849519144853)*sqrt(x[1])
        arg[3,4,0,2]=(0.101881023262)*sqrt(x[2])
        arg[3,4,1,0]=(-0.685891926653)*sqrt(x[0])
        arg[3,4,1,1]=(0.390887672448)*sqrt(x[1])
        arg[3,4,1,2]=(0.0911617386398)*sqrt(x[2])
        arg[3,4,2,0]=(-0.889618312513)*sqrt(x[0])
        arg[3,4,2,1]=(0.130310876981)*sqrt(x[1])
        arg[3,4,2,2]=(0.958562278204)*sqrt(x[2])
        ref=sqrt((12.9075635889)*3./2.+(5.72789470359)/0.5)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")

