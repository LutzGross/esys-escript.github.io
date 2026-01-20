
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

from test_util_spatial_functions2 import Test_Util_SpatialFunctions_noGradOnBoundary
from test_util_grad import Test_Util_Gradient


class Test_Util_SpatialFunctions(Test_Util_SpatialFunctions_noGradOnBoundary, Test_Util_Gradient):

   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onFunctionOnBoundary_fromData_ContinuousFunction(self):
      """
      tests divergence of Data on the FunctionOnBoundary

      assumptions: ContinuousFunction(self.domain) exists
                   self.domain supports div on FunctionOnBoundary
      """
      o=self.order
      dim=self.domain.getDim()
      w_ref=FunctionOnBoundary(self.domain)
      x_ref=w_ref.getX()
      w=ContinuousFunction(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(-0.39999188185)*x[0]**o+(-0.0985302665504)*x[0]+(0.103912148113)*x[1]**o+(-0.440893151166)*x[1]
        arg[1]=(-0.550867128413)*x[0]**o+(-0.141593933069)*x[0]+(0.70592012062)*x[1]**o+(-0.836308099686)*x[1]
        ref=o*(-0.39999188185)*x_ref[0]**(o-1)+o*(0.70592012062)*x_ref[1]**(o-1)+(-0.934838366236)
      else:
        arg[0]=(-0.00405079499937)*x[0]**o+(-0.282117848793)*x[0]+(0.979898720034)*x[1]**o+(-0.418106638625)*x[1]+(-0.443550810851)*x[2]**o+(-0.114976349698)*x[2]
        arg[1]=(0.194746221496)*x[0]**o+(0.324772666848)*x[0]+(0.887813387362)*x[1]**o+(-0.9362867149)*x[1]+(-0.837328978457)*x[2]**o+(0.666079514358)*x[2]
        arg[2]=(0.974797445328)*x[0]**o+(0.00365347225195)*x[0]+(-0.285413216102)*x[1]**o+(-0.253930177142)*x[1]+(0.306713249275)*x[2]**o+(-0.651133960743)*x[2]
        ref=o*(-0.00405079499937)*x_ref[0]**(o-1)+o*(0.887813387362)*x_ref[1]**(o-1)+o*(0.306713249275)*x_ref[2]**(o-1)+(-1.86953852444)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onFunctionOnBoundary_fromData_Solution(self):
      """
      tests divergence of Data on the FunctionOnBoundary

      assumptions: Solution(self.domain) exists
                   self.domain supports div on FunctionOnBoundary
      """
      o=self.order
      dim=self.domain.getDim()
      w_ref=FunctionOnBoundary(self.domain)
      x_ref=w_ref.getX()
      w=Solution(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(-0.756094856241)*x[0]**o+(-0.185186427186)*x[0]+(0.15709960228)*x[1]**o+(-0.408372399703)*x[1]
        arg[1]=(-0.680793761252)*x[0]**o+(0.815275471992)*x[0]+(-0.0564602888317)*x[1]**o+(0.542507598499)*x[1]
        ref=o*(-0.756094856241)*x_ref[0]**(o-1)+o*(-0.0564602888317)*x_ref[1]**(o-1)+(0.357321171313)
      else:
        arg[0]=(-0.934829971573)*x[0]**o+(-0.202157515712)*x[0]+(0.972095997183)*x[1]**o+(0.430890291958)*x[1]+(0.973333736054)*x[2]**o+(0.563067835398)*x[2]
        arg[1]=(-0.676905980848)*x[0]**o+(0.721377714968)*x[0]+(-0.594993636371)*x[1]**o+(0.678979489168)*x[1]+(-0.491382641255)*x[2]**o+(-0.237396572727)*x[2]
        arg[2]=(-0.429697616397)*x[0]**o+(-0.8280062193)*x[0]+(0.549112225165)*x[1]**o+(0.73490904828)*x[1]+(0.543450722632)*x[2]**o+(-0.910092410543)*x[2]
        ref=o*(-0.934829971573)*x_ref[0]**(o-1)+o*(-0.594993636371)*x_ref[1]**(o-1)+o*(0.543450722632)*x_ref[2]**(o-1)+(-0.433270437087)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onFunctionOnBoundary_fromData_ReducedSolution(self):
      """
      tests divergence of Data on the FunctionOnBoundary

      assumptions: ReducedSolution(self.domain) exists
                   self.domain supports div on FunctionOnBoundary
      """
      o=1
      dim=self.domain.getDim()
      w_ref=FunctionOnBoundary(self.domain)
      x_ref=w_ref.getX()
      w=ReducedSolution(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(-0.701055390328)*x[0]+(-0.873109677704)*x[1]
        arg[1]=(-0.459521943793)*x[0]+(-1.63064815751)*x[1]
        ref=(-2.33170354784)
      else:
        arg[0]=(-0.408202035271)*x[0]+(-0.114769362644)*x[1]+(1.5510990077)*x[2]
        arg[1]=(0.17481433607)*x[0]+(-0.466969145604)*x[1]+(1.25289823435)*x[2]
        arg[2]=(-0.215477376297)*x[0]+(0.437173510517)*x[1]+(1.00334359062)*x[2]
        ref=(0.128172409746)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onFunctionOnBoundary_fromData_ReducedContinuousFunction(self):
      """
      tests divergence of Data on the FunctionOnBoundary

      assumptions: ReducedContinuousFunction(self.domain) exists
                   self.domain supports div on FunctionOnBoundary
      """
      o=1
      dim=self.domain.getDim()
      w_ref=FunctionOnBoundary(self.domain)
      x_ref=w_ref.getX()
      w=ReducedContinuousFunction(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(-0.701055390328)*x[0]+(-0.873109677704)*x[1]
        arg[1]=(-0.459521943793)*x[0]+(-1.63064815751)*x[1]
        ref=(-2.33170354784)
      else:
        arg[0]=(-0.408202035271)*x[0]+(-0.114769362644)*x[1]+(1.5510990077)*x[2]
        arg[1]=(0.17481433607)*x[0]+(-0.466969145604)*x[1]+(1.25289823435)*x[2]
        arg[2]=(-0.215477376297)*x[0]+(0.437173510517)*x[1]+(1.00334359062)*x[2]
        ref=(0.128172409746)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onReducedFunctionOnBoundary_fromData_ContinuousFunction(self):
      """
      tests divergence of Data on the ReducedFunctionOnBoundary

      assumptions: ContinuousFunction(self.domain) exists
                   self.domain supports div on ReducedFunctionOnBoundary
      """
      o=self.order
      dim=self.domain.getDim()
      w_ref=ReducedFunctionOnBoundary(self.domain)
      x_ref=w_ref.getX()
      w=ContinuousFunction(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(-0.39999188185)*x[0]**o+(-0.0985302665504)*x[0]+(0.103912148113)*x[1]**o+(-0.440893151166)*x[1]
        arg[1]=(-0.550867128413)*x[0]**o+(-0.141593933069)*x[0]+(0.70592012062)*x[1]**o+(-0.836308099686)*x[1]
        ref=o*(-0.39999188185)*x_ref[0]**(o-1)+o*(0.70592012062)*x_ref[1]**(o-1)+(-0.934838366236)
      else:
        arg[0]=(-0.00405079499937)*x[0]**o+(-0.282117848793)*x[0]+(0.979898720034)*x[1]**o+(-0.418106638625)*x[1]+(-0.443550810851)*x[2]**o+(-0.114976349698)*x[2]
        arg[1]=(0.194746221496)*x[0]**o+(0.324772666848)*x[0]+(0.887813387362)*x[1]**o+(-0.9362867149)*x[1]+(-0.837328978457)*x[2]**o+(0.666079514358)*x[2]
        arg[2]=(0.974797445328)*x[0]**o+(0.00365347225195)*x[0]+(-0.285413216102)*x[1]**o+(-0.253930177142)*x[1]+(0.306713249275)*x[2]**o+(-0.651133960743)*x[2]
        ref=o*(-0.00405079499937)*x_ref[0]**(o-1)+o*(0.887813387362)*x_ref[1]**(o-1)+o*(0.306713249275)*x_ref[2]**(o-1)+(-1.86953852444)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onReducedFunctionOnBoundary_fromData_Solution(self):
      """
      tests divergence of Data on the ReducedFunctionOnBoundary

      assumptions: Solution(self.domain) exists
                   self.domain supports div on ReducedFunctionOnBoundary
      """
      o=self.order
      dim=self.domain.getDim()
      w_ref=ReducedFunctionOnBoundary(self.domain)
      x_ref=w_ref.getX()
      w=Solution(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(-0.756094856241)*x[0]**o+(-0.185186427186)*x[0]+(0.15709960228)*x[1]**o+(-0.408372399703)*x[1]
        arg[1]=(-0.680793761252)*x[0]**o+(0.815275471992)*x[0]+(-0.0564602888317)*x[1]**o+(0.542507598499)*x[1]
        ref=o*(-0.756094856241)*x_ref[0]**(o-1)+o*(-0.0564602888317)*x_ref[1]**(o-1)+(0.357321171313)
      else:
        arg[0]=(-0.934829971573)*x[0]**o+(-0.202157515712)*x[0]+(0.972095997183)*x[1]**o+(0.430890291958)*x[1]+(0.973333736054)*x[2]**o+(0.563067835398)*x[2]
        arg[1]=(-0.676905980848)*x[0]**o+(0.721377714968)*x[0]+(-0.594993636371)*x[1]**o+(0.678979489168)*x[1]+(-0.491382641255)*x[2]**o+(-0.237396572727)*x[2]
        arg[2]=(-0.429697616397)*x[0]**o+(-0.8280062193)*x[0]+(0.549112225165)*x[1]**o+(0.73490904828)*x[1]+(0.543450722632)*x[2]**o+(-0.910092410543)*x[2]
        ref=o*(-0.934829971573)*x_ref[0]**(o-1)+o*(-0.594993636371)*x_ref[1]**(o-1)+o*(0.543450722632)*x_ref[2]**(o-1)+(-0.433270437087)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onReducedFunctionOnBoundary_fromData_ReducedSolution(self):
      """
      tests divergence of Data on the ReducedFunctionOnBoundary

      assumptions: ReducedSolution(self.domain) exists
                   self.domain supports div on ReducedFunctionOnBoundary
      """
      o=1
      dim=self.domain.getDim()
      w_ref=ReducedFunctionOnBoundary(self.domain)
      x_ref=w_ref.getX()
      w=ReducedSolution(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(-0.701055390328)*x[0]+(-0.873109677704)*x[1]
        arg[1]=(-0.459521943793)*x[0]+(-1.63064815751)*x[1]
        ref=(-2.33170354784)
      else:
        arg[0]=(-0.408202035271)*x[0]+(-0.114769362644)*x[1]+(1.5510990077)*x[2]
        arg[1]=(0.17481433607)*x[0]+(-0.466969145604)*x[1]+(1.25289823435)*x[2]
        arg[2]=(-0.215477376297)*x[0]+(0.437173510517)*x[1]+(1.00334359062)*x[2]
        ref=(0.128172409746)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onReducedFunctionOnBoundary_fromData_ReducedContinuousFunction(self):
      """
      tests divergence of Data on the ReducedFunctionOnBoundary

      assumptions: ReducedContinuousFunction(self.domain) exists
                   self.domain supports div on ReducedFunctionOnBoundary
      """
      o=1
      dim=self.domain.getDim()
      w_ref=ReducedFunctionOnBoundary(self.domain)
      x_ref=w_ref.getX()
      w=ReducedContinuousFunction(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(-0.701055390328)*x[0]+(-0.873109677704)*x[1]
        arg[1]=(-0.459521943793)*x[0]+(-1.63064815751)*x[1]
        ref=(-2.33170354784)
      else:
        arg[0]=(-0.408202035271)*x[0]+(-0.114769362644)*x[1]+(1.5510990077)*x[2]
        arg[1]=(0.17481433607)*x[0]+(-0.466969145604)*x[1]+(1.25289823435)*x[2]
        arg[2]=(-0.215477376297)*x[0]+(0.437173510517)*x[1]+(1.00334359062)*x[2]
        ref=(0.128172409746)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onFunctionOnContactZero_fromData_ContinuousFunction(self):
      """
      tests divergence of Data on the FunctionOnContactZero

      assumptions: ContinuousFunction(self.domain) exists
                   self.domain supports div on FunctionOnContactZero
      """
      o=self.order
      dim=self.domain.getDim()
      w_ref=FunctionOnContactZero(self.domain)
      x_ref=w_ref.getX()
      w=ContinuousFunction(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(-0.168949755871)*x[0]**o+(-0.850028843808)*x[0]+(-0.120402543124)*x[1]**o+(-0.962848557803)*x[1]
        arg[1]=(-0.654809749448)*x[0]**o+(-0.86782287971)*x[0]+(0.063297622975)*x[1]**o+(0.400298481087)*x[1]
        ref=o*(-0.168949755871)*x_ref[0]**(o-1)+o*(0.063297622975)*x_ref[1]**(o-1)+(-0.449730362722)
      else:
        arg[0]=(-0.171176882625)*x[0]**o+(-0.983025148621)*x[0]+(-0.287555867477)*x[1]**o+(-0.689164237089)*x[1]+(0.700322207054)*x[2]**o+(0.224056812199)*x[2]
        arg[1]=(-0.112978574269)*x[0]**o+(-0.858210458939)*x[0]+(-0.9204250284)*x[1]**o+(0.123297295575)*x[1]+(-0.824406652345)*x[2]**o+(-0.514156554116)*x[2]
        arg[2]=(0.51215310895)*x[0]**o+(-0.484568237326)*x[0]+(0.324174170946)*x[1]**o+(0.00636990564332)*x[1]+(0.640205070676)*x[2]**o+(0.642111592847)*x[2]
        ref=o*(-0.171176882625)*x_ref[0]**(o-1)+o*(-0.9204250284)*x_ref[1]**(o-1)+o*(0.640205070676)*x_ref[2]**(o-1)+(-0.217616260198)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onFunctionOnContactZero_fromData_Solution(self):
      """
      tests divergence of Data on the FunctionOnContactZero

      assumptions: Solution(self.domain) exists
                   self.domain supports div on FunctionOnContactZero
      """
      o=self.order
      dim=self.domain.getDim()
      w_ref=FunctionOnContactZero(self.domain)
      x_ref=w_ref.getX()
      w=Solution(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(-0.878817592363)*x[0]**o+(-0.670917442592)*x[0]+(0.497874113437)*x[1]**o+(-0.923717382234)*x[1]
        arg[1]=(-0.26042852826)*x[0]**o+(0.0668495878627)*x[0]+(-0.24188250168)*x[1]**o+(-0.189346966858)*x[1]
        ref=o*(-0.878817592363)*x_ref[0]**(o-1)+o*(-0.24188250168)*x_ref[1]**(o-1)+(-0.86026440945)
      else:
        arg[0]=(-0.853829042745)*x[0]**o+(0.708253823804)*x[0]+(-0.420888906754)*x[1]**o+(0.94422817266)*x[1]+(-0.0454733411462)*x[2]**o+(-0.588636328569)*x[2]
        arg[1]=(0.502105411994)*x[0]**o+(0.618635348918)*x[0]+(-0.459898919811)*x[1]**o+(0.0206615239957)*x[1]+(-0.479378221929)*x[2]**o+(-0.0383028811998)*x[2]
        arg[2]=(-0.624273042239)*x[0]**o+(0.11218158851)*x[0]+(0.300065081794)*x[1]**o+(-0.823140829548)*x[1]+(-0.537442909756)*x[2]**o+(-0.537901446348)*x[2]
        ref=o*(-0.853829042745)*x_ref[0]**(o-1)+o*(-0.459898919811)*x_ref[1]**(o-1)+o*(-0.537442909756)*x_ref[2]**(o-1)+(0.191013901452)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onFunctionOnContactZero_fromData_ReducedSolution(self):
      """
      tests divergence of Data on the FunctionOnContactZero

      assumptions: ReducedSolution(self.domain) exists
                   self.domain supports div on FunctionOnContactZero
      """
      o=1
      dim=self.domain.getDim()
      w_ref=FunctionOnContactZero(self.domain)
      x_ref=w_ref.getX()
      w=ReducedSolution(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(1.40908912517)*x[0]+(0.980544777047)*x[1]
        arg[1]=(1.54735707911)*x[0]+(0.136617363775)*x[1]
        ref=(1.54570648895)
      else:
        arg[0]=(0.542899657951)*x[0]+(0.445471953165)*x[1]+(0.832668361385)*x[2]
        arg[1]=(-0.658976315278)*x[0]+(0.467903824571)*x[1]+(-1.00818943414)*x[2]
        arg[2]=(-1.10271841197)*x[0]+(1.48792855886)*x[1]+(-0.560275626738)*x[2]
        ref=(0.450527855784)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onFunctionOnContactZero_fromData_ReducedContinuousFunction(self):
      """
      tests divergence of Data on the FunctionOnContactZero

      assumptions: ReducedContinuousFunction(self.domain) exists
                   self.domain supports div on FunctionOnContactZero
      """
      o=1
      dim=self.domain.getDim()
      w_ref=FunctionOnContactZero(self.domain)
      x_ref=w_ref.getX()
      w=ReducedContinuousFunction(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(1.40908912517)*x[0]+(0.980544777047)*x[1]
        arg[1]=(1.54735707911)*x[0]+(0.136617363775)*x[1]
        ref=(1.54570648895)
      else:
        arg[0]=(0.542899657951)*x[0]+(0.445471953165)*x[1]+(0.832668361385)*x[2]
        arg[1]=(-0.658976315278)*x[0]+(0.467903824571)*x[1]+(-1.00818943414)*x[2]
        arg[2]=(-1.10271841197)*x[0]+(1.48792855886)*x[1]+(-0.560275626738)*x[2]
        ref=(0.450527855784)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onReducedFunctionOnContactZero_fromData_ContinuousFunction(self):
      """
      tests divergence of Data on the ReducedFunctionOnContactZero

      assumptions: ContinuousFunction(self.domain) exists
                   self.domain supports div on ReducedFunctionOnContactZero
      """
      o=self.order
      dim=self.domain.getDim()
      w_ref=ReducedFunctionOnContactZero(self.domain)
      x_ref=w_ref.getX()
      w=ContinuousFunction(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(-0.168949755871)*x[0]**o+(-0.850028843808)*x[0]+(-0.120402543124)*x[1]**o+(-0.962848557803)*x[1]
        arg[1]=(-0.654809749448)*x[0]**o+(-0.86782287971)*x[0]+(0.063297622975)*x[1]**o+(0.400298481087)*x[1]
        ref=o*(-0.168949755871)*x_ref[0]**(o-1)+o*(0.063297622975)*x_ref[1]**(o-1)+(-0.449730362722)
      else:
        arg[0]=(-0.171176882625)*x[0]**o+(-0.983025148621)*x[0]+(-0.287555867477)*x[1]**o+(-0.689164237089)*x[1]+(0.700322207054)*x[2]**o+(0.224056812199)*x[2]
        arg[1]=(-0.112978574269)*x[0]**o+(-0.858210458939)*x[0]+(-0.9204250284)*x[1]**o+(0.123297295575)*x[1]+(-0.824406652345)*x[2]**o+(-0.514156554116)*x[2]
        arg[2]=(0.51215310895)*x[0]**o+(-0.484568237326)*x[0]+(0.324174170946)*x[1]**o+(0.00636990564332)*x[1]+(0.640205070676)*x[2]**o+(0.642111592847)*x[2]
        ref=o*(-0.171176882625)*x_ref[0]**(o-1)+o*(-0.9204250284)*x_ref[1]**(o-1)+o*(0.640205070676)*x_ref[2]**(o-1)+(-0.217616260198)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onReducedFunctionOnContactZero_fromData_Solution(self):
      """
      tests divergence of Data on the ReducedFunctionOnContactZero

      assumptions: Solution(self.domain) exists
                   self.domain supports div on ReducedFunctionOnContactZero
      """
      o=self.order
      dim=self.domain.getDim()
      w_ref=ReducedFunctionOnContactZero(self.domain)
      x_ref=w_ref.getX()
      w=Solution(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(-0.878817592363)*x[0]**o+(-0.670917442592)*x[0]+(0.497874113437)*x[1]**o+(-0.923717382234)*x[1]
        arg[1]=(-0.26042852826)*x[0]**o+(0.0668495878627)*x[0]+(-0.24188250168)*x[1]**o+(-0.189346966858)*x[1]
        ref=o*(-0.878817592363)*x_ref[0]**(o-1)+o*(-0.24188250168)*x_ref[1]**(o-1)+(-0.86026440945)
      else:
        arg[0]=(-0.853829042745)*x[0]**o+(0.708253823804)*x[0]+(-0.420888906754)*x[1]**o+(0.94422817266)*x[1]+(-0.0454733411462)*x[2]**o+(-0.588636328569)*x[2]
        arg[1]=(0.502105411994)*x[0]**o+(0.618635348918)*x[0]+(-0.459898919811)*x[1]**o+(0.0206615239957)*x[1]+(-0.479378221929)*x[2]**o+(-0.0383028811998)*x[2]
        arg[2]=(-0.624273042239)*x[0]**o+(0.11218158851)*x[0]+(0.300065081794)*x[1]**o+(-0.823140829548)*x[1]+(-0.537442909756)*x[2]**o+(-0.537901446348)*x[2]
        ref=o*(-0.853829042745)*x_ref[0]**(o-1)+o*(-0.459898919811)*x_ref[1]**(o-1)+o*(-0.537442909756)*x_ref[2]**(o-1)+(0.191013901452)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onReducedFunctionOnContactZero_fromData_ReducedSolution(self):
      """
      tests divergence of Data on the ReducedFunctionOnContactZero

      assumptions: ReducedSolution(self.domain) exists
                   self.domain supports div on ReducedFunctionOnContactZero
      """
      o=1
      dim=self.domain.getDim()
      w_ref=ReducedFunctionOnContactZero(self.domain)
      x_ref=w_ref.getX()
      w=ReducedSolution(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(1.40908912517)*x[0]+(0.980544777047)*x[1]
        arg[1]=(1.54735707911)*x[0]+(0.136617363775)*x[1]
        ref=(1.54570648895)
      else:
        arg[0]=(0.542899657951)*x[0]+(0.445471953165)*x[1]+(0.832668361385)*x[2]
        arg[1]=(-0.658976315278)*x[0]+(0.467903824571)*x[1]+(-1.00818943414)*x[2]
        arg[2]=(-1.10271841197)*x[0]+(1.48792855886)*x[1]+(-0.560275626738)*x[2]
        ref=(0.450527855784)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onReducedFunctionOnContactZero_fromData_ReducedContinuousFunction(self):
      """
      tests divergence of Data on the ReducedFunctionOnContactZero

      assumptions: ReducedContinuousFunction(self.domain) exists
                   self.domain supports div on ReducedFunctionOnContactZero
      """
      o=1
      dim=self.domain.getDim()
      w_ref=ReducedFunctionOnContactZero(self.domain)
      x_ref=w_ref.getX()
      w=ReducedContinuousFunction(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(1.40908912517)*x[0]+(0.980544777047)*x[1]
        arg[1]=(1.54735707911)*x[0]+(0.136617363775)*x[1]
        ref=(1.54570648895)
      else:
        arg[0]=(0.542899657951)*x[0]+(0.445471953165)*x[1]+(0.832668361385)*x[2]
        arg[1]=(-0.658976315278)*x[0]+(0.467903824571)*x[1]+(-1.00818943414)*x[2]
        arg[2]=(-1.10271841197)*x[0]+(1.48792855886)*x[1]+(-0.560275626738)*x[2]
        ref=(0.450527855784)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onFunctionOnContactOne_fromData_ContinuousFunction(self):
      """
      tests divergence of Data on the FunctionOnContactOne

      assumptions: ContinuousFunction(self.domain) exists
                   self.domain supports div on FunctionOnContactOne
      """
      o=self.order
      dim=self.domain.getDim()
      w_ref=FunctionOnContactOne(self.domain)
      x_ref=w_ref.getX()
      w=ContinuousFunction(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(-0.506403770066)*x[0]**o+(0.32584268281)*x[0]+(-0.83529844868)*x[1]**o+(0.83697332673)*x[1]
        arg[1]=(-0.85254932286)*x[0]**o+(0.751060747065)*x[0]+(-0.672601723159)*x[1]**o+(0.33751172176)*x[1]
        ref=o*(-0.506403770066)*x_ref[0]**(o-1)+o*(-0.672601723159)*x_ref[1]**(o-1)+(0.66335440457)
      else:
        arg[0]=(-0.437811901761)*x[0]**o+(-0.135276240019)*x[0]+(0.683969739116)*x[1]**o+(-0.207547375093)*x[1]+(0.028594861031)*x[2]**o+(0.721966907461)*x[2]
        arg[1]=(0.315283027578)*x[0]**o+(0.724158692434)*x[0]+(-0.746372201791)*x[1]**o+(-0.370798791131)*x[1]+(-0.677893413238)*x[2]**o+(0.645317237187)*x[2]
        arg[2]=(0.0795681320657)*x[0]**o+(-0.165210743397)*x[0]+(0.866175461075)*x[1]**o+(0.10183311172)*x[1]+(-0.642389179022)*x[2]**o+(-0.917391184653)*x[2]
        ref=o*(-0.437811901761)*x_ref[0]**(o-1)+o*(-0.746372201791)*x_ref[1]**(o-1)+o*(-0.642389179022)*x_ref[2]**(o-1)+(-1.4234662158)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onFunctionOnContactOne_fromData_Solution(self):
      """
      tests divergence of Data on the FunctionOnContactOne

      assumptions: Solution(self.domain) exists
                   self.domain supports div on FunctionOnContactOne
      """
      o=self.order
      dim=self.domain.getDim()
      w_ref=FunctionOnContactOne(self.domain)
      x_ref=w_ref.getX()
      w=Solution(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(0.46385553428)*x[0]**o+(-0.91930502543)*x[0]+(0.316745110613)*x[1]**o+(-0.0996899204895)*x[1]
        arg[1]=(-0.473415736992)*x[0]**o+(0.514718374433)*x[0]+(-0.473429315265)*x[1]**o+(-0.194316615233)*x[1]
        ref=o*(0.46385553428)*x_ref[0]**(o-1)+o*(-0.473429315265)*x_ref[1]**(o-1)+(-1.11362164066)
      else:
        arg[0]=(0.740236725813)*x[0]**o+(0.887898498849)*x[0]+(-0.44803045798)*x[1]**o+(0.160930903094)*x[1]+(-0.530527194219)*x[2]**o+(0.0622347651074)*x[2]
        arg[1]=(0.36748785661)*x[0]**o+(-0.618052241987)*x[0]+(-0.57186325537)*x[1]**o+(-0.180576934525)*x[1]+(-0.381568503551)*x[2]**o+(0.0356882174068)*x[2]
        arg[2]=(0.070144257808)*x[0]**o+(0.823052625211)*x[0]+(-0.426522615959)*x[1]**o+(-0.294700653199)*x[1]+(0.209922433427)*x[2]**o+(-0.532056800218)*x[2]
        ref=o*(0.740236725813)*x_ref[0]**(o-1)+o*(-0.57186325537)*x_ref[1]**(o-1)+o*(0.209922433427)*x_ref[2]**(o-1)+(0.175264764106)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onFunctionOnContactOne_fromData_ReducedSolution(self):
      """
      tests divergence of Data on the FunctionOnContactOne

      assumptions: ReducedSolution(self.domain) exists
                   self.domain supports div on FunctionOnContactOne
      """
      o=1
      dim=self.domain.getDim()
      w_ref=FunctionOnContactOne(self.domain)
      x_ref=w_ref.getX()
      w=ReducedSolution(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(0.0337229572599)*x[0]+(0.820391566535)*x[1]
        arg[1]=(0.640567267185)*x[0]+(1.03537999206)*x[1]
        ref=(1.06910294932)
      else:
        arg[0]=(-0.0797501374913)*x[0]+(-0.798826730006)*x[1]+(-0.590979881432)*x[2]
        arg[1]=(-0.414547879032)*x[0]+(0.845538994985)*x[1]+(-1.21203164133)*x[2]
        arg[2]=(-0.438232052406)*x[0]+(-0.35977163998)*x[1]+(-1.01489636413)*x[2]
        ref=(-0.249107506635)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onFunctionOnContactOne_fromData_ReducedContinuousFunction(self):
      """
      tests divergence of Data on the FunctionOnContactOne

      assumptions: ReducedContinuousFunction(self.domain) exists
                   self.domain supports div on FunctionOnContactOne
      """
      o=1
      dim=self.domain.getDim()
      w_ref=FunctionOnContactOne(self.domain)
      x_ref=w_ref.getX()
      w=ReducedContinuousFunction(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(0.0337229572599)*x[0]+(0.820391566535)*x[1]
        arg[1]=(0.640567267185)*x[0]+(1.03537999206)*x[1]
        ref=(1.06910294932)
      else:
        arg[0]=(-0.0797501374913)*x[0]+(-0.798826730006)*x[1]+(-0.590979881432)*x[2]
        arg[1]=(-0.414547879032)*x[0]+(0.845538994985)*x[1]+(-1.21203164133)*x[2]
        arg[2]=(-0.438232052406)*x[0]+(-0.35977163998)*x[1]+(-1.01489636413)*x[2]
        ref=(-0.249107506635)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onReducedFunctionOnContactOne_fromData_ContinuousFunction(self):
      """
      tests divergence of Data on the ReducedFunctionOnContactOne

      assumptions: ContinuousFunction(self.domain) exists
                   self.domain supports div on ReducedFunctionOnContactOne
      """
      o=self.order
      dim=self.domain.getDim()
      w_ref=ReducedFunctionOnContactOne(self.domain)
      x_ref=w_ref.getX()
      w=ContinuousFunction(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(-0.506403770066)*x[0]**o+(0.32584268281)*x[0]+(-0.83529844868)*x[1]**o+(0.83697332673)*x[1]
        arg[1]=(-0.85254932286)*x[0]**o+(0.751060747065)*x[0]+(-0.672601723159)*x[1]**o+(0.33751172176)*x[1]
        ref=o*(-0.506403770066)*x_ref[0]**(o-1)+o*(-0.672601723159)*x_ref[1]**(o-1)+(0.66335440457)
      else:
        arg[0]=(-0.437811901761)*x[0]**o+(-0.135276240019)*x[0]+(0.683969739116)*x[1]**o+(-0.207547375093)*x[1]+(0.028594861031)*x[2]**o+(0.721966907461)*x[2]
        arg[1]=(0.315283027578)*x[0]**o+(0.724158692434)*x[0]+(-0.746372201791)*x[1]**o+(-0.370798791131)*x[1]+(-0.677893413238)*x[2]**o+(0.645317237187)*x[2]
        arg[2]=(0.0795681320657)*x[0]**o+(-0.165210743397)*x[0]+(0.866175461075)*x[1]**o+(0.10183311172)*x[1]+(-0.642389179022)*x[2]**o+(-0.917391184653)*x[2]
        ref=o*(-0.437811901761)*x_ref[0]**(o-1)+o*(-0.746372201791)*x_ref[1]**(o-1)+o*(-0.642389179022)*x_ref[2]**(o-1)+(-1.4234662158)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onReducedFunctionOnContactOne_fromData_Solution(self):
      """
      tests divergence of Data on the ReducedFunctionOnContactOne

      assumptions: Solution(self.domain) exists
                   self.domain supports div on ReducedFunctionOnContactOne
      """
      o=self.order
      dim=self.domain.getDim()
      w_ref=ReducedFunctionOnContactOne(self.domain)
      x_ref=w_ref.getX()
      w=Solution(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(0.46385553428)*x[0]**o+(-0.91930502543)*x[0]+(0.316745110613)*x[1]**o+(-0.0996899204895)*x[1]
        arg[1]=(-0.473415736992)*x[0]**o+(0.514718374433)*x[0]+(-0.473429315265)*x[1]**o+(-0.194316615233)*x[1]
        ref=o*(0.46385553428)*x_ref[0]**(o-1)+o*(-0.473429315265)*x_ref[1]**(o-1)+(-1.11362164066)
      else:
        arg[0]=(0.740236725813)*x[0]**o+(0.887898498849)*x[0]+(-0.44803045798)*x[1]**o+(0.160930903094)*x[1]+(-0.530527194219)*x[2]**o+(0.0622347651074)*x[2]
        arg[1]=(0.36748785661)*x[0]**o+(-0.618052241987)*x[0]+(-0.57186325537)*x[1]**o+(-0.180576934525)*x[1]+(-0.381568503551)*x[2]**o+(0.0356882174068)*x[2]
        arg[2]=(0.070144257808)*x[0]**o+(0.823052625211)*x[0]+(-0.426522615959)*x[1]**o+(-0.294700653199)*x[1]+(0.209922433427)*x[2]**o+(-0.532056800218)*x[2]
        ref=o*(0.740236725813)*x_ref[0]**(o-1)+o*(-0.57186325537)*x_ref[1]**(o-1)+o*(0.209922433427)*x_ref[2]**(o-1)+(0.175264764106)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onReducedFunctionOnContactOne_fromData_ReducedSolution(self):
      """
      tests divergence of Data on the ReducedFunctionOnContactOne

      assumptions: ReducedSolution(self.domain) exists
                   self.domain supports div on ReducedFunctionOnContactOne
      """
      o=1
      dim=self.domain.getDim()
      w_ref=ReducedFunctionOnContactOne(self.domain)
      x_ref=w_ref.getX()
      w=ReducedSolution(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(0.0337229572599)*x[0]+(0.820391566535)*x[1]
        arg[1]=(0.640567267185)*x[0]+(1.03537999206)*x[1]
        ref=(1.06910294932)
      else:
        arg[0]=(-0.0797501374913)*x[0]+(-0.798826730006)*x[1]+(-0.590979881432)*x[2]
        arg[1]=(-0.414547879032)*x[0]+(0.845538994985)*x[1]+(-1.21203164133)*x[2]
        arg[2]=(-0.438232052406)*x[0]+(-0.35977163998)*x[1]+(-1.01489636413)*x[2]
        ref=(-0.249107506635)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_div_onReducedFunctionOnContactOne_fromData_ReducedContinuousFunction(self):
      """
      tests divergence of Data on the ReducedFunctionOnContactOne

      assumptions: ReducedContinuousFunction(self.domain) exists
                   self.domain supports div on ReducedFunctionOnContactOne
      """
      o=1
      dim=self.domain.getDim()
      w_ref=ReducedFunctionOnContactOne(self.domain)
      x_ref=w_ref.getX()
      w=ReducedContinuousFunction(self.domain)
      x=w.getX()
      arg=Vector(0,w)
      if dim==2:
        arg[0]=(0.0337229572599)*x[0]+(0.820391566535)*x[1]
        arg[1]=(0.640567267185)*x[0]+(1.03537999206)*x[1]
        ref=(1.06910294932)
      else:
        arg[0]=(-0.0797501374913)*x[0]+(-0.798826730006)*x[1]+(-0.590979881432)*x[2]
        arg[1]=(-0.414547879032)*x[0]+(0.845538994985)*x[1]+(-1.21203164133)*x[2]
        arg[2]=(-0.438232052406)*x[0]+(-0.35977163998)*x[1]+(-1.01489636413)*x[2]
        ref=(-0.249107506635)
      res=div(arg,where=w_ref)
      self.assertTrue(isinstance(res,Data),"wrong type of result.")
      self.assertEqual(res.getShape(),(),"wrong shape of result.")
      self.assertEqual(res.getFunctionSpace(),w_ref,"wrong function space of result.")
      self.assertTrue(Lsup(res-ref)<=self.RES_TOL*Lsup(ref),"wrong result")

