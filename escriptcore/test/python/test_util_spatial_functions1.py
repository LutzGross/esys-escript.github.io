
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

from test_util_grad import Test_Util_Gradient_noBoundary
from test_util_integrals import Test_Util_Integration_noContact
from test_util_interpolation import Test_Util_Interpolation_noContact

class Test_Util_SpatialFunctions_noGradOnBoundary_noContact(Test_Util_Integration_noContact, Test_Util_Interpolation_noContact, Test_Util_Gradient_noBoundary):
   RES_TOL=1.e-8
   def test_x_ofDomain(self):
     """
     test getX() of the domain to be in the [0,1]^dim box
     """
     dim=self.domain.getDim()
     x=self.domain.getX()
     self.assertEqual(x.getShape(),(dim,),"wrong shape of result.")
     self.assertAlmostEqual(inf(x[0]),0.,int(-log10(self.RES_TOL)),"min x0 wrong")
     self.assertAlmostEqual(sup(x[0]),1.,int(-log10(self.RES_TOL)),"max x0 wrong")
     self.assertAlmostEqual(inf(x[1]),0.,int(-log10(self.RES_TOL)),"min x1 wrong")
     self.assertAlmostEqual(sup(x[1]),1.,int(-log10(self.RES_TOL)),"max x1 wrong")
     if dim>2:
       self.assertAlmostEqual(inf(x[2]),0.,int(-log10(self.RES_TOL)),"min x2 wrong")
       self.assertAlmostEqual(sup(x[2]),1.,int(-log10(self.RES_TOL)),"max x2 wrong")

   def test_SolutionOrder(self):
      """
      test the approximation order 
      """
      self.assertEqual(self.order, Solution(self.domain).getApproximationOrder(), "wrong order (Solution)")
      self.assertEqual(self.order, ContinuousFunction(self.domain).getApproximationOrder(), "wrong order (continuous function)")
      self.assertEqual(1, ReducedSolution(self.domain).getApproximationOrder(), "wrong order (ReducedSolution)")
      self.assertEqual(1, ReducedContinuousFunction(self.domain).getApproximationOrder(), "wrong order (Reduced continuous function)")
      for i in range(self.domain.getDim()):
         for k in range(Function(self.domain).getApproximationOrder()+1):
             self.assertAlmostEqual(integrate(Function(self.domain).getX()[i]**k),1./(k+1),8,"wrong integral (i=%s, order = %s)"%(i,k))
         for k in range(ReducedFunction(self.domain).getApproximationOrder()+1):
             self.assertAlmostEqual(integrate(ReducedFunction(self.domain).getX()[i]**k),1./(k+1),8,"wrong integral (i=%s, order = %s (reduced))"%(i,k))


   def test_normal_FunctionOnBoundary(self):
     """
     test getNormal() on boundary

     assumptions: FunctionOnBoundary(self.domain) exists
     """
     dim=self.domain.getDim()
     f=FunctionOnBoundary(self.domain)
     x=f.getX()
     ref=Vector(0.,what=f)
     if dim==3:
         ref.setTaggedValue(200,[0,0,1])
         ref.setTaggedValue(100,[0,0,-1])
         ref.setTaggedValue(20,[0,1,0])
         ref.setTaggedValue(10,[0,-1,0])
         ref.setTaggedValue(2,[1,0,0])
         ref.setTaggedValue(1,[-1,0,0])
     else:
         ref.setTaggedValue(2,[1,0])
         ref.setTaggedValue(1,[-1,0])
         ref.setTaggedValue(20, [0,1])
         ref.setTaggedValue(10, [0,-1])

     res=f.getNormal()
     self.assertEqual(res.getShape(),(dim,),"wrong shape of result.")
     self.assertEqual(res.getFunctionSpace(),f,"wrong functionspace of result.")
     self.assertLess(Lsup(ref-res), self.RES_TOL, "wrong result")

   def test_normal_ReducedFunctionOnBoundary(self):
     """
     test getNormal() on boundary

     assumptions: FunctionOnBoundary(self.domain) exists
     """
     dim=self.domain.getDim()
     f=ReducedFunctionOnBoundary(self.domain)
     x=f.getX()
     ref=Vector(0.,what=f)
     if dim==3:
         ref.setTaggedValue(200,[0,0,1])
         ref.setTaggedValue(100,[0,0,-1])
         ref.setTaggedValue(20,[0,1,0])
         ref.setTaggedValue(10,[0,-1,0])
         ref.setTaggedValue(2,[1,0,0])
         ref.setTaggedValue(1,[-1,0,0])
     else:
         ref.setTaggedValue(2,[1,0])
         ref.setTaggedValue(1,[-1,0])
         ref.setTaggedValue(20, [0,1])
         ref.setTaggedValue(10, [0,-1])

     res=f.getNormal()
     self.assertEqual(res.getShape(),(dim,),"wrong shape of result.")
     self.assertEqual(res.getFunctionSpace(),f,"wrong functionspace of result.")
     self.assertLess(Lsup(ref-res), self.RES_TOL, "wrong result")

   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunction_fromData_rank0(self):
      """
      tests L2-norm of Data on the Function

      assumptions: self.domain supports integration on Function
      """
      dim=self.domain.getDim()
      w=Function(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(),w)
        arg=(0.608797336225)*x[0]
        ref=sqrt((0.123544732198))

      else:
        arg=Data(0,(),w)
        arg=(0.136031275673)*x[0]
        ref=sqrt((0.00616816932037))

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunction_fromData_rank1(self):
      """
      tests L2-norm of Data on the Function

      assumptions: self.domain supports integration on Function
      """
      dim=self.domain.getDim()
      w=Function(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(2,),w)
        arg[0]=(-0.212143919436)*x[0]
        arg[1]=(-0.256194155686)*x[1]
        ref=sqrt((0.0368801626538))

      else:
        arg=Data(0,(3,),w)
        arg[0]=(0.0452831341416)*x[0]
        arg[1]=(-0.278640180656)*x[1]
        arg[2]=(-0.607035001062)*x[2]
        ref=sqrt((0.149394135009))

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunction_fromData_rank2(self):
      """
      tests L2-norm of Data on the Function

      assumptions: self.domain supports integration on Function
      """
      dim=self.domain.getDim()
      w=Function(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 2),w)
        arg[0,0]=(0.239448813076)*x[0]
        arg[0,1]=(-0.529349708753)*x[1]
        arg[1,0]=(-0.381557161859)*x[0]
        arg[1,1]=(0.731658534249)*x[1]
        arg[2,0]=(-0.813679062342)*x[0]
        arg[2,1]=(0.528100089704)*x[1]
        arg[3,0]=(-0.480867528161)*x[0]
        arg[3,1]=(-0.167862206972)*x[1]
        ref=sqrt((0.739610516051))

      else:
        arg=Data(0,(4, 3),w)
        arg[0,0]=(0.951209543612)*x[0]
        arg[0,1]=(0.735178735637)*x[1]
        arg[0,2]=(0.13074673272)*x[2]
        arg[1,0]=(0.412295676715)*x[0]
        arg[1,1]=(-0.657695950153)*x[1]
        arg[1,2]=(-0.900044734695)*x[2]
        arg[2,0]=(0.741773926224)*x[0]
        arg[2,1]=(0.0521828807406)*x[1]
        arg[2,2]=(0.797728501985)*x[2]
        arg[3,0]=(-0.61235554051)*x[0]
        arg[3,1]=(0.456652747412)*x[1]
        arg[3,2]=(-0.734303857319)*x[2]
        ref=sqrt((1.72901661926))

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunction_fromData_rank3(self):
      """
      tests L2-norm of Data on the Function

      assumptions: self.domain supports integration on Function
      """
      dim=self.domain.getDim()
      w=Function(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(6, 2, 2),w)
        arg[0,0,0]=(0.449174971953)*x[0]
        arg[0,0,1]=(-0.0109398763289)*x[1]
        arg[0,1,0]=(-0.202497187709)*x[0]
        arg[0,1,1]=(-0.12970879334)*x[1]
        arg[1,0,0]=(-0.138092481719)*x[0]
        arg[1,0,1]=(-0.528752200917)*x[1]
        arg[1,1,0]=(-0.605919441662)*x[0]
        arg[1,1,1]=(0.215615032334)*x[1]
        arg[2,0,0]=(-0.998734541972)*x[0]
        arg[2,0,1]=(0.725811901251)*x[1]
        arg[2,1,0]=(-0.966536503228)*x[0]
        arg[2,1,1]=(-0.528692217355)*x[1]
        arg[3,0,0]=(0.757633851466)*x[0]
        arg[3,0,1]=(-0.524660157377)*x[1]
        arg[3,1,0]=(0.983733431677)*x[0]
        arg[3,1,1]=(0.061279109546)*x[1]
        arg[4,0,0]=(0.85914215305)*x[0]
        arg[4,0,1]=(0.941714045112)*x[1]
        arg[4,1,0]=(0.172235529555)*x[0]
        arg[4,1,1]=(-0.108381454437)*x[1]
        arg[5,0,0]=(-0.736373697727)*x[0]
        arg[5,0,1]=(-0.599337929679)*x[1]
        arg[5,1,0]=(0.661072686392)*x[0]
        arg[5,1,1]=(-0.55107327409)*x[1]
        ref=sqrt((2.94641432714))

      else:
        arg=Data(0,(6, 2, 3),w)
        arg[0,0,0]=(0.69227064904)*x[0]
        arg[0,0,1]=(-0.968336177418)*x[1]
        arg[0,0,2]=(-0.634883146685)*x[2]
        arg[0,1,0]=(-0.12640661422)*x[0]
        arg[0,1,1]=(-0.637386589888)*x[1]
        arg[0,1,2]=(0.26060859356)*x[2]
        arg[1,0,0]=(-0.986864633297)*x[0]
        arg[1,0,1]=(-0.441589142379)*x[1]
        arg[1,0,2]=(-0.587865539582)*x[2]
        arg[1,1,0]=(0.596052465031)*x[0]
        arg[1,1,1]=(0.312732336652)*x[1]
        arg[1,1,2]=(-0.514423945092)*x[2]
        arg[2,0,0]=(-0.892391254794)*x[0]
        arg[2,0,1]=(0.377920185756)*x[1]
        arg[2,0,2]=(-0.120174597181)*x[2]
        arg[2,1,0]=(-0.469951576468)*x[0]
        arg[2,1,1]=(-0.788362249555)*x[1]
        arg[2,1,2]=(0.745625354986)*x[2]
        arg[3,0,0]=(0.542802498569)*x[0]
        arg[3,0,1]=(-0.814541028706)*x[1]
        arg[3,0,2]=(0.298410992196)*x[2]
        arg[3,1,0]=(0.981190341206)*x[0]
        arg[3,1,1]=(0.666421298608)*x[1]
        arg[3,1,2]=(-0.369751722626)*x[2]
        arg[4,0,0]=(-0.75379530597)*x[0]
        arg[4,0,1]=(0.283357267139)*x[1]
        arg[4,0,2]=(0.247787072861)*x[2]
        arg[4,1,0]=(0.301766692533)*x[0]
        arg[4,1,1]=(0.828183439224)*x[1]
        arg[4,1,2]=(-0.580824060547)*x[2]
        arg[5,0,0]=(0.637345610764)*x[0]
        arg[5,0,1]=(-0.234409115997)*x[1]
        arg[5,0,2]=(-0.192639300316)*x[2]
        arg[5,1,0]=(-0.62609237162)*x[0]
        arg[5,1,1]=(0.463404958552)*x[1]
        arg[5,1,2]=(-0.547814448738)*x[2]
        ref=sqrt((4.2381131862))

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunction_fromData_rank4(self):
      """
      tests L2-norm of Data on the Function

      assumptions: self.domain supports integration on Function
      """
      dim=self.domain.getDim()
      w=Function(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 5, 3, 2),w)
        arg[0,0,0,0]=(-0.232618585183)*x[0]
        arg[0,0,0,1]=(0.39796117869)*x[1]
        arg[0,0,1,0]=(-0.997336958262)*x[0]
        arg[0,0,1,1]=(-0.351780915076)*x[1]
        arg[0,0,2,0]=(-0.876764070136)*x[0]
        arg[0,0,2,1]=(0.808730805817)*x[1]
        arg[0,1,0,0]=(-0.197154744966)*x[0]
        arg[0,1,0,1]=(0.416246096086)*x[1]
        arg[0,1,1,0]=(0.708038457121)*x[0]
        arg[0,1,1,1]=(-0.00954021503188)*x[1]
        arg[0,1,2,0]=(-0.62630809425)*x[0]
        arg[0,1,2,1]=(0.430228727912)*x[1]
        arg[0,2,0,0]=(0.0306704747648)*x[0]
        arg[0,2,0,1]=(-0.913877199453)*x[1]
        arg[0,2,1,0]=(-0.697612800829)*x[0]
        arg[0,2,1,1]=(-0.17996376822)*x[1]
        arg[0,2,2,0]=(-0.304509578871)*x[0]
        arg[0,2,2,1]=(-0.610556755811)*x[1]
        arg[0,3,0,0]=(-0.452355972234)*x[0]
        arg[0,3,0,1]=(-0.368921242518)*x[1]
        arg[0,3,1,0]=(-0.478275554932)*x[0]
        arg[0,3,1,1]=(0.257178549127)*x[1]
        arg[0,3,2,0]=(0.530736487177)*x[0]
        arg[0,3,2,1]=(-0.567126272463)*x[1]
        arg[0,4,0,0]=(0.801519165938)*x[0]
        arg[0,4,0,1]=(-0.509816703951)*x[1]
        arg[0,4,1,0]=(-0.255412646934)*x[0]
        arg[0,4,1,1]=(0.437540101896)*x[1]
        arg[0,4,2,0]=(-0.815574969538)*x[0]
        arg[0,4,2,1]=(-0.94691547137)*x[1]
        arg[1,0,0,0]=(-0.732550722593)*x[0]
        arg[1,0,0,1]=(0.515752381704)*x[1]
        arg[1,0,1,0]=(-0.343590210899)*x[0]
        arg[1,0,1,1]=(-0.0601907964915)*x[1]
        arg[1,0,2,0]=(0.0199916154421)*x[0]
        arg[1,0,2,1]=(-0.136927227821)*x[1]
        arg[1,1,0,0]=(0.397994441702)*x[0]
        arg[1,1,0,1]=(0.953873148948)*x[1]
        arg[1,1,1,0]=(0.419416235967)*x[0]
        arg[1,1,1,1]=(0.700998577193)*x[1]
        arg[1,1,2,0]=(-0.497358799271)*x[0]
        arg[1,1,2,1]=(0.0851768858379)*x[1]
        arg[1,2,0,0]=(0.0936678875202)*x[0]
        arg[1,2,0,1]=(0.869883786896)*x[1]
        arg[1,2,1,0]=(0.582700123485)*x[0]
        arg[1,2,1,1]=(-0.433381106794)*x[1]
        arg[1,2,2,0]=(-0.527031777974)*x[0]
        arg[1,2,2,1]=(0.105105137652)*x[1]
        arg[1,3,0,0]=(-0.716750829134)*x[0]
        arg[1,3,0,1]=(0.774519209008)*x[1]
        arg[1,3,1,0]=(-0.568743372716)*x[0]
        arg[1,3,1,1]=(0.794732483944)*x[1]
        arg[1,3,2,0]=(0.246606002015)*x[0]
        arg[1,3,2,1]=(-0.988869494994)*x[1]
        arg[1,4,0,0]=(0.482379298083)*x[0]
        arg[1,4,0,1]=(-0.386268387903)*x[1]
        arg[1,4,1,0]=(0.137184889675)*x[0]
        arg[1,4,1,1]=(-0.140520035321)*x[1]
        arg[1,4,2,0]=(0.822755050415)*x[0]
        arg[1,4,2,1]=(-0.815562139522)*x[1]
        arg[2,0,0,0]=(-0.462891511962)*x[0]
        arg[2,0,0,1]=(-0.122643411631)*x[1]
        arg[2,0,1,0]=(-0.520861119962)*x[0]
        arg[2,0,1,1]=(-0.881189618018)*x[1]
        arg[2,0,2,0]=(-0.776157842774)*x[0]
        arg[2,0,2,1]=(-0.12354053207)*x[1]
        arg[2,1,0,0]=(0.395495230826)*x[0]
        arg[2,1,0,1]=(-0.388106659423)*x[1]
        arg[2,1,1,0]=(0.354250242834)*x[0]
        arg[2,1,1,1]=(-0.666514210192)*x[1]
        arg[2,1,2,0]=(0.951294655083)*x[0]
        arg[2,1,2,1]=(0.074024416386)*x[1]
        arg[2,2,0,0]=(0.335448485459)*x[0]
        arg[2,2,0,1]=(-0.40988282528)*x[1]
        arg[2,2,1,0]=(-0.805725968875)*x[0]
        arg[2,2,1,1]=(-0.949883082118)*x[1]
        arg[2,2,2,0]=(0.531549210683)*x[0]
        arg[2,2,2,1]=(-0.398401016682)*x[1]
        arg[2,3,0,0]=(-0.953963433205)*x[0]
        arg[2,3,0,1]=(0.643431126406)*x[1]
        arg[2,3,1,0]=(-0.167611998738)*x[0]
        arg[2,3,1,1]=(0.226130056552)*x[1]
        arg[2,3,2,0]=(0.0752687641131)*x[0]
        arg[2,3,2,1]=(-0.115742756362)*x[1]
        arg[2,4,0,0]=(0.579694491028)*x[0]
        arg[2,4,0,1]=(-0.112005738299)*x[1]
        arg[2,4,1,0]=(0.657291764224)*x[0]
        arg[2,4,1,1]=(0.62671154177)*x[1]
        arg[2,4,2,0]=(0.103695027944)*x[0]
        arg[2,4,2,1]=(0.462828491544)*x[1]
        arg[3,0,0,0]=(0.697692979998)*x[0]
        arg[3,0,0,1]=(-0.123481859619)*x[1]
        arg[3,0,1,0]=(-0.749745629459)*x[0]
        arg[3,0,1,1]=(-0.541969524069)*x[1]
        arg[3,0,2,0]=(0.819484470759)*x[0]
        arg[3,0,2,1]=(-0.860592326469)*x[1]
        arg[3,1,0,0]=(-0.716566084771)*x[0]
        arg[3,1,0,1]=(-0.949235434827)*x[1]
        arg[3,1,1,0]=(-0.826699498174)*x[0]
        arg[3,1,1,1]=(-0.138511521583)*x[1]
        arg[3,1,2,0]=(-0.951682890904)*x[0]
        arg[3,1,2,1]=(0.413293316925)*x[1]
        arg[3,2,0,0]=(0.909516836775)*x[0]
        arg[3,2,0,1]=(-0.919989721277)*x[1]
        arg[3,2,1,0]=(0.0994860369337)*x[0]
        arg[3,2,1,1]=(-0.933647246623)*x[1]
        arg[3,2,2,0]=(-0.759215183015)*x[0]
        arg[3,2,2,1]=(0.0975793309286)*x[1]
        arg[3,3,0,0]=(-0.130256739381)*x[0]
        arg[3,3,0,1]=(-0.582280862311)*x[1]
        arg[3,3,1,0]=(0.206970526192)*x[0]
        arg[3,3,1,1]=(-0.8678322258)*x[1]
        arg[3,3,2,0]=(0.133004501279)*x[0]
        arg[3,3,2,1]=(0.802921710935)*x[1]
        arg[3,4,0,0]=(-0.255737792764)*x[0]
        arg[3,4,0,1]=(-0.34168114937)*x[1]
        arg[3,4,1,0]=(-0.859309090399)*x[0]
        arg[3,4,1,1]=(0.245043986435)*x[1]
        arg[3,4,2,0]=(0.893062018695)*x[0]
        arg[3,4,2,1]=(0.709422742588)*x[1]
        ref=sqrt((13.7289280362))

      else:
        arg=Data(0,(4, 5, 3, 3),w)
        arg[0,0,0,0]=(0.0312828390439)*x[0]
        arg[0,0,0,1]=(-0.524970416212)*x[1]
        arg[0,0,0,2]=(0.561865217554)*x[2]
        arg[0,0,1,0]=(0.692457187384)*x[0]
        arg[0,0,1,1]=(0.946967182157)*x[1]
        arg[0,0,1,2]=(-0.863842279464)*x[2]
        arg[0,0,2,0]=(0.993922921598)*x[0]
        arg[0,0,2,1]=(0.322812768679)*x[1]
        arg[0,0,2,2]=(0.901876132204)*x[2]
        arg[0,1,0,0]=(0.967569979365)*x[0]
        arg[0,1,0,1]=(0.840979131355)*x[1]
        arg[0,1,0,2]=(0.0494811460856)*x[2]
        arg[0,1,1,0]=(0.315178456102)*x[0]
        arg[0,1,1,1]=(0.449848313024)*x[1]
        arg[0,1,1,2]=(0.765887852886)*x[2]
        arg[0,1,2,0]=(0.975541574352)*x[0]
        arg[0,1,2,1]=(-0.797851290751)*x[1]
        arg[0,1,2,2]=(0.628918775319)*x[2]
        arg[0,2,0,0]=(0.685635794312)*x[0]
        arg[0,2,0,1]=(0.10341799962)*x[1]
        arg[0,2,0,2]=(-0.964822756043)*x[2]
        arg[0,2,1,0]=(-0.56160368212)*x[0]
        arg[0,2,1,1]=(0.676344298102)*x[1]
        arg[0,2,1,2]=(-0.713924121843)*x[2]
        arg[0,2,2,0]=(-0.276655136263)*x[0]
        arg[0,2,2,1]=(0.336046973788)*x[1]
        arg[0,2,2,2]=(-0.68789392396)*x[2]
        arg[0,3,0,0]=(0.0172861311571)*x[0]
        arg[0,3,0,1]=(-0.301075956456)*x[1]
        arg[0,3,0,2]=(0.779442985415)*x[2]
        arg[0,3,1,0]=(-0.517629576558)*x[0]
        arg[0,3,1,1]=(0.584779586639)*x[1]
        arg[0,3,1,2]=(-0.53266435436)*x[2]
        arg[0,3,2,0]=(0.841533567102)*x[0]
        arg[0,3,2,1]=(0.0458746415489)*x[1]
        arg[0,3,2,2]=(0.921237870758)*x[2]
        arg[0,4,0,0]=(0.0548343238805)*x[0]
        arg[0,4,0,1]=(0.687022707412)*x[1]
        arg[0,4,0,2]=(-0.319803609795)*x[2]
        arg[0,4,1,0]=(0.409763007811)*x[0]
        arg[0,4,1,1]=(0.165501957435)*x[1]
        arg[0,4,1,2]=(0.116001692781)*x[2]
        arg[0,4,2,0]=(-0.515571394238)*x[0]
        arg[0,4,2,1]=(0.209467945147)*x[1]
        arg[0,4,2,2]=(-0.344827191247)*x[2]
        arg[1,0,0,0]=(0.57193838014)*x[0]
        arg[1,0,0,1]=(-0.0880683799076)*x[1]
        arg[1,0,0,2]=(0.956899617441)*x[2]
        arg[1,0,1,0]=(-0.783689636357)*x[0]
        arg[1,0,1,1]=(-0.25177506885)*x[1]
        arg[1,0,1,2]=(-0.97074584634)*x[2]
        arg[1,0,2,0]=(0.432543519806)*x[0]
        arg[1,0,2,1]=(0.481003021954)*x[1]
        arg[1,0,2,2]=(-0.0630751518268)*x[2]
        arg[1,1,0,0]=(-0.65152446796)*x[0]
        arg[1,1,0,1]=(-0.0323685084425)*x[1]
        arg[1,1,0,2]=(-0.508674033909)*x[2]
        arg[1,1,1,0]=(-0.533367818916)*x[0]
        arg[1,1,1,1]=(0.310738340288)*x[1]
        arg[1,1,1,2]=(0.694612234326)*x[2]
        arg[1,1,2,0]=(-0.622052473032)*x[0]
        arg[1,1,2,1]=(0.0498443793671)*x[1]
        arg[1,1,2,2]=(0.61023707512)*x[2]
        arg[1,2,0,0]=(0.0730267406859)*x[0]
        arg[1,2,0,1]=(0.146909334607)*x[1]
        arg[1,2,0,2]=(-0.641860284448)*x[2]
        arg[1,2,1,0]=(0.917976589737)*x[0]
        arg[1,2,1,1]=(0.50219672122)*x[1]
        arg[1,2,1,2]=(0.634559579812)*x[2]
        arg[1,2,2,0]=(0.0578772734534)*x[0]
        arg[1,2,2,1]=(0.288730973517)*x[1]
        arg[1,2,2,2]=(-0.0525978796154)*x[2]
        arg[1,3,0,0]=(-0.926152433388)*x[0]
        arg[1,3,0,1]=(0.0616647680855)*x[1]
        arg[1,3,0,2]=(-0.875889217846)*x[2]
        arg[1,3,1,0]=(-0.638931542845)*x[0]
        arg[1,3,1,1]=(0.708848122964)*x[1]
        arg[1,3,1,2]=(0.119066979792)*x[2]
        arg[1,3,2,0]=(0.853716218591)*x[0]
        arg[1,3,2,1]=(-0.92754322201)*x[1]
        arg[1,3,2,2]=(-0.671530626265)*x[2]
        arg[1,4,0,0]=(0.337424536231)*x[0]
        arg[1,4,0,1]=(0.335704451719)*x[1]
        arg[1,4,0,2]=(-0.484565969466)*x[2]
        arg[1,4,1,0]=(-0.855476192012)*x[0]
        arg[1,4,1,1]=(0.405674615553)*x[1]
        arg[1,4,1,2]=(0.728310771323)*x[2]
        arg[1,4,2,0]=(0.363651308265)*x[0]
        arg[1,4,2,1]=(0.174460594531)*x[1]
        arg[1,4,2,2]=(-0.0418244838617)*x[2]
        arg[2,0,0,0]=(-0.531341992511)*x[0]
        arg[2,0,0,1]=(0.584996796272)*x[1]
        arg[2,0,0,2]=(-0.752430968716)*x[2]
        arg[2,0,1,0]=(-0.341989849747)*x[0]
        arg[2,0,1,1]=(0.153572646953)*x[1]
        arg[2,0,1,2]=(-0.197130051737)*x[2]
        arg[2,0,2,0]=(-0.338082424082)*x[0]
        arg[2,0,2,1]=(0.000173657394772)*x[1]
        arg[2,0,2,2]=(0.365272907692)*x[2]
        arg[2,1,0,0]=(0.904304126564)*x[0]
        arg[2,1,0,1]=(0.161252368484)*x[1]
        arg[2,1,0,2]=(0.246854092422)*x[2]
        arg[2,1,1,0]=(-0.299880647529)*x[0]
        arg[2,1,1,1]=(-0.566917528608)*x[1]
        arg[2,1,1,2]=(0.243183337285)*x[2]
        arg[2,1,2,0]=(0.437406011474)*x[0]
        arg[2,1,2,1]=(0.727447394053)*x[1]
        arg[2,1,2,2]=(0.380752950664)*x[2]
        arg[2,2,0,0]=(0.172292846911)*x[0]
        arg[2,2,0,1]=(0.334201791643)*x[1]
        arg[2,2,0,2]=(0.739989926962)*x[2]
        arg[2,2,1,0]=(-0.0669843715042)*x[0]
        arg[2,2,1,1]=(-0.540497281635)*x[1]
        arg[2,2,1,2]=(-0.744217027088)*x[2]
        arg[2,2,2,0]=(-0.287295952259)*x[0]
        arg[2,2,2,1]=(-0.512411849183)*x[1]
        arg[2,2,2,2]=(0.953107417666)*x[2]
        arg[2,3,0,0]=(0.998168116695)*x[0]
        arg[2,3,0,1]=(0.960065646359)*x[1]
        arg[2,3,0,2]=(0.110048258832)*x[2]
        arg[2,3,1,0]=(-0.477271134724)*x[0]
        arg[2,3,1,1]=(0.707182612251)*x[1]
        arg[2,3,1,2]=(0.285500891755)*x[2]
        arg[2,3,2,0]=(-0.863497506661)*x[0]
        arg[2,3,2,1]=(-0.293917669879)*x[1]
        arg[2,3,2,2]=(-0.403384244295)*x[2]
        arg[2,4,0,0]=(0.848455277702)*x[0]
        arg[2,4,0,1]=(-0.530101455578)*x[1]
        arg[2,4,0,2]=(0.33887313048)*x[2]
        arg[2,4,1,0]=(-0.195313538124)*x[0]
        arg[2,4,1,1]=(-0.62754572008)*x[1]
        arg[2,4,1,2]=(-0.385132960582)*x[2]
        arg[2,4,2,0]=(0.240048012886)*x[0]
        arg[2,4,2,1]=(0.900766252969)*x[1]
        arg[2,4,2,2]=(0.669620533505)*x[2]
        arg[3,0,0,0]=(0.375766827301)*x[0]
        arg[3,0,0,1]=(0.705484960308)*x[1]
        arg[3,0,0,2]=(0.440931516034)*x[2]
        arg[3,0,1,0]=(-0.44724403177)*x[0]
        arg[3,0,1,1]=(-0.31558249626)*x[1]
        arg[3,0,1,2]=(-0.00419436365172)*x[2]
        arg[3,0,2,0]=(0.750599752032)*x[0]
        arg[3,0,2,1]=(0.367649951795)*x[1]
        arg[3,0,2,2]=(0.0488013073654)*x[2]
        arg[3,1,0,0]=(-0.992890068274)*x[0]
        arg[3,1,0,1]=(0.671447745511)*x[1]
        arg[3,1,0,2]=(0.85613331404)*x[2]
        arg[3,1,1,0]=(-0.46064764242)*x[0]
        arg[3,1,1,1]=(0.48138877715)*x[1]
        arg[3,1,1,2]=(0.396741761803)*x[2]
        arg[3,1,2,0]=(-0.879391967543)*x[0]
        arg[3,1,2,1]=(-0.44039462138)*x[1]
        arg[3,1,2,2]=(0.0330511573872)*x[2]
        arg[3,2,0,0]=(-0.367413701648)*x[0]
        arg[3,2,0,1]=(0.0359818324891)*x[1]
        arg[3,2,0,2]=(-0.307532667032)*x[2]
        arg[3,2,1,0]=(0.334663597166)*x[0]
        arg[3,2,1,1]=(0.541941978066)*x[1]
        arg[3,2,1,2]=(-0.609184079318)*x[2]
        arg[3,2,2,0]=(0.359349239826)*x[0]
        arg[3,2,2,1]=(0.0419272305685)*x[1]
        arg[3,2,2,2]=(0.557189794296)*x[2]
        arg[3,3,0,0]=(-0.85864165554)*x[0]
        arg[3,3,0,1]=(-0.185411404213)*x[1]
        arg[3,3,0,2]=(0.254294865253)*x[2]
        arg[3,3,1,0]=(0.870362177541)*x[0]
        arg[3,3,1,1]=(-0.439688612864)*x[1]
        arg[3,3,1,2]=(0.26006729357)*x[2]
        arg[3,3,2,0]=(-0.0724034754175)*x[0]
        arg[3,3,2,1]=(0.444871564246)*x[1]
        arg[3,3,2,2]=(0.485634530531)*x[2]
        arg[3,4,0,0]=(-0.744756961758)*x[0]
        arg[3,4,0,1]=(0.429761406102)*x[1]
        arg[3,4,0,2]=(-0.584963735834)*x[2]
        arg[3,4,1,0]=(0.684578379159)*x[0]
        arg[3,4,1,1]=(0.949460132601)*x[1]
        arg[3,4,1,2]=(-0.592179909559)*x[2]
        arg[3,4,2,0]=(0.707154437797)*x[0]
        arg[3,4,2,1]=(0.619200407063)*x[1]
        arg[3,4,2,2]=(-0.338547165)*x[2]
        ref=sqrt((19.2170638478))

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunction_fromData_rank0(self):
      """
      tests L2-norm of Data on the ReducedFunction

      assumptions: self.domain supports integration on ReducedFunction
      """
      dim=self.domain.getDim()
      w=ReducedFunction(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(),w)
        arg=1.*sqrt(x[0])
        ref=sqrt(0.5)

      else:
        arg=Data(0,(),w)
        arg=1.*sqrt(x[0])
        ref=sqrt(0.5)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunction_fromData_rank1(self):
      """
      tests L2-norm of Data on the ReducedFunction

      assumptions: self.domain supports integration on ReducedFunction
      """
      dim=self.domain.getDim()
      w=ReducedFunction(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(2,),w)
        arg[0]=1.*sqrt(x[0])
        arg[1]=2.*sqrt(x[1])
        ref=sqrt(2.5)

      else:
        arg=Data(0,(3,),w)
        arg[0]=1.*sqrt(x[0])
        arg[1]=2.*sqrt(x[1])
        arg[2]=3.*sqrt(x[2])
        ref=sqrt(7.)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunction_fromData_rank2(self):
      """
      tests L2-norm of Data on the ReducedFunction

      assumptions: self.domain supports integration on ReducedFunction
      """
      dim=self.domain.getDim()
      w=ReducedFunction(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 2),w)
        arg[0,0]=11.*sqrt(x[0])
        arg[0,1]=1.*sqrt(x[1])
        arg[1,0]=10.*sqrt(x[0])
        arg[1,1]=11.*sqrt(x[1])
        arg[2,0]=20.*sqrt(x[0])
        arg[2,1]=21.*sqrt(x[1])
        arg[3,0]=30.*sqrt(x[0])
        arg[3,1]=31.*sqrt(x[1])
        ref=sqrt(1522.5)

      else:
        arg=Data(0,(4, 3),w)
        arg[0,0]=11.*sqrt(x[0])
        arg[0,1]=1.*sqrt(x[1])
        arg[0,2]=2.*sqrt(x[2])
        arg[1,0]=10.*sqrt(x[0])
        arg[1,1]=11.*sqrt(x[1])
        arg[1,2]=12.*sqrt(x[2])
        arg[2,0]=20.*sqrt(x[0])
        arg[2,1]=21.*sqrt(x[1])
        arg[2,2]=22.*sqrt(x[2])
        arg[3,0]=30.*sqrt(x[0])
        arg[3,1]=31.*sqrt(x[1])
        arg[3,2]=32.*sqrt(x[2])
        ref=sqrt(2350.5)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunction_fromData_rank3(self):
      """
      tests L2-norm of Data on the ReducedFunction

      assumptions: self.domain supports integration on ReducedFunction
      """
      dim=self.domain.getDim()
      w=ReducedFunction(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(6, 2, 2),w)
        arg[0,0,0]=(0.449174971953)*sqrt(x[0])
        arg[0,0,1]=(-0.0109398763289)*sqrt(x[1])
        arg[0,1,0]=(-0.202497187709)*sqrt(x[0])
        arg[0,1,1]=(-0.12970879334)*sqrt(x[1])
        arg[1,0,0]=(-0.138092481719)*sqrt(x[0])
        arg[1,0,1]=(-0.528752200917)*sqrt(x[1])
        arg[1,1,0]=(-0.605919441662)*sqrt(x[0])
        arg[1,1,1]=(0.215615032334)*sqrt(x[1])
        arg[2,0,0]=(-0.998734541972)*sqrt(x[0])
        arg[2,0,1]=(0.725811901251)*sqrt(x[1])
        arg[2,1,0]=(-0.966536503228)*sqrt(x[0])
        arg[2,1,1]=(-0.528692217355)*sqrt(x[1])
        arg[3,0,0]=(0.757633851466)*sqrt(x[0])
        arg[3,0,1]=(-0.524660157377)*sqrt(x[1])
        arg[3,1,0]=(0.983733431677)*sqrt(x[0])
        arg[3,1,1]=(0.061279109546)*sqrt(x[1])
        arg[4,0,0]=(0.85914215305)*sqrt(x[0])
        arg[4,0,1]=(0.941714045112)*sqrt(x[1])
        arg[4,1,0]=(0.172235529555)*sqrt(x[0])
        arg[4,1,1]=(-0.108381454437)*sqrt(x[1])
        arg[5,0,0]=(-0.736373697727)*sqrt(x[0])
        arg[5,0,1]=(-0.599337929679)*sqrt(x[1])
        arg[5,1,0]=(0.661072686392)*sqrt(x[0])
        arg[5,1,1]=(-0.55107327409)*sqrt(x[1])
        ref=sqrt(4.4196214907099591)

      else:
        arg=Data(0,(6, 2, 3),w)
        arg[0,0,0]=(0.69227064904)*sqrt(x[0])
        arg[0,0,1]=(-0.968336177418)*sqrt(x[1])
        arg[0,0,2]=(-0.634883146685)*sqrt(x[2])
        arg[0,1,0]=(-0.12640661422)*sqrt(x[0])
        arg[0,1,1]=(-0.637386589888)*sqrt(x[1])
        arg[0,1,2]=(0.26060859356)*sqrt(x[2])
        arg[1,0,0]=(-0.986864633297)*sqrt(x[0])
        arg[1,0,1]=(-0.441589142379)*sqrt(x[1])
        arg[1,0,2]=(-0.587865539582)*sqrt(x[2])
        arg[1,1,0]=(0.596052465031)*sqrt(x[0])
        arg[1,1,1]=(0.312732336652)*sqrt(x[1])
        arg[1,1,2]=(-0.514423945092)*sqrt(x[2])
        arg[2,0,0]=(-0.892391254794)*sqrt(x[0])
        arg[2,0,1]=(0.377920185756)*sqrt(x[1])
        arg[2,0,2]=(-0.120174597181)*sqrt(x[2])
        arg[2,1,0]=(-0.469951576468)*sqrt(x[0])
        arg[2,1,1]=(-0.788362249555)*sqrt(x[1])
        arg[2,1,2]=(0.745625354986)*sqrt(x[2])
        arg[3,0,0]=(0.542802498569)*sqrt(x[0])
        arg[3,0,1]=(-0.814541028706)*sqrt(x[1])
        arg[3,0,2]=(0.298410992196)*sqrt(x[2])
        arg[3,1,0]=(0.981190341206)*sqrt(x[0])
        arg[3,1,1]=(0.666421298608)*sqrt(x[1])
        arg[3,1,2]=(-0.369751722626)*sqrt(x[2])
        arg[4,0,0]=(-0.75379530597)*sqrt(x[0])
        arg[4,0,1]=(0.283357267139)*sqrt(x[1])
        arg[4,0,2]=(0.247787072861)*sqrt(x[2])
        arg[4,1,0]=(0.301766692533)*sqrt(x[0])
        arg[4,1,1]=(0.828183439224)*sqrt(x[1])
        arg[4,1,2]=(-0.580824060547)*sqrt(x[2])
        arg[5,0,0]=(0.637345610764)*sqrt(x[0])
        arg[5,0,1]=(-0.234409115997)*sqrt(x[1])
        arg[5,0,2]=(-0.192639300316)*sqrt(x[2])
        arg[5,1,0]=(-0.62609237162)*sqrt(x[0])
        arg[5,1,1]=(0.463404958552)*sqrt(x[1])
        arg[5,1,2]=(-0.547814448738)*sqrt(x[2])
        ref=sqrt(6.3571697792950923)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunction_fromData_rank4(self):
      """
      tests L2-norm of Data on the ReducedFunction

      assumptions: self.domain supports integration on ReducedFunction
      """
      dim=self.domain.getDim()
      w=ReducedFunction(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 5, 3, 2),w)
        arg[0,0,0,0]=(-0.232618585183)*sqrt(x[0])
        arg[0,0,0,1]=(0.39796117869)*sqrt(x[1])
        arg[0,0,1,0]=(-0.997336958262)*sqrt(x[0])
        arg[0,0,1,1]=(-0.351780915076)*sqrt(x[1])
        arg[0,0,2,0]=(-0.876764070136)*sqrt(x[0])
        arg[0,0,2,1]=(0.808730805817)*sqrt(x[1])
        arg[0,1,0,0]=(-0.197154744966)*sqrt(x[0])
        arg[0,1,0,1]=(0.416246096086)*sqrt(x[1])
        arg[0,1,1,0]=(0.708038457121)*sqrt(x[0])
        arg[0,1,1,1]=(-0.00954021503188)*sqrt(x[1])
        arg[0,1,2,0]=(-0.62630809425)*sqrt(x[0])
        arg[0,1,2,1]=(0.430228727912)*sqrt(x[1])
        arg[0,2,0,0]=(0.0306704747648)*sqrt(x[0])
        arg[0,2,0,1]=(-0.913877199453)*sqrt(x[1])
        arg[0,2,1,0]=(-0.697612800829)*sqrt(x[0])
        arg[0,2,1,1]=(-0.17996376822)*sqrt(x[1])
        arg[0,2,2,0]=(-0.304509578871)*sqrt(x[0])
        arg[0,2,2,1]=(-0.610556755811)*sqrt(x[1])
        arg[0,3,0,0]=(-0.452355972234)*sqrt(x[0])
        arg[0,3,0,1]=(-0.368921242518)*sqrt(x[1])
        arg[0,3,1,0]=(-0.478275554932)*sqrt(x[0])
        arg[0,3,1,1]=(0.257178549127)*sqrt(x[1])
        arg[0,3,2,0]=(0.530736487177)*sqrt(x[0])
        arg[0,3,2,1]=(-0.567126272463)*sqrt(x[1])
        arg[0,4,0,0]=(0.801519165938)*sqrt(x[0])
        arg[0,4,0,1]=(-0.509816703951)*sqrt(x[1])
        arg[0,4,1,0]=(-0.255412646934)*sqrt(x[0])
        arg[0,4,1,1]=(0.437540101896)*sqrt(x[1])
        arg[0,4,2,0]=(-0.815574969538)*sqrt(x[0])
        arg[0,4,2,1]=(-0.94691547137)*sqrt(x[1])
        arg[1,0,0,0]=(-0.732550722593)*sqrt(x[0])
        arg[1,0,0,1]=(0.515752381704)*sqrt(x[1])
        arg[1,0,1,0]=(-0.343590210899)*sqrt(x[0])
        arg[1,0,1,1]=(-0.0601907964915)*sqrt(x[1])
        arg[1,0,2,0]=(0.0199916154421)*sqrt(x[0])
        arg[1,0,2,1]=(-0.136927227821)*sqrt(x[1])
        arg[1,1,0,0]=(0.397994441702)*sqrt(x[0])
        arg[1,1,0,1]=(0.953873148948)*sqrt(x[1])
        arg[1,1,1,0]=(0.419416235967)*sqrt(x[0])
        arg[1,1,1,1]=(0.700998577193)*sqrt(x[1])
        arg[1,1,2,0]=(-0.497358799271)*sqrt(x[0])
        arg[1,1,2,1]=(0.0851768858379)*sqrt(x[1])
        arg[1,2,0,0]=(0.0936678875202)*sqrt(x[0])
        arg[1,2,0,1]=(0.869883786896)*sqrt(x[1])
        arg[1,2,1,0]=(0.582700123485)*sqrt(x[0])
        arg[1,2,1,1]=(-0.433381106794)*sqrt(x[1])
        arg[1,2,2,0]=(-0.527031777974)*sqrt(x[0])
        arg[1,2,2,1]=(0.105105137652)*sqrt(x[1])
        arg[1,3,0,0]=(-0.716750829134)*sqrt(x[0])
        arg[1,3,0,1]=(0.774519209008)*sqrt(x[1])
        arg[1,3,1,0]=(-0.568743372716)*sqrt(x[0])
        arg[1,3,1,1]=(0.794732483944)*sqrt(x[1])
        arg[1,3,2,0]=(0.246606002015)*sqrt(x[0])
        arg[1,3,2,1]=(-0.988869494994)*sqrt(x[1])
        arg[1,4,0,0]=(0.482379298083)*sqrt(x[0])
        arg[1,4,0,1]=(-0.386268387903)*sqrt(x[1])
        arg[1,4,1,0]=(0.137184889675)*sqrt(x[0])
        arg[1,4,1,1]=(-0.140520035321)*sqrt(x[1])
        arg[1,4,2,0]=(0.822755050415)*sqrt(x[0])
        arg[1,4,2,1]=(-0.815562139522)*sqrt(x[1])
        arg[2,0,0,0]=(-0.462891511962)*sqrt(x[0])
        arg[2,0,0,1]=(-0.122643411631)*sqrt(x[1])
        arg[2,0,1,0]=(-0.520861119962)*sqrt(x[0])
        arg[2,0,1,1]=(-0.881189618018)*sqrt(x[1])
        arg[2,0,2,0]=(-0.776157842774)*sqrt(x[0])
        arg[2,0,2,1]=(-0.12354053207)*sqrt(x[1])
        arg[2,1,0,0]=(0.395495230826)*sqrt(x[0])
        arg[2,1,0,1]=(-0.388106659423)*sqrt(x[1])
        arg[2,1,1,0]=(0.354250242834)*sqrt(x[0])
        arg[2,1,1,1]=(-0.666514210192)*sqrt(x[1])
        arg[2,1,2,0]=(0.951294655083)*sqrt(x[0])
        arg[2,1,2,1]=(0.074024416386)*sqrt(x[1])
        arg[2,2,0,0]=(0.335448485459)*sqrt(x[0])
        arg[2,2,0,1]=(-0.40988282528)*sqrt(x[1])
        arg[2,2,1,0]=(-0.805725968875)*sqrt(x[0])
        arg[2,2,1,1]=(-0.949883082118)*sqrt(x[1])
        arg[2,2,2,0]=(0.531549210683)*sqrt(x[0])
        arg[2,2,2,1]=(-0.398401016682)*sqrt(x[1])
        arg[2,3,0,0]=(-0.953963433205)*sqrt(x[0])
        arg[2,3,0,1]=(0.643431126406)*sqrt(x[1])
        arg[2,3,1,0]=(-0.167611998738)*sqrt(x[0])
        arg[2,3,1,1]=(0.226130056552)*sqrt(x[1])
        arg[2,3,2,0]=(0.0752687641131)*sqrt(x[0])
        arg[2,3,2,1]=(-0.115742756362)*sqrt(x[1])
        arg[2,4,0,0]=(0.579694491028)*sqrt(x[0])
        arg[2,4,0,1]=(-0.112005738299)*sqrt(x[1])
        arg[2,4,1,0]=(0.657291764224)*sqrt(x[0])
        arg[2,4,1,1]=(0.62671154177)*sqrt(x[1])
        arg[2,4,2,0]=(0.103695027944)*sqrt(x[0])
        arg[2,4,2,1]=(0.462828491544)*sqrt(x[1])
        arg[3,0,0,0]=(0.697692979998)*sqrt(x[0])
        arg[3,0,0,1]=(-0.123481859619)*sqrt(x[1])
        arg[3,0,1,0]=(-0.749745629459)*sqrt(x[0])
        arg[3,0,1,1]=(-0.541969524069)*sqrt(x[1])
        arg[3,0,2,0]=(0.819484470759)*sqrt(x[0])
        arg[3,0,2,1]=(-0.860592326469)*sqrt(x[1])
        arg[3,1,0,0]=(-0.716566084771)*sqrt(x[0])
        arg[3,1,0,1]=(-0.949235434827)*sqrt(x[1])
        arg[3,1,1,0]=(-0.826699498174)*sqrt(x[0])
        arg[3,1,1,1]=(-0.138511521583)*sqrt(x[1])
        arg[3,1,2,0]=(-0.951682890904)*sqrt(x[0])
        arg[3,1,2,1]=(0.413293316925)*sqrt(x[1])
        arg[3,2,0,0]=(0.909516836775)*sqrt(x[0])
        arg[3,2,0,1]=(-0.919989721277)*sqrt(x[1])
        arg[3,2,1,0]=(0.0994860369337)*sqrt(x[0])
        arg[3,2,1,1]=(-0.933647246623)*sqrt(x[1])
        arg[3,2,2,0]=(-0.759215183015)*sqrt(x[0])
        arg[3,2,2,1]=(0.0975793309286)*sqrt(x[1])
        arg[3,3,0,0]=(-0.130256739381)*sqrt(x[0])
        arg[3,3,0,1]=(-0.582280862311)*sqrt(x[1])
        arg[3,3,1,0]=(0.206970526192)*sqrt(x[0])
        arg[3,3,1,1]=(-0.8678322258)*sqrt(x[1])
        arg[3,3,2,0]=(0.133004501279)*sqrt(x[0])
        arg[3,3,2,1]=(0.802921710935)*sqrt(x[1])
        arg[3,4,0,0]=(-0.255737792764)*sqrt(x[0])
        arg[3,4,0,1]=(-0.34168114937)*sqrt(x[1])
        arg[3,4,1,0]=(-0.859309090399)*sqrt(x[0])
        arg[3,4,1,1]=(0.245043986435)*sqrt(x[1])
        arg[3,4,2,0]=(0.893062018695)*sqrt(x[0])
        arg[3,4,2,1]=(0.709422742588)*sqrt(x[1])
        ref=sqrt(20.5933920543)
      else:
        arg=Data(0,(4, 5, 3, 3),w)
        arg[0,0,0,0]=(0.0312828390439)*sqrt(x[0])
        arg[0,0,0,1]=(-0.524970416212)*sqrt(x[1])
        arg[0,0,0,2]=(0.561865217554)*sqrt(x[2])
        arg[0,0,1,0]=(0.692457187384)*sqrt(x[0])
        arg[0,0,1,1]=(0.946967182157)*sqrt(x[1])
        arg[0,0,1,2]=(-0.863842279464)*sqrt(x[2])
        arg[0,0,2,0]=(0.993922921598)*sqrt(x[0])
        arg[0,0,2,1]=(0.322812768679)*sqrt(x[1])
        arg[0,0,2,2]=(0.901876132204)*sqrt(x[2])
        arg[0,1,0,0]=(0.967569979365)*sqrt(x[0])
        arg[0,1,0,1]=(0.840979131355)*sqrt(x[1])
        arg[0,1,0,2]=(0.0494811460856)*sqrt(x[2])
        arg[0,1,1,0]=(0.315178456102)*sqrt(x[0])
        arg[0,1,1,1]=(0.449848313024)*sqrt(x[1])
        arg[0,1,1,2]=(0.765887852886)*sqrt(x[2])
        arg[0,1,2,0]=(0.975541574352)*sqrt(x[0])
        arg[0,1,2,1]=(-0.797851290751)*sqrt(x[1])
        arg[0,1,2,2]=(0.628918775319)*sqrt(x[2])
        arg[0,2,0,0]=(0.685635794312)*sqrt(x[0])
        arg[0,2,0,1]=(0.10341799962)*sqrt(x[1])
        arg[0,2,0,2]=(-0.964822756043)*sqrt(x[2])
        arg[0,2,1,0]=(-0.56160368212)*sqrt(x[0])
        arg[0,2,1,1]=(0.676344298102)*sqrt(x[1])
        arg[0,2,1,2]=(-0.713924121843)*sqrt(x[2])
        arg[0,2,2,0]=(-0.276655136263)*sqrt(x[0])
        arg[0,2,2,1]=(0.336046973788)*sqrt(x[1])
        arg[0,2,2,2]=(-0.68789392396)*sqrt(x[2])
        arg[0,3,0,0]=(0.0172861311571)*sqrt(x[0])
        arg[0,3,0,1]=(-0.301075956456)*sqrt(x[1])
        arg[0,3,0,2]=(0.779442985415)*sqrt(x[2])
        arg[0,3,1,0]=(-0.517629576558)*sqrt(x[0])
        arg[0,3,1,1]=(0.584779586639)*sqrt(x[1])
        arg[0,3,1,2]=(-0.53266435436)*sqrt(x[2])
        arg[0,3,2,0]=(0.841533567102)*sqrt(x[0])
        arg[0,3,2,1]=(0.0458746415489)*sqrt(x[1])
        arg[0,3,2,2]=(0.921237870758)*sqrt(x[2])
        arg[0,4,0,0]=(0.0548343238805)*sqrt(x[0])
        arg[0,4,0,1]=(0.687022707412)*sqrt(x[1])
        arg[0,4,0,2]=(-0.319803609795)*sqrt(x[2])
        arg[0,4,1,0]=(0.409763007811)*sqrt(x[0])
        arg[0,4,1,1]=(0.165501957435)*sqrt(x[1])
        arg[0,4,1,2]=(0.116001692781)*sqrt(x[2])
        arg[0,4,2,0]=(-0.515571394238)*sqrt(x[0])
        arg[0,4,2,1]=(0.209467945147)*sqrt(x[1])
        arg[0,4,2,2]=(-0.344827191247)*sqrt(x[2])
        arg[1,0,0,0]=(0.57193838014)*sqrt(x[0])
        arg[1,0,0,1]=(-0.0880683799076)*sqrt(x[1])
        arg[1,0,0,2]=(0.956899617441)*sqrt(x[2])
        arg[1,0,1,0]=(-0.783689636357)*sqrt(x[0])
        arg[1,0,1,1]=(-0.25177506885)*sqrt(x[1])
        arg[1,0,1,2]=(-0.97074584634)*sqrt(x[2])
        arg[1,0,2,0]=(0.432543519806)*sqrt(x[0])
        arg[1,0,2,1]=(0.481003021954)*sqrt(x[1])
        arg[1,0,2,2]=(-0.0630751518268)*sqrt(x[2])
        arg[1,1,0,0]=(-0.65152446796)*sqrt(x[0])
        arg[1,1,0,1]=(-0.0323685084425)*sqrt(x[1])
        arg[1,1,0,2]=(-0.508674033909)*sqrt(x[2])
        arg[1,1,1,0]=(-0.533367818916)*sqrt(x[0])
        arg[1,1,1,1]=(0.310738340288)*sqrt(x[1])
        arg[1,1,1,2]=(0.694612234326)*sqrt(x[2])
        arg[1,1,2,0]=(-0.622052473032)*sqrt(x[0])
        arg[1,1,2,1]=(0.0498443793671)*sqrt(x[1])
        arg[1,1,2,2]=(0.61023707512)*sqrt(x[2])
        arg[1,2,0,0]=(0.0730267406859)*sqrt(x[0])
        arg[1,2,0,1]=(0.146909334607)*sqrt(x[1])
        arg[1,2,0,2]=(-0.641860284448)*sqrt(x[2])
        arg[1,2,1,0]=(0.917976589737)*sqrt(x[0])
        arg[1,2,1,1]=(0.50219672122)*sqrt(x[1])
        arg[1,2,1,2]=(0.634559579812)*sqrt(x[2])
        arg[1,2,2,0]=(0.0578772734534)*sqrt(x[0])
        arg[1,2,2,1]=(0.288730973517)*sqrt(x[1])
        arg[1,2,2,2]=(-0.0525978796154)*sqrt(x[2])
        arg[1,3,0,0]=(-0.926152433388)*sqrt(x[0])
        arg[1,3,0,1]=(0.0616647680855)*sqrt(x[1])
        arg[1,3,0,2]=(-0.875889217846)*sqrt(x[2])
        arg[1,3,1,0]=(-0.638931542845)*sqrt(x[0])
        arg[1,3,1,1]=(0.708848122964)*sqrt(x[1])
        arg[1,3,1,2]=(0.119066979792)*sqrt(x[2])
        arg[1,3,2,0]=(0.853716218591)*sqrt(x[0])
        arg[1,3,2,1]=(-0.92754322201)*sqrt(x[1])
        arg[1,3,2,2]=(-0.671530626265)*sqrt(x[2])
        arg[1,4,0,0]=(0.337424536231)*sqrt(x[0])
        arg[1,4,0,1]=(0.335704451719)*sqrt(x[1])
        arg[1,4,0,2]=(-0.484565969466)*sqrt(x[2])
        arg[1,4,1,0]=(-0.855476192012)*sqrt(x[0])
        arg[1,4,1,1]=(0.405674615553)*sqrt(x[1])
        arg[1,4,1,2]=(0.728310771323)*sqrt(x[2])
        arg[1,4,2,0]=(0.363651308265)*sqrt(x[0])
        arg[1,4,2,1]=(0.174460594531)*sqrt(x[1])
        arg[1,4,2,2]=(-0.0418244838617)*sqrt(x[2])
        arg[2,0,0,0]=(-0.531341992511)*sqrt(x[0])
        arg[2,0,0,1]=(0.584996796272)*sqrt(x[1])
        arg[2,0,0,2]=(-0.752430968716)*sqrt(x[2])
        arg[2,0,1,0]=(-0.341989849747)*sqrt(x[0])
        arg[2,0,1,1]=(0.153572646953)*sqrt(x[1])
        arg[2,0,1,2]=(-0.197130051737)*sqrt(x[2])
        arg[2,0,2,0]=(-0.338082424082)*sqrt(x[0])
        arg[2,0,2,1]=(0.000173657394772)*sqrt(x[1])
        arg[2,0,2,2]=(0.365272907692)*sqrt(x[2])
        arg[2,1,0,0]=(0.904304126564)*sqrt(x[0])
        arg[2,1,0,1]=(0.161252368484)*sqrt(x[1])
        arg[2,1,0,2]=(0.246854092422)*sqrt(x[2])
        arg[2,1,1,0]=(-0.299880647529)*sqrt(x[0])
        arg[2,1,1,1]=(-0.566917528608)*sqrt(x[1])
        arg[2,1,1,2]=(0.243183337285)*sqrt(x[2])
        arg[2,1,2,0]=(0.437406011474)*sqrt(x[0])
        arg[2,1,2,1]=(0.727447394053)*sqrt(x[1])
        arg[2,1,2,2]=(0.380752950664)*sqrt(x[2])
        arg[2,2,0,0]=(0.172292846911)*sqrt(x[0])
        arg[2,2,0,1]=(0.334201791643)*sqrt(x[1])
        arg[2,2,0,2]=(0.739989926962)*sqrt(x[2])
        arg[2,2,1,0]=(-0.0669843715042)*sqrt(x[0])
        arg[2,2,1,1]=(-0.540497281635)*sqrt(x[1])
        arg[2,2,1,2]=(-0.744217027088)*sqrt(x[2])
        arg[2,2,2,0]=(-0.287295952259)*sqrt(x[0])
        arg[2,2,2,1]=(-0.512411849183)*sqrt(x[1])
        arg[2,2,2,2]=(0.953107417666)*sqrt(x[2])
        arg[2,3,0,0]=(0.998168116695)*sqrt(x[0])
        arg[2,3,0,1]=(0.960065646359)*sqrt(x[1])
        arg[2,3,0,2]=(0.110048258832)*sqrt(x[2])
        arg[2,3,1,0]=(-0.477271134724)*sqrt(x[0])
        arg[2,3,1,1]=(0.707182612251)*sqrt(x[1])
        arg[2,3,1,2]=(0.285500891755)*sqrt(x[2])
        arg[2,3,2,0]=(-0.863497506661)*sqrt(x[0])
        arg[2,3,2,1]=(-0.293917669879)*sqrt(x[1])
        arg[2,3,2,2]=(-0.403384244295)*sqrt(x[2])
        arg[2,4,0,0]=(0.848455277702)*sqrt(x[0])
        arg[2,4,0,1]=(-0.530101455578)*sqrt(x[1])
        arg[2,4,0,2]=(0.33887313048)*sqrt(x[2])
        arg[2,4,1,0]=(-0.195313538124)*sqrt(x[0])
        arg[2,4,1,1]=(-0.62754572008)*sqrt(x[1])
        arg[2,4,1,2]=(-0.385132960582)*sqrt(x[2])
        arg[2,4,2,0]=(0.240048012886)*sqrt(x[0])
        arg[2,4,2,1]=(0.900766252969)*sqrt(x[1])
        arg[2,4,2,2]=(0.669620533505)*sqrt(x[2])
        arg[3,0,0,0]=(0.375766827301)*sqrt(x[0])
        arg[3,0,0,1]=(0.705484960308)*sqrt(x[1])
        arg[3,0,0,2]=(0.440931516034)*sqrt(x[2])
        arg[3,0,1,0]=(-0.44724403177)*sqrt(x[0])
        arg[3,0,1,1]=(-0.31558249626)*sqrt(x[1])
        arg[3,0,1,2]=(-0.00419436365172)*sqrt(x[2])
        arg[3,0,2,0]=(0.750599752032)*sqrt(x[0])
        arg[3,0,2,1]=(0.367649951795)*sqrt(x[1])
        arg[3,0,2,2]=(0.0488013073654)*sqrt(x[2])
        arg[3,1,0,0]=(-0.992890068274)*sqrt(x[0])
        arg[3,1,0,1]=(0.671447745511)*sqrt(x[1])
        arg[3,1,0,2]=(0.85613331404)*sqrt(x[2])
        arg[3,1,1,0]=(-0.46064764242)*sqrt(x[0])
        arg[3,1,1,1]=(0.48138877715)*sqrt(x[1])
        arg[3,1,1,2]=(0.396741761803)*sqrt(x[2])
        arg[3,1,2,0]=(-0.879391967543)*sqrt(x[0])
        arg[3,1,2,1]=(-0.44039462138)*sqrt(x[1])
        arg[3,1,2,2]=(0.0330511573872)*sqrt(x[2])
        arg[3,2,0,0]=(-0.367413701648)*sqrt(x[0])
        arg[3,2,0,1]=(0.0359818324891)*sqrt(x[1])
        arg[3,2,0,2]=(-0.307532667032)*sqrt(x[2])
        arg[3,2,1,0]=(0.334663597166)*sqrt(x[0])
        arg[3,2,1,1]=(0.541941978066)*sqrt(x[1])
        arg[3,2,1,2]=(-0.609184079318)*sqrt(x[2])
        arg[3,2,2,0]=(0.359349239826)*sqrt(x[0])
        arg[3,2,2,1]=(0.0419272305685)*sqrt(x[1])
        arg[3,2,2,2]=(0.557189794296)*sqrt(x[2])
        arg[3,3,0,0]=(-0.85864165554)*sqrt(x[0])
        arg[3,3,0,1]=(-0.185411404213)*sqrt(x[1])
        arg[3,3,0,2]=(0.254294865253)*sqrt(x[2])
        arg[3,3,1,0]=(0.870362177541)*sqrt(x[0])
        arg[3,3,1,1]=(-0.439688612864)*sqrt(x[1])
        arg[3,3,1,2]=(0.26006729357)*sqrt(x[2])
        arg[3,3,2,0]=(-0.0724034754175)*sqrt(x[0])
        arg[3,3,2,1]=(0.444871564246)*sqrt(x[1])
        arg[3,3,2,2]=(0.485634530531)*sqrt(x[2])
        arg[3,4,0,0]=(-0.744756961758)*sqrt(x[0])
        arg[3,4,0,1]=(0.429761406102)*sqrt(x[1])
        arg[3,4,0,2]=(-0.584963735834)*sqrt(x[2])
        arg[3,4,1,0]=(0.684578379159)*sqrt(x[0])
        arg[3,4,1,1]=(0.949460132601)*sqrt(x[1])
        arg[3,4,1,2]=(-0.592179909559)*sqrt(x[2])
        arg[3,4,2,0]=(0.707154437797)*sqrt(x[0])
        arg[3,4,2,1]=(0.619200407063)*sqrt(x[1])
        arg[3,4,2,2]=(-0.338547165)*sqrt(x[2])
        ref=sqrt(28.8255957718)
      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunctionOnBoundary_fromData_rank0(self):
      """
      tests L2-norm of Data on the FunctionOnBoundary

      assumptions: self.domain supports integration on FunctionOnBoundary
      """
      dim=self.domain.getDim()
      w=FunctionOnBoundary(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(),w)
        arg=(-0.245574919477)*x[0]
        ref=sqrt((0.0603070410759)*(2.*dim+1.)/3.)

      else:
        arg=Data(0,(),w)
        arg=(0.757324521515)*x[0]
        ref=sqrt((0.573540430888)*(2.*dim+1.)/3.)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunctionOnBoundary_fromData_rank1(self):
      """
      tests L2-norm of Data on the FunctionOnBoundary

      assumptions: self.domain supports integration on FunctionOnBoundary
      """
      dim=self.domain.getDim()
      w=FunctionOnBoundary(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(2,),w)
        arg[0]=(0.723421565407)*x[0]
        arg[1]=(-0.460477393103)*x[1]
        ref=sqrt((0.735378190855)*(2.*dim+1.)/3.)

      else:
        arg=Data(0,(3,),w)
        arg[0]=(-0.88528497163)*x[0]
        arg[1]=(-0.65510214636)*x[1]
        arg[2]=(0.399538866363)*x[2]
        ref=sqrt((1.37251960889)*(2.*dim+1.)/3.)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunctionOnBoundary_fromData_rank2(self):
      """
      tests L2-norm of Data on the FunctionOnBoundary

      assumptions: self.domain supports integration on FunctionOnBoundary
      """
      dim=self.domain.getDim()
      w=FunctionOnBoundary(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 2),w)
        arg[0,0]=(0.955620993904)*x[0]
        arg[0,1]=(-0.0987865813703)*x[1]
        arg[1,0]=(0.0288267231531)*x[0]
        arg[1,1]=(0.655440599879)*x[1]
        arg[2,0]=(0.685627284533)*x[0]
        arg[2,1]=(-0.989832824892)*x[1]
        arg[3,0]=(0.292184093194)*x[0]
        arg[3,1]=(0.149553857773)*x[1]
        ref=sqrt((2.91099532781)*(2.*dim+1.)/3.)

      else:
        arg=Data(0,(4, 3),w)
        arg[0,0]=(-0.325908541533)*x[0]
        arg[0,1]=(-0.992480479749)*x[1]
        arg[0,2]=(0.660360271799)*x[2]
        arg[1,0]=(0.173485908581)*x[0]
        arg[1,1]=(-0.328755199781)*x[1]
        arg[1,2]=(-0.943354674948)*x[2]
        arg[2,0]=(0.680713222646)*x[0]
        arg[2,1]=(-0.765971835693)*x[1]
        arg[2,2]=(0.0413284847528)*x[2]
        arg[3,0]=(0.990074004708)*x[0]
        arg[3,1]=(0.941801786766)*x[1]
        arg[3,2]=(0.886926192201)*x[2]
        ref=sqrt((6.26107155228)*(2.*dim+1.)/3.)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunctionOnBoundary_fromData_rank3(self):
      """
      tests L2-norm of Data on the FunctionOnBoundary

      assumptions: self.domain supports integration on FunctionOnBoundary
      """
      dim=self.domain.getDim()
      w=FunctionOnBoundary(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(6, 2, 2),w)
        arg[0,0,0]=(-0.0781611551598)*x[0]
        arg[0,0,1]=(0.291016249575)*x[1]
        arg[0,1,0]=(-0.107555233086)*x[0]
        arg[0,1,1]=(-0.559067108546)*x[1]
        arg[1,0,0]=(-0.0818406701266)*x[0]
        arg[1,0,1]=(-0.594866806483)*x[1]
        arg[1,1,0]=(-0.725814803863)*x[0]
        arg[1,1,1]=(0.59128101992)*x[1]
        arg[2,0,0]=(-0.15381555291)*x[0]
        arg[2,0,1]=(-0.679882948503)*x[1]
        arg[2,1,0]=(-0.58437193917)*x[0]
        arg[2,1,1]=(0.136304615849)*x[1]
        arg[3,0,0]=(0.0671365410096)*x[0]
        arg[3,0,1]=(-0.645687212187)*x[1]
        arg[3,1,0]=(-0.642492412392)*x[0]
        arg[3,1,1]=(-0.125760054735)*x[1]
        arg[4,0,0]=(0.731110824794)*x[0]
        arg[4,0,1]=(0.491668422979)*x[1]
        arg[4,1,0]=(-0.775841478292)*x[0]
        arg[4,1,1]=(0.728265567974)*x[1]
        arg[5,0,0]=(0.84511832373)*x[0]
        arg[5,0,1]=(-0.513796801068)*x[1]
        arg[5,1,0]=(0.113072243554)*x[0]
        arg[5,1,1]=(0.246630838744)*x[1]
        ref=sqrt((6.30829536252)*(2.*dim+1.)/3.)

      else:
        arg=Data(0,(6, 2, 3),w)
        arg[0,0,0]=(0.369748116859)*x[0]
        arg[0,0,1]=(-0.758056560031)*x[1]
        arg[0,0,2]=(-0.873984709951)*x[2]
        arg[0,1,0]=(0.311680165784)*x[0]
        arg[0,1,1]=(0.374400673651)*x[1]
        arg[0,1,2]=(0.712484217076)*x[2]
        arg[1,0,0]=(0.829379714484)*x[0]
        arg[1,0,1]=(-0.0551589596149)*x[1]
        arg[1,0,2]=(0.965672208426)*x[2]
        arg[1,1,0]=(-0.205044281547)*x[0]
        arg[1,1,1]=(0.238197452756)*x[1]
        arg[1,1,2]=(-0.33456139292)*x[2]
        arg[2,0,0]=(0.649928288926)*x[0]
        arg[2,0,1]=(-0.661384953389)*x[1]
        arg[2,0,2]=(-0.253241222975)*x[2]
        arg[2,1,0]=(-0.491716575992)*x[0]
        arg[2,1,1]=(-0.970872527468)*x[1]
        arg[2,1,2]=(0.222410198921)*x[2]
        arg[3,0,0]=(0.205752630262)*x[0]
        arg[3,0,1]=(0.864804362697)*x[1]
        arg[3,0,2]=(-0.417975564033)*x[2]
        arg[3,1,0]=(0.586425694033)*x[0]
        arg[3,1,1]=(0.952661122184)*x[1]
        arg[3,1,2]=(0.608680080453)*x[2]
        arg[4,0,0]=(0.625968903369)*x[0]
        arg[4,0,1]=(-0.573909405003)*x[1]
        arg[4,0,2]=(-0.762256394595)*x[2]
        arg[4,1,0]=(0.0710742394418)*x[0]
        arg[4,1,1]=(0.583378040574)*x[1]
        arg[4,1,2]=(0.719032893115)*x[2]
        arg[5,0,0]=(0.032173368884)*x[0]
        arg[5,0,1]=(-0.434042549492)*x[1]
        arg[5,0,2]=(0.363504765447)*x[2]
        arg[5,1,0]=(0.598817469198)*x[0]
        arg[5,1,1]=(-0.163967008775)*x[1]
        arg[5,1,2]=(0.546778730604)*x[2]
        ref=sqrt((11.9696343123)*(2.*dim+1.)/3.)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onFunctionOnBoundary_fromData_rank4(self):
      """
      tests L2-norm of Data on the FunctionOnBoundary

      assumptions: self.domain supports integration on FunctionOnBoundary
      """
      dim=self.domain.getDim()
      w=FunctionOnBoundary(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 5, 3, 2),w)
        arg[0,0,0,0]=(-0.273446520069)*x[0]
        arg[0,0,0,1]=(-0.913305910831)*x[1]
        arg[0,0,1,0]=(0.0745566747537)*x[0]
        arg[0,0,1,1]=(0.98803601919)*x[1]
        arg[0,0,2,0]=(-0.244120818875)*x[0]
        arg[0,0,2,1]=(0.247509644998)*x[1]
        arg[0,1,0,0]=(-0.548756417777)*x[0]
        arg[0,1,0,1]=(-0.354587911923)*x[1]
        arg[0,1,1,0]=(0.104268198867)*x[0]
        arg[0,1,1,1]=(-0.541700877072)*x[1]
        arg[0,1,2,0]=(-0.25900060689)*x[0]
        arg[0,1,2,1]=(-0.859660175231)*x[1]
        arg[0,2,0,0]=(0.195235560321)*x[0]
        arg[0,2,0,1]=(0.175518738589)*x[1]
        arg[0,2,1,0]=(-0.0638854232272)*x[0]
        arg[0,2,1,1]=(-0.586161016541)*x[1]
        arg[0,2,2,0]=(0.580258892247)*x[0]
        arg[0,2,2,1]=(-0.931927145435)*x[1]
        arg[0,3,0,0]=(-0.408209600298)*x[0]
        arg[0,3,0,1]=(-0.0344882667014)*x[1]
        arg[0,3,1,0]=(-0.131763534163)*x[0]
        arg[0,3,1,1]=(-0.787653739965)*x[1]
        arg[0,3,2,0]=(0.0910808104711)*x[0]
        arg[0,3,2,1]=(-0.280409023611)*x[1]
        arg[0,4,0,0]=(0.97745754012)*x[0]
        arg[0,4,0,1]=(-0.59829020936)*x[1]
        arg[0,4,1,0]=(0.890520260543)*x[0]
        arg[0,4,1,1]=(-0.600760090231)*x[1]
        arg[0,4,2,0]=(-0.897992297974)*x[0]
        arg[0,4,2,1]=(-0.841169898923)*x[1]
        arg[1,0,0,0]=(-0.249868816409)*x[0]
        arg[1,0,0,1]=(0.620375082228)*x[1]
        arg[1,0,1,0]=(-0.660480789306)*x[0]
        arg[1,0,1,1]=(-0.73638571806)*x[1]
        arg[1,0,2,0]=(0.339987316643)*x[0]
        arg[1,0,2,1]=(0.541112529894)*x[1]
        arg[1,1,0,0]=(-0.468808186705)*x[0]
        arg[1,1,0,1]=(-0.32919679792)*x[1]
        arg[1,1,1,0]=(0.917292803419)*x[0]
        arg[1,1,1,1]=(-0.834265058005)*x[1]
        arg[1,1,2,0]=(-0.247536849264)*x[0]
        arg[1,1,2,1]=(-0.197503469238)*x[1]
        arg[1,2,0,0]=(0.897591919909)*x[0]
        arg[1,2,0,1]=(-0.807446231234)*x[1]
        arg[1,2,1,0]=(-0.369878499382)*x[0]
        arg[1,2,1,1]=(0.985678692179)*x[1]
        arg[1,2,2,0]=(-0.709976427525)*x[0]
        arg[1,2,2,1]=(-0.368744647016)*x[1]
        arg[1,3,0,0]=(0.299630726462)*x[0]
        arg[1,3,0,1]=(-0.445295757899)*x[1]
        arg[1,3,1,0]=(-0.922386577254)*x[0]
        arg[1,3,1,1]=(0.234794853697)*x[1]
        arg[1,3,2,0]=(0.953377720197)*x[0]
        arg[1,3,2,1]=(0.409778183998)*x[1]
        arg[1,4,0,0]=(0.271967945488)*x[0]
        arg[1,4,0,1]=(-0.578629001202)*x[1]
        arg[1,4,1,0]=(0.210755226769)*x[0]
        arg[1,4,1,1]=(-0.0902751419945)*x[1]
        arg[1,4,2,0]=(0.70033387381)*x[0]
        arg[1,4,2,1]=(0.305733565661)*x[1]
        arg[2,0,0,0]=(-0.662483167298)*x[0]
        arg[2,0,0,1]=(0.585252048652)*x[1]
        arg[2,0,1,0]=(-0.398813785959)*x[0]
        arg[2,0,1,1]=(-0.797438697186)*x[1]
        arg[2,0,2,0]=(-0.508308971009)*x[0]
        arg[2,0,2,1]=(0.302249407524)*x[1]
        arg[2,1,0,0]=(0.208644491879)*x[0]
        arg[2,1,0,1]=(-0.604749055374)*x[1]
        arg[2,1,1,0]=(0.641654284594)*x[0]
        arg[2,1,1,1]=(0.456898593356)*x[1]
        arg[2,1,2,0]=(-0.398043867778)*x[0]
        arg[2,1,2,1]=(-0.0712344657587)*x[1]
        arg[2,2,0,0]=(0.0967860954865)*x[0]
        arg[2,2,0,1]=(0.520449905952)*x[1]
        arg[2,2,1,0]=(0.770198029595)*x[0]
        arg[2,2,1,1]=(-0.594004671621)*x[1]
        arg[2,2,2,0]=(-0.744571452885)*x[0]
        arg[2,2,2,1]=(0.544447367825)*x[1]
        arg[2,3,0,0]=(-0.137087966968)*x[0]
        arg[2,3,0,1]=(0.120672667497)*x[1]
        arg[2,3,1,0]=(0.204800088057)*x[0]
        arg[2,3,1,1]=(0.626526346076)*x[1]
        arg[2,3,2,0]=(-0.696480393227)*x[0]
        arg[2,3,2,1]=(0.188533741996)*x[1]
        arg[2,4,0,0]=(-0.403523821067)*x[0]
        arg[2,4,0,1]=(-0.428048989483)*x[1]
        arg[2,4,1,0]=(-0.244186366584)*x[0]
        arg[2,4,1,1]=(0.00866444909003)*x[1]
        arg[2,4,2,0]=(-0.445991308853)*x[0]
        arg[2,4,2,1]=(-0.899951068935)*x[1]
        arg[3,0,0,0]=(0.609340085418)*x[0]
        arg[3,0,0,1]=(0.878750391425)*x[1]
        arg[3,0,1,0]=(0.258064654464)*x[0]
        arg[3,0,1,1]=(-0.482402612985)*x[1]
        arg[3,0,2,0]=(0.943732283389)*x[0]
        arg[3,0,2,1]=(0.65514211843)*x[1]
        arg[3,1,0,0]=(-0.894551979619)*x[0]
        arg[3,1,0,1]=(0.220116541042)*x[1]
        arg[3,1,1,0]=(0.386887699577)*x[0]
        arg[3,1,1,1]=(-0.422560075108)*x[1]
        arg[3,1,2,0]=(0.00387273783493)*x[0]
        arg[3,1,2,1]=(0.465673505613)*x[1]
        arg[3,2,0,0]=(0.987383428982)*x[0]
        arg[3,2,0,1]=(0.376320055964)*x[1]
        arg[3,2,1,0]=(-0.463778689128)*x[0]
        arg[3,2,1,1]=(0.179816566227)*x[1]
        arg[3,2,2,0]=(0.961522856801)*x[0]
        arg[3,2,2,1]=(-0.257779627946)*x[1]
        arg[3,3,0,0]=(0.748886531458)*x[0]
        arg[3,3,0,1]=(-0.257342282566)*x[1]
        arg[3,3,1,0]=(0.377494024401)*x[0]
        arg[3,3,1,1]=(-0.334588346017)*x[1]
        arg[3,3,2,0]=(0.502495149189)*x[0]
        arg[3,3,2,1]=(0.534612429702)*x[1]
        arg[3,4,0,0]=(-0.308551337355)*x[0]
        arg[3,4,0,1]=(-0.471825826745)*x[1]
        arg[3,4,1,0]=(-0.262606531584)*x[0]
        arg[3,4,1,1]=(0.766089616367)*x[1]
        arg[3,4,2,0]=(-0.136526755642)*x[0]
        arg[3,4,2,1]=(0.675111459363)*x[1]
        ref=sqrt((37.3453550914)*(2.*dim+1.)/3.)

      else:
        arg=Data(0,(4, 5, 3, 3),w)
        arg[0,0,0,0]=(0.532946231146)*x[0]
        arg[0,0,0,1]=(0.269364089513)*x[1]
        arg[0,0,0,2]=(-0.207412457081)*x[2]
        arg[0,0,1,0]=(-0.843104704858)*x[0]
        arg[0,0,1,1]=(0.0416216508473)*x[1]
        arg[0,0,1,2]=(-0.836074693662)*x[2]
        arg[0,0,2,0]=(0.943609268731)*x[0]
        arg[0,0,2,1]=(0.0154543737816)*x[1]
        arg[0,0,2,2]=(0.0726268788381)*x[2]
        arg[0,1,0,0]=(0.108422740078)*x[0]
        arg[0,1,0,1]=(-0.296667916638)*x[1]
        arg[0,1,0,2]=(-0.769732600535)*x[2]
        arg[0,1,1,0]=(-0.428575493834)*x[0]
        arg[0,1,1,1]=(0.421245456722)*x[1]
        arg[0,1,1,2]=(-0.588277652277)*x[2]
        arg[0,1,2,0]=(0.145294576795)*x[0]
        arg[0,1,2,1]=(0.323206623794)*x[1]
        arg[0,1,2,2]=(0.788115602892)*x[2]
        arg[0,2,0,0]=(-0.227877282292)*x[0]
        arg[0,2,0,1]=(-0.630647460719)*x[1]
        arg[0,2,0,2]=(0.58754882135)*x[2]
        arg[0,2,1,0]=(0.347191113403)*x[0]
        arg[0,2,1,1]=(0.464093634725)*x[1]
        arg[0,2,1,2]=(-0.0412800774497)*x[2]
        arg[0,2,2,0]=(0.223364317185)*x[0]
        arg[0,2,2,1]=(0.257201130157)*x[1]
        arg[0,2,2,2]=(0.063203467463)*x[2]
        arg[0,3,0,0]=(-0.723240451643)*x[0]
        arg[0,3,0,1]=(-0.862468295097)*x[1]
        arg[0,3,0,2]=(-0.149283247587)*x[2]
        arg[0,3,1,0]=(0.15680097839)*x[0]
        arg[0,3,1,1]=(0.421563637547)*x[1]
        arg[0,3,1,2]=(0.111549188549)*x[2]
        arg[0,3,2,0]=(-0.272783329363)*x[0]
        arg[0,3,2,1]=(-0.420352789853)*x[1]
        arg[0,3,2,2]=(0.570865117722)*x[2]
        arg[0,4,0,0]=(-0.321910078414)*x[0]
        arg[0,4,0,1]=(0.988695599439)*x[1]
        arg[0,4,0,2]=(0.920200893398)*x[2]
        arg[0,4,1,0]=(0.0260910072651)*x[0]
        arg[0,4,1,1]=(0.460012578184)*x[1]
        arg[0,4,1,2]=(0.848099524112)*x[2]
        arg[0,4,2,0]=(0.242157803251)*x[0]
        arg[0,4,2,1]=(0.394528777004)*x[1]
        arg[0,4,2,2]=(0.562996837311)*x[2]
        arg[1,0,0,0]=(0.459886225958)*x[0]
        arg[1,0,0,1]=(-0.721868942003)*x[1]
        arg[1,0,0,2]=(0.432203082994)*x[2]
        arg[1,0,1,0]=(0.409831045482)*x[0]
        arg[1,0,1,1]=(-0.481677513473)*x[1]
        arg[1,0,1,2]=(0.439387853437)*x[2]
        arg[1,0,2,0]=(0.261583198434)*x[0]
        arg[1,0,2,1]=(0.290993423577)*x[1]
        arg[1,0,2,2]=(0.477993114134)*x[2]
        arg[1,1,0,0]=(0.586344598248)*x[0]
        arg[1,1,0,1]=(-0.105390792831)*x[1]
        arg[1,1,0,2]=(0.335990751314)*x[2]
        arg[1,1,1,0]=(-0.191500562856)*x[0]
        arg[1,1,1,1]=(0.244514598216)*x[1]
        arg[1,1,1,2]=(-0.804402720669)*x[2]
        arg[1,1,2,0]=(-0.455225710648)*x[0]
        arg[1,1,2,1]=(-0.505052700585)*x[1]
        arg[1,1,2,2]=(-0.0240295199362)*x[2]
        arg[1,2,0,0]=(-0.718487964893)*x[0]
        arg[1,2,0,1]=(-0.0899522570462)*x[1]
        arg[1,2,0,2]=(-0.293353754696)*x[2]
        arg[1,2,1,0]=(-0.180013826342)*x[0]
        arg[1,2,1,1]=(0.793689231922)*x[1]
        arg[1,2,1,2]=(0.673066555571)*x[2]
        arg[1,2,2,0]=(0.705362155032)*x[0]
        arg[1,2,2,1]=(0.54476742883)*x[1]
        arg[1,2,2,2]=(-0.331195064878)*x[2]
        arg[1,3,0,0]=(-0.360927441647)*x[0]
        arg[1,3,0,1]=(0.230772030282)*x[1]
        arg[1,3,0,2]=(0.912342489431)*x[2]
        arg[1,3,1,0]=(-0.817510690014)*x[0]
        arg[1,3,1,1]=(0.397583721353)*x[1]
        arg[1,3,1,2]=(-0.982551067917)*x[2]
        arg[1,3,2,0]=(0.86380240427)*x[0]
        arg[1,3,2,1]=(-0.415018976841)*x[1]
        arg[1,3,2,2]=(0.271582572267)*x[2]
        arg[1,4,0,0]=(0.252845347406)*x[0]
        arg[1,4,0,1]=(0.687786802906)*x[1]
        arg[1,4,0,2]=(0.465501171342)*x[2]
        arg[1,4,1,0]=(-0.613703721675)*x[0]
        arg[1,4,1,1]=(-0.110297640533)*x[1]
        arg[1,4,1,2]=(-0.836768056501)*x[2]
        arg[1,4,2,0]=(-0.0400898232224)*x[0]
        arg[1,4,2,1]=(0.0358172759009)*x[1]
        arg[1,4,2,2]=(-0.335751455408)*x[2]
        arg[2,0,0,0]=(-0.309992915015)*x[0]
        arg[2,0,0,1]=(-0.721404217867)*x[1]
        arg[2,0,0,2]=(-0.548000635629)*x[2]
        arg[2,0,1,0]=(0.651175831531)*x[0]
        arg[2,0,1,1]=(0.158960783491)*x[1]
        arg[2,0,1,2]=(-0.310676926155)*x[2]
        arg[2,0,2,0]=(-0.122289734411)*x[0]
        arg[2,0,2,1]=(-0.252405938421)*x[1]
        arg[2,0,2,2]=(-0.938280244213)*x[2]
        arg[2,1,0,0]=(0.559495801686)*x[0]
        arg[2,1,0,1]=(-0.547182622716)*x[1]
        arg[2,1,0,2]=(0.397441517898)*x[2]
        arg[2,1,1,0]=(-0.406112472071)*x[0]
        arg[2,1,1,1]=(0.355063810677)*x[1]
        arg[2,1,1,2]=(0.760400203215)*x[2]
        arg[2,1,2,0]=(0.992201320481)*x[0]
        arg[2,1,2,1]=(0.0580660882576)*x[1]
        arg[2,1,2,2]=(-0.643170879939)*x[2]
        arg[2,2,0,0]=(-0.280644461832)*x[0]
        arg[2,2,0,1]=(-0.0467430285531)*x[1]
        arg[2,2,0,2]=(0.314050593255)*x[2]
        arg[2,2,1,0]=(-0.230032618609)*x[0]
        arg[2,2,1,1]=(0.0996058698273)*x[1]
        arg[2,2,1,2]=(-0.0270266073208)*x[2]
        arg[2,2,2,0]=(0.767914132956)*x[0]
        arg[2,2,2,1]=(0.496930363612)*x[1]
        arg[2,2,2,2]=(-0.599525033616)*x[2]
        arg[2,3,0,0]=(-0.326433376073)*x[0]
        arg[2,3,0,1]=(-0.0366374501025)*x[1]
        arg[2,3,0,2]=(0.22555705749)*x[2]
        arg[2,3,1,0]=(-0.162548813895)*x[0]
        arg[2,3,1,1]=(-0.110074212194)*x[1]
        arg[2,3,1,2]=(-0.143600895553)*x[2]
        arg[2,3,2,0]=(0.771148880174)*x[0]
        arg[2,3,2,1]=(0.112528116552)*x[1]
        arg[2,3,2,2]=(-0.955735294341)*x[2]
        arg[2,4,0,0]=(-0.968392951034)*x[0]
        arg[2,4,0,1]=(-0.36901708507)*x[1]
        arg[2,4,0,2]=(0.283692515492)*x[2]
        arg[2,4,1,0]=(0.997238032837)*x[0]
        arg[2,4,1,1]=(-0.625794124653)*x[1]
        arg[2,4,1,2]=(0.533386027556)*x[2]
        arg[2,4,2,0]=(0.977311695557)*x[0]
        arg[2,4,2,1]=(0.693009976689)*x[1]
        arg[2,4,2,2]=(0.711179347652)*x[2]
        arg[3,0,0,0]=(-0.155585788931)*x[0]
        arg[3,0,0,1]=(0.0228078851234)*x[1]
        arg[3,0,0,2]=(0.510104938032)*x[2]
        arg[3,0,1,0]=(0.74865995369)*x[0]
        arg[3,0,1,1]=(0.672153736284)*x[1]
        arg[3,0,1,2]=(0.588012355098)*x[2]
        arg[3,0,2,0]=(-0.924508475715)*x[0]
        arg[3,0,2,1]=(-0.392784674758)*x[1]
        arg[3,0,2,2]=(-0.36371454642)*x[2]
        arg[3,1,0,0]=(-0.709783490337)*x[0]
        arg[3,1,0,1]=(0.844136172222)*x[1]
        arg[3,1,0,2]=(0.621011730043)*x[2]
        arg[3,1,1,0]=(0.428807337181)*x[0]
        arg[3,1,1,1]=(0.126300214574)*x[1]
        arg[3,1,1,2]=(0.795972806221)*x[2]
        arg[3,1,2,0]=(-0.252334324004)*x[0]
        arg[3,1,2,1]=(-0.722829467938)*x[1]
        arg[3,1,2,2]=(-0.551540062366)*x[2]
        arg[3,2,0,0]=(-0.134668475963)*x[0]
        arg[3,2,0,1]=(-0.598747540536)*x[1]
        arg[3,2,0,2]=(0.426422436624)*x[2]
        arg[3,2,1,0]=(-0.363050323762)*x[0]
        arg[3,2,1,1]=(0.980891457977)*x[1]
        arg[3,2,1,2]=(0.162831912555)*x[2]
        arg[3,2,2,0]=(-0.126505493475)*x[0]
        arg[3,2,2,1]=(-0.578567864811)*x[1]
        arg[3,2,2,2]=(-0.509843129095)*x[2]
        arg[3,3,0,0]=(-0.446171262265)*x[0]
        arg[3,3,0,1]=(-0.715175197494)*x[1]
        arg[3,3,0,2]=(-0.881016888806)*x[2]
        arg[3,3,1,0]=(-0.942020866327)*x[0]
        arg[3,3,1,1]=(0.156434646828)*x[1]
        arg[3,3,1,2]=(0.523624761583)*x[2]
        arg[3,3,2,0]=(-0.683550923926)*x[0]
        arg[3,3,2,1]=(0.857075218033)*x[1]
        arg[3,3,2,2]=(0.297672594023)*x[2]
        arg[3,4,0,0]=(0.74317121113)*x[0]
        arg[3,4,0,1]=(0.076464540756)*x[1]
        arg[3,4,0,2]=(0.781965468281)*x[2]
        arg[3,4,1,0]=(0.417750169098)*x[0]
        arg[3,4,1,1]=(0.82275428729)*x[1]
        arg[3,4,1,2]=(0.919072321093)*x[2]
        arg[3,4,2,0]=(-0.0246706472217)*x[0]
        arg[3,4,2,1]=(0.179863245513)*x[1]
        arg[3,4,2,2]=(0.539115287766)*x[2]
        ref=sqrt((52.492676775)*(2.*dim+1.)/3.)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunctionOnBoundary_fromData_rank0(self):
      """
      tests L2-norm of Data on the ReducedFunctionOnBoundary

      assumptions: self.domain supports integration on ReducedFunctionOnBoundary
      """
      dim=self.domain.getDim()
      w=ReducedFunctionOnBoundary(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(),w)
        arg=(-0.245574919477)*sqrt(x[0])
        ref=sqrt((0.0603070410759)*dim)

      else:
        arg=Data(0,(),w)
        arg=(0.757324521515)*sqrt(x[0])
        ref=sqrt((0.573540430888)*dim)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunctionOnBoundary_fromData_rank1(self):
      """
      tests L2-norm of Data on the ReducedFunctionOnBoundary

      assumptions: self.domain supports integration on ReducedFunctionOnBoundary
      """
      dim=self.domain.getDim()
      w=ReducedFunctionOnBoundary(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(2,),w)
        arg[0]=(0.723421565407)*sqrt(x[0])
        arg[1]=(-0.460477393103)*sqrt(x[1])
        ref=sqrt((0.735378190855)*dim)

      else:
        arg=Data(0,(3,),w)
        arg[0]=(-0.88528497163)*sqrt(x[0])
        arg[1]=(-0.65510214636)*sqrt(x[1])
        arg[2]=(0.399538866363)*sqrt(x[2])
        ref=sqrt((1.37251960889)*dim)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunctionOnBoundary_fromData_rank2(self):
      """
      tests L2-norm of Data on the ReducedFunctionOnBoundary

      assumptions: self.domain supports integration on ReducedFunctionOnBoundary
      """
      dim=self.domain.getDim()
      w=ReducedFunctionOnBoundary(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 2),w)
        arg[0,0]=(0.955620993904)*sqrt(x[0])
        arg[0,1]=(-0.0987865813703)*sqrt(x[1])
        arg[1,0]=(0.0288267231531)*sqrt(x[0])
        arg[1,1]=(0.655440599879)*sqrt(x[1])
        arg[2,0]=(0.685627284533)*sqrt(x[0])
        arg[2,1]=(-0.989832824892)*sqrt(x[1])
        arg[3,0]=(0.292184093194)*sqrt(x[0])
        arg[3,1]=(0.149553857773)*sqrt(x[1])
        ref=sqrt((2.91099532781)*dim)

      else:
        arg=Data(0,(4, 3),w)
        arg[0,0]=(-0.325908541533)*sqrt(x[0])
        arg[0,1]=(-0.992480479749)*sqrt(x[1])
        arg[0,2]=(0.660360271799)*sqrt(x[2])
        arg[1,0]=(0.173485908581)*sqrt(x[0])
        arg[1,1]=(-0.328755199781)*sqrt(x[1])
        arg[1,2]=(-0.943354674948)*sqrt(x[2])
        arg[2,0]=(0.680713222646)*sqrt(x[0])
        arg[2,1]=(-0.765971835693)*sqrt(x[1])
        arg[2,2]=(0.0413284847528)*sqrt(x[2])
        arg[3,0]=(0.990074004708)*sqrt(x[0])
        arg[3,1]=(0.941801786766)*sqrt(x[1])
        arg[3,2]=(0.886926192201)*sqrt(x[2])
        ref=sqrt((6.26107155228)*dim)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunctionOnBoundary_fromData_rank3(self):
      """
      tests L2-norm of Data on the ReducedFunctionOnBoundary

      assumptions: self.domain supports integration on ReducedFunctionOnBoundary
      """
      dim=self.domain.getDim()
      w=ReducedFunctionOnBoundary(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(6, 2, 2),w)
        arg[0,0,0]=(-0.0781611551598)*sqrt(x[0])
        arg[0,0,1]=(0.291016249575)*sqrt(x[1])
        arg[0,1,0]=(-0.107555233086)*sqrt(x[0])
        arg[0,1,1]=(-0.559067108546)*sqrt(x[1])
        arg[1,0,0]=(-0.0818406701266)*sqrt(x[0])
        arg[1,0,1]=(-0.594866806483)*sqrt(x[1])
        arg[1,1,0]=(-0.725814803863)*sqrt(x[0])
        arg[1,1,1]=(0.59128101992)*sqrt(x[1])
        arg[2,0,0]=(-0.15381555291)*sqrt(x[0])
        arg[2,0,1]=(-0.679882948503)*sqrt(x[1])
        arg[2,1,0]=(-0.58437193917)*sqrt(x[0])
        arg[2,1,1]=(0.136304615849)*sqrt(x[1])
        arg[3,0,0]=(0.0671365410096)*sqrt(x[0])
        arg[3,0,1]=(-0.645687212187)*sqrt(x[1])
        arg[3,1,0]=(-0.642492412392)*sqrt(x[0])
        arg[3,1,1]=(-0.125760054735)*sqrt(x[1])
        arg[4,0,0]=(0.731110824794)*sqrt(x[0])
        arg[4,0,1]=(0.491668422979)*sqrt(x[1])
        arg[4,1,0]=(-0.775841478292)*sqrt(x[0])
        arg[4,1,1]=(0.728265567974)*sqrt(x[1])
        arg[5,0,0]=(0.84511832373)*sqrt(x[0])
        arg[5,0,1]=(-0.513796801068)*sqrt(x[1])
        arg[5,1,0]=(0.113072243554)*sqrt(x[0])
        arg[5,1,1]=(0.246630838744)*sqrt(x[1])
        ref=sqrt((6.30829536252)*dim)

      else:
        arg=Data(0,(6, 2, 3),w)
        arg[0,0,0]=(0.369748116859)*sqrt(x[0])
        arg[0,0,1]=(-0.758056560031)*sqrt(x[1])
        arg[0,0,2]=(-0.873984709951)*sqrt(x[2])
        arg[0,1,0]=(0.311680165784)*sqrt(x[0])
        arg[0,1,1]=(0.374400673651)*sqrt(x[1])
        arg[0,1,2]=(0.712484217076)*sqrt(x[2])
        arg[1,0,0]=(0.829379714484)*sqrt(x[0])
        arg[1,0,1]=(-0.0551589596149)*sqrt(x[1])
        arg[1,0,2]=(0.965672208426)*sqrt(x[2])
        arg[1,1,0]=(-0.205044281547)*sqrt(x[0])
        arg[1,1,1]=(0.238197452756)*sqrt(x[1])
        arg[1,1,2]=(-0.33456139292)*sqrt(x[2])
        arg[2,0,0]=(0.649928288926)*sqrt(x[0])
        arg[2,0,1]=(-0.661384953389)*sqrt(x[1])
        arg[2,0,2]=(-0.253241222975)*sqrt(x[2])
        arg[2,1,0]=(-0.491716575992)*sqrt(x[0])
        arg[2,1,1]=(-0.970872527468)*sqrt(x[1])
        arg[2,1,2]=(0.222410198921)*sqrt(x[2])
        arg[3,0,0]=(0.205752630262)*sqrt(x[0])
        arg[3,0,1]=(0.864804362697)*sqrt(x[1])
        arg[3,0,2]=(-0.417975564033)*sqrt(x[2])
        arg[3,1,0]=(0.586425694033)*sqrt(x[0])
        arg[3,1,1]=(0.952661122184)*sqrt(x[1])
        arg[3,1,2]=(0.608680080453)*sqrt(x[2])
        arg[4,0,0]=(0.625968903369)*sqrt(x[0])
        arg[4,0,1]=(-0.573909405003)*sqrt(x[1])
        arg[4,0,2]=(-0.762256394595)*sqrt(x[2])
        arg[4,1,0]=(0.0710742394418)*sqrt(x[0])
        arg[4,1,1]=(0.583378040574)*sqrt(x[1])
        arg[4,1,2]=(0.719032893115)*sqrt(x[2])
        arg[5,0,0]=(0.032173368884)*sqrt(x[0])
        arg[5,0,1]=(-0.434042549492)*sqrt(x[1])
        arg[5,0,2]=(0.363504765447)*sqrt(x[2])
        arg[5,1,0]=(0.598817469198)*sqrt(x[0])
        arg[5,1,1]=(-0.163967008775)*sqrt(x[1])
        arg[5,1,2]=(0.546778730604)*sqrt(x[2])
        ref=sqrt((11.9696343123)*dim)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_L2_onReducedFunctionOnBoundary_fromData_rank4(self):
      """
      tests L2-norm of Data on the ReducedFunctionOnBoundary

      assumptions: self.domain supports integration on ReducedFunctionOnBoundary
      """
      dim=self.domain.getDim()
      w=ReducedFunctionOnBoundary(self.domain)
      x=w.getX()
      if dim==2:
        arg=Data(0,(4, 5, 3, 2),w)
        arg[0,0,0,0]=(-0.273446520069)*sqrt(x[0])
        arg[0,0,0,1]=(-0.913305910831)*sqrt(x[1])
        arg[0,0,1,0]=(0.0745566747537)*sqrt(x[0])
        arg[0,0,1,1]=(0.98803601919)*sqrt(x[1])
        arg[0,0,2,0]=(-0.244120818875)*sqrt(x[0])
        arg[0,0,2,1]=(0.247509644998)*sqrt(x[1])
        arg[0,1,0,0]=(-0.548756417777)*sqrt(x[0])
        arg[0,1,0,1]=(-0.354587911923)*sqrt(x[1])
        arg[0,1,1,0]=(0.104268198867)*sqrt(x[0])
        arg[0,1,1,1]=(-0.541700877072)*sqrt(x[1])
        arg[0,1,2,0]=(-0.25900060689)*sqrt(x[0])
        arg[0,1,2,1]=(-0.859660175231)*sqrt(x[1])
        arg[0,2,0,0]=(0.195235560321)*sqrt(x[0])
        arg[0,2,0,1]=(0.175518738589)*sqrt(x[1])
        arg[0,2,1,0]=(-0.0638854232272)*sqrt(x[0])
        arg[0,2,1,1]=(-0.586161016541)*sqrt(x[1])
        arg[0,2,2,0]=(0.580258892247)*sqrt(x[0])
        arg[0,2,2,1]=(-0.931927145435)*sqrt(x[1])
        arg[0,3,0,0]=(-0.408209600298)*sqrt(x[0])
        arg[0,3,0,1]=(-0.0344882667014)*sqrt(x[1])
        arg[0,3,1,0]=(-0.131763534163)*sqrt(x[0])
        arg[0,3,1,1]=(-0.787653739965)*sqrt(x[1])
        arg[0,3,2,0]=(0.0910808104711)*sqrt(x[0])
        arg[0,3,2,1]=(-0.280409023611)*sqrt(x[1])
        arg[0,4,0,0]=(0.97745754012)*sqrt(x[0])
        arg[0,4,0,1]=(-0.59829020936)*sqrt(x[1])
        arg[0,4,1,0]=(0.890520260543)*sqrt(x[0])
        arg[0,4,1,1]=(-0.600760090231)*sqrt(x[1])
        arg[0,4,2,0]=(-0.897992297974)*sqrt(x[0])
        arg[0,4,2,1]=(-0.841169898923)*sqrt(x[1])
        arg[1,0,0,0]=(-0.249868816409)*sqrt(x[0])
        arg[1,0,0,1]=(0.620375082228)*sqrt(x[1])
        arg[1,0,1,0]=(-0.660480789306)*sqrt(x[0])
        arg[1,0,1,1]=(-0.73638571806)*sqrt(x[1])
        arg[1,0,2,0]=(0.339987316643)*sqrt(x[0])
        arg[1,0,2,1]=(0.541112529894)*sqrt(x[1])
        arg[1,1,0,0]=(-0.468808186705)*sqrt(x[0])
        arg[1,1,0,1]=(-0.32919679792)*sqrt(x[1])
        arg[1,1,1,0]=(0.917292803419)*sqrt(x[0])
        arg[1,1,1,1]=(-0.834265058005)*sqrt(x[1])
        arg[1,1,2,0]=(-0.247536849264)*sqrt(x[0])
        arg[1,1,2,1]=(-0.197503469238)*sqrt(x[1])
        arg[1,2,0,0]=(0.897591919909)*sqrt(x[0])
        arg[1,2,0,1]=(-0.807446231234)*sqrt(x[1])
        arg[1,2,1,0]=(-0.369878499382)*sqrt(x[0])
        arg[1,2,1,1]=(0.985678692179)*sqrt(x[1])
        arg[1,2,2,0]=(-0.709976427525)*sqrt(x[0])
        arg[1,2,2,1]=(-0.368744647016)*sqrt(x[1])
        arg[1,3,0,0]=(0.299630726462)*sqrt(x[0])
        arg[1,3,0,1]=(-0.445295757899)*sqrt(x[1])
        arg[1,3,1,0]=(-0.922386577254)*sqrt(x[0])
        arg[1,3,1,1]=(0.234794853697)*sqrt(x[1])
        arg[1,3,2,0]=(0.953377720197)*sqrt(x[0])
        arg[1,3,2,1]=(0.409778183998)*sqrt(x[1])
        arg[1,4,0,0]=(0.271967945488)*sqrt(x[0])
        arg[1,4,0,1]=(-0.578629001202)*sqrt(x[1])
        arg[1,4,1,0]=(0.210755226769)*sqrt(x[0])
        arg[1,4,1,1]=(-0.0902751419945)*sqrt(x[1])
        arg[1,4,2,0]=(0.70033387381)*sqrt(x[0])
        arg[1,4,2,1]=(0.305733565661)*sqrt(x[1])
        arg[2,0,0,0]=(-0.662483167298)*sqrt(x[0])
        arg[2,0,0,1]=(0.585252048652)*sqrt(x[1])
        arg[2,0,1,0]=(-0.398813785959)*sqrt(x[0])
        arg[2,0,1,1]=(-0.797438697186)*sqrt(x[1])
        arg[2,0,2,0]=(-0.508308971009)*sqrt(x[0])
        arg[2,0,2,1]=(0.302249407524)*sqrt(x[1])
        arg[2,1,0,0]=(0.208644491879)*sqrt(x[0])
        arg[2,1,0,1]=(-0.604749055374)*sqrt(x[1])
        arg[2,1,1,0]=(0.641654284594)*sqrt(x[0])
        arg[2,1,1,1]=(0.456898593356)*sqrt(x[1])
        arg[2,1,2,0]=(-0.398043867778)*sqrt(x[0])
        arg[2,1,2,1]=(-0.0712344657587)*sqrt(x[1])
        arg[2,2,0,0]=(0.0967860954865)*sqrt(x[0])
        arg[2,2,0,1]=(0.520449905952)*sqrt(x[1])
        arg[2,2,1,0]=(0.770198029595)*sqrt(x[0])
        arg[2,2,1,1]=(-0.594004671621)*sqrt(x[1])
        arg[2,2,2,0]=(-0.744571452885)*sqrt(x[0])
        arg[2,2,2,1]=(0.544447367825)*sqrt(x[1])
        arg[2,3,0,0]=(-0.137087966968)*sqrt(x[0])
        arg[2,3,0,1]=(0.120672667497)*sqrt(x[1])
        arg[2,3,1,0]=(0.204800088057)*sqrt(x[0])
        arg[2,3,1,1]=(0.626526346076)*sqrt(x[1])
        arg[2,3,2,0]=(-0.696480393227)*sqrt(x[0])
        arg[2,3,2,1]=(0.188533741996)*sqrt(x[1])
        arg[2,4,0,0]=(-0.403523821067)*sqrt(x[0])
        arg[2,4,0,1]=(-0.428048989483)*sqrt(x[1])
        arg[2,4,1,0]=(-0.244186366584)*sqrt(x[0])
        arg[2,4,1,1]=(0.00866444909003)*sqrt(x[1])
        arg[2,4,2,0]=(-0.445991308853)*sqrt(x[0])
        arg[2,4,2,1]=(-0.899951068935)*sqrt(x[1])
        arg[3,0,0,0]=(0.609340085418)*sqrt(x[0])
        arg[3,0,0,1]=(0.878750391425)*sqrt(x[1])
        arg[3,0,1,0]=(0.258064654464)*sqrt(x[0])
        arg[3,0,1,1]=(-0.482402612985)*sqrt(x[1])
        arg[3,0,2,0]=(0.943732283389)*sqrt(x[0])
        arg[3,0,2,1]=(0.65514211843)*sqrt(x[1])
        arg[3,1,0,0]=(-0.894551979619)*sqrt(x[0])
        arg[3,1,0,1]=(0.220116541042)*sqrt(x[1])
        arg[3,1,1,0]=(0.386887699577)*sqrt(x[0])
        arg[3,1,1,1]=(-0.422560075108)*sqrt(x[1])
        arg[3,1,2,0]=(0.00387273783493)*sqrt(x[0])
        arg[3,1,2,1]=(0.465673505613)*sqrt(x[1])
        arg[3,2,0,0]=(0.987383428982)*sqrt(x[0])
        arg[3,2,0,1]=(0.376320055964)*sqrt(x[1])
        arg[3,2,1,0]=(-0.463778689128)*sqrt(x[0])
        arg[3,2,1,1]=(0.179816566227)*sqrt(x[1])
        arg[3,2,2,0]=(0.961522856801)*sqrt(x[0])
        arg[3,2,2,1]=(-0.257779627946)*sqrt(x[1])
        arg[3,3,0,0]=(0.748886531458)*sqrt(x[0])
        arg[3,3,0,1]=(-0.257342282566)*sqrt(x[1])
        arg[3,3,1,0]=(0.377494024401)*sqrt(x[0])
        arg[3,3,1,1]=(-0.334588346017)*sqrt(x[1])
        arg[3,3,2,0]=(0.502495149189)*sqrt(x[0])
        arg[3,3,2,1]=(0.534612429702)*sqrt(x[1])
        arg[3,4,0,0]=(-0.308551337355)*sqrt(x[0])
        arg[3,4,0,1]=(-0.471825826745)*sqrt(x[1])
        arg[3,4,1,0]=(-0.262606531584)*sqrt(x[0])
        arg[3,4,1,1]=(0.766089616367)*sqrt(x[1])
        arg[3,4,2,0]=(-0.136526755642)*sqrt(x[0])
        arg[3,4,2,1]=(0.675111459363)*sqrt(x[1])
        ref=sqrt((37.3453550914)*dim)

      else:
        arg=Data(0,(4, 5, 3, 3),w)
        arg[0,0,0,0]=(0.532946231146)*sqrt(x[0])
        arg[0,0,0,1]=(0.269364089513)*sqrt(x[1])
        arg[0,0,0,2]=(-0.207412457081)*sqrt(x[2])
        arg[0,0,1,0]=(-0.843104704858)*sqrt(x[0])
        arg[0,0,1,1]=(0.0416216508473)*sqrt(x[1])
        arg[0,0,1,2]=(-0.836074693662)*sqrt(x[2])
        arg[0,0,2,0]=(0.943609268731)*sqrt(x[0])
        arg[0,0,2,1]=(0.0154543737816)*sqrt(x[1])
        arg[0,0,2,2]=(0.0726268788381)*sqrt(x[2])
        arg[0,1,0,0]=(0.108422740078)*sqrt(x[0])
        arg[0,1,0,1]=(-0.296667916638)*sqrt(x[1])
        arg[0,1,0,2]=(-0.769732600535)*sqrt(x[2])
        arg[0,1,1,0]=(-0.428575493834)*sqrt(x[0])
        arg[0,1,1,1]=(0.421245456722)*sqrt(x[1])
        arg[0,1,1,2]=(-0.588277652277)*sqrt(x[2])
        arg[0,1,2,0]=(0.145294576795)*sqrt(x[0])
        arg[0,1,2,1]=(0.323206623794)*sqrt(x[1])
        arg[0,1,2,2]=(0.788115602892)*sqrt(x[2])
        arg[0,2,0,0]=(-0.227877282292)*sqrt(x[0])
        arg[0,2,0,1]=(-0.630647460719)*sqrt(x[1])
        arg[0,2,0,2]=(0.58754882135)*sqrt(x[2])
        arg[0,2,1,0]=(0.347191113403)*sqrt(x[0])
        arg[0,2,1,1]=(0.464093634725)*sqrt(x[1])
        arg[0,2,1,2]=(-0.0412800774497)*sqrt(x[2])
        arg[0,2,2,0]=(0.223364317185)*sqrt(x[0])
        arg[0,2,2,1]=(0.257201130157)*sqrt(x[1])
        arg[0,2,2,2]=(0.063203467463)*sqrt(x[2])
        arg[0,3,0,0]=(-0.723240451643)*sqrt(x[0])
        arg[0,3,0,1]=(-0.862468295097)*sqrt(x[1])
        arg[0,3,0,2]=(-0.149283247587)*sqrt(x[2])
        arg[0,3,1,0]=(0.15680097839)*sqrt(x[0])
        arg[0,3,1,1]=(0.421563637547)*sqrt(x[1])
        arg[0,3,1,2]=(0.111549188549)*sqrt(x[2])
        arg[0,3,2,0]=(-0.272783329363)*sqrt(x[0])
        arg[0,3,2,1]=(-0.420352789853)*sqrt(x[1])
        arg[0,3,2,2]=(0.570865117722)*sqrt(x[2])
        arg[0,4,0,0]=(-0.321910078414)*sqrt(x[0])
        arg[0,4,0,1]=(0.988695599439)*sqrt(x[1])
        arg[0,4,0,2]=(0.920200893398)*sqrt(x[2])
        arg[0,4,1,0]=(0.0260910072651)*sqrt(x[0])
        arg[0,4,1,1]=(0.460012578184)*sqrt(x[1])
        arg[0,4,1,2]=(0.848099524112)*sqrt(x[2])
        arg[0,4,2,0]=(0.242157803251)*sqrt(x[0])
        arg[0,4,2,1]=(0.394528777004)*sqrt(x[1])
        arg[0,4,2,2]=(0.562996837311)*sqrt(x[2])
        arg[1,0,0,0]=(0.459886225958)*sqrt(x[0])
        arg[1,0,0,1]=(-0.721868942003)*sqrt(x[1])
        arg[1,0,0,2]=(0.432203082994)*sqrt(x[2])
        arg[1,0,1,0]=(0.409831045482)*sqrt(x[0])
        arg[1,0,1,1]=(-0.481677513473)*sqrt(x[1])
        arg[1,0,1,2]=(0.439387853437)*sqrt(x[2])
        arg[1,0,2,0]=(0.261583198434)*sqrt(x[0])
        arg[1,0,2,1]=(0.290993423577)*sqrt(x[1])
        arg[1,0,2,2]=(0.477993114134)*sqrt(x[2])
        arg[1,1,0,0]=(0.586344598248)*sqrt(x[0])
        arg[1,1,0,1]=(-0.105390792831)*sqrt(x[1])
        arg[1,1,0,2]=(0.335990751314)*sqrt(x[2])
        arg[1,1,1,0]=(-0.191500562856)*sqrt(x[0])
        arg[1,1,1,1]=(0.244514598216)*sqrt(x[1])
        arg[1,1,1,2]=(-0.804402720669)*sqrt(x[2])
        arg[1,1,2,0]=(-0.455225710648)*sqrt(x[0])
        arg[1,1,2,1]=(-0.505052700585)*sqrt(x[1])
        arg[1,1,2,2]=(-0.0240295199362)*sqrt(x[2])
        arg[1,2,0,0]=(-0.718487964893)*sqrt(x[0])
        arg[1,2,0,1]=(-0.0899522570462)*sqrt(x[1])
        arg[1,2,0,2]=(-0.293353754696)*sqrt(x[2])
        arg[1,2,1,0]=(-0.180013826342)*sqrt(x[0])
        arg[1,2,1,1]=(0.793689231922)*sqrt(x[1])
        arg[1,2,1,2]=(0.673066555571)*sqrt(x[2])
        arg[1,2,2,0]=(0.705362155032)*sqrt(x[0])
        arg[1,2,2,1]=(0.54476742883)*sqrt(x[1])
        arg[1,2,2,2]=(-0.331195064878)*sqrt(x[2])
        arg[1,3,0,0]=(-0.360927441647)*sqrt(x[0])
        arg[1,3,0,1]=(0.230772030282)*sqrt(x[1])
        arg[1,3,0,2]=(0.912342489431)*sqrt(x[2])
        arg[1,3,1,0]=(-0.817510690014)*sqrt(x[0])
        arg[1,3,1,1]=(0.397583721353)*sqrt(x[1])
        arg[1,3,1,2]=(-0.982551067917)*sqrt(x[2])
        arg[1,3,2,0]=(0.86380240427)*sqrt(x[0])
        arg[1,3,2,1]=(-0.415018976841)*sqrt(x[1])
        arg[1,3,2,2]=(0.271582572267)*sqrt(x[2])
        arg[1,4,0,0]=(0.252845347406)*sqrt(x[0])
        arg[1,4,0,1]=(0.687786802906)*sqrt(x[1])
        arg[1,4,0,2]=(0.465501171342)*sqrt(x[2])
        arg[1,4,1,0]=(-0.613703721675)*sqrt(x[0])
        arg[1,4,1,1]=(-0.110297640533)*sqrt(x[1])
        arg[1,4,1,2]=(-0.836768056501)*sqrt(x[2])
        arg[1,4,2,0]=(-0.0400898232224)*sqrt(x[0])
        arg[1,4,2,1]=(0.0358172759009)*sqrt(x[1])
        arg[1,4,2,2]=(-0.335751455408)*sqrt(x[2])
        arg[2,0,0,0]=(-0.309992915015)*sqrt(x[0])
        arg[2,0,0,1]=(-0.721404217867)*sqrt(x[1])
        arg[2,0,0,2]=(-0.548000635629)*sqrt(x[2])
        arg[2,0,1,0]=(0.651175831531)*sqrt(x[0])
        arg[2,0,1,1]=(0.158960783491)*sqrt(x[1])
        arg[2,0,1,2]=(-0.310676926155)*sqrt(x[2])
        arg[2,0,2,0]=(-0.122289734411)*sqrt(x[0])
        arg[2,0,2,1]=(-0.252405938421)*sqrt(x[1])
        arg[2,0,2,2]=(-0.938280244213)*sqrt(x[2])
        arg[2,1,0,0]=(0.559495801686)*sqrt(x[0])
        arg[2,1,0,1]=(-0.547182622716)*sqrt(x[1])
        arg[2,1,0,2]=(0.397441517898)*sqrt(x[2])
        arg[2,1,1,0]=(-0.406112472071)*sqrt(x[0])
        arg[2,1,1,1]=(0.355063810677)*sqrt(x[1])
        arg[2,1,1,2]=(0.760400203215)*sqrt(x[2])
        arg[2,1,2,0]=(0.992201320481)*sqrt(x[0])
        arg[2,1,2,1]=(0.0580660882576)*sqrt(x[1])
        arg[2,1,2,2]=(-0.643170879939)*sqrt(x[2])
        arg[2,2,0,0]=(-0.280644461832)*sqrt(x[0])
        arg[2,2,0,1]=(-0.0467430285531)*sqrt(x[1])
        arg[2,2,0,2]=(0.314050593255)*sqrt(x[2])
        arg[2,2,1,0]=(-0.230032618609)*sqrt(x[0])
        arg[2,2,1,1]=(0.0996058698273)*sqrt(x[1])
        arg[2,2,1,2]=(-0.0270266073208)*sqrt(x[2])
        arg[2,2,2,0]=(0.767914132956)*sqrt(x[0])
        arg[2,2,2,1]=(0.496930363612)*sqrt(x[1])
        arg[2,2,2,2]=(-0.599525033616)*sqrt(x[2])
        arg[2,3,0,0]=(-0.326433376073)*sqrt(x[0])
        arg[2,3,0,1]=(-0.0366374501025)*sqrt(x[1])
        arg[2,3,0,2]=(0.22555705749)*sqrt(x[2])
        arg[2,3,1,0]=(-0.162548813895)*sqrt(x[0])
        arg[2,3,1,1]=(-0.110074212194)*sqrt(x[1])
        arg[2,3,1,2]=(-0.143600895553)*sqrt(x[2])
        arg[2,3,2,0]=(0.771148880174)*sqrt(x[0])
        arg[2,3,2,1]=(0.112528116552)*sqrt(x[1])
        arg[2,3,2,2]=(-0.955735294341)*sqrt(x[2])
        arg[2,4,0,0]=(-0.968392951034)*sqrt(x[0])
        arg[2,4,0,1]=(-0.36901708507)*sqrt(x[1])
        arg[2,4,0,2]=(0.283692515492)*sqrt(x[2])
        arg[2,4,1,0]=(0.997238032837)*sqrt(x[0])
        arg[2,4,1,1]=(-0.625794124653)*sqrt(x[1])
        arg[2,4,1,2]=(0.533386027556)*sqrt(x[2])
        arg[2,4,2,0]=(0.977311695557)*sqrt(x[0])
        arg[2,4,2,1]=(0.693009976689)*sqrt(x[1])
        arg[2,4,2,2]=(0.711179347652)*sqrt(x[2])
        arg[3,0,0,0]=(-0.155585788931)*sqrt(x[0])
        arg[3,0,0,1]=(0.0228078851234)*sqrt(x[1])
        arg[3,0,0,2]=(0.510104938032)*sqrt(x[2])
        arg[3,0,1,0]=(0.74865995369)*sqrt(x[0])
        arg[3,0,1,1]=(0.672153736284)*sqrt(x[1])
        arg[3,0,1,2]=(0.588012355098)*sqrt(x[2])
        arg[3,0,2,0]=(-0.924508475715)*sqrt(x[0])
        arg[3,0,2,1]=(-0.392784674758)*sqrt(x[1])
        arg[3,0,2,2]=(-0.36371454642)*sqrt(x[2])
        arg[3,1,0,0]=(-0.709783490337)*sqrt(x[0])
        arg[3,1,0,1]=(0.844136172222)*sqrt(x[1])
        arg[3,1,0,2]=(0.621011730043)*sqrt(x[2])
        arg[3,1,1,0]=(0.428807337181)*sqrt(x[0])
        arg[3,1,1,1]=(0.126300214574)*sqrt(x[1])
        arg[3,1,1,2]=(0.795972806221)*sqrt(x[2])
        arg[3,1,2,0]=(-0.252334324004)*sqrt(x[0])
        arg[3,1,2,1]=(-0.722829467938)*sqrt(x[1])
        arg[3,1,2,2]=(-0.551540062366)*sqrt(x[2])
        arg[3,2,0,0]=(-0.134668475963)*sqrt(x[0])
        arg[3,2,0,1]=(-0.598747540536)*sqrt(x[1])
        arg[3,2,0,2]=(0.426422436624)*sqrt(x[2])
        arg[3,2,1,0]=(-0.363050323762)*sqrt(x[0])
        arg[3,2,1,1]=(0.980891457977)*sqrt(x[1])
        arg[3,2,1,2]=(0.162831912555)*sqrt(x[2])
        arg[3,2,2,0]=(-0.126505493475)*sqrt(x[0])
        arg[3,2,2,1]=(-0.578567864811)*sqrt(x[1])
        arg[3,2,2,2]=(-0.509843129095)*sqrt(x[2])
        arg[3,3,0,0]=(-0.446171262265)*sqrt(x[0])
        arg[3,3,0,1]=(-0.715175197494)*sqrt(x[1])
        arg[3,3,0,2]=(-0.881016888806)*sqrt(x[2])
        arg[3,3,1,0]=(-0.942020866327)*sqrt(x[0])
        arg[3,3,1,1]=(0.156434646828)*sqrt(x[1])
        arg[3,3,1,2]=(0.523624761583)*sqrt(x[2])
        arg[3,3,2,0]=(-0.683550923926)*sqrt(x[0])
        arg[3,3,2,1]=(0.857075218033)*sqrt(x[1])
        arg[3,3,2,2]=(0.297672594023)*sqrt(x[2])
        arg[3,4,0,0]=(0.74317121113)*sqrt(x[0])
        arg[3,4,0,1]=(0.076464540756)*sqrt(x[1])
        arg[3,4,0,2]=(0.781965468281)*sqrt(x[2])
        arg[3,4,1,0]=(0.417750169098)*sqrt(x[0])
        arg[3,4,1,1]=(0.82275428729)*sqrt(x[1])
        arg[3,4,1,2]=(0.919072321093)*sqrt(x[2])
        arg[3,4,2,0]=(-0.0246706472217)*sqrt(x[0])
        arg[3,4,2,1]=(0.179863245513)*sqrt(x[1])
        arg[3,4,2,2]=(0.539115287766)*sqrt(x[2])
        ref=sqrt((52.492676775)*dim)

      res=L2(arg)
      self.assertTrue(isinstance(res,float),"wrong type of result.")
      self.assertAlmostEqual(res,ref,int(-log10(self.RES_TOL)),"wrong result")

