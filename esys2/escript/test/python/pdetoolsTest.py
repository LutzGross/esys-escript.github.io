# $Id$

from esys.escript import *
from esys.escript.pdetools import Locator,Projector
from esys.finley import Rectangle

def testLocator(domain):
      """runs a few test of the Locator"""

      error_max=0.
      error_text=""
      x=domain.getX()
      l=Locator(domain,[1.,1.])
      print l.getFunctionSpace()
      print l.getId()

      l=Locator(ContinuousFunction(domain),[1.,1.])
      print l.getFunctionSpace()
      print l.getId()

      text="l.getX()"
      err=Lsup(l.getX()-numarray.ones((2,)))
      print text," error = ",err
      if err>=error_max: error_max,error_text=err,text

      text="l(x)"
      err=Lsup(l(x)-numarray.ones((2,)))
      print text," error = ",err
      if err>=error_max: error_max,error_text=err,text

      text="l(x[0]+x[1])"
      err=abs(l(x[0]+x[1])-2.)
      print text," error = ",err
      if err>=error_max: error_max,error_text=err,text

      return error_max,error_text
      
def testProjector(domain):
      """runs a few test of the Projector factory and returns the largest error plus a description of the test this error occured"""
      error_max=0.
      error_text=""
      x=ContinuousFunction(domain).getX()
      for f in [True,False]:
         p=Projector(domain,reduce=False,fast=f)
         for r in range(5):
            text="range %s , fast=%s"%(r,f)
            if r==0:
               td_ref=x[0]
            elif r==1:
               td_ref=x
            elif r==2:
               td_ref=[[11.,12.],[21,22.]]*(x[0]+x[1])
            elif r==3:
               td_ref=[[[111.,112.],[121,122.]],[[211.,212.],[221,222.]]]*(x[0]+x[1])
            elif r==3:
               td_ref=[[[[1111.,1112.],[1121,1122.]],[[1211.,1212.],[1221,1222.]]],[[[2111.,2112.],[2121,2122.]],[[2211.,2212.],[2221,2222.]]]]*(x[0]+x[1])
            td=p(td_ref.interpolate(Function(domain)))
            er=Lsup(td-td_ref)/Lsup(td_ref)
            print text," error = ",er
            if er>error_max:
                 error_max=er
                 error_text=text
      return error_max,error_text


txt=testLocator(Rectangle(56,61))
print "test Locator: ",txt

txt=testProjector(Rectangle(56,61))
print "test Projector: ",txt



