# $Id$
"""tests a variety of functions in connection with contact elements"""
from escript import *
from linearPDEs import LinearPDE
import finley
import math
numElements=4
numEquations=2

seed=1
rand_factor=sqrt(2.)
def randomNum():
     global seed
     seed+=1
     s=rand_factor*seed
     return s-int(s)


def mkMesh(dim,order,numElm,onElem):
    if dim==2:
         ms1=finley.Rectangle(numElm,numElm,order,l1=0.5,useElementsOnFace=onElem)
         ms2=finley.Rectangle(numElm,numElm,order,l1=0.5,useElementsOnFace=onElem)
         ms2.setX(ms2.getX()+[0,0.5])
    else:
         ms1=finley.Brick(numElm,numElm,numElm,order,l2=0.5,useElementsOnFace=onElem)
         ms2=finley.Brick(numElm,numElm,numElm,order,l2=0.5,useElementsOnFace=onElem)
         ms2.setX(ms2.getX()+[0,0,0.5])
    return finley.JoinFaces([ms1,ms2])


def mkCharateristicFunction(msh):
      """returns a scalar function on nodes which is -1 below and 1 above the fault"""
      options = {
          "verbose" : True,
          "reordering" : NO_REORDERING,
          "tolerance" : 1.E-13,
          "final_residual" : 0.,
          "iterative" : True,
          "iterative_method" : BICGSTAB,
          "preconditioner" : JACOBI,
          "iter_max" :  4000,
          "iter" : 0,
          "drop_tolerance" : 1.10,
          "drop_storage" : 2.
      }
      e=Function(msh)
      d=msh.getDim()
      x=e.getX()[d-1]
      mypde=LinearPDE(D=1,Y=(x-0.5).whatPositive())
      return 2*mypde.getSolution(**options)-1

max_error_text=""
max_error=0.
for dim in [2,3]:
    for order in [1,2]:
        for onElem in [False,True]:
            if onElem:
                  onElem_text=",elements on face"
            else:
                  onElem_text=""
            case="Contact: %dD, order %d%s"%(dim,order,onElem_text)
            print case
            msh=mkMesh(dim,order,numElements,onElem)
            c0=FunctionOnContact0(msh)
            c1=FunctionOnContact1(msh)
            n=ContinuousFunction(msh)
            #
            # check the normal on the fault 0:
            #
            refNormal=Vector(0,what=c0)
            if dim==3:
               refNormal.setTaggedValue(100,[0,0,1])
            else:
               refNormal.setTaggedValue(10,[0,1])
            error_norm=Lsup(c0.getNormal()-refNormal)
            text=" %s : error normals 0= %e"%(case,error_norm)
            print "@@@ %s"%text
            if error_norm>max_error: max_error_text,max_error=text,error_norm
            #
            # check the normal on the fault 1:
            #
            refNormal=Vector(0,what=c1)
            if dim==3:
               refNormal.setTaggedValue(100,[0,0,-1])
            else:
               refNormal.setTaggedValue(10,[0,-1])
            error_norm=Lsup(c1.getNormal()-refNormal)
            text=" %s : error normals 1= %e"%(case,error_norm)
            print "@@@ %s"%text
            if error_norm>max_error: max_error_text,max_error=text,error_norm
            #
            #  integration on 0:
            #
            g=c0.getX()[dim-1]**2
            error_norm=abs(g.integrate()-0.25)
            text=" %s : error integrate 0= %e"%(case,error_norm)
            print "@@@ %s"%text
            if error_norm>max_error: max_error_text,max_error=text,error_norm
            #
            #  integration on 1:
            #
            g=c1.getX()[0]**2
            error_norm=abs(g.integrate()-1./3.)
            text=" %s : error integrate 1= %e"%(case,error_norm)
            print "@@@ %s"%text
            if error_norm>max_error: max_error_text,max_error=text,error_norm
            #
            #   gradient:
            #
            if onElem:
               x=n.getX()
               char=(x[dim-1]-0.5).whereNegative()
               foo=2*x[dim-1]*char+(10*x[dim-1]-4)*(1-char)
 
               error_norm=Lsup(foo.grad(c0)[dim-1]-2)
               text=" %s : error gradient on 0= %e"%(case,error_norm)
               print "@@@ %s"%text
               if error_norm>max_error: max_error_text,max_error=text,error_norm

               error_norm=Lsup(foo.grad(c1)[dim-1]-10)
               text=" %s : error gradient on 1= %e"%(case,error_norm)
               print "@@@ %s"%text
               if error_norm>max_error: max_error_text,max_error=text,error_norm

               char=mkCharateristicFunction(msh)
               for red in [False,True]:
                  if red: 
                     case_red=",reduced"
                  else:
                     case_red=""
                  for ne in range(1,numEquations+1):
                     u_ex=Data(0,shape=(ne,),what=n)
                     for i in range(ne):
                        u_ex[i]=char*(i+1)*n.getX()[0]
                     u0=u_ex.interpolate(c0)
                     u1=u_ex.interpolate(c1)
                     if ne==1:
                       d=randomNum()
                       y=d*(u1-u0)
                     else:
                       d=Id(ne)
                       for i in range(ne):
                         for j in range(ne):
                           d[i,j]=randomNum()
                       y=Data(0,shape=(ne,),what=c0)
                       for i in range(ne):
                         for j in range(ne):
                           y[i]+=d[i,j]*(u1[j]-u0[j])
                     mypde=linearPDE(d_contact=d,y_contact=y)
                     mypde.setReducedOrderForEquationsTo(red)
                     error_norm=Lsup(mypde.getResidual(u_ex))/Lsup(mypde.getRightHandSide())
                     text=" %s%s, assembling %d equations: error = %e"%(case,case_red,ne,error_norm)
                     print "@@@ %s"%text
                     if error_norm>max_error: max_error_text,max_error=text,error_norm


print "maximal error for ",max_error_text

