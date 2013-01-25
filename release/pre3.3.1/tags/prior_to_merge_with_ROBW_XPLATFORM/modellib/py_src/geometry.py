# $Id$

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""


from esys.escript import *
from esys.escript.modelframe import Model,ParameterSet
from esys import finley

class FinleyReader(ParameterSet):
       """
       Generates a mesh over a rectangular domain finley.

       @ivar filename:
       @ivar intergrationOrder
       @ivar domain: 
       """
       def __init__(self,debug=False):
           super(FinleyReader,self).__init__(debug=debug)
           ParameterSet.__init__(self,debug=debug)
           self.declareParameter(source="none",\
                                 integrationOrder=-1)
           self._domain=None

       def domain(self):
          if self._domain==None:
              self._domain=finley.ReadMesh(self.source,self.integrationOrder) 
              self.trace("mesh read from %s"%self.source)           
          return self._domain
                       
class RectangularDomain(ParameterSet):
       """
       Generates a mesh over a rectangular domain finley.

       @ivar dim:
       @ivar l:
       @ivar n:
       @ivar order:
       @ivar periodic:
       @ivar intergration order:
       @ivar domain: 
       """
       def __init__(self,debug=False):
           super(RectangularDomain,self).__init__(debug=debug)
           self.declareParameter(dim=2,\
                                 l=[1.,1.,1.],\
                                 n=[10,10,10], \
				 order=1,\
                                 periodic=[False,False,False],\
                                 integrationOrder=-1)
           self._domain=None

       def domain(self):
          if self._domain==None:
             if self.dim==2:
                self._domain=finley.Rectangle(n0=self.n[0],\
                                             n1=self.n[1],\
                                             l0=self.l[0],\
                                             l1=self.l[1],\
                                             order=self.order, \
                                             periodic0=self.periodic[0], \
                                             periodic1=self.periodic[1], \
                                             integrationOrder=self.integrationOrder)
             else:
                self._domain=finley.Brick(n0=self.n[0],\
                                         n1=self.n[1],\
                                         n2=self.n[2],\
                                         l0=self.l[0],\
                                         l1=self.l[1],\
                                         l2=self.l[2],\
                                         order=self.order, \
                                         periodic0=self.periodic[0], \
                                         periodic1=self.periodic[1], \
                                         periodic2=self.periodic[2], \
                                         integrationOrder=self.integrationOrder)

          return self._domain

class ConstrainValue(Model):
       """
       selects values for a given distribution to be used as a constrain. the location of the
       constrain are he faces of a rectangular domain. This Model is typically used in 
       time dependend problems to fix the values in a given initial condition.
       """
       def __init__(self,debug=False):
           Model.__init__(self,debug=debug)
           self.declareParameter(domain=None, \
                                 value=0,  \
                                 top=True,  \
                                 bottom=True,\
                                 front=False, \
                                 back=False,\
                                 left=False,\
                                 right=False,\
                                 constrain_value = None,  \
                                 location_constrained_value=None)
       def doInitialization(self):
           """
           initialize time stepping
           """
           tol=1.e-8
           x=self.domain.getX()
           d=self.domain.getDim()
           self.location_constrained_value=0
           x0=x[0]
           mx=sup(x0)
           mn=inf(x0)
           if self.left: 
               self.location_constrained_value=self.location_constrained_value+whereZero(x0-mn,tol*(mx-mn))
           if self.right: 
               self.location_constrained_value=self.location_constrained_value+whereZero(x0-mx,tol*(mx-mn))
           x0=x[d-1]
           mx=sup(x0)
           mn=inf(x0)
           if self.bottom: 
               self.location_constrained_value=self.location_constrained_value+whereZero(x0-mn,tol*(mx-mn))
           if self.top: 
               self.location_constrained_value=self.location_constrained_value+whereZero(x0-mx,tol*(mx-mn))
           if d>2:
              x0=x[1]
              mx=sup(x0)
              mn=inf(x0)
              if self.front: 
                 self.location_constrained_value=self.location_constrained_value+whereZero(x0-mn,tol*(mx-mn))
              if self.back: 
                 self.location_constrained_value=self.location_constrained_value+whereZero(x0-mx,tol*(mx-mn))           
           self.constrain_value=self.value*self.location_constrained_value
           
class ScalarConstrainer(ParameterSet):
     """
     Creates a characteristic function for the location of constraints 
     for a scalar value.

     In the case that the spatial dimension is two, the arguments front 
     and back are ignored.

     @ivar domain (in): rectangular domain
     @ivar left (in): True to set a constraint at the left face of the 
               domain (x[0]=min x[0]), default is False
     @ivar right (in): True to set a constraint at the left face of the 
               domain (x[0]=max x[0]), default is False
     @ivar top (in): True to set a constraint at the left face of the 
               domain (x[1]=min x[1]), default is False
     @ivar bottom (in): True to set a constraint at the left face of the 
               domain (x[1]=max x[1]), default is False
     @ivar front (in): True to set a constraint at the left face of the 
               domain (x[2]=min x[2]), default is False
     @ivar back (in): True to set a constraint at the left face of the 
               domain (x[2]=max x[2]), default is False
     @ivar location_of_constraint (out): object that defines the location 
               of the constraints.
     """
     def __init__(self,debug=False):
           ParameterSet.__init__(self,debug=debug)
           self.declareParameter(domain=None, \
                                 left=False, \
                                 right=False, \
                                 top=False, \
                                 bottom=False, \
                                 front=False, \
                                 back=False)
           self._location_of_constraint=None

     def location_of_constraint(self):
          """
          Returns the mask of the location of constraint.
          """
          if self._location_of_constraint==None:
             x=self.domain.getX()
             self._location_of_constraint=Scalar(0,x.getFunctionSpace())
             if self.domain.getDim()==3:
                if self.left: self._location_of_constraint+=whereZero(x[0]-inf(x[0]))
                if self.right: self._location_of_constraint+=whereZero(x[0]-sup(x[0]))
                if self.front: self._location_of_constraint+=whereZero(x[1]-inf(x[1]))
                if self.back: self._location_of_constraint+=whereZero(x[1]-sup(x[1]))
                if self.bottom: self._location_of_constraint+=whereZero(x[2]-inf(x[2]))
                if self.top: self._location_of_constraint+=whereZero(x[2]-sup(x[2]))
             else:
                if self.left: self._location_of_constraint+=whereZero(x[0]-inf(x[0]))
                if self.right: self._location_of_constraint+=whereZero(x[0]-sup(x[0]))
                if self.bottom: self._location_of_constraint+=whereZero(x[1]-inf(x[1]))
                if self.top: self._location_of_constraint+=whereZero(x[1]-sup(x[1]))
          return self._location_of_constraint

class VectorConstrainer(ParameterSet):
      """
      Creates a characteristic function for the location of constraints 
      for a scalar value.

      @ivar domain (in): rectangular domain
      @ivar left (in): list of three boolean. left[i]==True sets a 
                constraint for the i-th component at the left
                face of the domain (x[0]=min x[0]), 
                default is [False,False,False]
      @ivar right (in): list of three boolean. left[i]==True sets a 
                constraint for the i-th component at the right
                face of the domain (x[0]=max x[0]), 
                default is [False,False,False]
      @ivar top (in): list of three boolean. left[i]==True sets a 
                constraint for the i-th component at the top
                face of the domain (x[1]=min x[1]), 
                default is [False,False,False]
      @ivar bottom (in): list of three boolean. left[i]==True sets a 
                constraint for the i-th component at the bottom
                face of the domain (x[1]=min x[1]), 
                default is [False,False,False]
      @ivar front (in): list of three boolean. left[i]==True sets a 
                constraint for the i-th component at the front
                face of the domain (x[2]=min x[2]), 
                default is [False,False,False]
      @ivar back (in): list of three boolean. left[i]==True sets a 
                constraint for the i-th component at the back
                face of the domain (x[2]=max x[2]), 
                default is [False,False,False]
      @ivar location_of_constraint (callable): object that defines the location of the constraints for each vector component.

      In the case that the spatial dimension is two, thh arguments front and
      back as well as the third component of each argument is ignored.
      """
      def __init__(self,debug=False):
           ParameterSet.__init__(self,debug=debug)
           self.declareParameter(domain=None, \
                                 left=[0,0,0],  \
                                 right=[0,0,0],  \
                                 top=[0,0,0],  \
                                 bottom=[0,0,0],  \
                                 front=[0,0,0],
                                 back=[0,0,0])
           self._location_of_constraint=None
      def location_of_constraint(self):
          """
          Returns the mask of the location of constraint.
          """
          if self._location_of_constraint==None:
             x=self.domain.getX()
             self._location_of_constraint=Vector(0,x.getFunctionSpace())
             if self.domain.getDim()==3:
                left_mask=whereZero(x[0]-inf(x[0]))
                if self.left[0]: self._location_of_constraint+=left_mask*[1.,0.,0.]
                if self.left[1]: self._location_of_constraint+=left_mask*[0.,1.,0.]
                if self.left[2]: self._location_of_constraint+=left_mask*[0.,0.,1.]
                right_mask=whereZero(x[0]-sup(x[0]))
                if self.right[0]: self._location_of_constraint+=right_mask*[1.,0.,0.]
                if self.right[1]: self._location_of_constraint+=right_mask*[0.,1.,0.]
                if self.right[2]: self._location_of_constraint+=right_mask*[0.,0.,1.]
                front_mask=whereZero(x[1]-inf(x[1]))
                if self.front[0]: self._location_of_constraint+=front_mask*[1.,0.,0.]
                if self.front[1]: self._location_of_constraint+=front_mask*[0.,1.,0.]
                if self.front[2]: self._location_of_constraint+=front_mask*[0.,0.,1.]
                back_mask=whereZero(x[1]-sup(x[1]))
                if self.back[0]: self._location_of_constraint+=back_mask*[1.,0.,0.]
                if self.back[1]: self._location_of_constraint+=back_mask*[0.,1.,0.]
                if self.back[2]: self._location_of_constraint+=back_mask*[0.,0.,1.]
                bottom_mask=whereZero(x[2]-inf(x[2]))
                if self.bottom[0]: self._location_of_constraint+=bottom_mask*[1.,0.,0.]
                if self.bottom[1]: self._location_of_constraint+=bottom_mask*[0.,1.,0.]
                if self.bottom[2]: self._location_of_constraint+=bottom_mask*[0.,0.,1.]
                top_mask=whereZero(x[2]-sup(x[2]))
                if self.top[0]: self._location_of_constraint+=top_mask*[1.,0.,0.]
                if self.top[1]: self._location_of_constraint+=top_mask*[0.,1.,0.]
                if self.top[2]: self._location_of_constraint+=top_mask*[0.,0.,1.]
             else:
                left_mask=whereZero(x[0]-inf(x[0]))
                if self.left[0]: self._location_of_constraint+=left_mask*[1.,0.]
                if self.left[1]: self._location_of_constraint+=left_mask*[0.,1.]
                right_mask=whereZero(x[0]-sup(x[0]))
                if self.right[0]: self._location_of_constraint+=right_mask*[1.,0.]
                if self.right[1]: self._location_of_constraint+=right_mask*[0.,1.]
                bottom_mask=whereZero(x[1]-inf(x[1]))
                if self.bottom[0]: self._location_of_constraint+=bottom_mask*[1.,0.]
                if self.bottom[1]: self._location_of_constraint+=bottom_mask*[0.,1.]
                top_mask=whereZero(x[1]-sup(x[1]))
                if self.top[0]: self._location_of_constraint+=top_mask*[1.,0.]
                if self.top[1]: self._location_of_constraint+=top_mask*[0.,1.]
          return self._location_of_constraint

# vim: expandtab shiftwidth=4:
