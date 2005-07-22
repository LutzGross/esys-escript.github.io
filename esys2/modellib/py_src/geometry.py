# $Id$


from esys.modelframe import Model
from esys.escript import *
import esys.finley as finley

class RectangularDomain(Model):
       """   """
       def __init__(self,debug=False):
           Model.__init__(self,debug=debug)
           self.declareParameter(domain=None, dim=2,\
                                 l=[1.,1.,1.],\
                                 n=[10,10,10], \
				 order=1,\
                                 periodic=[False,False,False],\
                                 integrationOrder=-1)
       def doInitialization(self,t):
          if self.dim==2:
             self.domain=finley.Rectangle(n0=self.n[0],\
                                          n1=self.n[1],\
                                          l0=self.l[0],\
                                          l1=self.l[1],\
                                          order=self.order, \
                                          periodic0=self.periodic[0], \
                                          periodic1=self.periodic[1], \
                                          integrationOrder=self.integrationOrder)
          else:
             self.domain=finley.Brick(n0=self.n[0],\
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

class ScalarConstrainer(Model):
     """@brief creates a characteristic function for the location of constraints for a scalar value

              @param domain (in) - rectangular domain
              @param left (in)  - True to set a constraint at the left face of the domain (x[0]=min x[0]), default is False
              @param right (in) - True to set a constraint at the left face of the domain (x[0]=max x[0]), default is False
              @param top (in)  - True to set a constraint at the left face of the domain (x[1]=min x[1]), default is False
              @param bottom (in) - True to set a constraint at the left face of the domain (x[1]=max x[1]), default is False
              @param front (in)  - True to set a constraint at the left face of the domain (x[2]=min x[2]), default is False
              @param back (in) - True to set a constraint at the left face of the domain (x[2]=max x[2]), default is False
              @param location_of_constraint (out) - object that defines the location of the constraints.

       In the case that the spatial dimension is two, teh arguments front and back are ignored

     """
     def __init__(self,debug=False):
           Model.__init__(self,debug=debug)
           self.declareParameter(domain=None, \
                                 left=False,
                                 right=False,
                                 top=False,
                                 bottom=False,
                                 front=False,
                                 back=False,
                                 location_of_constraint=Data())
     def doInitialization(self,t):
          x=self.domain.getX()
          self.location_of_constraint=Scalar(0,x.getFunctionSpace())
          if self.domain.getDim()==3:
             if self.left: self.location_of_constraint+=(x[0]-inf(x[0])).whereZero()
             if self.right: self.location_of_constraint+=(x[0]-sup(x[0])).whereZero()
             if self.front: self.location_of_constraint+=(x[1]-inf(x[1])).whereZero()
             if self.back: self.location_of_constraint+=(x[1]-sup(x[1])).whereZero()
             if self.bottom: self.location_of_constraint+=(x[2]-inf(x[2])).whereZero()
             if self.top: self.location_of_constraint+=(x[2]-sup(x[2])).whereZero()
          else:
             if self.left: self.location_of_constraint+=(x[0]-inf(x[0])).whereZero()
             if self.right: self.location_of_constraint+=(x[0]-sup(x[0])).whereZero()
             if self.bottom: self.location_of_constraint+=(x[1]-inf(x[1])).whereZero()
             if self.top: self.location_of_constraint+=(x[1]-sup(x[1])).whereZero()

class VectorConstrainer(Model):
      """@brief creates a characteristic function for the location of constraints for a scalar value

              @param domain (in) - rectangular domain
              @param left (in)  - list of three boolean. left[i]==True sets a constraint for the i-th component at the left
                                  face of the domain (x[0]=min x[0]), default is [False,False,False]
              @param right (in) - list of three boolean. left[i]==True sets a constraint for the i-th component at the right
                                  face of the domain (x[0]=max x[0]), default is [False,False,False]
              @param top (in)  - list of three boolean. left[i]==True sets a constraint for the i-th component at the top
                                  face of the domain (x[1]=min x[1]), default is [False,False,False]
              @param bottom (in) - list of three boolean. left[i]==True sets a constraint for the i-th component at the bottom
                                  face of the domain (x[1]=min x[1]), default is [False,False,False]
              @param front (in)  - list of three boolean. left[i]==True sets a constraint for the i-th component at the front
                                  face of the domain (x[2]=min x[2]), default is [False,False,False]
              @param back (in) - list of three boolean. left[i]==True sets a constraint for the i-th component at the back
                                  face of the domain (x[2]=max x[2]), default is [False,False,False]
              @param location_of_constraint (out) - object that defines the location of the constraints for each vector component.

       In the case that the spatial dimension is two, thh arguments front and back as well as the third component of each argument is ignored.

      """
      def __init__(self,debug=False):
           Model.__init__(self,debug=debug)
           self.declareParameter(domain=None, \
                                 left=[0,0,0],
                                 right=[0,0,0],
                                 top=[0,0,0],
                                 bottom=[0,0,0],
                                 front=[0,0,0],
                                 back=[0,0,0],
                                 location_of_constraint=Data())
      def doInitialization(self,t):
          x=self.domain.getX()
          self.location_of_constraint=Vector(0,x.getFunctionSpace())
          if self.domain.getDim()==3:
             left_mask=(x[0]-inf(x[0])).whereZero()
             if self.left[0]: self.location_of_constraint+=left_mask*[1.,0.,0.]
             if self.left[1]: self.location_of_constraint+=left_mask*[0.,1.,0.]
             if self.left[2]: self.location_of_constraint+=left_mask*[0.,0.,1.]
             right_mask=(x[0]-inf(x[0])).whereZero()
             if self.right[0]: self.location_of_constraint+=right_mask*[1.,0.,0.]
             if self.right[1]: self.location_of_constraint+=right_mask*[0.,1.,0.]
             if self.right[2]: self.location_of_constraint+=right_mask*[0.,0.,1.]
             front_mask=(x[1]-inf(x[1])).whereZero()
             if self.front[0]: self.location_of_constraint+=front_mask*[1.,0.,0.]
             if self.front[1]: self.location_of_constraint+=front_mask*[0.,1.,0.]
             if self.front[2]: self.location_of_constraint+=front_mask*[0.,0.,1.]
             back_mask=(x[1]-sup(x[1])).whereZero()
             if self.back[0]: self.location_of_constraint+=back_mask*[1.,0.,0.]
             if self.back[1]: self.location_of_constraint+=back_mask*[0.,1.,0.]
             if self.back[2]: self.location_of_constraint+=back_mask*[0.,0.,1.]
             bottom_mask=(x[2]-inf(x[2])).whereZero()
             if self.bottom[0]: self.location_of_constraint+=bottom_mask*[1.,0.,0.]
             if self.bottom[1]: self.location_of_constraint+=bottom_mask*[0.,1.,0.]
             if self.bottom[2]: self.location_of_constraint+=bottom_mask*[0.,0.,1.]
             top_mask=(x[2]-sup(x[2])).whereZero()
             if self.top[0]: self.location_of_constraint+=top_mask*[1.,0.,0.]
             if self.top[1]: self.location_of_constraint+=top_mask*[0.,1.,0.]
             if self.top[2]: self.location_of_constraint+=top_mask*[0.,0.,1.]
          else:
             left_mask=(x[0]-inf(x[0])).whereZero()
             if self.left[0]: self.location_of_constraint+=left_mask*[1.,0.]
             if self.left[1]: self.location_of_constraint+=left_mask*[0.,1.]
             right_mask=(x[0]-inf(x[0])).whereZero()
             if self.right[0]: self.location_of_constraint+=right_mask*[1.,0.]
             if self.right[1]: self.location_of_constraint+=right_mask*[0.,1.]
             bottom_mask=(x[1]-inf(x[1])).whereZero()
             if self.bottom[0]: self.location_of_constraint+=bottom_mask*[1.,0.]
             if self.bottom[1]: self.location_of_constraint+=bottom_mask*[0.,1.]
             top_mask=(x[1]-sup(x[1])).whereZero()
             if self.top[0]: self.location_of_constraint+=top_mask*[1.,0.]
             if self.top[1]: self.location_of_constraint+=top_mask*[0.,1.]
