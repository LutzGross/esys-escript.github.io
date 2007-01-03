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
       reads finley mesh file.

       @ivar source: file name of the finley input file
       @type source: C{str}
       @ivar intergrationOrder: integration order, default -1 (in).
       @type intergrationOrder: C{int}
       """
       def __init__(self,**kwargs):
          """
          initializes the object
          """
          super(FinleyReader,self).__init__(**kwargs)
          self.declareParameter(source="none",
                                 integrationOrder=-1)
          self.__domain=None

       def domain(self):
          """
          returns the domain

          @return: the domain
          @rtype: L{Domain}
          """
          if self.__domain == None:
             self.__domain=finley.ReadMesh(self.source.getLocalFileName(),self.integrationOrder) 
             self.trace("mesh read from %s"%self.source)           
          return self.__domain
          
                       
class RectangularDomain(ParameterSet):
       """
       Generates a mesh over a rectangular domain finley.

       @ivar dim: spatial dimension, default =2 (in).
       @type dim: spatial dimension
       @ivar l: spatial lengths, default [1.,1.,1.] (in).
       @type l: C{list} of C{floats}s
       @ivar n: number of elements, default [10,10,10] (in).
       @type n: C{list} of C{int}s
       @ivar order: element order, default 1 (in).
       @type order: C{int}
       @ivar periodic: flags for periodicity, default [False,False,False] (in).
       @type periodic: C{list} of C{bool}s
       @ivar intergrationOrder: integration order, default -1 (in).
       @type intergrationOrder: C{int}
       """
       def __init__(self,**kwargs):
           """
           initializes the object
           """
           super(RectangularDomain,self).__init__(**kwargs)
           self.declareParameter(dim=2,\
                                 l=[1.,1.,1.],\
                                 n=[10,10,10], \
				 order=1,\
                                 periodic=[False,False,False],
                                 integrationOrder=-1)
           self.__domain=None

       def domain(self):
           """
           returns the domain

           @return: the domain
           @rtype: L{Domain}
           """
           if self.__domain==None:
              if self.dim==2:
                   self.__domain=finley.Rectangle(n0=self.n[0],\
                                                n1=self.n[1],\
                                                l0=self.l[0],\
                                                l1=self.l[1],\
                                                order=self.order, \
                                                periodic0=self.periodic[0], \
                                                periodic1=self.periodic[1], \
                                                integrationOrder=self.integrationOrder)
              else:
                   self.__domain=finley.Brick(n0=self.n[0],\
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
           return self.__domain

class UpdateGeometry(Model):
      """
      applies a displacement field to a domain
      
      @ivar displacement: displacements applied to the original mesh coordinates (in).
      @type displacement: L{escript.Vector}
      @ivar domain: domain
      @type domain: L{escript.Domain}
      """
      def __init__(self,**kwargs):
           """
           set-up the object
           """
           super(UpdateGeometry, self).__init__(**kwargs)
           self.declareParameter(domain=None,\
                                 displacement=None)


      def doInitialization(self):
         """
         initialize model
         """
         self.__x=self.domain.getX()
         self.__reset=True
         
      def doStepPreprocessing(self,dt):
         """
         applies the current L{displacement} to mesh nodes if required.
         """
         if self.__reset:
            self.trace("mesh nodes updated.")
            self.domain.setX(self.__x+self.displacement)
         self.__reset=False

      def doStep(self,dt):
         """
         applies the current L{displacement} to mesh nodes. 
         """
         self.trace("mesh nodes updated.")
         self.domain.setX(self.__x+self.displacement)
         self.__reset=True

      def doStepPostprocessing(self,dt):
         """
         marks nodes as beeing updated.
         """
         self.__reset=False

class ScalarConstrainerOverBox(Model):
      """
      Creates a characteristic function for the location of constraints 
      for a scalar value and selects the value from an initial value 
      ate these locations.

      In the case that the spatial dimension is two, the arguments front and back are ignored.

      @ivar domain: domain (in).
      @ivar left:  True to set a constraint at the left face of the domain (x[0]=min x[0]), default False (in).
      @ivar right: True to set a constraint at the left face of the domain (x[0]=max x[0]), default False (in).
      @ivar top: True to set a constraint at the left face of the domain (x[1]=min x[1]), default False (in).
      @ivar bottom: True to set a constraint at the left face of the domain (x[1]=max x[1]), default False (in).
      @ivar front: True to set a constraint at the left face of the domain (x[2]=min x[2]), default False (in).
      @ivar back: True to set a constraint at the left face of the domain (x[2]=max x[2]), default False (in).
      @ivar tol: absolute tolerance for "x=max x" condition, default 1.e-8 (in).
      """
      def __init__(self,**kwargs):
           super(ScalarConstrainerOverBox, self).__init__(**kwargs)
           self.declareParameter(domain=None, \
                                 value=None,  \
                                 left=False, \
                                 right=False, \
                                 top=False, \
                                 bottom=False, \
                                 front=False, \
                                 back=False, \
                                 tol=1.e-8)
           self.__value_of_constraint = None
           self.__location_of_constraint=None
      def location_of_constraint(self):
          """
          return the values used to constrain a solution

          @return: the mask marking the locations of the constraints
          @rtype: L{escript.Scalar}
          """
          if self.__location_of_constraint == None: self.__setOutput()
          return self.__location_of_constraint
         
      def value_of_constraint(self):
          """
          return the values used to constrain a solution

          @return: values to be used at the locations of the constraints. If 
                  L{value} is not given C{None} is rerturned.
          @rtype: L{escript.Scalar}
          """
          if self.__location_of_constraint == None: self.__setOutput()
          return self.__value_of_constraint
         
      def __setOutput(self):
          x=self.domain.getX()
          self.__location_of_constraint=Scalar(0,x.getFunctionSpace())
          if self.domain.getDim()==3:
                x0,x1,x2=x[0],x[1],x[2]
                if self.left: self.__location_of_constraint+=whereZero(x0-inf(x0),self.tol)
                if self.right: self.__location_of_constraint+=whereZero(x0-sup(x0),self.tol)
                if self.front: self.__location_of_constraint+=whereZero(x1-inf(x1),self.tol)
                if self.back: self.__location_of_constraint+=whereZero(x1-sup(x1),self.tol)
                if self.bottom: self.__location_of_constraint+=whereZero(x2-inf(x2),self.tol)
                if self.top: self.__location_of_constraint+=whereZero(x2-sup(x2),self.tol)
          else:
                x0,x1=x[0],x[1]
                if self.left: self.__location_of_constraint+=whereZero(x0-inf(x0),self.tol)
                if self.right: self.__location_of_constraint+=whereZero(x0-sup(x0),self.tol)
                if self.bottom: self.__location_of_constraint+=whereZero(x1-inf(x1),self.tol)
                if self.top: self.__location_of_constraint+=whereZero(x1-sup(x1),self.tol)
          if self.value:
              self.__value_of_constraint=self.__location_of_constraint*self.value

class VectorConstrainerOverBox(Model):
      """
      Creates a characteristic function for the location of constraints vector value.
      In the case that the spatial dimension is two, the arguments front and
      back as well as the third component of each argument is ignored.

      @ivar domain: domain
      @ivar left: list of three boolean. left[i]==True sets a constraint for the i-th component at the left face of the domain (x[0]=min x[0]),
                       default [False,False,False] (in).
      @ivar right: list of three boolean. left[i]==True sets a constraint for the i-th component at the right face of the domain (x[0]=max x[0]), 
                default [False,False,False] (in).
      @ivar top: list of three boolean. left[i]==True sets a constraint for the i-th component at the top face of the domain (x[1]=min x[1]), 
                default [False,False,False] (in).
      @ivar bottom: list of three boolean. left[i]==True sets a constraint for the i-th component at the bottom face of the domain (x[1]=min x[1]), 
                default [False,False,False] (in).
      @ivar front: list of three boolean. left[i]==True sets a constraint for the i-th component at the front face of the domain (x[2]=min x[2]), 
                default [False,False,False] (in).
      @ivar back: list of three boolean. left[i]==True sets a constraint for the i-th component at the back face of the domain (x[2]=max x[2]), 
                default [False,False,False] (in).
      @ivar tol: absolute tolerance for "x=max x" condition, default 1.e-8 (in).
      """
      def __init__(self, **kwargs):
           super(VectorConstrainerOverBox, self).__init__(**kwargs)
           self.declareParameter(domain=None, \
                                 value=None,  \
                                 left=[0,0,0],  \
                                 right=[0,0,0],  \
                                 top=[0,0,0],  \
                                 bottom=[0,0,0],  \
                                 front=[0,0,0], \
                                 back=[0,0,0], \
                                 tol=1.e-8)
           self.__value_of_constraint = None
           self.__location_of_constraint=None

      def location_of_constraint(self):
          """
          return the values used to constrain a solution

          @return: the mask marking the locations of the constraints
          @rtype: L{escript.Vector}
          """
          if self.__location_of_constraint == None: self.__setOutput()
          return self.__location_of_constraint
         
      def value_of_constraint(self):
          """
          return the values used to constrain a solution

          @return: values to be used at the locations of the constraints. If 
                  L{value} is not given C{None} is rerturned.
          @rtype: L{escript.Vector}
          """
          if self.__location_of_constraint == None: self.__setOutput()
          return self.__value_of_constraint
         
      def __setOutput(self):
          x=self.domain.getX()
          self.__location_of_constraint=Vector(0,x.getFunctionSpace())
          if self.domain.getDim()==3:
             x0,x1,x2=x[0],x[1],x[2]
             left_mask=whereZero(x0-inf(x0),self.tol)
             if self.left[0]: self.__location_of_constraint+=left_mask*[1.,0.,0.]
             if self.left[1]: self.__location_of_constraint+=left_mask*[0.,1.,0.]
             if self.left[2]: self.__location_of_constraint+=left_mask*[0.,0.,1.]
             right_mask=whereZero(x0-sup(x0),self.tol)
             if self.right[0]: self.__location_of_constraint+=right_mask*[1.,0.,0.]
             if self.right[1]: self.__location_of_constraint+=right_mask*[0.,1.,0.]
             if self.right[2]: self.__location_of_constraint+=right_mask*[0.,0.,1.]
             front_mask=whereZero(x1-inf(x1),self.tol)
             if self.front[0]: self.__location_of_constraint+=front_mask*[1.,0.,0.]
             if self.front[1]: self.__location_of_constraint+=front_mask*[0.,1.,0.]
             if self.front[2]: self.__location_of_constraint+=front_mask*[0.,0.,1.]
             back_mask=whereZero(x1-sup(x1),self.tol)
             if self.back[0]: self.__location_of_constraint+=back_mask*[1.,0.,0.]
             if self.back[1]: self.__location_of_constraint+=back_mask*[0.,1.,0.]
             if self.back[2]: self.__location_of_constraint+=back_mask*[0.,0.,1.]
             bottom_mask=whereZero(x2-inf(x2),self.tol)
             if self.bottom[0]: self.__location_of_constraint+=bottom_mask*[1.,0.,0.]
             if self.bottom[1]: self.__location_of_constraint+=bottom_mask*[0.,1.,0.]
             if self.bottom[2]: self.__location_of_constraint+=bottom_mask*[0.,0.,1.]
             top_mask=whereZero(x2-sup(x2),self.tol)
             if self.top[0]: self.__location_of_constraint+=top_mask*[1.,0.,0.]
             if self.top[1]: self.__location_of_constraint+=top_mask*[0.,1.,0.]
             if self.top[2]: self.__location_of_constraint+=top_mask*[0.,0.,1.]
             if self.value:
                self.__value_of_constraint=self.__location_of_constraint*self.value
          else:
             x0,x1=x[0],x[1]
             left_mask=whereZero(x0-inf(x0),self.tol)
             if self.left[0]: self.__location_of_constraint+=left_mask*[1.,0.]
             if self.left[1]: self.__location_of_constraint+=left_mask*[0.,1.]
             right_mask=whereZero(x0-sup(x0),self.tol)
             if self.right[0]: self.__location_of_constraint+=right_mask*[1.,0.]
             if self.right[1]: self.__location_of_constraint+=right_mask*[0.,1.]
             bottom_mask=whereZero(x1-inf(x1),self.tol)
             if self.bottom[0]: self.__location_of_constraint+=bottom_mask*[1.,0.]
             if self.bottom[1]: self.__location_of_constraint+=bottom_mask*[0.,1.]
             top_mask=whereZero(x1-sup(x1),self.tol)
             if self.top[0]: self.__location_of_constraint+=top_mask*[1.,0.]
             if self.top[1]: self.__location_of_constraint+=top_mask*[0.,1.]
             if self.value:
                self.__value_of_constraint=self.__location_of_constraint*self.value[:2]

# vim: expandtab shiftwidth=4:
