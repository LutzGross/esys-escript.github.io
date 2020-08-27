
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

from __future__ import division, print_function

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import esys.escript.modelframe as esmf
import esys.escript as es

try:
    import esys.finley
except ImportError:
    raise ImportError("Finley module not available")

class DomainReader(esmf.ParameterSet):
#class DudleyReader(ParameterSet):
       """
       """
       def __init__(self,domainmodule=None, **kwargs):
          """
          initializes the object
          """
          super(DomainReader,self).__init__(**kwargs)
          self.declareParameter(source="none",
                                dim=None,
                                optimizeLabeling=True,
                                reducedIntegrationOrder=-1,
                                integrationOrder=-1)
          if domainmodule==None:
              domainmodule=esys.finley
          self.__domainModule=domainmodule
          self.__domain=None


       def domain(self):
          """
          returns the domain

          :return: the domain
          :rtype: `Domain`
          """
          if self.__domain is None:
             if  self.source.fileformat == "fly":
                self.__domain=self.__domainModule.ReadMesh(self.source.getLocalFileName(),self.integrationOrder)
             elif self.source.fileformat == "gmsh":
                if self.dim==None:
                   dim=3
                else:
                   dim=self.dim
                self.__domain=self.__domainModule.ReadGmsh(self.source.getLocalFileName(),dim,self.integrationOrder,self.reducedIntegrationOrder, self.optimizeLabeling)
             else:
                raise TypeError("unknown mesh file format %s."%self.source.fileformat)
             self.trace("mesh read from %s in %s format."%(self.source.getLocalFileName(), self.source.fileformat))           
          return self.__domain
                
class FinleyReader(DomainReader):
        def __init__(self, **kw):
            DomainReader.__init__(self, esys.finley, **kw)

class RectangularDomain(esmf.ParameterSet):
       """
       Generates a mesh over a rectangular domain.

       :ivar dim: spatial dimension, default =2 (in).
       :type dim: spatial dimension
       :ivar l: spatial lengths, default [1.,1.,1.] (in).
       :type l: ``list`` of ``float``
       :ivar n: number of elements, default [10,10,10] (in).
       :type n: ``list`` of ``int``
       :ivar order: element order, default 1 (in).
       :type order: ``int``
       :ivar periodic: flags for periodicity, default [False,False,False] (in).
       :type periodic: ``list`` of ``bool``
       :ivar intergrationOrder: integration order, default -1 (in).
       :type intergrationOrder: ``int``
       """
       def __init__(self,domainmodule=None,**kwargs):
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
           self.__domainModule=domainmodule
           if self.__domainModule==None:
                self.__domainModule=esys.finley

       def domain(self):
           """
           returns the domain

           :return: the domain
           :rtype: `Domain`
           """
           if self.__domain==None:
              if self.dim==2:
                   self.__domain=self.__domainModule.Rectangle(n0=self.n[0],\
                                                n1=self.n[2],\
                                                l0=self.l[0],\
                                                l1=self.l[2],\
                                                order=self.order, \
                                                periodic0=self.periodic[0], \
                                                periodic1=self.periodic[2], \
                                                integrationOrder=self.integrationOrder)
              else:
                   self.__domain=self__domainModule.Brick(n0=self.n[0],\
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

class UpdateGeometry(esmf.Model):
      """
      applies a displacement field to a domain
      
      :note: Instance variable displacement - displacements applied to the original mesh coordinates (in).
      :note: Instance variable displacement - `escript.Vector`
      :note: Instance variable domain - domain
      :note: Instance variable domain - `escript.Domain`
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
         applies the current `displacement` to mesh nodes if required.
         """
         if self.__reset:
            self.trace("mesh nodes updated.")
            self.domain.setX(self.__x+self.displacement)
         self.__reset=False

      def doStep(self,dt):
         """
         applies the current `displacement` to mesh nodes. 
         """
         self.trace("mesh nodes updated.")
         self.domain.setX(self.__x+self.displacement)
         self.__reset=True

      def doStepPostprocessing(self,dt):
         """
         marks nodes as beeing updated.
         """
         self.__reset=False

class ConstrainerOverBox(esmf.Model):
      """
      Creates a characteristic function for the location of constraints 
      for all components of a value and selects the value from an initial value 
      ate these locations.

      In the case that the spatial dimension is two, the arguments front and back are ignored.

      :note: Instance variable - domain (in).
      :note: Instance variable left -  True to set a constraint at the left face of the domain (x[0]=min x[0]), default False (in).
      :note: Instance variable right - True to set a constraint at the left face of the domain (x[0]=max x[0]), default False (in).
      :note: Instance variable top - True to set a constraint at the left face of the domain (x[1]=min x[1]), default False (in).
      :note: Instance variable bottom - True to set a constraint at the left face of the domain (x[1]=max x[1]), default False (in).
      :note: Instance variable front - True to set a constraint at the left face of the domain (x[2]=min x[2]), default False (in).
      :note: Instance variable back - True to set a constraint at the left face of the domain (x[2]=max x[2]), default False (in).
      :note: Instance variable tol - absolute tolerance for "x=max x" condition, default 1.e-8 (in).
      """
      def __init__(self,**kwargs):
           super(ConstrainerOverBox, self).__init__(**kwargs)
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

          :return: the mask marking the locations of the constraints
          :rtype: `escript.Scalar`
          """
          if self.__location_of_constraint is None: self.__setOutput()
          return self.__location_of_constraint
         
      def value_of_constraint(self):
          """
          return the values used to constrain a solution

          :return: values to be used at the locations of the constraints. If 
                  ``value`` is not given ``None`` is rerturned.
          :rtype: `escript.Scalar`
          """
          if self.__location_of_constraint is None: self.__setOutput()
          return self.__value_of_constraint
         
      def __setOutput(self):
          if self.__location_of_constraint is None:
             x=self.domain.getX()
             val=self.value
             if isinstance(val, int) or isinstance(val, float):
                shape=()
             elif isinstance(val, list) or isinstance(val, tuple) :
                shape=(len(val),)
             elif isinstance(val, numpy.ndarray):
                 shape=val.shape
             elif val is None:
                  shape=()
             else: 
                 shape=val.getShape()
             self.__location_of_constraint=Data(0,shape,x.getFunctionSpace())
             if self.domain.getDim()==3:
                   x0,x1,x2=x[0],x[1],x[2]
                   if self.left: self.__location_of_constraint+=es.whereZero(x0-es.inf(x0),self.tol)
                   if self.right: self.__location_of_constraint+=es.whereZero(x0-es.sup(x0),self.tol)
                   if self.front: self.__location_of_constraint+=es.whereZero(x1-es.inf(x1),self.tol)
                   if self.back: self.__location_of_constraint+=es.whereZero(x1-es.sup(x1),self.tol)
                   if self.bottom: self.__location_of_constraint+=es.whereZero(x2-es.inf(x2),self.tol)
                   if self.top: self.__location_of_constraint+=es.whereZero(x2-es.sup(x2),self.tol)
             else:
                   x0,x1=x[0],x[1]
                   if self.left: self.__location_of_constraint+=es.whereZero(x0-es.inf(x0),self.tol)
                   if self.right: self.__location_of_constraint+=es.whereZero(x0-es.sup(x0),self.tol)
                   if self.bottom: self.__location_of_constraint+=es.whereZero(x1-es.inf(x1),self.tol)
                   if self.top: self.__location_of_constraint+=es.whereZero(x1-es.sup(x1),self.tol)
             if not self.value is None:
                   self.__value_of_constraint=self.__location_of_constraint*self.value
class ScalarConstrainerOverBox(esmf.Model):
      """
      Creates a characteristic function for the location of constraints 
      for a scalar value and selects the value from an initial value 
      ate these locations.

      In the case that the spatial dimension is two, the arguments front and back are ignored.

      :note: Instance variable domain - domain (in).
      :note: Instance variable left -  True to set a constraint at the left face of the domain (x[0]=min x[0]), default False (in).
      :note: Instance variable right - True to set a constraint at the left face of the domain (x[0]=max x[0]), default False (in).
      :note: Instance variable top - True to set a constraint at the left face of the domain (x[1]=min x[1]), default False (in).
      :note: Instance variable bottom - True to set a constraint at the left face of the domain (x[1]=max x[1]), default False (in).
      :note: Instance variable front - True to set a constraint at the left face of the domain (x[2]=min x[2]), default False (in).
      :note: Instance variable back - True to set a constraint at the left face of the domain (x[2]=max x[2]), default False (in).
      :note: Instance variable tol - absolute tolerance for "x=max x" condition, default 1.e-8 (in).
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

          :return: the mask marking the locations of the constraints
          :rtype: `escript.Scalar`
          """
          if self.__location_of_constraint is None: self.__setOutput()
          return self.__location_of_constraint
         
      def value_of_constraint(self):
          """
          return the values used to constrain a solution

          :return: values to be used at the locations of the constraints. If 
                  ``value`` is not given ``None`` is rerturned.
          :rtype: `escript.Scalar`
          """
          if self.__location_of_constraint is None: self.__setOutput()
          return self.__value_of_constraint
         
      def __setOutput(self):
          x=self.domain.getX()
          self.__location_of_constraint=es.Scalar(0,x.getFunctionSpace())
          if self.domain.getDim()==3:
                x0,x1,x2=x[0],x[1],x[2]
                d=max(es.sup(x0)-es.inf(x0), sup(x1)-es.inf(x1), sup(x2)-es.inf(x2))
                if self.left: self.__location_of_constraint+=es.whereZero(x0-es.inf(x0),self.tol*d)
                if self.right: self.__location_of_constraint+=es.whereZero(x0-es.sup(x0),self.tol*d)
                if self.front: self.__location_of_constraint+=es.whereZero(x1-es.inf(x1),self.tol*d)
                if self.back: self.__location_of_constraint+=es.whereZero(x1-es.sup(x1),self.tol*d)
                if self.bottom: self.__location_of_constraint+=es.whereZero(x2-es.inf(x2),self.tol*d)
                if self.top: self.__location_of_constraint+=es.whereZero(x2-es.sup(x2),self.tol*d)
          else:
                x0,x1=x[0],x[1]
                d=max(es.sup(x0)-es.inf(x0), es.sup(x1)-es.inf(x1))
                if self.left: self.__location_of_constraint+=es.whereZero(x0-es.inf(x0),self.tol*d)
                if self.right: self.__location_of_constraint+=es.whereZero(x0-es.sup(x0),self.tol*d)
                if self.bottom: self.__location_of_constraint+=es.whereZero(x1-es.inf(x1),self.tol*d)
                if self.top: self.__location_of_constraint+=es.whereZero(x1-es.sup(x1),self.tol*d)
          if not self.value is None:
              self.__value_of_constraint=self.__location_of_constraint*self.value

class VectorConstrainerOverBox(esmf.Model):
      """
      Creates a characteristic function for the location of constraints vector value.
      In the case that the spatial dimension is two, the arguments front and
      back as well as the third component of each argument is ignored.

      :note: Instance variable domain
      :note: Instance variable left - list of three boolean. left[i]==True sets a constraint for the i-th component at the left face of the domain (x[0]=min x[0]),
                       default [False,False,False] (in).
      :note: Instance variable right - list of three boolean. left[i]==True sets a constraint for the i-th component at the right face of the domain (x[0]=max x[0]), 
                default [False,False,False] (in).
      :note: Instance variable top - list of three boolean. left[i]==True sets a constraint for the i-th component at the top face of the domain (x[1]=min x[1]), 
                default [False,False,False] (in).
      :note: Instance variable bottom - list of three boolean. left[i]==True sets a constraint for the i-th component at the bottom face of the domain (x[1]=min x[1]), 
                default [False,False,False] (in).
      :note: Instance variable front - list of three boolean. left[i]==True sets a constraint for the i-th component at the front face of the domain (x[2]=min x[2]), 
                default [False,False,False] (in).
      :note: Instance variable back - list of three boolean. left[i]==True sets a constraint for the i-th component at the back face of the domain (x[2]=max x[2]), 
                default [False,False,False] (in).
      :note: Instance variable tol - absolute tolerance for "x=max x" condition, default 1.e-8 (in).
      """
      def __init__(self, **kwargs):
           super(VectorConstrainerOverBox, self).__init__(**kwargs)
           self.declareParameter(domain=None, \
                                 value=None,  \
                                 left=[False ,False ,False ],  \
                                 right=[False ,False ,False ],  \
                                 top=[False ,False ,False ],  \
                                 bottom=[False ,False ,False ],  \
                                 front=[False ,False ,False ], \
                                 back=[False ,False ,False ], \
                                 tol=1.e-8)
           self.__value_of_constraint = None
           self.__location_of_constraint=None

      def location_of_constraint(self):
          """
          return the values used to constrain a solution

          :return: the mask marking the locations of the constraints
          :rtype: `escript.Vector`
          """
          if self.__location_of_constraint is None: self.__setOutput()
          return self.__location_of_constraint
         
      def value_of_constraint(self):
          """
          return the values used to constrain a solution

          :return: values to be used at the locations of the constraints. If 
                  ``value`` is not given ``None`` is rerturned.
          :rtype: `escript.Vector`
          """
          if self.__location_of_constraint is None: self.__setOutput()
          return self.__value_of_constraint
         
      def __setOutput(self):
          x=self.domain.getX()
          self.__location_of_constraint=es.Vector(0,x.getFunctionSpace())
          if self.domain.getDim()==3:
             x0,x1,x2=x[0],x[1],x[2]
             d=max(es.sup(x0)-es.inf(x0), es.sup(x1)-es.inf(x1), es.sup(x2)-es.inf(x2))
             left_mask=es.whereZero(x0-es.inf(x0),self.tol*d)
             if self.left[0]: self.__location_of_constraint+=left_mask*[1.,0.,0.]
             if self.left[1]: self.__location_of_constraint+=left_mask*[0.,1.,0.]
             if self.left[2]: self.__location_of_constraint+=left_mask*[0.,0.,1.]
             right_mask=es.whereZero(x0-es.sup(x0),self.tol*d)
             if self.right[0]: self.__location_of_constraint+=right_mask*[1.,0.,0.]
             if self.right[1]: self.__location_of_constraint+=right_mask*[0.,1.,0.]
             if self.right[2]: self.__location_of_constraint+=right_mask*[0.,0.,1.]
             front_mask=es.whereZero(x1-es.inf(x1),self.tol*d)
             if self.front[0]: self.__location_of_constraint+=front_mask*[1.,0.,0.]
             if self.front[1]: self.__location_of_constraint+=front_mask*[0.,1.,0.]
             if self.front[2]: self.__location_of_constraint+=front_mask*[0.,0.,1.]
             back_mask=es.whereZero(x1-es.sup(x1),self.tol*d)
             if self.back[0]: self.__location_of_constraint+=back_mask*[1.,0.,0.]
             if self.back[1]: self.__location_of_constraint+=back_mask*[0.,1.,0.]
             if self.back[2]: self.__location_of_constraint+=back_mask*[0.,0.,1.]
             bottom_mask=es.whereZero(x2-es.inf(x2),self.tol*d)
             if self.bottom[0]: self.__location_of_constraint+=bottom_mask*[1.,0.,0.]
             if self.bottom[1]: self.__location_of_constraint+=bottom_mask*[0.,1.,0.]
             if self.bottom[2]: self.__location_of_constraint+=bottom_mask*[0.,0.,1.]
             top_mask=es.whereZero(x2-es.sup(x2),self.tol*d)
             if self.top[0]: self.__location_of_constraint+=top_mask*[1.,0.,0.]
             if self.top[1]: self.__location_of_constraint+=top_mask*[0.,1.,0.]
             if self.top[2]: self.__location_of_constraint+=top_mask*[0.,0.,1.]
             if not self.value is None:
                self.__value_of_constraint=self.__location_of_constraint*self.value
          else:
             x0,x1=x[0],x[1]
             d=max(es.sup(x0)-es.inf(x0), es.sup(x1)-es.inf(x1))
             left_mask=es.whereZero(x0-es.inf(x0),self.tol*d)
             if self.left[0]: self.__location_of_constraint+=left_mask*[1.,0.]
             if self.left[1]: self.__location_of_constraint+=left_mask*[0.,1.]
             right_mask=es.whereZero(x0-es.sup(x0),self.tol*d)
             if self.right[0]: self.__location_of_constraint+=right_mask*[1.,0.]
             if self.right[1]: self.__location_of_constraint+=right_mask*[0.,1.]
             bottom_mask=es.whereZero(x1-es.inf(x1),self.tol*d)
             if self.bottom[0]: self.__location_of_constraint+=bottom_mask*[1.,0.]
             if self.bottom[1]: self.__location_of_constraint+=bottom_mask*[0.,1.]
             top_mask=es.whereZero(x1-es.sup(x1),self.tol*d)
             if self.top[0]: self.__location_of_constraint+=top_mask*[1.,0.]
             if self.top[1]: self.__location_of_constraint+=top_mask*[0.,1.]
             if not self.value is None:
                self.__value_of_constraint=self.__location_of_constraint*self.value[:2]

class ConstrainerAtBoxVertex(esmf.Model):
      """
      Creates a characteristic function for the location of constraints 
      for all components of a value and selects the value from an initial value 
      ate these locations.

      In the case that the spatial dimension is two, the arguments front and back are ignored.
      
      :note: Instance variable domain
      :note: Instance variable tol - absolute tolerance for "x=left, front, bottom vertex" condition, default 1.e-8 (in).
      """
      def __init__(self,**kwargs):
           super(ConstrainerAtBoxVertex, self).__init__(**kwargs)
           self.declareParameter(domain=None, \
                                 value=None,  \
                                 tol=1.e-8)
           self.__value_of_constraint = None
           self.__location_of_constraint=None
      def location_of_constraint(self):
          """
          return the values used to constrain a solution

          :return: the mask marking the locations of the constraints
          :rtype: `escript.Scalar`
          """
          if self.__location_of_constraint is None: self.__setOutput()
          return self.__location_of_constraint
         
      def value_of_constraint(self):
          """
          return the values used to constrain a solution

          :return: values to be used at the locations of the constraints. If 
                  ``value`` is not given ``None`` is rerturned.
          :rtype: `escript.Scalar`
          """
          if self.__location_of_constraint is None: self.__setOutput()
          return self.__value_of_constraint
         
      def __setOutput(self):
          if self.__location_of_constraint is None:
             x=self.domain.getX()
             val=self.value
             if isinstance(val, int) or isinstance(val, float):
                shape=()
             elif isinstance(val, list) or isinstance(val, tuple) :
                shape=(len(val),)
             elif isinstance(val, numpy.ndarray):
                 shape=val.shape
             elif val is None:
                  shape=()
             else: 
                 shape=val.getShape()
             if self.domain.getDim()==3:
                   vertex=[es.inf(x[0]),es.inf(x[1]),es.inf(x[2])]
             else:
                   vertex=[es.inf(x[0]),es.inf(x[1])]
             self.__location_of_constraint=es.whereZero(es.length(x-vertex),self.tol)*numpy.ones(shape)
             if not self.value is None:
                   self.__value_of_constraint=self.__location_of_constraint*self.value

class ScalarConstrainerAtBoxVertex(esmf.Model):
      """
      Creates a characteristic function for the location of constraints 
      for a scalar value and selects the value from an initial value 
      ate these locations.

      In the case that the spatial dimension is two, the arguments front and back are ignored.

      :note: Instance variable domain
      :note: Instance variable tol - absolute tolerance for "x=left, front, bottom vertex" condition, default 1.e-8 (in).
      """
      def __init__(self,**kwargs):
           super(ScalarConstrainerAtBoxVertex, self).__init__(**kwargs)
           self.declareParameter(domain=None, \
                                 value=None,  \
                                 tol=1.e-8)
           self.__value_of_constraint = None
           self.__location_of_constraint=None
      def location_of_constraint(self):
          """
          return the values used to constrain a solution

          :return: the mask marking the locations of the constraints
          :rtype: `escript.Scalar`
          """
          if self.__location_of_constraint is None: self.__setOutput()
          return self.__location_of_constraint
         
      def value_of_constraint(self):
          """
          return the values used to constrain a solution

          :return: values to be used at the locations of the constraints. If 
                  ``value`` is not given ``None`` is rerturned.
          :rtype: `escript.Scalar`
          """
          if self.__location_of_constraint is None: self.__setOutput()
          return self.__value_of_constraint
         
      def __setOutput(self):
          x=self.domain.getX()
          self.__location_of_constraint=es.Scalar(0,x.getFunctionSpace())
          if self.domain.getDim()==3:
                   vertex=[es.inf(x[0]),es.inf(x[1]),es.inf(x[2])]
          else:
                 vertex=[es.inf(x[0]),es.inf(x[1])]
          self.__location_of_constraint=es.whereZero(es.length(x-vertex),self.tol)
          if not self.value is None:
              self.__value_of_constraint=self.__location_of_constraint*self.value

class VectorConstrainerAtBoxVertex(esmf.Model):
      """
      Creates a characteristic function for the location of constraints vector value.
      In the case that the spatial dimension is two, the arguments front and
      back as well as the third component of each argument is ignored.

      :note: Instance variable domain
      :note: Instance variable comp_mask - list of three boolean. comp_mask[i]==True sets a constraint for the i-th component at the left, front, bottom vertex, default [False,False,False] (in).
      :note: Instance variable tol - absolute tolerance for "x=left, front, bottom vertex" condition, default 1.e-8 (in).
      """
      def __init__(self, **kwargs):
           super(VectorConstrainerAtBoxVertex, self).__init__(**kwargs)
           self.declareParameter(domain=None, \
                                 value=None,  \
                                 comp_mask=[False, False, False], 
                                 tol=1.e-8)
           self.__value_of_constraint = None
           self.__location_of_constraint=None

      def location_of_constraint(self):
          """
          return the values used to constrain a solution

          :return: the mask marking the locations of the constraints
          :rtype: `escript.Vector`
          """
          if self.__location_of_constraint is None: self.__setOutput()
          return self.__location_of_constraint
         
      def value_of_constraint(self):
          """
          return the values used to constrain a solution

          :return: values to be used at the locations of the constraints. If 
                  ``value`` is not given ``None`` is rerturned.
          :rtype: `escript.Vector`
          """
          if self.__location_of_constraint is None: self.__setOutput()
          return self.__value_of_constraint
         
      def __setOutput(self):
          x=self.domain.getX()
          self.__location_of_constraint=es.Vector(0,x.getFunctionSpace())
          if self.domain.getDim()==3:
             vertex=[es.inf(x[0]),es.inf(x[1]),es.inf(x[2])]
             msk=numpy.zeros((3,))
             if self.comp_mask[0]: msk[0]=1
             if self.comp_mask[1]: msk[1]=1
             if self.comp_mask[2]: msk[2]=1
          else:
             vertex=[es.inf(x[0]),es.inf(x[1])]
             msk=numpy.zeros((2,))
             if self.comp_mask[0]: msk[0]=1
             if self.comp_mask[1]: msk[1]=1
          self.__location_of_constraint=es.whereZero(es.length(x-vertex),self.tol)*numpy.ones(shape)
          if not self.value is None:
                self.__value_of_constraint=self.__location_of_constraint*self.value

