########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
import numpy
import math

class FaultSystem:
  """
  The FaultSystem class defines a system of faults in the Earth crust. 
 
  A fault system is defined by set of faults index by a tag. Each fault is defined as a polygon describing the top of
  the fault and - in case of a 3D model - a polygon describing the bottom of the fault. This class provides a mechanism
  to parametrise each fault with the domain [0,w0_max] x [0, w1_max]  (to [0,w0_max] in the 2D case).
  """
  NOTAG="__NOTAG__"
  def __init__(self,dim=3):
    """
    Sets up the fault system

    :param dim: spatial dimension
    :type dim: ``int`` of value 2 or 3
    """
    if not (dim == 2 or dim == 3):
       raise ValueError,"only dimension2 2 and 3 are supported."
    self.__dim=dim
    self.__faults={}
    self.__length={}
    self.__depth={}
    self.__offsets={}
    self.__w1_max={}
    self.__w0_max={}
    self.__center=None
    self.__orientation = None

  def getDim(self):
     """
     returns the spatial dimension
     :rtype: ``int``
     """
     return self.__dim
  def getLength(self,tag=None):
     """
     returns the unrolled length of fault ``tag``
     :rtype: ``float``
     """
     if tag==None: tag=self.NOTAG
     return self.__length[tag]

  def getDepth(self,tag=None):
     """
     returns the medium depth of fault ``tag``
     :rtype: ``float``
     """
     if tag==None: tag=self.NOTAG
     return self.__depth[tag]

  def getTags(self):
     """
     returns a list of the tags used by the fault system
     :rtype: ``list``
     """
     return self.__faults.keys()

  def getW0Range(self,tag=None):
     """
     returns the range of the parameterization in ``w0``
     :rtype: two ``float``
     """
     return self.getW0Offsets(tag)[0], self.getW0Offsets(tag)[-1]

  def getW1Range(self,tag=None):
     """
     returns the range of the parameterization in ``w1``
     :rtype: two ``float``
     """
     if tag==None: tag=self.NOTAG
     return -self.__w1_max[tag],0

  def getW0Offsets(self, tag=None):
     """
     returns the offsets for the parametrization of fault ``tag``. 

     :return: the offsets in the parametrization
     :rtype: ``list`` of ``float``
     """
     if tag==None: tag=self.NOTAG
     return self.__offsets[tag]


  def getFaultSegments(self, tag=None):
     """
     returns the polygons used to describe fault tagged by ``tag`` 
      
     :param tag: the tag of the fault
     :type tag: ``float`` or ``str``
     :return: the list of vertices defining the top and the list of the vertices defining the bottom of the fault indexed by ``tag``. The coordinates are `numpy.ndarray'.
     """
     if tag==None: tag=self.NOTAG
     return self.__faults[tag]

  def getCenterOnSurface(self):
      """
      returns the center point of the fault system at the surface 
      :rtype: `numpy.ndarray`
      """
      if self.__center == None:
        self.__center=numpy.zeros((3,),numpy.float)
        counter=0
        for t in self.getTags():
            for s in self.getFaultSegments(t)[0]:
                self.__center[:2]+=s[:2]
                counter+=1
        self.__center/=counter
      return self.__center[:self.getDim()]

  def getOrientationOnSurface(self):
      """
      returns the orientation of the fault system in RAD on the surface around the fault system center
      :rtype: ``float``
      """
      if self.__orientation == None:
          center=self.getCenterOnSurface()
          covariant=numpy.zeros((2,2))
          for t in self.getTags():
              for s in self.getFaultSegments(t)[0]:
                covariant[0,0]+=(center[0]-s[0])**2
                covariant[0,1]+=(center[1]-s[1])*(center[0]-s[0])
                covariant[1,1]+=(center[1]-s[1])**2
                covariant[1,0]+=(center[1]-s[1])*(center[0]-s[0])
          e, V=numpy.linalg.eigh(covariant)
          if e[0]>e[1]:
             d=V[:,1]
          else:
             d=V[:,0]
          if abs(d[0])>0.:
             self.__orientation=atan(d[1]/d[0])
          else:
             self.__orientation=math.pi/2
      return self.__orientation 
  def transform(self, rot=0, shift=numpy.zeros((3,))):
     """
     applies a shift and a consecutive rotation in the x0x1 plane.
    
     :param rot: rotation angle in RAD
     :type rot: ``float``
     :param shift: shift vector to be applied before rotation
     :type shift: `numpy.ndarray` of size 2 or 3
     """
     if self.getDim() == 2:
        mat=numpy.array([[cos(rot), -sin(rot)], [sin(rot), cos(rot)] ])
     else:
        mat=numpy.array([[cos(rot), -sin(rot),0.], [sin(rot), cos(rot),0.], [0.,0.,1.] ])
     for t in self.getTags():
           top=[]
           for s in self.getFaultSegments(t)[0]: top.append(numpy.dot(mat,(s+shift[:self.getDim()])))
           if self.getDim() == 3:
              bottom=[]
              for s in self.getFaultSegments(t)[1]: bottom.append(numpy.dot(mat,(s+shift[:self.getDim()])))
           else:
              bottom=None
           self.addFault(top=top, tag=t, bottom=bottom, w0_offsets=self.getW0Offsets(t), w1_max=-self.getW1Range(t)[0])

  def addFault(self, top, tag=None, bottom=None, w0_offsets=None, w1_max=None):
     """
     adds a new fault to the fault system. The fault is named by ``tag``.

     a fault is subdivided into segments where segment ``i`` is represented by the four corner points 
     ``v0``=``top[i]``, ``v1``=``bottom[i]``, ``v2``=``bottom[i+1]`` and ``v3``==``top[i+1]``. The segment is parametrized 
     by ``w0`` and ``w1`` with ``w0_offsets[i]<=w0<=w0_offsets[i+1]`` and ``-w1_max<=w1<=0``. Moreover 
   
     - ``(w0,w1)=(w0_offsets[i]  ,       0)->v0``
     - ``(w0,w1)=(w0_offsets[i]  , -w1_max)->v1``
     - ``(w0,w1)=(w0_offsets[i+1], -w1_max)->v2``
     - ``(w0,w1)=(w0_offsets[i+1],       0)->v3``

     If no ``w0_offsets`` is given, 
  
     - ``w0_offsets[0]=0``
     - ``w0_offsets[i]=w0_offsets[i-1]+length(v0-v4)``

     If no ``w1_max`` is given, the average fault depth (=length(v0-v1)) is used.

     :param top: list of coordinates. top[i] and top[i+1] are the top corners of the i-th segment of the fault
     :type top: ``list`` of objects which can be converted to `numpy.ndarray` vector of length 3 (or 2 in the 2D case)
     :param tag: the tag of the fault. If fault ``tag`` already exists it is overwritten.
     :type tag: ``float`` or ``str``
     :param bottom: list of coordinates. bottom[i] and bottom[i+1] are the bottom corners of the i-th segment of the fault
     :type bottom: ``list`` of objects which can be converted to `numpy.ndarray` vector of length 3 (or 2 in the 2D case)
     :param w0_offsets: ``w0_offsets[i]`` defines the offset of the segment ``i`` in the fault to be used in the parametrization of the fault. If not present the cumulative length of the fault segments is used. 
     :type w0_offsets: ``list`` of ``float`` or ``None``
     :param w1_max: the maximum value used for parametrization of the fault in the depth direction. If not present the mean depth of the fault segments is used.
     :type w1_max: ``float``
     :note: In the three dimensional case the lists ``bottom`` and ``top`` must have the same length.
     """
     if tag==None: 
         tag=self.NOTAG
     else:
         if self.NOTAG in self.getTags():
              raise ValueError,'Attempt to add a fault with no tag to a set of existing faults'
     if len(top)<2:
         raise Value,"at least two vertices must be given."
     if self.getDim()==2 and not bottom==None:
           raise ValueError,'Spatial dimension two does not support bottom polygon for faults.'
     if not bottom == None:
        if len(top) != len(bottom):
           raise ValueError,'length of top and bottom polygon must be identical.'
     if w0_offsets != None:
       if len(w0_offsets) != len(top):
          raise ValueError,'expected length of w0_offsets is %s'%(len(top))
     self.__center=None
     self.__orientation = None
     # 
     #   translate to numpy
     #
     if self.getDim()==2:
        self.__faults[tag]=([ numpy.array(t[:2],numpy.double) for t in top ],)
     else:
        self.__faults[tag]=( [ numpy.array(t[:3],numpy.double) for t in top ], [ numpy.array(b[:3],numpy.double) for b in bottom ])
     faults=self.__faults[tag]
     len_faults=len(faults[0])
     # 
     #    calculate the depth
     #
     if self.getDim()==2:
        self.__depth[tag]=0.
     else:
        self.__depth[tag]=float(sum([numpy.linalg.norm(faults[1][i]-faults[0][i]) for i in xrange(len_faults)])/len_faults)
     if w1_max==None or self.getDim()==2: w1_max=self.__depth[tag]
     self.__length[tag]=float(sum([numpy.linalg.norm(faults[0][i]-faults[0][i-1]) for i in xrange(1,len_faults)]))
     self.__w1_max[tag]=w1_max
     #
     # 
     #
     if w0_offsets!=None:
        self.__offsets[tag]=w0_offsets
     else:
        self.__offsets[tag]=[0.]
        for  i in xrange(1,len_faults):
            self.__offsets[tag].append(self.__offsets[tag][-1]+float(numpy.linalg.norm(faults[0][i]-faults[0][i-1])))
     w0_max=max(self.__offsets[tag])
     self.__w0_max[tag]=w0_max

  def getMaxValue(self,f, tol=sqrt(EPSILON)):
     """
     returns the maximum value of ``f``, the fault and the location on the fault in fault coordinates.

     :param f: a distribution of values 
     :type f: `escript.Data`
     :param tol: relative tolerance used to decide if point is on fault 
     :type f: ``tol``
     :return: the maximum value of across all faults in the fault system, the fault tag the maximum is taken, and the coordinates of the location in the coordinates of the fault. The returned fault tag is ``None`` if no points on the fault are available
     """
     ref=-Lsup(f)*2
     f_max=ref
     t_max=None
     p_max=None
     x=f.getFunctionSpace().getX()
     for t in self.getTags():
        p,m=self.getParametrization(x,tag=t, tol=tol)
        loc=((m*f)+(1.-m)*ref).maxGlobalDataPoint()
        f_t=f.getTupleForGlobalDataPoint(*loc)[0]
        if f_t>f_max:
           f_max=f_t
           t_max=t
           p_max=p.getTupleForGlobalDataPoint(*loc)[0]

     return f_max, t_max, p_max

  def getParametrization(self,x,tag=None, tol=sqrt(EPSILON), outsider=None):
    """
    returns the parametrization of the fault ``tag`` in the fault system. In fact the values of the parametrization for at given coordinates ``x`` is returned. In addition to the value of the parametrization a mask is returned indicating if the given location is on the fault with given tolerance ``tol``.

    Typical usage of the this method is

    dom=Domain(..)
    x=dom.getX()
    fs=FaultSystem()
    fs.addFault(tag=3,...)
    p, m=fs.getParametrization(x, outsider=0,tag=3)
    saveDataCSV('x.csv',p=p, x=x, mask=m)

    to create a file with the coordinates of the points in ``x`` which are on the fault (as ``mask=m``) together with their location ``p`` in the fault coordinate system.

    :param x: location(s)
    :type x: `escript.Data` object or `numpy.ndarray`
    :param tag: the tage of the fault
    :param tol: relative tolerance to check if location is on fault.
    :type tol: ``float``
    :param outsider: value used for parametrization values outside the fault. If not present an appropriate value is choosen.
    :type outsider: ``float``
    :return: the coordinates ``x`` in the coordinate system of the fault and a mask indicating coordinates in the fault by 1 (0 elsewhere)
    :rtype: `escript.Data` object or `numpy.ndarray`
    """
    offsets=self.getW0Offsets(tag)
    w0_range=self.getW0Range(tag)[1]-self.getW0Range(tag)[0]
    if outsider == None:
       outsider=min(self.getW0Range(tag)[0],self.getW0Range(tag)[1])-abs(w0_range)/sqrt(EPSILON)
        
    p=x[0]*0 + outsider
    updated=x[0]*0 

    if self.getDim()==2:
        #
        #
        s=self.getFaultSegments(tag)[0]
        for i in xrange(1,len(s)):
           d=s[i]-s[i-1]
           h=x-s[i-1]
           h_l=length(h)
           d_l=length(d)
           if not d_l>0:
              raise ValueError,"segement %s in fault %s has zero length."%(i,tag)
           t=inner(h,d)/d_l**2
           t=t*whereNonPositive(t-1.)*whereNonNegative(t)
           r=length(h-t*d)/maximum(h_l,d_l)
           m=whereNonPositive(r-tol)*(1.-updated)
           p=(1.-m)*p+m*(offsets[i-1]+(offsets[i]-offsets[i-1])*t)
           updated=wherePositive(updated+m)
    else:
       raise ValueError,"3D is not supported yet."
    
    return p, updated
 
  def getSideAndDistance(self,x,tag=None):
    """
    returns the side and the distance at ``x`` from the fault ``tag``. 

    :param x: location(s)
    :type x: `escript.Data` object or `numpy.ndarray`
    :param tag: the tage of the fault
    :return: the side of ``x`` (positive means to the right of the fault, negative to the left) and the distance to the fault. Note that a value zero for the side means that that the side is undefined.
    """
    d=None
    side=None
    if self.getDim()==2:
        mat=numpy.array([[0., 1.], [-1., 0.] ])
        s=self.getFaultSegments(tag)[0]
        for i in xrange(1,len(s)):
           q=(s[i]-s[i-1])
           h=x-s[i-1]
           q_l=length(q)
           if not q_l>0:
              raise ValueError,"segement %s in fault %s has zero length."%(i,tag)
           qt=matrixmult(mat,q)   # orthogonal direction
           t=inner(q,h)/q_l**2
           t=maximum(minimum(t,1,),0.)
           p=h-t*q
           dist=length(p)
           lside=sign(inner(p,qt))
           if d == None:
               d=dist
               side=lside
           else:
               m=whereNegative(d-dist)
               m2=wherePositive(whereZero(abs(lside))+m)
               d=dist*(1-m)+d*m
               side=lside*(1-m2)+side*m2
    else:
       raise ValueError,"3D is not supported yet."
    return side, d

def patchMap(v,v1,v2,v3,v4,x_offset,x_length,depth):
    """
    maps a patch with corners v1,v2,v3,v4 onto 
    the rectangle [x_offset,0] x [x_offset+x_length,-depth]
    """
    #
    #   inspect the mapping v-> w=(w0,w1) with 0<=s<=1 0<=t<=1
    #   then we apply x_offset+s*x_length, -depth*t to map to the final domain
    #
    #   we need to find s and t such that
    #
    #   v=v1*(1-t)*(1-s)+v2*t*(1-s)+v3*t*s+v4*(1-t)*s
    #    =v1*(1-t-s+t*s)+v2*(t-t*s)+v3*t*s+v4*(s-t*s)
    #    =v1+t*(v2-v1)+s*(v4-v1)+s*t*(v1-v2+v3-v4)
    #   
    #   with d4=(v2-v1)-dot(v2-v1,v4-v1)/dot(v4-v1,v4-v1)*(v4-v1) => dot(d4,v2-v1)=0
    #   with d2=(v4-v1)-dot(v2-v1,v4-v1)/dot(v2-v1,v2-v1)*(v2-v1) => dot(d2,v4-v1)=0
    # 
    #   dot(v-v1,d4) = s*dot(v4-v1,d4)+s*t*dot(v3-v4,d4)   (a)
    #   dot(v-v1,d2) = t*dot(v2-v1,d2)+s*t*dot(v3-v2,d2)   (b)
    # 
    #   in addition we set n=cross(d2,d4) => dot(v2-v1,n)=dot(v4-v1,n)=0
    # 
    #   dot(v-v1,n)=s*t*dot(v3-v1,n)
    # 
    #   => s*t=dot(v-v1,n)/dot(v3-v1,n) if dot(v3-v1,n)!=0
    #      otherwise on can set s*t=0 to solve (a) and (b)
    #     
    h=v-v1
    h2=(v2-v1)
    h3=(v3-v1)
    h4=(v4-v1)
    d4=h4-numpy.dot(h2,h4)/numpy.dot(h2,h2)*h2
    d2=h2-numpy.dot(h2,h4)/numpy.dot(h4,h4)*h4
    n=numpy.cross(d4,d2)
    h31n=numpy.dot(n,v3-v1) 
    if abs(h31n)>0:
        st=n/h31n
    else:
        st=0*n
    s=numpy.dot(h,d4-st*numpy.dot(v3-v4,d4))/numpy.dot(v4-v1,d4)
    t=numpy.dot(h,d2-st*numpy.dot(v3-v2,d2))/numpy.dot(v2-v1,d2)
    return x_offset+s*x_length, -depth*t

