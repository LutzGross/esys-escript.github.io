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
 
  A fault system is defined by set of faults index by a tag. Each fault is defined by a starting point V0 and a list of 
  strikes ``strike`` and length ``l``. The strikes and the length is used to define a polyline with points ``V[i]`` such that

  - ``V[0]=V0``
  - ``V[i]=V[i]+ls[i]*array(cos(strike[i]),sin(strike[i]),0)''

  So ``strike`` defines the angle between the direction of the fault segment and the x0 axis. ''ls[i]``==0 is allowed.

  In case of a 3D model a fault plane is defined through a dip and depth. 

  The class provides a mechanism to parametrise each fault with the domain [0,w0_max] x [0, w1_max]  (to [0,w0_max] in the 2D case).
  """
  NOTAG="__NOTAG__"
  MIN_DEPTH_ANGLE=0.1
  def __init__(self,dim=3):
    """
    Sets up the fault system

    :param dim: spatial dimension
    :type dim: ``int`` of value 2 or 3
    """
    if not (dim == 2 or dim == 3):
       raise ValueError,"only dimension2 2 and 3 are supported."
    self.__dim=dim
    self.__top={}
    self.__ls={}
    self.__strikes={}
    self.__strike_vectors={}
    self.__medDepth={}
    self.__total_length={}
    if dim ==2:
       self.__depths=None
       self.__depth_vectors=None
       self.__dips=None
       self.__bottom=None
       self.__normals=None
    else:
       self.__depths={}
       self.__depth_vectors={}
       self.__dips={}
       self.__bottom={}
       self.__normals={}
    self.__offsets={}
    self.__w1_max={}
    self.__w0_max={}
    self.__center=None
    self.__orientation = None
  def getStart(self,tag=None):
     """
     returns the starting point of fault ``tag``
     :rtype: `numpy.ndarray'.
     """
     return self.getTopPolyline(tag)[0]

  def getTags(self):
     """
     returns a list of the tags used by the fault system
     :rtype: ``list``
     """
     return self.__top.keys()
  def getDim(self):
     """
     returns the spatial dimension
     :rtype: ``int``
     """
     return self.__dim

  def getTopPolyline(self, tag=None):
     """
     returns the polyline used to describe fault tagged by ``tag`` 
      
     :param tag: the tag of the fault
     :type tag: ``float`` or ``str``
     :return: the list of vertices defining the top of the fault.  The coordinates are `numpy.ndarray'.
     """
     if tag==None: tag=self.NOTAG
     return self.__top[tag]
  def getStrikes(self, tag=None):
     """
     returns the strike of the segements in fault ``tag``
     :rtype: ``list`` of ``float``
     """
     if tag==None: tag=self.NOTAG
     return self.__strikes[tag]
  def getStrikeVectors(self, tag=None):
     """
     returns the strike vectors of fault ``tag``
     :rtype: ``list`` of `numpy.ndarray'.
     """
     if tag==None: tag=self.NOTAG
     return self.__strike_vectors[tag]
  def getLengths(self, tag=None):
     """
     returns the lengths of segments in fault ``tag``
     :rtype: ``list`` of `float``
     """
     if tag==None: tag=self.NOTAG
     return self.__ls[tag]

  def getTotalLength(self, tag=None):
     """
     returns the total unrolled length of fault ``tag``
     :rtype: `float``
     """
     if tag==None: tag=self.NOTAG
     return self.__total_length[tag]

  def getMediumDepth(self,tag=None):
     """
     returns the medium depth of fault ``tag``
     :rtype: ``float``
     """
     if tag==None: tag=self.NOTAG
     return self.__medDepth[tag]

  def getDips(self, tag=None):
     """
     returns the list of the dips of the segements in fault ``tag``
     :param tag: the tag of the fault
     :type tag: ``float`` or ``str``
     :return: the list of segment dips. In the 2D case None is returned.
     """
     if tag==None: tag=self.NOTAG
     if self.getDim()==3:
         return self.__dips[tag]
     else:
         return None

  def getBottomPolyline(self, tag=None):
     """
     returns the list of the vertices defining the bottom of the fault ``tag``
     :param tag: the tag of the fault
     :type tag: ``float`` or ``str``
     :return: the list of vertices. In the 2D case None is returned.
     """
     if tag==None: tag=self.NOTAG
     if self.getDim()==3:
         return self.__bottom[tag]
     else:
         return None

  def getSegmentNormals(self, tag=None):
     """
     returns the list of the normals of the segments in fault ``tag``
     :param tag: the tag of the fault
     :type tag: ``float`` or ``str``
     :return: the list of vectors normal to the segments. In the 2D case None is returned.
     """
     if tag==None: tag=self.NOTAG
     if self.getDim()==3:
         return self.__normals[tag]
     else:
         return None

  def getDepthVectors(self, tag=None):
     """
     returns the list of the depth vector at top vertices in fault ``tag``.
     :param tag: the tag of the fault
     :type tag: ``float`` or ``str``
     :return: the list of segment depths. In the 2D case None is returned.
     """
     if tag==None: tag=self.NOTAG
     if self.getDim()==3:
         return self.__depth_vectors[tag]
     else:
         return None
  def getDepths(self, tag=None):
     """
     returns the list of the depths of the segements in fault ``tag``.
     :param tag: the tag of the fault
     :type tag: ``float`` or ``str``
     :return: the list of segment depths. In the 2D case None is returned.
     """
     if tag==None: tag=self.NOTAG
     if self.getDim()==3:
         return self.__depths[tag]
     else:
         return None

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


  def getCenterOnSurface(self):
      """
      returns the center point of the fault system at the surface 
      :rtype: `numpy.ndarray`
      """
      if self.__center == None:
        self.__center=numpy.zeros((3,),numpy.float)
        counter=0
        for t in self.getTags():
            for s in self.getTopPolyline(t):
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
              for s in self.getTopPolyline(t):
                covariant[0,0]+=(center[0]-s[0])**2
                covariant[0,1]+=(center[1]-s[1])*(center[0]-s[0])
                covariant[1,1]+=(center[1]-s[1])**2
                covariant[1,0]+=(center[1]-s[1])*(center[0]-s[0])
          e, V=numpy.linalg.eigh(covariant)
          if e[0]>e[1]:
             d=V[:,0]
          else:
             d=V[:,1]
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
         strikes=[ s+ rot for s in self.getStrikes(t) ]
         V0=self.getStart(t)

         self.addFault(strikes = [ s+ rot for s in self.getStrikes(t) ], \
                       ls = self.getLengths(t), \
                       V0=numpy.dot(mat,self.getStart(t)+shift), \
                       tag =t, \
                       dips=self.getDips(t),\
                       depths=self.getDepths(t), \
                       w0_offsets=self.getW0Offsets(t), \
                       w1_max=-self.getW1Range(t)[0]) 

  def addFault(self, strikes, ls, V0=[0.,0.,0.],tag=None, dips=None, depths= None, w0_offsets=None, w1_max=None):
     """
     adds a new fault to the fault system. The fault is named by ``tag``. 

     The fault is defined by a starting point V0 and a list of ``strikes`` and length ``ls``. The strikes and the length 
     is used to define a polyline with points ``V[i]`` such that

     - ``V[0]=V0``
     - ``V[i]=V[i]+ls[i]*array(cos(strikes[i]),sin(strikes[i]),0)''

     So ``strikes`` defines the angle between the direction of the fault segment and the x0 axis. In 3D ''ls[i]``==0 is allowed.

     In case of a 3D model a fault plane is defined through a dip ``dips`` and depth ``depths``. 
     From the dip and the depth the polyline ``bottom`` of the bottom of the fault is computed.


     Each segment in the fault is decribed by the for vertices ``v0=top[i]``, ``v1``==``top[i+1]``, ``v2=bottom[i]`` and ``v3=bottom[i+1]`` 
     The segment is parametrized by ``w0`` and ``w1`` with ``w0_offsets[i]<=w0<=w0_offsets[i+1]`` and ``-w1_max<=w1<=0``. Moreover 
   
     - ``(w0,w1)=(w0_offsets[i]  ,       0)->v0``
     - ``(w0,w1)=(w0_offsets[i+1],       0)->v1``
     - ``(w0,w1)=(w0_offsets[i]  , -w1_max)->v2``
     - ``(w0,w1)=(w0_offsets[i+1], -w1_max)->v3``

     If no ``w0_offsets`` is given, 
  
     - ``w0_offsets[0]=0``
     - ``w0_offsets[i]=w0_offsets[i-1]+L[i]``

     where ``L[i]`` is the length of the segments on the top in 2D and in the middle of the segment in 3D. 

     If no ``w1_max`` is given, the average fault depth is used.


     :param strikes: list of strikes. This is the angle of the fault segment direction with x0 axis. Right hand rule applies.
     :type strikes: ``list`` of ``float``
     :param ls: list of fault lengths. In the case of a 3D fault a segment may have length 0. 
     :type ls: ``list`` of ``float``
     :param V0: start point of the fault
     :type V0: ``list`` or `numpy.ndarray' with 2 or 3 components. ``V0[2]`` must be zero.
     :param tag: the tag of the fault. If fault ``tag`` already exists it is overwritten.
     :type tag: ``float`` or ``str``
     :param dips: list of dip angles. Right hand rule around strike direction applies.
     :type dips: ``list`` of ``float``
     :param depths: list of segment depth. Value mut be positive in the 3D case.
     :type depths: ``list`` of ``float``
     :param w0_offsets: ``w0_offsets[i]`` defines the offset of the segment ``i`` in the fault to be used in the parametrization of the fault. If not present the cumulative length of the fault segments is used. 
     :type w0_offsets: ``list`` of ``float`` or ``None``
     :param w1_max: the maximum value used for parametrization of the fault in the depth direction. If not present the mean depth of the fault segments is used.
     :type w1_max: ``float``
     :note: In the three dimensional case the lists ``dip`` and ``top`` must have the same length.
     """
     if tag==None: 
         tag=self.NOTAG
     else:
         if self.NOTAG in self.getTags():
              raise ValueError,'Attempt to add a fault with no tag to a set of existing faults'
     if not isinstance(strikes, list): strikes=[strikes, ]
     n_segs=len(strikes)
     if not isinstance(ls, list): ls=[ ls for i in xrange(n_segs) ]
     if not n_segs==len(ls):
         raise ValueError,"number of strike direction and length must match."
     if len(V0)>2:
          if abs(V0[2])>0: raise Value,"start point needs to be surface (3rd component ==0)"
     if self.getDim()==2 and not  (dips==None and depths == None) :
           raise ValueError,'Spatial dimension two does not support dip and depth for faults.'
     if not dips == None:
        if not isinstance(dips, list): dips=[dips for i in xrange(n_segs) ]
        if n_segs != len(dips):
           raise ValueError,'length of dips must be one less then the length of top.'
     if not depths == None:
        if not isinstance(depths, list): depths=[depths for i in xrange(n_segs+1) ]
        if n_segs+1 != len(depths):
           raise ValueError,'length of depths must be one less then the length of top.'
     if w0_offsets != None:
       if len(w0_offsets) != n_segs+1:
          raise ValueError,'expected length of w0_offsets is %s'%(n_segs)
     self.__center=None
     self.__orientation = None
     #
     #  in the 2D case we don't allow zero length:
     #
     if self.getDim() == 2:
        for l in ls:
            if l<=0: raise ValueError,"length must be positive"
     else:
        for l in ls:
            if l<0: raise ValueError,"length must be non-negative"
        for i in xrange(n_segs+1):
           if depths[i]<0: raise ValueError,"negative depth."
     # 
     #   translate start point to numarray
     #
     V0= numpy.array(V0[:self.getDim()],numpy.double)
     #
     #  set strike vectors:
     #
     strike_vectors=[]
     top_polyline=[V0]   
     total_length=0
     for i in xrange(n_segs):
         v=numpy.zeros((self.getDim(),))
         v[0]=cos(strikes[i])
         v[1]=sin(strikes[i])
         strike_vectors.append(v)
         top_polyline.append(top_polyline[-1]+ls[i]*v)
         total_length+=ls[i]
     #
     #    normal and depth direction
     #
     if self.getDim()==3:
        normals=[]
        for i in xrange(n_segs):
           normals.append(numpy.array([sin(dips[i])*strike_vectors[i][1], -sin(dips[i])*strike_vectors[i][0], cos(dips[i])]) )
  
        d=numpy.cross(strike_vectors[0],normals[0])
        if d[2]>0:
             f=-1
        else:
             f=1
        depth_vectors=[f*depths[0]*d/numpy.linalg.norm(d) ]
        for i in xrange(1,n_segs):
            d=-numpy.cross(normals[i-1],normals[i])
            d_l=numpy.linalg.norm(d)
            if d_l<=0:
                 d=numpy.cross(strike_vectors[i],normals[i])
                 d_l=numpy.linalg.norm(d)
            else:
                 for L in [ strike_vectors[i], strike_vectors[i-1]]:
                    if numpy.linalg.norm(numpy.cross(L,d)) <= self.MIN_DEPTH_ANGLE * numpy.linalg.norm(L) * d_l:
                         raise ValueError,"%s-th depth vector %s too flat."%(i, d)
            if d[2]>0:
                f=-1
            else:
                f=1
            depth_vectors.append(f*d*depths[i]/d_l)
        d=numpy.cross(strike_vectors[n_segs-1],normals[n_segs-1])
        if d[2]>0:
             f=-1
        else:
             f=1
        depth_vectors.append(f*depths[n_segs]*d/numpy.linalg.norm(d))
        bottom_polyline=[ top_polyline[i]+depth_vectors[i] for i in xrange(n_segs+1) ]
     #
     #   calculate offsets if required:
     #
     if w0_offsets==None:
        w0_offsets=[0.] 
        for  i in xrange(n_segs):
            if self.getDim()==3:
               w0_offsets.append(w0_offsets[-1]+(float(numpy.linalg.norm(bottom_polyline[i+1]-bottom_polyline[i]))+ls[i])/2.)
            else:
               w0_offsets.append(w0_offsets[-1]+ls[i])
     w0_max=max(w0_offsets)
     if self.getDim()==3:
        self.__normals[tag]=normals
        self.__depth_vectors[tag]=depth_vectors
        self.__depths[tag]=depths
        self.__dips[tag]=dips
        self.__bottom[tag]=bottom_polyline
     self.__ls[tag]=ls
     self.__strikes[tag]=strikes
     self.__strike_vectors[tag]=strike_vectors
     self.__top[tag]=top_polyline
     self.__total_length[tag]=total_length
     self.__offsets[tag]=w0_offsets

     if self.getDim()==2:
        self.__medDepth[tag]=0.
     else:
        self.__medDepth[tag]=sum([ numpy.linalg.norm(v) for v in depth_vectors])/len(depth_vectors)
     if w1_max==None or self.getDim()==2: w1_max=self.__medDepth[tag]
     self.__w0_max[tag]=w0_max
     self.__w1_max[tag]=w1_max

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
    :param tag: the tag of the fault
    :param tol: relative tolerance to check if location is on fault.
    :type tol: ``float``
    :param outsider: value used for parametrization values outside the fault. If not present an appropriate value is choosen.
    :type outsider: ``float``
    :return: the coordinates ``x`` in the coordinate system of the fault and a mask indicating coordinates in the fault by 1 (0 elsewhere)
    :rtype: `escript.Data` object or `numpy.ndarray`
    """
    offsets=self.getW0Offsets(tag)
    w1_range=self.getW1Range(tag)
    w0_range=self.getW0Range(tag)[1]-self.getW0Range(tag)[0]
    if outsider == None:
       outsider=min(self.getW0Range(tag)[0],self.getW0Range(tag)[1])-abs(w0_range)/sqrt(EPSILON)
        
    if isinstance(x,list): x=numpy.array(x, numpy.double)
    updated=x[0]*0 

    if self.getDim()==2:
        #
        #
        p=x[0]*0 + outsider
        top=self.getTopPolyline(tag)
        for i in xrange(1,len(top)):
           d=top[i]-top[i-1]
           h=x-top[i-1]
           h_l=length(h)
           d_l=length(d)
           s=inner(h,d)/d_l**2
           s=s*whereNonPositive(s-1.-tol)*whereNonNegative(s+tol)
           m=whereNonPositive(length(h-s*d)-tol*maximum(h_l,d_l))*(1.-updated)
           p=(1.-m)*p+m*(offsets[i-1]+(offsets[i]-offsets[i-1])*s)
           updated=wherePositive(updated+m)
    else:
        p=x[:2]*0 + outsider
        top=self.getTopPolyline(tag)
        bottom=self.getBottomPolyline(tag)
        n=self.getSegmentNormals(tag)
        for i in xrange(len(top)-1):
            h=x-top[i]
            R=top[i+1]-top[i]
            r=bottom[i+1]-bottom[i]
            D0=bottom[i]-top[i]
            D1=bottom[i+1]-top[i+1]
            s_upper=matrix_mult(numpy.linalg.pinv(numpy.vstack((R,D1)).T),h)
            s_lower=matrix_mult(numpy.linalg.pinv(numpy.vstack((r,D0)).T),h)
            m_ul=wherePositive(s_upper[0]-s_upper[1])
            s=s_upper*m_ul+s_lower*(1-m_ul)
            s0=s[0]
            s1=s[1]
            m=whereNonNegative(s0+tol)*whereNonPositive(s0-1.-tol)*whereNonNegative(s1+tol)*whereNonPositive(s1-1.-tol)
            s0=s0*m
            s1=s1*m
            atol=tol*maximum(length(h),length(top[i]-bottom[i+1]))
            m=whereNonPositive(length(h-s0*R-s1*D1)*m_ul+(1-m_ul)*length(h-s0*r-s1*D0)-atol)
            p[0]=(1.-m)*p[0]+m*(offsets[i]+(offsets[i+1]-offsets[i])*s0)
            p[1]=(1.-m)*p[1]+m*(w1_range[1]+(w1_range[0]-w1_range[1])*s1)
            updated=wherePositive(updated+m)
    
    return p, updated
 
  def getSideAndDistance(self,x,tag=None):
    """
    returns the side and the distance at ``x`` from the fault ``tag``. 

    :param x: location(s)
    :type x: `escript.Data` object or `numpy.ndarray`
    :param tag: the tag of the fault
    :return: the side of ``x`` (positive means to the right of the fault, negative to the left) and the distance to the fault. Note that a value zero for the side means that that the side is undefined.
    """
    d=None
    side=None
    if self.getDim()==2:
        mat=numpy.array([[0., 1.], [-1., 0.] ])
        s=self.getTopPolyline(tag)
        for i in xrange(1,len(s)):
           q=(s[i]-s[i-1])
           h=x-s[i-1]
           q_l=length(q)
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
        ns=self.getSegmentNormals(tag)
        top=self.getTopPolyline(tag)
        bottom=self.getBottomPolyline(tag)
        for i in xrange(len(top)-1):
            h=x-top[i]
            R=top[i+1]-top[i]
            r=bottom[i+1]-bottom[i]
            D0=bottom[i]-top[i]
            D1=bottom[i+1]-top[i+1]
            s_upper=matrix_mult(numpy.linalg.pinv(numpy.vstack((R,D1)).T),h)
            s_lower=matrix_mult(numpy.linalg.pinv(numpy.vstack((r,D0)).T),h)
            m_ul=wherePositive(s_upper[0]-s_upper[1])
            s=s_upper*m_ul+s_lower*(1-m_ul)
            s=maximum(minimum(s,1.),0)
            p=h-s[0]*(R*m_ul+(1-m_ul)*r)-s[1]*(D1*m_ul+(1-m_ul)*D0)
            dist=length(p)
            lside=sign(inner(p,ns[i]))
            if d == None:
               d=dist
               side=lside
            else:
               m=whereNegative(d-dist)
               m2=wherePositive(whereZero(abs(lside))+m)
               d=dist*(1-m)+d*m
               side=lside*(1-m2)+side*m2

    return side, d

