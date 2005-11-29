# $Id:$


import os
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.escript.modelframe import Model
import numarray
import math
import urllib

EPS=1.e-5

#====================================================================================================================================================
class GridData:
      """
      this object is used to store data on grid.
      it will be replaced by Bruce at a later stage.

      data[i,j] are data are x=j*s[0]+o[0] and y=i*s[1]+o[1]

      for 0<=j<n[0] and 0<=i<n[1]
      """
      def __init__(self,s,o,n,data):
         self.s=tuple(s)
         self.o=tuple(o)
         self.n=tuple(n)
         self.data=data

      def getOrigin(self):
          return self.o
      def getSpacing(self):
          return self.s
      def getData(self):
          return self.data

      def interpolate(self,x):
          if hasattr(x,"convertToNumArray"):
             x_array=x.convertToNumArray()
             return_data_object=True
          else:
             x_array=numarray.array(x)
             return_data_object=False
          data=numarray.zeros(x_array.shape[0],numarray.Float)
          ox,oy=self.getOrigin()
          dx,dy=self.getSpacing()
          data_array=self.getData()
          i_dx=1
          i_dy=1
          for i in range(x_array.shape[0]):
             x0=x_array[i,0]-ox
             x1=x_array[i,1]-oy
             j0=min(max(int(x0/dx),0),data_array.shape[1]-1-i_dy)
             j1=min(max(int(x1/dy),0),data_array.shape[0]-1-i_dx)
             f01=(x0-j0*dx)/dx
             f00=1.-f01
             f11=(x1-j1*dy)/dy
             f10=1.-f11
             H00=data_array[j1,j0]
             H01=data_array[j1,j0+i_dx]
             H11=data_array[j1+i_dy,j0+i_dx]
             H10=data_array[j1+i_dy,j0]
             data[i]=(H00*f00+H01*f01)*f10+(H10*f00+H11*f01)*f11
          if return_data_object:
             out=Scalar(0,x.getFunctionSpace())
             out.fillFromNumArray(data)
             return out
          else:
             return data

      def getVTK(self):
          pass


class PointOnEarthSurface:
   """
   coordinates of a  point on the surface of the Earth
   """
   def __init__(self,phi=0,theta=0):
       self.phi=phi
       self.theta=theta
   def __str__(self):
       return "(%s,%s)"%(self.phi,self.theta)

   def __sub__(self,other):
       return self.dist(other)

   def split(self,p,t):
      phi=self.phi+t*(p.phi-self.phi)
      theta=self.theta+t*(p.theta-self.theta)
      return PointOnEarthSurface(phi,theta)

   def midPoint(self,other):
      return PointOnEarthSurface((self.phi+other.phi)/2,(self.theta+other.theta)/2)

   def dist(self,other):
      return math.sqrt((self.phi-other.phi)**2+(self.theta-other.theta)**2)

class RegionOnEarthSurface:
     """
     defines an region by a south,east and north,west corner
     """
     RADIUS=6378.8e3
     GRAVITY=9.81
     WEST=0
     NORTH=1
     EAST=2
     SOUTH=3
     def __init__(self,south_east=PointOnEarthSurface(0.,0.),north_west=PointOnEarthSurface(10.,10.),resolution=1.):
          if resolution<=0:
              raise ValueError,"resolution must be positive"
          if  south_east.phi>=north_west.phi:
              raise ValueError,"east south corner must be further east than west north corner"
          if  south_east.theta>=north_west.theta:
              raise ValueError,"east south corner must be further down than west north corner"
          if  north_west.theta-south_east.theta < resolution/2:
              raise ValueError,"latitude length of region must be 2*bigger then resolution"
          if  north_west.phi-south_east.phi < resolution/2:
              raise ValueError,"longitude length of region must be 2*bigger then resolution"

          self.south_east=south_east
          self.north_west=north_west
          self.resolution=resolution

     def __str__(self):
       return "RegionOnEarthSurface between %s and %s"%(str(self.south_east),str(self.north_west))

     def isOnFace(self,p):
         return self.isOnThisFace(p,self.SOUTH) or self.isOnThisFace(p,self.NORTH) or self.isOnThisFace(p,self.WEST) or self.isOnThisFace(p,self.EAST)
     def isOnThisFace(self,p,face):
        if face==self.WEST:
           return self.south_east.phi==p.phi
        if face==self.EAST:
         return p.phi==self.north_west.phi
        if face==self.SOUTH:
         return self.south_east.theta==p.theta
        if face==self.NORTH:
           return p.theta==self.north_west.theta

     def isBoxVertex(self,p):
         return ( self.isOnThisFace(p,self.WEST) and  self.isOnThisFace(p,self.SOUTH) ) or \
                ( self.isOnThisFace(p,self.WEST) and  self.isOnThisFace(p,self.NORTH) ) or \
                ( self.isOnThisFace(p,self.EAST) and  self.isOnThisFace(p,self.NORTH) ) or \
                ( self.isOnThisFace(p,self.EAST) and  self.isOnThisFace(p,self.SOUTH) )


     def getFace(self,p):
        # order is critical
        if self.south_east.phi==p.phi: return self.WEST
        if p.theta==self.north_west.theta: return self.NORTH
        if p.phi==self.north_west.phi: return self.EAST
        if self.south_east.theta==p.theta: return self.SOUTH

     def comparePointsOnFace(self,p0,p1):
         f0=self.getFace(p0)
         f1=self.getFace(p1)

         if f0<f1:
            return -1
         elif f0>f1:
            return 1
         else:
            if f0 == self.WEST:
                if p0.theta<p1.theta:
                   return -1
                elif p0.theta>p1.theta:
                   return 1
                else:
                   return 0
            elif f0 == self.EAST:
                if p0.theta<p1.theta:
                   return 1
                elif p0.theta>p1.theta:
                   return -1
                else:
                   return 0
            elif f0 == self.NORTH:
                if p0.phi<p1.phi:
                   return -1
                elif p0.phi>p1.phi:
                   return 1
                else:
                   return 0
            else:
                if p0.phi<p1.phi:
                   return 1
                elif p0.phi>p1.phi:
                   return -1
                else:
                   return 0


     def isInRegion(self,p):
         return  self.south_east.phi<=p.phi \
             and p.phi<=self.north_west.phi \
             and self.south_east.theta<=p.theta \
             and p.theta<=self.north_west.theta

     def cutSegment(self,p0,p1):
            t=None
            p=None
            tmp=self.interseptOfSegment(p0,p1,d=0,v=self.south_east.phi)
            if not tmp==None:
                p_tmp=p0.split(p1,tmp)
                if self.isInRegion(p_tmp) and t<tmp:
                   t=tmp
                   p=p_tmp

            tmp=self.interseptOfSegment(p0,p1,d=0,v=self.north_west.phi)
            if not tmp==None:
                p_tmp=p0.split(p1,tmp)
                if self.isInRegion(p_tmp) and t<tmp:
                   t=tmp
                   p=p_tmp

            tmp=self.interseptOfSegment(p0,p1,d=1,v=self.south_east.theta)
            if not tmp==None:
                p_tmp=p0.split(p1,tmp)
                if self.isInRegion(p_tmp) and t<tmp:
                   t=tmp
                   p=p_tmp

            tmp=self.interseptOfSegment(p0,p1,d=1,v=self.north_west.theta)
            if not tmp==None:
                p_tmp=p0.split(p1,tmp)
                if self.isInRegion(p_tmp) and t<tmp:
                   t=tmp
                   p=p_tmp
            return p

     def interseptOfSegment(self,p0,p1,d=0,v=0.):
         """
         find t isn [0,1] such that p0+t*(p0-p1) cuts x[d]=v. if it does nit exist None is returned.
         """
         if d==0:
           a=p0.phi
           b=p1.phi
         else:
           a=p0.theta
           b=p1.theta
         if b==a:
           if a==v:
              t=0
           else:
              t=None
         else:
           t=(v-a)/(b-a)
           if not (0<=t and t<=1): t=None
         return t

class Polyline:
     """
     defines set of segments through a list of coordinates
     """
     def __init__(self,list_of_coordinates=[],name="none"):
         c=[]
         if len(list_of_coordinates)>0:
             for i in range(len(list_of_coordinates)-1):
                if list_of_coordinates[i]-list_of_coordinates[i+1]>0: c.append(list_of_coordinates[i])
             c.append(list_of_coordinates[-1])
         self.list_of_coordinates=c
         self.name=name
     def getDiameter(self):
         out=0.
         for p in self.list_of_coordinates:
            for q in self.list_of_coordinates:
               out=max(out,p-q)
         return out
     def isLoop(self):
           if len(self)>0:
             return  not self[0]-self[-1]>EPS
           else:
             return True

     def insert(self,index,coordinate):
          """
          insert point before index
          """
          if self.list_of_coordinates[index]-coordinate<EPS:
               return index
          elif self.list_of_coordinates[index-1]-coordinate<EPS:
               return index-1
          else:
             self.list_of_coordinates.insert(index,coordinate)
             return index

     def split(self,index):
         """
         splits the polyline at point index
         """
         return Polyline(self.list_of_coordinates[:index],self.name),Polyline(self.list_of_coordinates[index:],self.name)

     def __getitem__(self,s):
           return self.list_of_coordinates.__getitem__(s)
     def __iter__(self):
          return iter(self.list_of_coordinates)

     def __str__(self):
         if self.isLoop():
            out="loop["
         else:
            out="["
         for i in self:
            if out[-1]=="[":
                out+="%s"%str(i)
            else:
                out+=",%s"%str(i)
         return out+"]"

     def __len__(self):
         return len(self.list_of_coordinates)

     def orientation(self):
         """
         returns the orientation of the polyline
         """
         if not self.isLoop():
              raise TypeError,"polyline is not a loop"

         integ=0.
         for i in range(len(self.list_of_coordinates)-1):
             p0=self.list_of_coordinates[i]
             p1=self.list_of_coordinates[i+1]
             integ+=(p1.theta-p0.theta)*(p1.phi+p0.phi)-(p1.phi-p0.phi)*(p1.theta+p0.theta)
         if integ>=0.:
            return 1.
         else:
            return -1.

     def givePositiveOrientation(self):
         if self.orientation()<0: self.list_of_coordinates.reverse()

class Coastline:
     """
     defines a coast line by a Polyline within a RegionOnEarthSurface
     """
     def __init__(self,region,name="none"):
         self.region=region
         self.name=name
         self.polylines=[]
     def __str__(self):
         out=self.name+" in "+str(self.region)
         for pl in self.polylines:
             out+="\n"+str(pl)
         return out
     def makeTriangulation(self,south_east_is_water=True,south_west_is_water=True,north_west_is_water=True,north_east_is_water=True):
         self.clean()
         vertices=[]
         segments=[]
         holes=[]
         vertices_on_face=[]
         for pl in self.polylines:
             if pl.getDiameter()>self.region.resolution:
                short_pl=[pl[0]]
                for i in range(1,len(pl)):
                      if pl[i]-short_pl[-1]>self.region.resolution/5:
                         short_pl.append(pl[i])
                      elif i==len(pl)-1:
                         short_pl[-1]=pl[i]
                if pl.isLoop():
                   if len(short_pl)>3:
                       offset=len(vertices)
                       v_tmp=[short_pl[0]]
                       s_tmp=[]
                       for i in range(1,len(short_pl)):
                          if i==len(short_pl)-1:
                            s_tmp.append((len(v_tmp)-1+offset,offset))
                          else:
                            s_tmp.append((len(v_tmp)-1+offset,len(v_tmp)+offset))
                            v_tmp.append(short_pl[i])
                       vertices+=v_tmp
                       segments+=s_tmp
                       # find a point in the loop:
                       d_phi=v_tmp[1].phi-v_tmp[0].phi
                       d_theta=v_tmp[1].theta-v_tmp[0].theta
                       l=math.sqrt(d_phi**2+d_theta**2)
                       mid=v_tmp[0].midPoint(v_tmp[1])
                       n_phi=-d_theta/l
                       n_theta=d_phi/l
                       s=self.region.resolution
                       for seg in s_tmp:
                         p0=vertices[seg[0]]
                         p1=vertices[seg[1]]
                         a_phi=p1.phi-p0.phi
                         a_theta=p1.theta-p0.theta
                         d=a_theta*n_phi-a_phi*n_theta
                         if abs(d)>0.:
                            t=((mid.theta-p0.theta)*n_phi-(mid.phi-p0.phi)*n_theta)/d
                            if 0<=t and t<=1:
                                s_tmp=((mid.theta-p0.theta)*a_phi-(mid.phi-p0.phi)*a_theta)/d
                                if s_tmp>EPS: s=min(s,s_tmp)
                       h=PointOnEarthSurface(mid.phi+s/2*n_phi,mid.theta+s/2*n_theta)
                       holes.append(h)
                else:
                   if len(short_pl)>1:
                       if self.region.isOnFace(short_pl[0]): vertices_on_face.append(short_pl[0])
                       if self.region.isOnFace(short_pl[-1]): vertices_on_face.append(short_pl[-1])
                       vertices.append(short_pl[0])
                       for i in range(1,len(short_pl)):
                          segments.append((len(vertices)-1,len(vertices)))
                          vertices.append(short_pl[i])
         # put into the bounding box:
         new_vertices=[]
         if south_east_is_water:
            new_vertices.append(PointOnEarthSurface(self.region.south_east.phi,self.region.south_east.theta))
         if north_east_is_water:
            new_vertices.append(PointOnEarthSurface(self.region.south_east.phi,self.region.north_west.theta))
         if north_west_is_water:
            new_vertices.append(PointOnEarthSurface(self.region.north_west.phi,self.region.north_west.theta))
         if south_west_is_water:
            new_vertices.append(PointOnEarthSurface(self.region.north_west.phi,self.region.south_east.theta))

         # add new vertices if they don't exist yet
         for q in new_vertices:
             for p2 in vertices_on_face:
                 if p2-q<EPS:
                     q=None
                     raise ValueError,"coast line crosses boundary box vertex. This case is currenrly not supported."
             if not q==None:
                vertices.append(q)
                vertices_on_face.append(q)
         vertices_on_face.sort(self.region.comparePointsOnFace)
         index=0
         walking_on_water=south_east_is_water
         l=len(vertices_on_face)
         while index<l:
             p1=vertices_on_face[(index+1)%l]
             p0=vertices_on_face[index]
             if walking_on_water:
                 segments.append((vertices.index(p0),vertices.index(p1)))
                 walking_on_water=False
             else:
                 if self.region.isBoxVertex(p0):
                     segments.append((vertices.index(p0),vertices.index(p1)))
                 else:
                     walking_on_water=True
             index+=1
         return EarthTriangulation(vertices,segments,holes,self.region.resolution)

     def clean(self):
        """
        cleans up the coast line by joining polylines to loops or connecting faces of the region
        """
        # find a poylines that are linked
        while True:
            k0=None
            for pl in self.polylines:
              if not pl.isLoop():
                 for k in [0,-1]:
                     for ql in self.polylines:
                        if not (ql==pl or ql.isLoop()):
                           for k2 in [0,-1]:
                               if ql[k2]-pl[k]<EPS:
                                    pl0=pl
                                    pl1=ql
                                    k0=k
                                    k1=k2
                                    break
                        if not k0==None: break # ql
                     if not k0==None: break # k
              if not k0==None: break # pl

            if k0==None:
                break
            else:
                self.polylines.remove(pl0)
                self.polylines.remove(pl1)
                pl0c=pl0.list_of_coordinates
                pl1c=pl1.list_of_coordinates
                if k0==0: pl0c.reverse()
                if k1==-1: pl1c.reverse()
                pl=Polyline(pl0c+pl1c[1:],pl0.name+" + "+pl1.name)
                self.append(pl)

        # find a polyline that is not a loop and has an end or start point not on the face of the region:
        while True:
            pl=None
            k=None
            for pl2 in self.polylines:
              if not pl2.isLoop():
                 pl=pl2
                 if not self.region.isOnFace(pl[0]): k=0
                 if not self.region.isOnFace(pl[-1]): k=-1
                 if not k==None: break
            if k==None: break
            self.polylines.remove(pl)
            d_min=50000.
            k_min=None
            for pl2 in self.polylines:
               if not pl2.isLoop():
                  for k2 in [0,-1]:
                      if not self.region.isOnFace(pl2[k2]):
                        d2=pl2[k2]-pl[k]
                        if d2<d_min:
                            d_min=d2
                            pl_min=pl2
                            k_min=k2
            if k_min==None:
                 raise ValueError,"cannot link coastline %s to any other coastline."%pl.name
            plc=pl.list_of_coordinates
            plc_min=pl_min.list_of_coordinates
            if k==0: plc.reverse()
            if k_min==-1: plc_min.reverse()
            if d_min<EPS:
               new_pl=Polyline(plc+plc_min[1:],pl.name+" + "+pl_min.name)
            else:
               new_pl=Polyline(plc+plc_min,pl.name+" + "+pl_min.name)
            self.polylines.remove(pl_min)
            self.append(new_pl)
        # give positive orientation to loops:
        for pl in self.polylines:
             if pl.isLoop(): pl.givePositiveOrientation()

     def append(self,polyline=Polyline()):
        """append a polyline """
        if len(polyline)>1:
           pl=[]
           outside_region=None
           for i in range(len(polyline)):
              if not self.region.isInRegion(polyline[i]):
                  outside_region=i
                  break
#              pl.append(self.region.nudgeToFace(polyline[i]))
              pl.append(polyline[i])
           if not outside_region==None:
             if outside_region==0:
                for i in range(outside_region+1,len(polyline)):
                    if self.region.isInRegion(polyline[i]):
                       polyline.insert(i,self.region.cutSegment(polyline[i-1],\
                                                                polyline[i]))
                       pl1=polyline.split(i)[1]
                       self.append(pl1)
                       break
             else:
                # split polyline in two part first one is fully within the region the other starts with
                # point outside the region
                c=self.region.cutSegment(polyline[outside_region-1],polyline[outside_region])
                i=polyline.insert(outside_region,c)
                pl0,pl1=polyline.split(i+1)
                self.append(pl0)
                self.append(pl1)
           else:
             if len(pl)>1:
                  pply= Polyline(pl,polyline.name)
                  self.polylines.append(pply)

class EarthTriangulation:
      GENERATOR="triangle -pqa%g %s"

      def __init__(self,vertices=[],segments=[],holes=[],resolution=1.):
           self.fn=os.tempnam()
           #   write triangle input file
           poly_file=self.fn+".poly"
           f=open(poly_file,"w")
           f.writelines("%d %d %d %d\n"%(len(vertices),2,0,0))
           for i in range(len(vertices)): f.writelines("%d %e %e\n"%(i,vertices[i].phi,vertices[i].theta))
           f.writelines("%d %d\n"%(len(segments),0))
           for i in range(len(segments)): f.writelines("%d %d %d\n"%(i,segments[i][0],segments[i][1]))
           f.writelines("%d\n"%(len(holes)))
           for i in range(len(holes)):  f.writelines("%d %e %e\n"%(i,holes[i].phi,holes[i].theta))
           f.close()
           # start mesh generator:
           os.system(self.GENERATOR%(resolution**2,poly_file))
           # read mesh file:
           self.node_coordinates=[]
           self.node_tags=[]
           self.node_ids=[]
           self.triangles_nodes=[]
           self.triangles_id=[]
           node_file=open("%s.1.node"%self.fn,"r")
           nodes=node_file.readline().strip().split()
           nn=int(nodes[0])
           for i in range(nn):
               nodes=node_file.readline().strip().split()
               self.node_coordinates.append((float(nodes[1]),float(nodes[2])))
               self.node_tags.append(int(nodes[3]))
               self.node_ids.append(int(nodes[0]))
           node_file.close()
           ele_file=open("%s.1.ele"%self.fn,"r")
           elem=ele_file.readline().strip().split()
           ne=int(elem[0])
           for i in range(ne):
               elem=ele_file.readline().strip().split()
               self.triangles_id.append(int(elem[0]))
               self.triangles_nodes.append((int(elem[1]),int(elem[2]),int(elem[3])))

           ele_file.close()
           # os.remove("%s.1.node"%self.fn)
           # os.remove("%s.1.ele"%self.fn)
           # os.remove("%s.ploy"%self.fn)

      def getFinleyDomain(self):
           from esys.finley import ReadMesh
           finley_file=open("%s.msh"%self.fn,"w")
           finley_file.writelines("%s\n2D-Nodes %d\n"%(self.fn,len(self.node_coordinates)))
           for i in range(len(self.node_coordinates)):
             finley_file.writelines("%s %s %s %e %e\n"%(self.node_ids[i],self.node_ids[i],self.node_tags[i],\
                                                        self.node_coordinates[i][0],self.node_coordinates[i][1]))

           finley_file.writelines("Tri3 %d\n"%len(self.triangles_nodes))
           for i in range(len(self.triangles_nodes)):
              finley_file.writelines("%s 0 %s %s %s\n"%(self.triangles_id[i], \
                                                        self.triangles_nodes[i][0], \
                                                        self.triangles_nodes[i][1], \
                                                        self.triangles_nodes[i][2]))
           finley_file.writelines("Line2 %d\n"%0)
           finley_file.writelines("Line2_Contact %d\n"%0)
           finley_file.writelines("Point1 %d\n"%0)
           finley_file.close()
           # get the mesh
           out=ReadMesh("%s.msh"%self.fn)
           # os.remove("%s.msh"%self.fn)
           return out


#=================================
#   Model interfaces:
#================================
class OceanRegionCollector(Model):
      """


      """
      def __init__(self,debug=False):
           Model.__init__(self,debug=debug)
           self.declareParameter(coastline_source="http://jamboree.esscc.uq.edu.au/cgi-bin/doreen/datafile.txt?west=%%west%%&east=%%east%%&south=%%south%%&north=%%north%%&resolution=%%resolution%%&range360=%%range360%%&river=false&city=false&citytyp=wcity.dat&volcano=false&hotspot=false&shoreline=true&state=false&WSMdata=false&plate=false",
                                 bathymetry_source="http://jamboree.esscc.uq.edu.au/cgi-bin/doreen/grdfile.xyz?west=%%west%%&east=%%east%%&south=%%south%%&north=%%north%%&resolution=%%resolution%%",
                                 resolution=1.,
                                 south=0.,
                                 north=10.,
                                 east=0.,
                                 west=20.,
                                 range360=True,
                                 coastline_stream=None,
                                 bathymetry_stream=None)


      def doInitialization(self):
           """
           Initializes the ocean region
           """
           c=self.__mergeParameters(self.coastline_source)
           b=self.__mergeParameters(self.bathymetry_source)
           self.coastline_stream=urllib.urlopen(c)
           self.bathymetry_stream=urllib.urlopen(b)

      def __mergeParameters(self,txt):
           return txt.replace("%%west%%",str(self.west))\
                     .replace("%%east%%",str(self.east))\
                     .replace("%%south%%",str(self.south))\
                     .replace("%%north%%",str(self.north))\
                     .replace("%%resolution%%",str(self.resolution)) \
                     .replace("%%range360%%",str(self.range360).lower())

class Bathymetry(Model):
       """
       generates the bathymetry data within a region on the earth
       """
       def __init__(self,debug=False):
           Model.__init__(self,debug=debug)
           self.declareParameter(source="none",
                                 bathymetry=1.)

       def doInitialization(self):
           """
           Initializes the
           """
           if hasattr(self.source,"readline"):
               f=self.source
           else:
               f=open(filename,"r")
           x_grd_list=[]
           y_grd_list=[]
           data_grd_list=[]
           line=f.readline().strip()
           while line!="":
               v=line.split()
               x_grd_list.append(float(v[0]))
               y_grd_list.append(float(v[1]))
               data_grd_list.append(float(v[2]))
               line=f.readline().strip()
           self.trace("%s data have been read from %s."%(len(data_grd_list),self.source))
           data_grd=numarray.array(data_grd_list)
           x_grd=numarray.array(x_grd_list)
           y_grd=numarray.array(y_grd_list)
           if len(x_grd)<2:
               raise ValueError,"%s: data base is two small"%str(self)
           ox=x_grd[0]
           oy=y_grd[0]
           diam=max(abs(x_grd[len(x_grd)-1]-ox),abs(y_grd[len(y_grd)-1]-oy))
           dx=x_grd[1]-ox
           nx=1
           nx=1
           while abs(x_grd[nx]-ox)>1.e-10*diam:
               nx+=1
           dy=y_grd[nx]-oy
           ny=len(x_grd)/nx
           data_grd.resize((ny,nx))
           self.bathymetry=GridData(s=[dx,dy],o=[ox,oy],n=[nx,ny],data=data_grd)
           self.trace("%s x %s grid with %s x %s spacing."%(nx,ny,dx,dy))


class OceanRegion(Model):
       """
       generates the ocean region with a coast line and a bathymetry

       """
       def __init__(self,debug=False):
           Model.__init__(self,debug=debug)
           self.declareParameter(domain=None, \
                                 resolution=1.,
                                 south=0.,
                                 north=10.,
                                 east=0.,
                                 west=20.,
                                 bathymetry=None,
                                 bathymetry_data=None,
                                 coastline=None,
                                 source="none")

       def doInitialization(self):
           """
           Initializes the ocean region
           """
           if hasattr(self.source,"readline"):
               f=self.source
               data_name=f.geturl()
           else:
               f=open(self.source,"r")
               data_name=self.source

           segs=[]
           name=""
           line=f.readline()
           while not line=="":
              line=line.strip()
              if line[:7]=="#region":
                 data=line[line.index("[")+1:line.index("]")].split(",")
                 reg=RegionOnEarthSurface(PointOnEarthSurface(self.south,self.west),PointOnEarthSurface(self.north,self.east),self.resolution)
                 self.trace(str(reg))
                 self.coastline=Coastline(region=reg,name=data_name)
              elif line.find("Shore Bin")>-1:
                   self.coastline.append(Polyline(segs,name))
                   name=line[2:]
                   segs=[]
              if not (line=="" or line[0]=="#" or line[0]==">") :
                  x=line.split()
                  segs.append(PointOnEarthSurface(float(x[0]),float(x[1])))
              line=f.readline()
           self.coastline.append(Polyline(segs,name))
           d=self.bathymetry_data.interpolate([[self.east,self.south],[self.west,self.south],[self.east,self.north],[self.west,self.north]])
           self.domain=self.coastline.makeTriangulation(south_east_is_water=d[0]<0,south_west_is_water=d[1]<0,north_east_is_water=d[2]<0,north_west_is_water=d[3]<0).getFinleyDomain()
           self.bathymetry=maximum(-self.bathymetry_data.interpolate(Function(self.domain).getX()),0.)


class TsunamiSource(Model):
       """
       defines a wave in Gaussean form between start and end.
       """
       GAMMA=0.05
       def __init__(self,debug=False):
           Model.__init__(self,debug=debug)
           self.declareParameter(domain=None, 
                                 start_lat=-10.,
                                 start_long=110.,
                                 end_lat=-12.,
                                 end_long=120.,
                                 width=5.,
                                 decay_zone=0.1,
                                 amplitude=1.,
                                 wave_height=1.)

       def doInitialization(self):
           """
           set initial wave form
           """
           beta=math.sqrt(-math.log(self.GAMMA))/self.decay_zone
           x=self.domain.getX()
           x0=x[0]
           x1=x[1]
           mid0=(self.start_long+self.end_long)/2
           mid1=(self.start_lat+self.end_lat)/2
           dist=math.sqrt((mid0-self.end_long)**2+(mid1-self.end_lat)**2)
           a=(self.start_long-mid0)/dist
           b=(self.end_long-mid1)/dist
           self.trace("source length = %s"%(dist*2.))

           x0_hat= a*(x0-mid0)+b*(x1-mid1)
           x1_hat=-b*(x0-mid0)+a*(x1-mid1)
           # x1-direction
           s=abs(x1_hat)-self.width
           m=s.whereNegative()
           f1=(1.-m)*exp(-(s*beta)**2)+m
           # x0-direction
           s=abs(x0_hat)-dist
           m=s.whereNegative()
           f0=(1.-m)*exp(-(s*beta)**2)+m
           self.wave_height=f1*f0*self.amplitude

#====================================================================================================================================================
class TsunamiInDeepWater(Model):
       """
       Runs the deep water tsunami model based on a simplfied form of the shallow water equation.


       M{d^2 h/dt^2 =div(c grad(h)) }

       where h is the wave height above sea level, and with H the depth of the water level and g is the gravity constant
       c=sqrt(g*H).

       The simulation uses the Verlet scheme.

       """
       def __init__(self,debug=False):
           Model.__init__(self,debug=debug)
           self.declareParameter(domain=None, \
                                 wave_height=1.,
                                 wave_velocity=0.,
                                 initial_time_step=None,
                                 bathymetry=1.,
                                 safety_factor=0.1)

       def doInitialization(self):
           """
           Initializes the time integartion scheme
           """
           self.__pde=LinearPDE(self.domain)
           self.__pde.setSolverMethod(self.__pde.LUMPING)
           self.__pde.setValue(D=1.)
           self.__c2=RegionOnEarthSurface.GRAVITY*self.bathymetry/RegionOnEarthSurface.RADIUS**2
           c_max=math.sqrt(Lsup(self.__c2))
           self.__dt=self.safety_factor*inf(self.domain.getSize()/(sqrt(self.__c2)+EPS*c_max))
           if self.initial_time_step==None: self.initial_time_step=self.__dt
           self.trace("maximum wave velocity %s m/sec^2"%c_max)
           self.trace("Time step size is %s sec"%self.__dt)


       def getSafeTimeStepSize(self,dt):
           """
           returns new step size

           @param dt: last time step size used
           @type dt: C{float}
           @return: time step size that can savely be used
           @rtype: C{float}
           """
           return self.__dt

       def doStepPostprocessing(self,dt):
           """
           perform the time step using the Valet scheme

           @param dt: time step size to be used
           @type dt: C{float}
           """
           self.__pde.setValue(X=-self.__c2*grad(self.wave_height))
           new_height=self.wave_height+dt*self.wave_velocity+dt*(self.initial_time_step+dt)/2*self.__pde.getSolution()
           self.wave_velocity=(new_height-self.wave_height)/dt
           self.wave_height=new_height
           self.initial_time_step=dt
           self.trace("Wave height range is %e %e"%(inf(self.wave_height),sup(self.wave_height)))

class SurfMovie(Model):
       """
       movie from a wave propagation on the sea

       @ivar time: current time
       @ivar bathymetry: scalar data set
       @ivar wave_height: vector data set
       @ivar filename: name of the movie file
       """
       def __init__(self,debug=False):
           Model.__init__(self,debug=debug)
           self.declareParameter(bathymetry=1.,
                                 wave_height=1.,
                                 coastline=None,
                                 t=0.,
                                 dt=1.,
                                 south=2.,
                                 north=5.,
                                 east=3.,
                                 west=15.,
                                 filename="movie.mpg")

       def doInitialization(self):
          """
          Initializes the time integartion scheme
          """
          self.__fn=os.tempnam()+".xml"
          self.__frame_name=os.tempnam()
          self.__next_t=self.dt
          # self.coastline.getVTK()
          # self.bathymetry.getVTK()
          # wndow(south,west,north,east)

       def doStepPostprocessing(self, dt):
        """
        Does any necessary postprocessing after each step

        @param dt:
        """
        if self.t>=self.__next_t:
             print self.t,"write ",Lsup(self.wave_height)
             saveVTK(self.__fn,h=self.wave_height)
             # vtkobj=...
             # save(self.__frame_name)
             self.__next_t+=self.dt

       def getSafeTimeStepSize(self,dt):
           """
           returns new step size

           @param dt: last time step size used
           @type dt: C{float}
           @return: time step size that can savely be used
           @rtype: C{float}
           """
           return self.__next_t-self.t

       def doFinalization(self):
          """
          Finalises the visualisation.  For instance, makes a movie of the
          image files.
          """
          # make the movie into self.filename
          pass


if __name__=="__main__":
   from esys.escript.modelframe import Link,Simulation
   from esys.modellib.input import Sequencer

   oc=OceanRegionCollector()
   oc.resolution=0.5
   oc.south=-45.5
   oc.north=-5.4
   oc.east=161.1
   oc.west=108.2
   oc.range360=True


   b=Bathymetry()
   b.source=Link(oc,"bathymetry_stream")

   oreg=OceanRegion()
   oreg.source=Link(oc,"coastline_stream")
   oreg.resolution=Link(oc,"resolution")
   oreg.south=Link(oc,"south")
   oreg.north=Link(oc,"north")
   oreg.east=Link(oc,"east")
   oreg.west=Link(oc,"west")
   oreg.bathymetry_data=Link(b,"bathymetry")

   src=TsunamiSource()
   src.domain=Link(oreg,"domain")
   src.start_lat=-10.
   src.start_long=110.
   src.end_lat=-12.
   src.end_long=120.
   src.width=0.1
   src.decay_zone=0.01
   src.amplitude=1.

   ts=TsunamiInDeepWater()
   ts.domain=Link(oreg,"domain")
   ts.wave_height=Link(src,"wave_height")
   ts.wave_velocity=0.
   ts.bathymetry=Link(oreg,"bathymetry")

   sq=Sequencer()
   sq.t_end=1000.

   sm=SurfMovie()
   sm.bathymetry=Link(b,"bathymetry")
   sm.wave_height=Link(ts,"wave_height")
   sm.coastline=Link(oreg,"coastline")
   sm.t=Link(sq,"t")
   sm.dt=100.
   sm.filename="movie.mpg"
   
   s=Simulation([sq,oc,b,oreg,src,ts,sm])
   s.writeXML()
   s.run()
