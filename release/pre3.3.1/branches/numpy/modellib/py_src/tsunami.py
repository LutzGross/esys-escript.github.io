
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import os, sys
import vtk
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.escript.modelframe import Model
import numpy
import math
from tempfile import mkstemp

WORKDIR="./work/"
EPS=1.e-8

#=============================================================================

class GridData:
    """
    This object is used to store data on a grid.
    It will be replaced by Bruce at a later stage.

    data[i,j] are x=j*s[0]+o[0] and y=i*s[1]+o[1]

    for 0<=j<n[0] and 0<=i<n[1]
    """
    def __init__(self, s, o, n, data):
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
        if hasattr(x, "getNumberOfDataPoints"):
            x_shape0 = x.getNumberOfDataPoints()
            return_data_object = True
        else:
            x_array = numpy.array(x)
            x_shape0 = x_array.shape[0]
            return_data_object = False

        data=numpy.zeros(x_shape0, numpy.float_)
        ox,oy = self.getOrigin()
        dx,dy = self.getSpacing()
        data_array = self.getData()
        i_dx = 1
        i_dy = 1
        for i in range(x_shape0):
            if return_data_object:
                x_array = x.getTupleForDataPoint(i)
                x_long = x_array[0]-ox
                x_lat = x_array[1]-oy
            else:
                x_long = x_array[i,0]-ox
                x_lat = x_array[i,1]-oy
            j0 = min(max(int(x_long/dx),0), data_array.shape[1]-1-i_dy)
            j1 = min(max(int(x_lat/dy),0), data_array.shape[0]-1-i_dx)
            f01 = (x_long-j0*dx)/dx
            f00 = 1.-f01
            f11 = (x_lat-j1*dy)/dy
            f10 = 1.-f11
            H00 = data_array[j1,j0]
            H01 = data_array[j1,j0+i_dx]
            H11 = data_array[j1+i_dy,j0+i_dx]
            H10 = data_array[j1+i_dy,j0]
            data[i] = (H00*f00+H01*f01)*f10 + (H10*f00+H11*f01)*f11

        if return_data_object:
            dataObj = Scalar(0, x.getFunctionSpace(), True)
            for i in range(x_shape0):
                dataObj.setValueOfDataPoint(i, data[i])
            return dataObj
        else:
            return data


class PointOnEarthSurface:
    """
    Coordinates of a point on the surface of the Earth
    """
    def __init__(self, long=0, lat=0):
        self.long=long
        self.lat=lat

    def __str__(self):
        return "(%s,%s)"%(self.long,self.lat)

    def __sub__(self, other):
        return self.dist(other)

    def split(self,p,t):
        return PointOnEarthSurface(long=self.long+t*(p.long-self.long),
                                   lat=self.lat+t*(p.lat-self.lat) )

    def midPoint(self,other):
        return PointOnEarthSurface(long=(self.long+other.long)/2,
                                   lat=(self.lat+other.lat)/2 )

    def dist(self,other):
        return math.sqrt((self.long-other.long)**2+(self.lat-other.lat)**2)


class RegionOnEarthSurface:
    """
    Defines a region by a south-east and north-west corner
    """
    RADIUS=6378.8e3
    GRAVITY=9.81
    WEST=0
    NORTH=1
    EAST=2
    SOUTH=3
    def __init__(self, west_south=PointOnEarthSurface(), east_north=PointOnEarthSurface(), resolution=1.):
        if resolution <= 0:
            raise ValueError, "resolution must be positive"
        if west_south.long >= east_north.long:
            raise ValueError, "south-west corner must be west of north-east corner"
        if west_south.lat >= east_north.lat:
            raise ValueError, "south-east corner must be south of north-west corner"
        if east_north.lat-west_south.lat < resolution/2:
            raise ValueError, "latitude length of region must be at least 2*larger than the resolution"
        if  east_north.long-west_south.long < resolution/2:
            raise ValueError, "longitude length of region must be at least 2*larger than the resolution"

        self.west_south = west_south
        self.east_north = east_north
        self.resolution = resolution

    def __str__(self):
        return "RegionOnEarthSurface between %s and %s" % (str(self.west_south), str(self.east_north))

    def isOnFace(self, p):
        return self.isOnThisFace(p,self.SOUTH) or \
               self.isOnThisFace(p,self.NORTH) or \
               self.isOnThisFace(p,self.WEST) or \
               self.isOnThisFace(p,self.EAST)

    def isOnThisFace(self, p, face):
        if face==self.WEST:
            return self.west_south.long==p.long
        if face==self.EAST:
            return p.long==self.east_north.long
        if face==self.SOUTH:
            return self.west_south.lat==p.lat
        if face==self.NORTH:
            return p.lat==self.east_north.lat

    def isBoxVertex(self, p):
        return ( self.isOnThisFace(p,self.WEST) and self.isOnThisFace(p,self.SOUTH) ) or \
               ( self.isOnThisFace(p,self.WEST) and self.isOnThisFace(p,self.NORTH) ) or \
               ( self.isOnThisFace(p,self.EAST) and self.isOnThisFace(p,self.NORTH) ) or \
               ( self.isOnThisFace(p,self.EAST) and self.isOnThisFace(p,self.SOUTH) )

    def getFace(self, p):
        # order is critical
        if self.west_south.long == p.long: return self.WEST
        if p.lat == self.east_north.lat: return self.NORTH
        if p.long == self.east_north.long: return self.EAST
        if self.west_south.lat == p.lat: return self.SOUTH

    def comparePointsOnFace(self, p0, p1):
        f0 = self.getFace(p0)
        f1 = self.getFace(p1)

        if f0 < f1:
            return -1
        elif f0 > f1:
            return 1
        else:
            if f0 == self.WEST:
                if p0.lat < p1.lat:
                    return -1
                elif p0.lat > p1.lat:
                    return 1
                else:
                    return 0
            elif f0 == self.EAST:
                if p0.lat < p1.lat:
                    return 1
                elif p0.lat > p1.lat:
                    return -1
                else:
                    return 0
            elif f0 == self.NORTH:
                if p0.long < p1.long:
                    return -1
                elif p0.long > p1.long:
                    return 1
                else:
                    return 0
            else:
                if p0.long < p1.long:
                    return 1
                elif p0.long > p1.long:
                    return -1
                else:
                    return 0

    def isInRegion(self, p):
        return self.west_south.long <= p.long \
            and p.long <= self.east_north.long \
            and self.west_south.lat <= p.lat \
            and p.lat <= self.east_north.lat

    def cutSegment(self, p0, p1):
        t = None
        p = None
        tmp = self.interceptOfSegment(p0, p1, d=0, v=self.west_south.long)
        if not tmp == None:
            p_tmp = p0.split(p1, tmp)
            if self.isInRegion(p_tmp) and t < tmp:
                t = tmp
                p = p_tmp

        tmp = self.interceptOfSegment(p0, p1, d=0, v=self.east_north.long)
        if not tmp == None:
            p_tmp = p0.split(p1, tmp)
            if self.isInRegion(p_tmp) and t < tmp:
                t = tmp
                p = p_tmp

        tmp = self.interceptOfSegment(p0, p1, d=1, v=self.west_south.lat)
        if not tmp == None:
            p_tmp = p0.split(p1, tmp)
            if self.isInRegion(p_tmp) and t < tmp:
                t = tmp
                p = p_tmp

        tmp = self.interceptOfSegment(p0, p1, d=1, v=self.east_north.lat)
        if not tmp == None:
            p_tmp = p0.split(p1, tmp)
            if self.isInRegion(p_tmp) and t < tmp:
                t = tmp
                p = p_tmp
        return p

    def interceptOfSegment(self, p0, p1, d=0, v=0.):
        """
        Finds t in [0,1] such that p0+t*(p0-p1) cuts x[d]=v. If it does not
        exist None is returned.
        """
        if d == 0:
            a = p0.long
            b = p1.long
        else:
            a = p0.lat
            b = p1.lat
        if b == a:
            if a == v:
                t = 0
            else:
                t = None
        else:
            t = (v-a)/(b-a)
            if not (0<=t and t<=1):
                t = None
        return t


class Polyline:
    """
    Defines a set of segments through a list of coordinates
    """
    def __init__(self, list_of_coordinates=[], name="none"):
        c=[]
        if len(list_of_coordinates)>0:
            for i in range(len(list_of_coordinates)-1):
                if list_of_coordinates[i]-list_of_coordinates[i+1] > 0:
                    c.append(list_of_coordinates[i])
            c.append(list_of_coordinates[-1])
        self.list_of_coordinates = c
        self.name = name

    def getDiameter(self):
        out = 0.
        for p in self.list_of_coordinates:
            for q in self.list_of_coordinates:
                out = max(out, p-q)
        return out

    def isLoop(self):
        if len(self) > 0:
            return not self[0]-self[-1]>EPS
        else:
            return True

    def insert(self, index, coordinate):
        """
        Inserts point before index
        """
        if self.list_of_coordinates[index]-coordinate < EPS:
            return index
        elif self.list_of_coordinates[index-1]-coordinate < EPS:
            return index-1
        else:
            self.list_of_coordinates.insert(index, coordinate)
            return index

    def split(self, index):
        """
        Splits the polyline at point index
        """
        return Polyline(self.list_of_coordinates[:index], self.name), \
                Polyline(self.list_of_coordinates[index:], self.name)

    def __getitem__(self, s):
        return self.list_of_coordinates.__getitem__(s)

    def __iter__(self):
        return iter(self.list_of_coordinates)

    def __str__(self):
        if self.isLoop():
            out = "loop["
        else:
            out = "["
        for i in self:
            if out[-1] == "[":
                out += "%s" % str(i)
            else:
                out += ",%s" % str(i)
        return out+"]"

    def __len__(self):
        return len(self.list_of_coordinates)

    def orientation(self):
        """
        Returns the orientation of the polyline
        """
        if not self.isLoop():
            raise TypeError, "polyline is not a loop"

        integ = 0.
        for i in range(len(self.list_of_coordinates)-1):
            p0 = self.list_of_coordinates[i]
            p1 = self.list_of_coordinates[i+1]
            integ += (p1.lat-p0.lat)*(p1.long+p0.long)-(p1.long-p0.long)*(p1.lat+p0.lat)
        if integ >= 0.:
            return 1.
        else:
            return -1.

    def givePositiveOrientation(self):
        if self.orientation() < 0:
            self.list_of_coordinates.reverse()


class Coastline:
    """
    Defines a coast line by a Polyline within a RegionOnEarthSurface
    """
    def __init__(self, region, name="none"):
        self.region = region
        self.name = name
        self.polylines = []

    def __str__(self):
        out = self.name + " in " + str(self.region)
        for pl in self.polylines:
            out += "\n" + str(pl)
        return out

    def makeTriangulation(self, west_south_is_water=True, east_south_is_water=True, west_north_is_water=True, east_north_is_water=True):
        self.clean()
        vertices=[]
        segments=[]
        holes=[]
        vertices_on_face=[]
        for pl in self.polylines:
            if pl.getDiameter() > self.region.resolution:
                short_pl = [pl[0]]
                for i in range(1, len(pl)):
                    if pl[i]-short_pl[-1] > 0*EPS+self.region.resolution/10:
                        short_pl.append(pl[i])
                    elif i == len(pl)-1:
                        short_pl[-1] = pl[i]
                if pl.isLoop():
                    if len(short_pl) > 3:
                        offset = len(vertices)
                        v_tmp = [short_pl[0]]
                        s_tmp = []
                        for i in range(1, len(short_pl)):
                            if i == len(short_pl)-1:
                                s_tmp.append((len(v_tmp)-1+offset, offset))
                            else:
                                s_tmp.append((len(v_tmp)-1+offset, len(v_tmp)+offset))
                                v_tmp.append(short_pl[i])
                        vertices += v_tmp
                        segments += s_tmp
                        # find a point in the loop:
                        d_long = v_tmp[1].long-v_tmp[0].long
                        d_lat = v_tmp[1].lat-v_tmp[0].lat
                        l = math.sqrt(d_long**2+d_lat**2)
                        mid = v_tmp[0].midPoint(v_tmp[1])
                        n_long = -d_lat/l
                        n_lat = d_long/l
                        s = self.region.resolution
                        for seg in s_tmp:
                            p0 = vertices[seg[0]]
                            p1 = vertices[seg[1]]
                            a_long = p1.long-p0.long
                            a_lat = p1.lat-p0.lat
                            d = a_lat*n_long-a_long*n_lat
                            if abs(d) > 0.:
                                t = ((mid.lat-p0.lat)*n_long-(mid.long-p0.long)*n_lat)/d
                                if 0<=t and t<=1:
                                    s_tmp = ((mid.lat-p0.lat)*a_long-(mid.long-p0.long)*a_lat)/d
                                    if s_tmp > EPS:
                                        s = min(s,s_tmp)
                        h = PointOnEarthSurface(long=mid.long+s/2*n_long,lat=mid.lat+s/2*n_lat)
                        holes.append(h)
                else:
                    if len(short_pl) > 1:
                        if self.region.isOnFace(short_pl[0]):
                            vertices_on_face.append(short_pl[0])
                        if self.region.isOnFace(short_pl[-1]):
                            vertices_on_face.append(short_pl[-1])
                        vertices.append(short_pl[0])
                        for i in range(1, len(short_pl)):
                            segments.append((len(vertices)-1, len(vertices)))
                            vertices.append(short_pl[i])
        # put into the bounding box:
        new_vertices = []
        if west_south_is_water:
            new_vertices.append(PointOnEarthSurface(long=self.region.west_south.long, lat=self.region.west_south.lat))
        if east_south_is_water:
            new_vertices.append(PointOnEarthSurface(long=self.region.east_north.long, lat=self.region.west_south.lat))
        if west_north_is_water:
            new_vertices.append(PointOnEarthSurface(long=self.region.west_south.long, lat=self.region.east_north.lat))
        if east_north_is_water:
            new_vertices.append(PointOnEarthSurface(long=self.region.east_north.long, lat=self.region.east_north.lat))

        # add new vertices if they don't exist yet
        for q in new_vertices:
            for p2 in vertices_on_face:
                if p2-q < EPS:
                    q = None
                    raise ValueError, "coast line crosses boundary box vertex. This case is currently not supported."
            if not q == None:
                vertices.append(q)
                vertices_on_face.append(q)
        vertices_on_face.sort(self.region.comparePointsOnFace)
        index = 0
        walking_on_water = True
        l = len(vertices_on_face)
        while index < l:
            p1 = vertices_on_face[(index+1)%l]
            p0 = vertices_on_face[index]
            if walking_on_water:
                segments.append((vertices.index(p0), vertices.index(p1)))
                walking_on_water = False
            else:
                if self.region.isBoxVertex(p0):
                    segments.append((vertices.index(p0), vertices.index(p1)))
                else:
                    walking_on_water = True
            index += 1
        return EarthTriangulation(vertices, segments, holes, self.region.resolution)

    def clean(self):
        """
        Cleans up the coastline by joining polylines to loops or connecting
        faces of the region
        """
        # find poylines that are linked
        while True:
            k0 = None
            for pl in self.polylines:
                if not pl.isLoop():
                    for k in [0,-1]:
                        for ql in self.polylines:
                            if not (ql==pl or ql.isLoop()):
                                for k2 in [0,-1]:
                                    if ql[k2]-pl[k] < EPS:
                                        pl0 = pl
                                        pl1 = ql
                                        k0 = k
                                        k1 = k2
                                        break
                            if not k0==None: break # ql
                        if not k0==None: break # k
                if not k0==None: break # pl

            if k0 == None:
                break
            else:
                self.polylines.remove(pl0)
                self.polylines.remove(pl1)
                pl0c = pl0.list_of_coordinates
                pl1c = pl1.list_of_coordinates
                if k0 == 0: pl0c.reverse()
                if k1 == -1: pl1c.reverse()
                pl = Polyline(pl0c+pl1c[1:], pl0.name+" + "+pl1.name)
                self.append(pl)

        # find a polyline that is not a loop and has an end or start point
        # not on the face of the region:
        while True:
            pl = None
            k = None
            for pl2 in self.polylines:
                if not pl2.isLoop():
                    pl = pl2
                    if not self.region.isOnFace(pl[0]): k=0
                    if not self.region.isOnFace(pl[-1]): k=-1
                    if not k == None: break
            if k == None: break
            self.polylines.remove(pl)
            d_min = 50000.
            k_min = None
            for pl2 in self.polylines:
                if not pl2.isLoop():
                    for k2 in [0,-1]:
                        if not self.region.isOnFace(pl2[k2]):
                            d2 = pl2[k2]-pl[k]
                            if d2 < d_min:
                                d_min = d2
                                pl_min = pl2
                                k_min = k2
            if k_min == None:
                raise ValueError, "cannot link coastline %s to any other coastline." % pl.name
            plc = pl.list_of_coordinates
            plc_min = pl_min.list_of_coordinates
            if k == 0: plc.reverse()
            if k_min == -1: plc_min.reverse()
            if d_min < EPS:
                new_pl = Polyline(plc+plc_min[1:], pl.name+" + "+pl_min.name)
            else:
                new_pl = Polyline(plc+plc_min, pl.name+" + "+pl_min.name)
            self.polylines.remove(pl_min)
            self.append(new_pl)
        # give positive orientation to loops:
        for pl in self.polylines:
            if pl.isLoop(): pl.givePositiveOrientation()

    def append(self, polyline=Polyline()):
        """
        Appends a polyline
        """
        if len(polyline) > 1:
            pl = []
            outside_region = None
            for i in range(len(polyline)):
                if not self.region.isInRegion(polyline[i]):
                    outside_region = i
                    break
                #pl.append(self.region.nudgeToFace(polyline[i]))
                pl.append(polyline[i])
            if not outside_region == None:
                if outside_region == 0:
                    for i in range(outside_region+1, len(polyline)):
                        if self.region.isInRegion(polyline[i]):
                            polyline.insert(i,self.region.cutSegment( \
                                    polyline[i-1], polyline[i]))
                            pl1 = polyline.split(i)[1]
                            self.append(pl1)
                            break
                else:
                    # split polyline in two parts: first is fully within the
                    # region, the other starts with point outside the region
                    c = self.region.cutSegment(polyline[outside_region-1], polyline[outside_region])
                    i = polyline.insert(outside_region, c)
                    pl0,pl1 = polyline.split(i+1)
                    self.append(pl0)
                    self.append(pl1)
            else:
                if len(pl) > 1:
                    pply = Polyline(pl,polyline.name)
                    self.polylines.append(pply)


class EarthTriangulation:
    """
    Generates earth mesh and triangulates it
    """
    GENERATOR = "triangle -pqa%g %s"

    def __init__(self,vertices=[], segments=[], holes=[], resolution=1.):
        (fd, self.fn) = mkstemp(suffix=".poly")
        self.fn = self.fn.replace(".poly", "")
        f=os.fdopen(fd, "w")
        # write triangle input file
        f.writelines("%d %d %d %d\n" % (len(vertices),2,0,0))
        for i in range(len(vertices)):
            f.writelines("%d %e %e\n" % (i,vertices[i].long,vertices[i].lat))
        f.writelines("%d %d\n"%(len(segments),0))
        for i in range(len(segments)):
            f.writelines("%d %d %d\n" % (i,segments[i][0],segments[i][1]))
        f.writelines("%d\n"%(len(holes)))
        for i in range(len(holes)):
            f.writelines("%d %e %e\n"%(i,holes[i].long,holes[i].lat))
        f.close()

        # create mesh using generator tool
        os.system(self.GENERATOR % (resolution**2, self.fn))

        # read the resulting mesh file
        self.node_coordinates=[]
        self.node_tags=[]
        self.node_ids=[]
        self.triangles_nodes=[]
        self.triangles_id=[]
        node_file=open("%s.1.node" % self.fn, "r")
        nodes=node_file.readline().strip().split()
        nn=int(nodes[0])
        for i in range(nn):
            nodes=node_file.readline().strip().split()
            self.node_coordinates.append((float(nodes[1]),float(nodes[2])))
            self.node_tags.append(int(nodes[3]))
            self.node_ids.append(int(nodes[0]))
        node_file.close()
        ele_file=open("%s.1.ele" % self.fn, "r")
        elem=ele_file.readline().strip().split()
        ne=int(elem[0])
        for i in range(ne):
            elem=ele_file.readline().strip().split()
            self.triangles_id.append(int(elem[0]))
            self.triangles_nodes.append((int(elem[1]),int(elem[2]),int(elem[3])))
        ele_file.close()

        # clean up
        os.remove("%s.poly" % self.fn)
        os.remove("%s.1.poly" % self.fn)
        os.remove("%s.1.node" % self.fn)
        os.remove("%s.1.ele" % self.fn)

    def getFinleyDomain(self):
        from esys.finley import ReadMesh
        finley_file = open("%s.msh" % self.fn, "w")
        finley_file.writelines("%s\n2D-Nodes %d\n" % (self.fn,
                               len(self.node_coordinates)))
        for i in range(len(self.node_coordinates)):
            finley_file.writelines("%s %s %s %e %e\n" % (self.node_ids[i],
                                   self.node_ids[i],self.node_tags[i],
                                   self.node_coordinates[i][0],
                                   self.node_coordinates[i][1]))

        finley_file.writelines("Tri3 %d\n" % len(self.triangles_nodes))
        for i in range(len(self.triangles_nodes)):
            finley_file.writelines("%s 0 %s %s %s\n" % (
                                   self.triangles_id[i],
                                   self.triangles_nodes[i][0],
                                   self.triangles_nodes[i][1],
                                   self.triangles_nodes[i][2]))
        finley_file.writelines("Line2 0\n")
        finley_file.writelines("Line2_Contact 0\n")
        finley_file.writelines("Point1 0\n")
        finley_file.close()
        # read the mesh with finley
        out=ReadMesh("%s.msh" % self.fn)
        os.remove("%s.msh" % self.fn)
        return out


#=================================
#        Model interfaces        #
#=================================
class OceanRegionCollector(Model):
    """
    Opens input streams for coastline and bathymetry data (generated by GMT)
    """

    def __init__(self,**kwargs):
        Model.__init__(self,**kwargs)
        self.declareParameter(coastline_source = "coastline_default.dat",
                              bathymetry_source = "bathymetry_default.dat",
                              resolution = 1.,
                              south = 0.,
                              north = 10.,
                              east = 0.,
                              west = 20.,
                              range360 = True,
                              coastline_stream = None,
                              bathymetry_stream = None)

    def doInitialization(self):
        """
        Initializes the ocean region by reading coordinate files
        """
        # We presume that %%'s in the source name requires calling a command
        # otherwise it is a filename containing the required data
        if self.coastline_source.find("%%") != -1:
            c=self.__mergeParameters(self.coastline_source)
            self.coastline_stream=os.popen(c).read()
        else:
            self.coastline_stream=self.coastline_source

        if self.bathymetry_source.find("%%") != -1:
            b=self.__mergeParameters(self.bathymetry_source)
            self.bathymetry_stream=os.popen(b).read()
        else:
            self.bathymetry_stream=self.bathymetry_source

    def __mergeParameters(self,txt):
        return txt.replace("%%west%%", str(self.west)) \
                  .replace("%%east%%", str(self.east)) \
                  .replace("%%south%%", str(self.south)) \
                  .replace("%%north%%", str(self.north)) \
                  .replace("%%resolution%%", str(self.resolution)) \
                  .replace("%%range360%%", str(self.range360).lower())


#=============================================================================
class Bathymetry(Model):
    """
    Generates the bathymetry data within a region on the earth
    """
    def __init__(self,**kwargs):
        Model.__init__(self,**kwargs)
        self.declareParameter(source = "none",
                              bathymetry = 1.)

    def doInitialization(self):
        """
        Initializes the bathymetry grid
        """

        print "Intializing bathymetry..."
        if hasattr(self.source,"readline"):
            f=self.source
        else:
            f=open(self.source,"r")
        x_grd_list=[]
        y_grd_list=[]
        data_grd_list=[]
        line=f.readline().strip()
        while line != "":
            v=line.split()
            x_grd_list.append(float(v[0]))
            y_grd_list.append(float(v[1]))
            data_grd_list.append(float(v[2]))
            line=f.readline().strip()
        self.trace("%s lines have been read from %s." % (len(data_grd_list), self.source))
        data_grd=numpy.array(data_grd_list)
        x_grd=numpy.array(x_grd_list)
        y_grd=numpy.array(y_grd_list)
        if len(x_grd)<2:
            raise ValueError,"%s: data base is too small"%str(self)
        ox=x_grd[0]
        oy=y_grd[0]
        diam=max(abs(x_grd[len(x_grd)-1]-ox),abs(y_grd[len(y_grd)-1]-oy))
        dx=x_grd[1]-ox
        nx=1
        nx=1
        while abs(x_grd[nx]-ox)>EPS*diam:
            nx+=1
        dy=y_grd[nx]-oy
        ny=len(x_grd)/nx
        data_grd.resize((ny,nx))
        self.bathymetry=GridData(s=[dx,dy],o=[ox,oy],n=[nx,ny],data=data_grd)
        self.trace("%sx%s grid with %sx%s spacing." % (nx,ny,dx,dy))


#=============================================================================
class OceanRegion(Model):
    """
    Generates the ocean region with a coast line and a bathymetry
    """
    def __init__(self,**kwargs):
        Model.__init__(self,**kwargs)
        self.declareParameter(domain = None,
                              resolution = 1.,
                              south = 0.,
                              north = 10.,
                              east = 0.,
                              west = 20.,
                              bathymetry = None,
                              bathymetry_data = None,
                              coastline = None,
                              source = "none")

    def doInitialization(self):
        """
        Initializes the ocean region
        """

        # create output directory if necessary
        if not os.path.isdir(WORKDIR): os.mkdir(WORKDIR)

        print "Initializing ocean region..."
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
                reg=RegionOnEarthSurface(PointOnEarthSurface(lat=self.south,long=self.west),PointOnEarthSurface(lat=self.north,long=self.east),self.resolution)
                self.coastline=Coastline(region=reg,name=data_name)
            elif line.find("Shore Bin")>-1:
                self.coastline.append(Polyline(segs,name))
                name=line[2:]
                segs=[]
            if not (line=="" or line[0]=="#" or line[0]==">") :
                x=line.split()
                segs.append(PointOnEarthSurface(long=float(x[0]),lat=float(x[1])))
            line=f.readline()
        self.coastline.append(Polyline(segs,name))
        d=self.bathymetry_data.interpolate([[self.east,self.south],
                                            [self.west,self.south],
                                            [self.east,self.north],
                                            [self.west,self.north]])
        self.domain = self.coastline.makeTriangulation(
                                east_south_is_water=d[0]<=0,
                                west_south_is_water=d[1]<=0,
                                east_north_is_water=d[2]<=0,
                                west_north_is_water=d[3]<=0
                      ).getFinleyDomain()
        self.domain.dump(os.path.join(WORKDIR, "coastline.nc"))
        self.bathymetry = maximum(-self.bathymetry_data.interpolate(Function(self.domain).getX()),0.)


#=============================================================================
class TsunamiSource(Model):
    """
    Defines a wave in Gaussian form between start and end.
    """
    GAMMA=0.05
    def __init__(self,**kwargs):
        Model.__init__(self,**kwargs)
        self.declareParameter(domain = None, 
                              start_lat = -10.,
                              start_long = 110.,
                              end_lat = -12.,
                              end_long = 120.,
                              width = 5.,
                              decay_zone = 0.1,
                              amplitude = 1.,
                              wave_height = 1.)

    def doInitialization(self):
        """
        Sets the initial wave form
        """
        beta=math.sqrt(-math.log(self.GAMMA))/self.decay_zone
        x=self.domain.getX()
        x_long = x[0]
        x_lat = x[1]
        mid_long = (self.start_long+self.end_long)/2
        mid_lat = (self.start_lat+self.end_lat)/2
        dist = math.sqrt((mid_long-self.end_long)**2+(mid_lat-self.end_lat)**2)
        a = (self.start_long-mid_long)/dist
        b = (self.start_lat-mid_lat)/dist
        self.trace("source length = %s" % (dist*2.))
        x_long_hat = a*(x_long-mid_long)+b*(x_lat-mid_lat)
        x_lat_hat = -b*(x_long-mid_long)+a*(x_lat-mid_lat)
        # x_lat-direction
        s = abs(x_lat_hat)-self.width
        m = whereNegative(s)
        f1 = (1.-m)*exp(-(s*beta)**2)+m

        # x_long-direction
        s = abs(x_long_hat)-dist
        m = whereNegative(s)
        f0 = (1.-m)*exp(-(s*beta)**2)+m
        self.wave_height = f1*f0*self.amplitude


#=============================================================================
class TsunamiInDeepWater(Model):
    """
    Runs the deep water tsunami model based on a simplified form of the
    shallow water equation.

    M{d^2 h/dt^2 =div(c grad(h)) }

    where h is the wave height above sea level, and c=sqrt(g*H),
    with H - depth of the water level, g - gravity constant

    The simulation uses the Verlet scheme.
    """
    def __init__(self,**kwargs):
        Model.__init__(self,**kwargs)
        self.declareParameter(domain = None,
                              wave_height = 1.,
                              wave_velocity = 0.,
                              initial_time_step = None,
                              bathymetry = 1.,
                              safety_factor = 1.)

    def doInitialization(self):
        """
        Initializes the time integration scheme
        """

        self.__pde = LinearPDE(self.domain)
        self.__pde.setSolverMethod(self.__pde.LUMPING)
        self.__pde.setValue(D=1.)
        self.__c2 = RegionOnEarthSurface.GRAVITY*self.bathymetry/(RegionOnEarthSurface.RADIUS*2*numpy.pi/360.)**2
        c_max = math.sqrt(Lsup(self.__c2))
        self.__dt = self.safety_factor*inf(self.domain.getSize()/(sqrt(self.__c2)+EPS*c_max))
        if self.initial_time_step==None:
            self.initial_time_step=self.__dt
        self.trace("Maximum wave velocity %s m/sec" % c_max)
        self.trace("Time step size is %s sec" % self.__dt)

    def getSafeTimeStepSize(self, dt):
        """
        Returns new step size

        @param dt: last time step size used
        @type dt: C{float}
        @return: time step size that can safely be used
        @rtype: C{float}
        """
        return self.__dt

    def doStepPostprocessing(self,dt):
        """
        Performs the time step using the Verlet scheme

        @param dt: time step size to be used
        @type dt: C{float}
        """
        self.__pde.setValue(X=-self.__c2*grad(self.wave_height))

        new_height = self.wave_height+dt*self.wave_velocity+dt*(self.initial_time_step+dt)/2*self.__pde.getSolution()

        self.wave_velocity = (new_height-self.wave_height)/dt
        self.wave_height = new_height
        self.initial_time_step = dt
        self.trace("Wave height range is %e %e" % (inf(self.wave_height), sup(self.wave_height)))


#=============================================================================
class SurfMovie(Model):
    """
    movie from a wave propagation on the sea

    @ivar time: current time
    @ivar bathymetry: scalar data set
    @ivar wave_height: vector data set
    @ivar filename: name of the movie file
    """
    def __init__(self,**kwargs):
        Model.__init__(self,**kwargs)
        self.declareParameter(bathymetry = 1.,
                              wave_height = 1.,
                              coastline = None,
                              t = 0.,
                              dt = 1.,
                              south = 2.,
                              north = 5.,
                              east = 3.,
                              west = 15.,
                              max_height = 1.,
                              filename = "tsunami.mpg")

    def doInitialization(self):
        """
        Initializes the time integration scheme
        """
        self.__frame_name="tsunami"
        self.__next_t=self.dt
        self.__fileIndex=0

        self.firstFrame = True
        self.imageFiles = []

        factGraphics = vtk.vtkGraphicsFactory()
        factGraphics.SetUseMesaClasses(1)
        factImage = vtk.vtkImagingFactory()
        factImage.SetUseMesaClasses(1)

        print "Preparing bathymetry data..."
        # get depth (z), origin, dx and dy of the bathymetry
        bathZData = self.bathymetry.getData()
        bathOrigin = self.bathymetry.getOrigin()
        bathSpacing = self.bathymetry.getSpacing()

        # now construct the x and y data from the spacing and the origin
        numXPoints = bathZData.shape[1]
        numYPoints = bathZData.shape[0]
        numPoints = numXPoints*numYPoints

        bathXData = numpy.zeros(numXPoints, numpy.float_)
        bathYData = numpy.zeros(numYPoints, numpy.float_)
        for i in range(numXPoints):
            bathXData[i] = bathOrigin[0] + i*bathSpacing[0]

        for i in range(numYPoints):
            bathYData[i] = bathOrigin[1] + i*bathSpacing[1]

        # prepare the bathymetry data set
        bathPoints = vtk.vtkPoints()
        bathPoints.SetNumberOfPoints(numPoints)
        bathData = vtk.vtkFloatArray()
        bathData.SetNumberOfComponents(1)
        bathData.SetNumberOfTuples(numPoints)

        # add the points and data values
        count = 0
        for i in range(numXPoints):
            for j in range(numYPoints):
                bathPoints.InsertPoint(count, bathXData[i], bathYData[j], 0.0)
                bathData.InsertTuple1(count, bathZData[j,i])
                count += 1
        
        # create an unstructured grid for the bathymetry
        bathGrid = vtk.vtkUnstructuredGrid()
        bathGrid.SetPoints(bathPoints)
        bathGrid.GetPointData().SetScalars(bathData)

        # do a delaunay on the grid
        bathDelaunay = vtk.vtkDelaunay2D()
        bathDelaunay.SetInput(bathGrid)
        bathDelaunay.SetTolerance(0.001)
        bathDelaunay.SetAlpha(2.5)

        # save static bathymetry file
        writer = vtk.vtkPolyDataWriter()
        writer.SetInput(bathDelaunay.GetOutput())
        writer.SetScalarsName("bathymetry")
        writer.SetFileName(os.path.join(WORKDIR, "bathymetry.vtk"))
        writer.Write()

        # create the bathymetry colourmap
        data = []
        data.append([-8000, 0,   0,   0])
        data.append([-7000, 0,   5,   25])
        data.append([-6000, 0,   10,  50])
        data.append([-5000, 0,   80,  125])
        data.append([-4000, 0,   150, 200])
        data.append([-3000, 86,  197, 184])
        data.append([-2000, 172, 245, 168])
        data.append([-1000, 211, 250, 211])
        data.append([0,     16,  48,  16])
        data.append([300,   70,  150, 50])
        data.append([500,   146, 126, 60])
        data.append([1000,  198, 178, 80])
        data.append([1250,  250, 230, 100])
        data.append([1500,  250, 234, 126])
        data.append([1750,  252, 238, 152])
        data.append([2000,  252, 243, 177])
        data.append([2250,  253, 249, 216])
        data.append([2500,  255, 255, 255])

        # the amount to scale the data by
        scale = 255.0
        numColours = len(data)

        # convert the colourmap into something vtk is more happy with
        height = numpy.zeros(numColours, numpy.float_)
        red = numpy.zeros(numColours, numpy.float_)
        green = numpy.zeros(numColours, numpy.float_)
        blue = numpy.zeros(numColours, numpy.float_)
        for i in range(numColours):
            (h, r, g, b) = data[i]
            height[i] = float(h)
            red[i] = float(r)/scale
            green[i] = float(g)/scale
            blue[i] = float(b)/scale

        transFunc = vtk.vtkColorTransferFunction()
        transFunc.SetColorSpaceToRGB()
        for i in range(numColours):
            transFunc.AddRGBPoint(height[i], red[i], green[i], blue[i])
            h = height[i]
            while i < numColours-1 and h < height[i+1]:
                h += 1
                transFunc.AddRGBPoint(h, red[i], green[i], blue[i])

        # set up the mapper
        bathMapper = vtk.vtkDataSetMapper()
        bathMapper.SetInput(bathDelaunay.GetOutput())
        bathMapper.SetLookupTable(transFunc)

        # set up the actor
        self.bathActor = vtk.vtkActor()
        self.bathActor.SetMapper(bathMapper)

        print "Done preparing bathymetry data"

        ### prepare and add the coastline ###

        print "Preparing coastline data..."
        # create the coastline grid
        coastGrid = vtk.vtkUnstructuredGrid()
        coastPoints = vtk.vtkPoints()
        totalCoastPoints = 0
        for polyline in self.coastline.polylines:
            numPoints = len(polyline)
            coastLine = vtk.vtkPolyLine()
            coastLine.GetPointIds().SetNumberOfIds(numPoints)
            j = 0
            for point in polyline:
                coastLine.GetPointIds().SetId(j, j+totalCoastPoints)
                coastPoints.InsertNextPoint(point.long, point.lat, 0.0)
                j += 1
            coastGrid.InsertNextCell(coastLine.GetCellType(),
                    coastLine.GetPointIds())
            totalCoastPoints += numPoints

        coastGrid.SetPoints(coastPoints)
        coastMapper = vtk.vtkDataSetMapper()
        coastMapper.SetInput(coastGrid)

        # create an actor for the coastline
        self.coastActor = vtk.vtkActor()
        self.coastActor.SetMapper(coastMapper)
        self.coastActor.GetProperty().SetColor(0,0,0)

        print "Done preparing coastline data"

        # set up the lookup table for the wave data
        refLut = vtk.vtkLookupTable()
        self.lutTrans = vtk.vtkLookupTable()
        refLut.Build()
        alpha = 0.7   # alpha channel value
        for i in range(256):
            (r,g,b,a) = refLut.GetTableValue(255-i)
            if g == 1.0 and b < 0.5 and r < 0.5:
                self.lutTrans.SetTableValue(i, r, g, b, 0.0)
            else:
                self.lutTrans.SetTableValue(i, r, g-0.2, b, alpha)

        # create an actor for the wave
        self.waveActor = vtk.vtkActor()

        # calculate window size
        dataWidth = max(bathXData) - min(bathXData)
        dataHeight = max(bathYData) - min(bathYData)
        self._winHeight = 600
        self._winWidth = int(dataWidth*float(self._winHeight)/dataHeight)

    def saveImage(self, prefix):
        """
        Saves current state as PNM image by saving a VTK file with the
        wave height first and reading it back afterwards.
        """

        # do stuff here that should only be done once
        if self.firstFrame:
            ### set up the renderer and the render window ###
            self.ren = vtk.vtkRenderer()
            self.renWin = vtk.vtkRenderWindow()
            self.renWin.AddRenderer(self.ren)
            self.renWin.SetSize(self._winWidth, self._winHeight)
            # add actors to the renderer
            self.ren.AddActor(self.bathActor)
            self.ren.AddActor(self.coastActor)
            self.ren.AddActor(self.waveActor)
            self.ren.GetActiveCamera().Zoom(2.0)

            # use interactor...
            self.iren = vtk.vtkRenderWindowInteractor()
            self.iren.SetRenderWindow(self.renWin)
            self.iren.Initialize()
            self.iren.Start()

            # ...or make sure rendering is offscreen
            #self.renWin.OffScreenRenderingOn()

            self.firstFrame = False

        vtuFile = prefix+".vtu"
        imgFile = prefix+".pnm"

        # save the VTK file first
        print "Writing", vtuFile
        saveVTK(vtuFile, h=self.wave_height)

        # make a reader for the data
        waveReader = vtk.vtkXMLUnstructuredGridReader()
        waveReader.SetFileName(vtuFile)
        waveReader.Update()

        waveGrid = waveReader.GetOutput()
        waveGrid.Update()

        waveMapper = vtk.vtkDataSetMapper()
        waveMapper.SetInput(waveGrid)
        waveMapper.SetLookupTable(self.lutTrans)
        waveMapper.SetScalarRange(-self.max_height*0.3, self.max_height*0.3)

        self.waveActor.SetMapper(waveMapper)

        # now render the window and save the image file
        print "Rendering to", imgFile
        self.renWin.Render()

        # convert the render window to an image
        win2imgFilter = vtk.vtkWindowToImageFilter()
        win2imgFilter.SetInput(self.renWin)

        # save the image to file
        outWriter = vtk.vtkPNMWriter()
        outWriter.SetInput(win2imgFilter.GetOutput())
        outWriter.SetFileName(imgFile)
        outWriter.Write()
        self.imageFiles.append(imgFile)

    def doStepPostprocessing(self, dt):
        """
        Does any necessary postprocessing after each step

        @param dt:
        """
        if self.t >= self.__next_t:
            prefix = os.path.join(WORKDIR, \
                        "%s%04d" % (self.__frame_name, self.__fileIndex))

            zMin = inf(self.wave_height)
            zMax = sup(self.wave_height)
            print "T = %7.1f: Wave height range=%f...%f" % (self.t, zMin, zMax)

            #self.saveImage(prefix)
            self.wave_height.dump(os.path.join(WORKDIR, \
                    "waveheight%04d.nc" % self.__fileIndex))

            self.__next_t += self.dt
            self.__fileIndex += 1

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

        if len(self.imageFiles) == 0:
            return

        # set up the movie parameters
        paramsFileString = "REFERENCE_FRAME DECODED\n"
        paramsFileString += "FRAME_RATE 24\n"
        paramsFileString += "OUTPUT %s\n" % self.filename
        paramsFileString += "PATTERN IBBPBBPBB\n"
        paramsFileString += "FORCE_ENCODE_LAST_FRAME\n"
        paramsFileString += "GOP_SIZE 20\n"
        paramsFileString += "BSEARCH_ALG CROSS2\n"
        paramsFileString += "PSEARCH_ALG TWOLEVEL\n"
        paramsFileString += "IQSCALE 10\n"
        paramsFileString += "PQSCALE 11\n"
        paramsFileString += "BQSCALE 16\n"
        paramsFileString += "RANGE 8\n"
        paramsFileString += "SLICES_PER_FRAME 1\n"
        paramsFileString += "BASE_FILE_FORMAT PNM\n"
        paramsFileString += "INPUT_DIR .\n"
        paramsFileString += "INPUT_CONVERT *\n"
        paramsFileString += "INPUT\n"
        for fname in self.imageFiles:
            paramsFileString += "%s\n" % fname
        paramsFileString += "END_INPUT\n"
        paramsFileString += "PIXEL HALF\n"
        paramsFileString += "ASPECT_RATIO 1\n"

        # write the parameter file
        fp = open("%s.params" % self.filename, "w")
        fp.write(paramsFileString)
        fp.close()

        print "Performing conversion to mpeg"
        convertString = "ppmtompeg %s.params" % self.filename
        # write the movie into self.filename
        result = os.system(convertString)
        if result != 0:
            print "An error occurred in mpeg conversion"

        # now clean up the image files
        if result == 0:
            print "Removing temporary image files"
            os.unlink("%s.params" % self.filename)
            for fname in self.imageFiles:
                os.unlink(fname)


def main():
    from esys.escript.modelframe import Link,Simulation
    from esys.modellib.input import Sequencer

    oc = OceanRegionCollector()
    oc.range360 = True
    oc.west = 102.9
    oc.east = 232.6
    oc.south = -71.3
    oc.north = 26.7
    oc.resolution = 0.25

    # Uncomment and use the following if GMT is available and the scripts work
    #oc.coastline_source="tsunami_create_coast.py west=%%west%% east=%%east%% south=%%south%% north=%%north%% resolution=%%resolution%% range360=%%range360%% river=false city=false citytyp=wcity.dat volcano=false hotspot=false shoreline=true state=false WSMdata=false plate=false earthquake=false"
    #oc.bathymetry_source = "tsunami_create_bath.py west=%%west%% east=%%east%% south=%%south%% north=%%north%% resolution=%%resolution%%"

    # ...otherwise use pre-created data files (range parameters are ignored)
    oc.coastline_source = "coastline_default.dat"
    oc.bathymetry_source = "bathymetry_default.dat"

    b = Bathymetry()
    b.source = Link(oc,"bathymetry_stream")

    oreg = OceanRegion()
    oreg.source = Link(oc,"coastline_stream")
    oreg.resolution = Link(oc, "resolution")
    oreg.south = Link(oc, "south")
    oreg.north = Link(oc, "north")
    oreg.east = Link(oc, "east")
    oreg.west = Link(oc, "west")
    oreg.bathymetry_data = Link(b, "bathymetry")

    src = TsunamiSource()
    src.domain = Link(oreg, "domain")
    src.decay_zone = 0.01
    src.end_long = 185.
    src.end_lat = -37.
    src.width = 0.5
    src.start_long = 174.
    src.start_lat = -15.
    src.amplitude = 3

    ts = TsunamiInDeepWater()
    ts.domain = Link(oreg, "domain")
    ts.wave_height = Link(src, "wave_height")
    ts.wave_velocity = 0.
    ts.bathymetry = Link(oreg, "bathymetry")

    sq = Sequencer()
    sq.t_end = 30000.

    sm = SurfMovie()
    sm.bathymetry = Link(b, "bathymetry")
    sm.wave_height = Link(ts, "wave_height")
    sm.coastline = Link(oreg, "coastline")
    sm.t = Link(sq,"t")
    sm.filename = "tsunami.mpg"
    sm.north = 8.7
    sm.west = 138.9
    sm.dt = 50.
    sm.east = 196.6
    sm.south = -53.3
    sm.max_height = Link(src, "amplitude")
   
    s = Simulation([sq,oc,b,oreg,src,ts,sm])
    #s.writeXML()
    s.run()

if __name__=="__main__":
    #from esys.modellib import tsunami
    main()

# vim: expandtab shiftwidth=4:
