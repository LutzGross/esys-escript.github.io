"""
mining activities
                                                                                                                                                                                                     
@var __author__: name of author
@var __licence__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Lutz Gross, l.gross@uq.edu.au"
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision$"
__date__="$Date$"

from xml.dom import minidom
from esys.pycad import *
from math import pi, cos, sin, sqrt
from esys.pycad.gmsh import Design
#======== that should be in a different file =============================
class LocationOnEarth(object):
    AGD84="AGD84"
    AGD84_semi_major_axis=6378160.00
    AGD84_inverse_flatening=298.2500
    DEG_UNIT="D.MM"
    LENGTH_UNIT="M"
    def __init__(self, longitude=0., latitude=0., altitude=0.):
         """
         """
         self.__longitude=longitude
         self.__latitude=latitude
         self.__altitude=altitude
    def getLLA(self):
        return (self.__longitude,self.__latitude,self.__altitude)

    def getRelativeXYZ(self,center):
       """
       return the location in global coordinates relative to the center of the earth/
       default reference grid is AGD84.

       @note: formula for AGD84 from http://www.ga.gov.au/geodesy/datums/transxyz.xls
       """
       C_LONG,C_LAT,C_ALT=center.getLLA()
       C_LONG=unitConverter(C_LONG, center.DEG_UNIT,"RAD")
       C_LAT=unitConverter(C_LAT, center.DEG_UNIT,"RAD")
       I7=unitConverter(self.__latitude,self.DEG_UNIT,"RAD")
       I8=unitConverter(self.__longitude,self.DEG_UNIT,"RAD")
       C_ALT=unitConverter(C_ALT, center.LENGTH_UNIT,self.LENGTH_UNIT)
       X=self.AGD84_semi_major_axis*sin(I7-C_LAT)
       Y=self.AGD84_semi_major_axis*sin(I8-C_LONG)
       Z=(self.__altitude-C_ALT)
       return numarray.array([X,Y,Z])
  

    def getXYZ(self,center=None,reference_grid=None):
       """
       return the location in global coordinates relative to the center of the earth/
       default reference grid is AGD84.

       @note: formula for AGD84 from http://www.ga.gov.au/geodesy/datums/transxyz.xls
       """
       if reference_grid==None:
          I7=unitConverter(self.__latitude,self.DEG_UNIT,"RAD")
          I8=unitConverter(self.__longitude,self.DEG_UNIT,"RAD")
          D15=1./self.AGD84_inverse_flatening
          D16=2*D15-D15*D15
          D17=self.AGD84_semi_major_axis/sqrt(1-D16*sin(I7)*sin(I7))
          X=(D17        +self.__altitude)*cos(I7)*cos(I8)
          Y=(D17        +self.__altitude)*cos(I7)*sin(I8)
          Z=((1-D16)*D17+self.__altitude)*sin(I7)
       else:
          raise ValueError("reference  grid %s is not supported."%reference_grid)
       return numarray.array([X,Y,Z])

def unitConverter(val,val_unit,out_unit):
    VAL_UNIT=val_unit.upper()
    OUT_UNIT=out_unit.upper()
    out=None
    if OUT_UNIT == "D.MM":
       if VAL_UNIT == "D.MM":
          out=float(val)
    elif OUT_UNIT == "RAD":
       if VAL_UNIT == "D.MM":
         out=float(val)/180.*pi
    elif OUT_UNIT == "M":
       if VAL_UNIT == "M":
          out=float(val)
       elif VAL_UNIT == "KM":
          out=float(val)*1000.
    elif OUT_UNIT == "KG":
        if VAL_UNIT == "MT":
            out=float(val)*1e9
    if out == None:
          raise ValueError("cannot convert from unit %s to %s."%(val_unit,out_unit))
    return out
#======== end of snip ====================================================
class BoxInTheCrust(object):
     """
     defines a box in the outer crust
     """
     def __init__(self, longitude_length, latitude_length, depth, location=LocationOnEarth(), name="unknown", tilt=0.):
        self.__longitude_length=longitude_length
        self.__latitude_length=latitude_length
        self.__name=name
        self.__location=location
        self.__tilt=tilt
        self.__depth=depth
        self.__centerxyz=self.__location.getXYZ()
     def getDiameter(self):
         return sqrt(self.__depth**2+self.__longitude_length**2+self.__latitude_length**2)
     def getCenterXYZ(self):
         return self.__centerxyz
     def getCenter(self):
         return self.__location
     def getName(self):
         return self.__name
     def getLongitudeLength(self):
            return self.__longitude_length
     def getLatitudeLength(self):
            return self.__latitude_length
     def getDepth(self):
            return self.__depth
     def getTotalDepth(self):
            return self.__depth-self.__location.getLLA()[2]
         
     def createSurfaceLoop(self):
        """
        this creates a L{SurfaceLoop} describing the region
        """
        p1=Point(-self.__longitude_length/2,-self.__latitude_length/2,-self.__depth)
        p2=Point( self.__longitude_length/2, self.__latitude_length/2,0.)
        out=Brick(p1,p2)
        out*=Rotatation(axis=[0,0,1], angle=self.__tilt)
        return out
       

class MiningArea(BoxInTheCrust):
     """
     defines minuing activites in an area
     """
     def __init__(self, longitude_length, latitude_length, depth, location=LocationOnEarth(), name="unknown", tilt=0., mines=[]):
        super(MiningArea, self).__init__(longitude_length, latitude_length, depth, location, name, tilt)
        self.__mines=mines
     def getStartOfRecords(self):
         return min([ m.getStartOfRecords() for m in self.__mines ])
     def getMassChanges(self,t):
         out={}
         for m in self.__mines: out[m.getName()]=m.getMassChanges(t)
         return out
     def getNextTimeMarker(self,t):
         return min([ m.getNextTimeMarker(t) for m in self.__mines ])
          
     def fillDesign(self,design):
         """
         this puts the mining area and all the mines into the given design
         """
         # define a bounding box around the mines:
         X_MIN=0
         X_MAX=0
         Y_MIN=0
         Y_MAX=0
         DEPTH_MAX=0.
         for m in self.__mines:
               long=m.getLongitudeLength()
               lat=m.getLatitudeLength()
               pos=m.getCenter().getRelativeXYZ(self.getCenter())
               X_MIN=min(X_MIN,pos[0]-lat/2.)
               X_MAX=max(X_MAX,pos[0]+lat/2.)
               Y_MIN=min(Y_MIN,pos[1]-long/2.)
               Y_MAX=max(Y_MAX,pos[1]+long/2.)
               DEPTH_MAX=max(DEPTH_MAX,m.getDepth()-pos[2])
         lx=(X_MAX-X_MIN)/2*1.6
         cx=(X_MIN+X_MAX)/2
         X_MAX=cx+lx
         X_MIN=cx-lx
         ly=(Y_MAX-Y_MIN)/2*1.6
         cy=(Y_MIN+Y_MAX)/2
         Y_MAX=cy+ly
         Y_MIN=cy-ly
         DEPTH_MAX*=2.
         dx=[X_MAX-X_MIN, Y_MAX-Y_MIN, DEPTH_MAX]
         bb_p000=Point(X_MIN,Y_MIN,-DEPTH_MAX)
         bb_p100=bb_p000+[dx[0],0.,0.]
         bb_p010=bb_p000+[0.,dx[1],0.]
         bb_p110=bb_p000+[dx[0],dx[1],0.]
         bb_p001=bb_p000+[0.,0.,dx[2]]
         bb_p101=bb_p000+[dx[0],0.,dx[2]]
         bb_p011=bb_p000+[0.,dx[1],dx[2]]
         bb_p111=bb_p000+[dx[0],dx[1],dx[2]]
         bb_l10=Line(bb_p000,bb_p100)
         bb_l20=Line(bb_p100,bb_p110)
         bb_l30=Line(bb_p110,bb_p010)
         bb_l40=Line(bb_p010,bb_p000)
         bb_l11=Line(bb_p000,bb_p001)
         bb_l21=Line(bb_p100,bb_p101)
         bb_l31=Line(bb_p110,bb_p111)
         bb_l41=Line(bb_p010,bb_p011)
         bb_l12=Line(bb_p001,bb_p101)
         bb_l22=Line(bb_p101,bb_p111)
         bb_l32=Line(bb_p111,bb_p011)
         bb_l42=Line(bb_p011,bb_p001)
         bb_bottom=PlaneSurface(CurveLoop(-bb_l10,-bb_l40,-bb_l30,-bb_l20))
         bb_top=PlaneSurface(CurveLoop(bb_l12,bb_l22,bb_l32,bb_l42))
         bb_front=PlaneSurface(CurveLoop(-bb_l11,bb_l10,bb_l21,-bb_l12))
         bb_back=PlaneSurface(CurveLoop(bb_l30,bb_l41,-bb_l32,-bb_l31))
         bb_left=PlaneSurface(CurveLoop(bb_l11,-bb_l42,-bb_l41,bb_l40))
         bb_right=PlaneSurface(CurveLoop(-bb_l21,bb_l20,bb_l31,-bb_l22))
         bb=SurfaceLoop(bb_bottom,bb_top,bb_front,bb_back,bb_left,bb_right)
         bb.setLocalScale(min(X_MAX-X_MIN, Y_MAX-Y_MIN, DEPTH_MAX)/design.getElementSize())
         # put all mines in to the domain:
         mboxes=[]
         props=[]
         for m in self.__mines:
               mb=m.createSurfaceLoop()
               mb.setLocalScale(m.getDiameter()/design.getElementSize())
               mb+=m.getCenter().getRelativeXYZ(self.getCenter())
               mboxes.append(-mb)
               props.append(PropertySet(m.getName(),Volume(mb)))
         design.addItems(*tuple(props))
         long=self.getLongitudeLength()
         lat=self.getLatitudeLength()
         dep=self.getDepth()
         p000=Point(-long/2,-lat/2,-dep)
         p100=p000+[long,0.,0.]
         p010=p000+[0.,lat,0.]
         p110=p000+[long,lat,0.]
         p001=p000+[0.,0.,dep]
         p101=p000+[long,0.,dep]
         p011=p000+[0.,lat,dep]
         p111=p000+[long,lat,dep]
         l10=Line(p000,p100)
         l20=Line(p100,p110)
         l30=Line(p110,p010)
         l40=Line(p010,p000)
         l11=Line(p000,p001)
         l21=Line(p100,p101)
         l31=Line(p110,p111)
         l41=Line(p010,p011)
         l12=Line(p001,p101)
         l22=Line(p101,p111)
         l32=Line(p111,p011)
         l42=Line(p011,p001)
         bottom=PlaneSurface(CurveLoop(-l10,-l40,-l30,-l20))
         top=PlaneSurface(CurveLoop(l12,l22,l32,l42),holes=[CurveLoop(bb_l12,bb_l22,bb_l32,bb_l42)])
         front=PlaneSurface(CurveLoop(-l11,l10,l21,-l12))
         back=PlaneSurface(CurveLoop(l30,l41,-l32,-l31))
         left=PlaneSurface(CurveLoop(l11,-l42,-l41,l40))
         right=PlaneSurface(CurveLoop(-l21,l20,l31,-l22))
         dom=SurfaceLoop(bottom,top,front,back,left,right,-bb_bottom,-bb_front,-bb_back,-bb_left,-bb_right)
         design.addItems(PropertySet("matrix",Volume(bb, holes= mboxes), Volume(dom)))
         return

class Mine(BoxInTheCrust):
     """
     defines a mine
     """
     def __init__(self, longitude_length, latitude_length, depth, location=LocationOnEarth(), name="unknown", tilt=0.):
        super(Mine, self).__init__(longitude_length, latitude_length, depth, location, name, tilt)
        self.__record=[]

     def addRecord(self, material, year, extraction):
        TYPE=material.upper()
        for y, e in self.__record:
           if y == year:
               if e.has_key(TYPE):
                  e[TYPE]+=extraction
               else:
                  e[TYPE]=extraction
               return
           if y>year:
               self.__record.insert(self.__record.index((y,e)), { TYPE : extraction} )
               return
        self.__record.append((year, { TYPE : extraction} ))
        return
     def getStartOfRecords(self):
          if len(self.__record)>0:
             return self.__record[0][0]
          else: 
             raise ValueError("empty record of %s mine."%self.getName())
     def getMassChanges(self,t):
         m0=None
         t0=self.getStartOfRecords()
         if t<=t0:
            return 0.
         for y, e in self.__record:
             m=sum(e.values())
             if t<=y:
                  if m0==None:
                     return m
                  else:
                     return (m-m0)/(y-t0)*(t-t0)+m0
             else:
               t0,m0=y,m
         return m0
     def getNextTimeMarker(self,t):
         for y, e in self.__record:
           if y>t: return y
         else:
           return 999999999


def _parse(root):
     """
     parses child roots on root
     """
     if isinstance(root, minidom.Element):
        if root.tagName == 'MiningArea':
                  NAME="unknown"
                  LOC=LocationOnEarth()
                  TILT=0.
                  LONGLENGTH=0.
                  LATLENGTH=0.
                  DEPTH=0.
                  MINES=[]
                  for node in root.childNodes:
                       if isinstance(node, minidom.Element):
                          if node.tagName == 'Tilt': TILT=_parse(node)
                          if node.tagName == 'Location': LOC=_parse(node)
                          if node.tagName == 'Name': NAME=_parse(node)
                          if node.tagName == 'Mine': MINES.append(_parse(node))
                          if node.tagName == 'Depth': DEPTH=_parse(node)
                          if node.tagName == 'LongitudeLength': LONGLENGTH=_parse(node)
                          if node.tagName == 'LatitudeLength': LATLENGTH=_parse(node)
                  return MiningArea(location=LOC, name=NAME, longitude_length=LONGLENGTH, latitude_length=LATLENGTH, tilt=TILT, mines=MINES, depth=DEPTH)
        elif root.tagName == 'Mine':
                 NAME="unknown"
                 LOC=LocationOnEarth()
                 TILT=0.
                 DEPTH=0.
                 LONGLENGTH=0.
                 LATLENGTH=0.
                 RECORD=[]
                 for node in root.childNodes:
                       if isinstance(node, minidom.Element):
                          if node.tagName == 'Name': NAME=_parse(node)
                          if node.tagName == 'Location': LOC=_parse(node)
                          if node.tagName == 'Tilt': TILT=_parse(node)
                          if node.tagName == 'Depth': DEPTH=_parse(node)
                          if node.tagName == 'LongitudeLength': LONGLENGTH=_parse(node)
                          if node.tagName == 'LatitudeLength': LATLENGTH=_parse(node)
                          if node.tagName == 'Record': RECORD.append(_parse(node))
                 m=Mine(location=LOC, tilt=TILT,name=NAME, longitude_length=LONGLENGTH, latitude_length=LATLENGTH, depth=DEPTH)
                 for r in RECORD:
                    m.addRecord(material=r[0],year=r[1],extraction=r[2])
                 return m
        elif root.tagName == 'LocationOnEarth':
                    long=0.
                    lat=0.
                    alt=0.
                    for node in  root.childNodes:
                       if isinstance(node, minidom.Element):
                          if node.tagName == 'Longitude': long=_parse(node)
                          if node.tagName == 'Latitude': lat=_parse(node)
                          if node.tagName == 'Altitude': alt=_parse(node)
                    return LocationOnEarth(longitude=long, latitude=lat, altitude=alt)
        elif root.tagName == 'Record':
                    year=0
                    mat="unknown"
                    ext=0.
                    for node in  root.childNodes:
                       if isinstance(node, minidom.Element):
                          if node.tagName == 'Year': year=_parse(node)
                          if node.tagName == 'Material': mat=_parse(node)
                          if node.tagName == 'Extraction': ext=_parse(node)
                    return (mat,year,ext)
        elif root.tagName == 'SouthWest':
                  return _parse(root.getElementsByTagName('LocationOnEarth')[0])
        elif root.tagName == 'NorthEast':
                  return _parse(root.getElementsByTagName('LocationOnEarth')[0])
        elif root.tagName == 'Longitude':
                   if root.hasAttribute("unit"):
                       unit=root.getAttribute("unit")
                   else:
                       unit="D.MM"
                   for node in root.childNodes:
                      if isinstance(node, minidom.Text):
                        return unitConverter(node.nodeValue.strip(), unit, LocationOnEarth.DEG_UNIT)
                   return 0.
        elif root.tagName== 'Latitude':
                   if root.hasAttribute("unit"):
                       unit=root.getAttribute("unit")
                   else:
                       unit="D.MM"
                   for node in root.childNodes:
                      if isinstance(node, minidom.Text):
                        return unitConverter(node.nodeValue.strip(), unit, LocationOnEarth.DEG_UNIT)
                   return 0.
        elif root.tagName == 'Altitude':
                   if root.hasAttribute("unit"):
                       unit=root.getAttribute("unit")
                   else:
                       unit="M"
                   for node in root.childNodes:
                      if isinstance(node, minidom.Text):
                        return unitConverter(node.nodeValue.strip(), unit, LocationOnEarth.LENGTH_UNIT)
                   return 0.
        elif root.tagName == 'LongitudeLength':
                   if root.hasAttribute("unit"):
                       unit=root.getAttribute("unit")
                   else:
                       unit="M"
                   for node in root.childNodes:
                      if isinstance(node, minidom.Text):
                        return unitConverter(node.nodeValue.strip(), unit, LocationOnEarth.LENGTH_UNIT)
                   return 0.
        elif root.tagName == 'LatitudeLength':
                   if root.hasAttribute("unit"):
                       unit=root.getAttribute("unit")
                   else:
                       unit="M"
                   for node in root.childNodes:
                      if isinstance(node, minidom.Text):
                        return unitConverter(node.nodeValue.strip(), unit, LocationOnEarth.LENGTH_UNIT)
                   return 0.
        elif root.tagName == 'Depth':
                   if root.hasAttribute("unit"):
                       unit=root.getAttribute("unit")
                   else:
                       unit="M"
                   for node in root.childNodes:
                      if isinstance(node, minidom.Text):
                        return unitConverter(node.nodeValue.strip(), unit, LocationOnEarth.LENGTH_UNIT)
                   return 0.
        elif root.tagName == 'LongitudeLength':
                   if root.hasAttribute("unit"):
                       unit=root.getAttribute("unit")
                   else:
                       unit="M"
                   for node in root.childNodes:
                      if isinstance(node, minidom.Text):
                        return unitConverter(node.nodeValue.strip(), unit, LocationOnEarth.LENGTH_UNIT)
                   return 0.
        elif root.tagName == 'Tilt':
                   if root.hasAttribute("unit"):
                       unit=root.getAttribute("unit")
                   else:
                       unit="D.MM"
                   for node in root.childNodes:
                      if isinstance(node, minidom.Text):
                        return unitConverter(node.nodeValue.strip(), unit, LocationOnEarth.DEG_UNIT)
                   return 0.
        elif root.tagName == 'Name':
                   for node in root.childNodes:
                      if isinstance(node, minidom.Text):
                        return node.nodeValue.strip()
                   return "unknown"
        elif root.tagName == 'Year':
                  for node in root.childNodes:
                      if isinstance(node, minidom.Text):
                          return int(node.nodeValue.strip())
                  return 0
        elif root.tagName == 'Material':
                   for node in root.childNodes:
                      if isinstance(node, minidom.Text):
                        return node.nodeValue.strip()
                   return "unknown"
        elif root.tagName == 'Extraction':
                   if root.hasAttribute("unit"):
                       unit=root.getAttribute("unit")
                   else:
                       unit="MT"
                   for node in root.childNodes:
                      if isinstance(node, minidom.Text):
                        return unitConverter(node.nodeValue.strip(), unit, "kg")
                   return 0.
     for node in root.childNodes:
           if isinstance(node, minidom.Element): return _parse(node)
   
    
def parse(xml):
     if isinstance(xml,str):
          dom=minidom.parseString(xml)
     else:
          dom=minidom.parse(xml)
     root=dom.getElementsByTagName('ESys')[0]
     return _parse(root)

     
if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage = "usage: %prog [options] filename")
    parser.add_option("-g", "--geo", dest="geofile", type="string", action = "store", default=None,
                        help="geometry file (output)") 
    parser.add_option("-m", "--msh", dest="mshfile", type="string", action = "store", default=None,
                        help="mesh file (output)")
    parser.add_option("-t", "--tag", dest="tagfile", type="string", action = "store", default=None,
                        help="tags file (output)")
    parser.add_option("-s", "--size", dest="rel_size", type="float", action = "store", default=0.2,
                        help="relative mesh size")
    (options, args) = parser.parse_args()
    if not len(args)==1:
        raise parser.error("input file missing.")
    FILE=args[0]
    mine=parse(open(FILE,'r'))
    dsgn=Design(element_size=mine.getDiameter()*options.rel_size)
    if not options.geofile == None:
        dsgn.setScriptFileName(options.geofile)
    if not options.mshfile == None:
        dsgn.setMeshFileName(options.mshfile)
    mine.fillDesign(dsgn)
    print dsgn.getCommandString()
    print "mesh in gmsh format is written to ",dsgn.getMeshHandler()
    if not options.tagfile == None:
        dsgn.getTagMap().writeXML(open(options.tagfile,"w"))
        print "tag map written to %s."%options.tagfile
