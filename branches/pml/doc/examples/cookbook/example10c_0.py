##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################
from __future__ import division, print_function

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.pycad import *
from esys.pycad.gmsh import Design
from esys.escript.unitsSI import *
from math import tan

try:
    # This imports the rectangle domain function 
    from esys.finley import MakeDomain
    HAVE_FINLEY = True
except ImportError:
    print("Finley module not available")
    HAVE_FINLEY = False
#the folder to put our outputs in, leave blank "" for script path 
save_path= os.path.join("data","example10")
#
#   input data:
#
width=5*km
length=5*km
depth=2*km
fault_w_front=200*m
fault_w_back=200*m
fault_mid_front=2*km
fault_mid_back=4*km
fault_dip_front=30*DEG
fault_dip_back=30*DEG

layers_left_at_front=[ [ 'limestone' ,  1100*m ] , [  'xx' , 150*m ] ,['shale',300*m],  [ 'limestone' ] ]
layers_right_at_front=[ [ 'limestone' , 500*m ], [ 'xx' , 150*m ],['shale',300*m], ['limestone' ] ]
slop=2*DEG
#
# ====================================================================================
#
def setOffsetsFromThickness(layers):
   l=0
   offsets=[]
   for n in layers:
       if  len(n)<2:
          d=depth-l
          if depth-l <=0:
              raise "Layer structure is too thick. Increase depth > %e"%l
       else:
          d=n[1]
       l+=d
       offsets.append([ n[0], -(l-d), d, -l ])
   return offsets
def setVerticalPositionsFromDip(layers, dip):
   offsets=[]
   s=1/tan(90*DEG-dip)
   for n in layers:
        name=n[0]
        top=n[1]
        d=n[2]
        bot=n[3]
        offsets.append([ name, top, s*top, d, bot, s*bot ] )
   return offsets

def setBackLayers(layers, slop, dip):
      if layers[0][3]-tan(slop)*length <=0:
          raise ValueError("negative thickness in %s at back"%(layers[0][0],))
      out=[[ layers[0][0], layers[0][3]-tan(slop)*length ]]
      for i in range(1,len(layers)-1): out.append([layers[i][0], layers[i][3] ])
      out.append([ layers[-1][0] ])
      out=setOffsetsFromThickness(out)
      return setVerticalPositionsFromDip(out, dip)

def getCutLine(p0,layers,offset=True):
    out=[ ]
    p=p0
    for n in layers:
       if offset:
          p1,p = p, Point(p0.getCoordinates()[0]+n[5],p0.getCoordinates()[1],p0.getCoordinates()[2]+n[4])
       else:
          p1,p = p, Point(p0.getCoordinates()[0],p0.getCoordinates()[1],p0.getCoordinates()[2]+n[4])
       out.append(Line(p1,p))
    return out

def getBackLine(front, back):
    out = [ Line(front[0].getStartPoint(), back[0].getStartPoint()) ]
    for i in range(len(front)):
        out.append(Line(front[i].getEndPoint(), back[i].getEndPoint()))
    return out

def getCrossLine(left, right):
    return getBackLine(left,right)     

def addVolume(front_left, back_left, front_right, back_right, PS, FF, map, filter_left=False):
  
    front_to_back_left = getBackLine(front_left , back_left)
    front_to_back_right = getBackLine(front_right , back_right)
    front_left_to_right= getCrossLine(front_left, front_right)
    back_left_to_right= getCrossLine(back_left, back_right)

     #if filter_left:
      #  out1=front_to_back_left[0]
      #  out2=front_to_back_left[-1]
     #else:
    out1=front_to_back_right[0]
    out2=front_to_back_right[-1]

    topface=PlaneSurface(CurveLoop(front_left_to_right[0], front_to_back_right[0], -back_left_to_right[0], -front_to_back_left[0]))
    for i in range(len(map)):
        name=map[i][0]
        face2=PlaneSurface(CurveLoop(front_left_to_right[i+1], front_to_back_right[i+1], -back_left_to_right[i+1], -front_to_back_left[i+1]))
        face_front=PlaneSurface(CurveLoop(front_left_to_right[i+1], -front_right[i], -front_left_to_right[i],  front_left[i]))
        face_back=PlaneSurface(CurveLoop(back_left_to_right[i+1], -back_right[i], -back_left_to_right[i],  back_left[i]))
        face_left=PlaneSurface(CurveLoop(front_left[i], front_to_back_left[i+1], -back_left[i], -front_to_back_left[i]))
        face_right=PlaneSurface(CurveLoop(front_right[i], front_to_back_right[i+1], -back_right[i], -front_to_back_right[i]))
        
        v=Volume(SurfaceLoop(topface,-face2, face_front, -face_back, -face_left, face_right))
        print(v)
        topface=face2
        if filter_left:
            FF.append(face_right)
        else:
            FF.append(-face_right)
        if name not in PS: PS[name]=PropertySet(name)
        PS[name].addItem(v)
    return out1, out2, PS, FF
    
    
       
if HAVE_FINLEY:
    layers_left_at_edge_at_front=setOffsetsFromThickness(layers_left_at_front)
    layers_right_at_edge_at_front=setOffsetsFromThickness(layers_right_at_front)

    layers_left_at_front=setVerticalPositionsFromDip(layers_left_at_edge_at_front, fault_dip_front)
    layers_right_at_front=setVerticalPositionsFromDip(layers_right_at_edge_at_front, fault_dip_front)


    layers_left_at_back=setBackLayers(layers_left_at_front, slop, fault_dip_back)
    layers_right_at_back=setBackLayers(layers_right_at_front, slop, fault_dip_back)


    left_front_edge=getCutLine(Point(0.,0.,0.), layers_left_at_front, offset=False)
    left_front_fault=getCutLine(Point(fault_mid_front-fault_w_front/2,0.,0.), layers_left_at_front, offset=True)

    right_front_fault=getCutLine(Point(fault_mid_front+fault_w_front/2,0.,0.), layers_right_at_front, offset=True)
    right_front_edge=getCutLine(Point(width,0.,0.), layers_right_at_front, offset=False)

    left_back_edge=getCutLine(Point(0.,length,0.), layers_left_at_back, offset=False)
    left_back_fault=getCutLine(Point(fault_mid_back-fault_w_back/2,length,0.), layers_left_at_back, offset=True)

    right_back_fault=getCutLine(Point(fault_mid_back+fault_w_back/2,length,0.), layers_right_at_back, offset=True)
    right_back_edge=getCutLine(Point(width,length,0.), layers_right_at_back, offset=False)

    PS={}
    FF=[]
    front_to_back_left_top, front_to_back_left_bot, PS, FF=addVolume(left_front_edge, left_back_edge, left_front_fault, left_back_fault, PS, FF, layers_left_at_front, filter_left=False)
    front_to_back_right_top, front_to_back_right_bot, PS, FF=addVolume(right_front_edge, right_back_edge, right_front_fault, right_back_fault, PS, FF, layers_right_at_front, filter_left=True)

    fault_line_top_front=Line(front_to_back_left_top.getStartPoint(), front_to_back_right_top.getStartPoint())
    fault_line_bot_front=Line(front_to_back_left_bot.getStartPoint(), front_to_back_right_bot.getStartPoint())
    fault_line_top_back=Line(front_to_back_left_top.getEndPoint(), front_to_back_right_top.getEndPoint())
    fault_line_bot_back=Line(front_to_back_left_bot.getEndPoint(), front_to_back_right_bot.getEndPoint())

    FF.append(PlaneSurface(CurveLoop(front_to_back_left_top,fault_line_top_back,-front_to_back_right_top,-fault_line_top_front)))
    FF.append(-PlaneSurface(CurveLoop(front_to_back_left_bot,fault_line_bot_back,-front_to_back_right_bot,-fault_line_bot_front)))

    FF.append(PlaneSurface(CurveLoop(*tuple([ -fault_line_top_front,fault_line_bot_front ]+left_front_fault+[ -l for l in right_front_fault ]))))
    FF.append(-PlaneSurface(CurveLoop(*tuple([ -fault_line_top_back,fault_line_bot_back ]+left_back_fault+[ -l for l in right_back_fault ]))))


    # war 120
    des=Design(dim=3, order=1, element_size = 400*m, keep_files=True)
    des.addItems(*tuple(PS.values()))
    des.addItems(PropertySet("fault",Volume(SurfaceLoop( *tuple(FF)))))
    des.setMeshFileName(os.path.join(save_path,"fault.msh"))
    dom=MakeDomain(des)
    dom.write(os.path.join(save_path,"fault.fly"))

