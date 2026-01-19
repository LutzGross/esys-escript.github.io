#*******************************************************
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
#******************************************************/
#
#  this code generates the ReferenceElement_InfoList in ReferenceElements.cpp
#

GEOBASE = {
"Point": (1, 1, "Point", [0] ),
"Line":  (2, 2, "Point", [0,2, 2,1]),
"Tri":(4, 3,   "Line",  [0, 3, 5,  5, 4, 2,  3, 1, 4,  4, 5, 3]),
"Rec" :(4,  4, "Line", [0, 4, 8, 7,  7, 8, 6, 3,   4, 1, 5, 8,  8, 5, 2, 6] ),
"Tet": (8, 4, "Tri",  [0, 4, 6, 7,  4, 1, 5, 8,  6, 5, 2, 9,  7, 8, 9, 3,  4, 5, 6, 8,  5, 9, 8, 6,  9, 7, 6, 8,  7, 4, 6, 8] ),
"Hex": (8, 8, "Rec", [0, 8, 20, 11, 12, 21, 26, 24,  8, 1, 9, 20, 21, 13, 22, 26,  11, 20, 10, 3, 24, 26, 23, 15,  20, 9, 2, 10, 26, 22, 14, 23,  12, 21, 26, 24, 4, 16, 25, 19,  21, 13, 22, 26, 16, 5, 17, 25,  24, 26, 23, 15, 19, 25, 18, 7,  26, 22, 14, 23, 25, 17, 6, 18] )

}

RELEVANTGEO =  {
  "Point1" : [0],
  "Line2" : [0],
  "Line3" : [0],
  "Line4" : [0],
  "Tri3" : [0,1],
  "Tri6" : [0,1,3],
  "Tri9" : [0,1,3,4], 
  "Tri10" : [0,1,3,4],
  "Rec4" : [0,1], 
  "Rec8" : [0,1,4], 
  "Rec9" : [0,1,4],
  "Rec12" : [0,1,4,5], 
  "Rec16" : [0,1,4,5], 
  "Tet4" : [0,1,2], 
  "Tet10" : [0,1,2,4,5,6], 
  "Tet16" : [0,1,2,4,5,6,7,8,9],
  "Hex8" : [0,1,2,3], 
  "Hex20" : [0,1,2,3,8,9,10,11],
  "Hex27" : [0,1,2,3,8,9,10,11,20],
  "Hex32" : [0,1,2,3,8,9,10,11,12,13,14,15]
}

FACENODES= {
  "Line2" : [0], 
  "Line3" : [0], 
  "Line4" : [0], 
  "Tri3" : [0,1],
  "Tri6" : [0,1,3], 
  "Tri9" : [0,1,3,4], 
  "Tri10" : [0,1,3,4], 
  "Rec4" : [0,1], 
  "Rec8" : [0,1,4], 
  "Rec9" : [0,1,4], 
  "Rec12" : [0,1,4,5], 
  "Rec16" : [0,1,4,5], 
  "Tet4" : [0,1,2,3] ,
  "Tet10" : [0,1,2,4,5,6],
  "Tet16" : [0,1,2,4,5,6,7,8,9], 
  "Hex8" : [0,1,2,3], 
  "Hex20" : [0,1,2,3,8,9,10,11], 
  "Hex27" : [0,1,2,3,8,9,10,11,20],
  "Hex32" : [0,1,2,3,8,9,10,11,12,13,14,15]
} 

SHIFTNODES= {
  "Point1" : ( (0,),(-1,)),
  "Line2" : ( (1,0),(-1,)),
  "Line3" : ( (1,0,2),(-1,)),
  "Line4" : ( (1,0,3,2),(-1,)),
  "Tri3" : ( (1,2,0),(0,2,1)),
  "Tri6" : ( (1,2,0,4,5,3),(0,2,1,5,4,3)),
  "Tri9" : ( (1,2,0,5,6,7,8,3,4),(0,2,1,8,7,6,5,4,3)),
  "Tri10" : ( (1,2,0,5,6,7,8,3,4,9),(0,2,1,8,7,6,5,4,3,9)),
  "Rec4" : ( (1,2,3,0),(0,3,2,1)),
  "Rec8" : ( (1,2,3,0,5,6,7,4),(0,3,2,1,7,6,5,4)),
  "Rec9" : ( (1,2,3,0,5,6,7,4,8),(0,3,2,1,7,6,5,4,8)),
  "Rec12" : ( (1,2,3,0,6,7,8,9,10,11,4,5),(0,3,2,1,11,10,9,8,7,6,5,4)) ,
  "Rec16" : ( (1,2,3,0,6,7,8,9,10,11,4,5,13,14,15,12),(0,3,2,1,11,10,9,8,7,6,5,4,12,15,14,13)),
  "Tet4" : ( (-1,),(-1,)),
  "Tet10" : ( (-1,),(-1,)),
  "Tet16" : ( (-1,),(-1,)),
  "Hex8" : ( (-1,),(-1,)),
  "Hex20" : ( (-1,),(-1,)),
  "Hex27" : ( (-1,),(-1,)),
  "Hex32" : ( (-1,),(-1,)),
  "Line2Face" : ( (0,1,2),(-1,)),
  "Line3Face" : ( (0,1,2),(-1,)),
  "Line4Face" : ( (0,1,2),(-1,)),
  "Tri3Face" : ( (1,0,2), (-1,)),
  "Tri6Face" : ( (1,0,2,3,5,4),(-1,)),
  "Tri9Face" : ( (1,0,2,4,3,7,8,6,5),(-1,)),
  "Tri10Face" : ( (1,0,2,4,3,7,8,6,5,9),(-1,)),
  "Rec4Face" : ( (1,0,3,2),(-1,)),
  "Rec8Face" : ( (1,0,3,2,4,7,6,5),(-1,)),
  "Rec9Face" : ( (1,0,3,2,4,7,6,5,8),(-1,)),
  "Rec12Face" : ( (1,0,3,2,5,4,11,10,9,8,7,6),(-1,)),
  "Rec16Face" : ( (1,0,3,2,5,4,11,10,9,8,7,6,13,12,15,14),(-1,)),
  "Tet4Face" : ( (1,2,0,3),(0,2,1,3)) ,
  "Tet10Face" : ( (1,2,0,3,5,6,4,8,9,7),(0,2,1,3,6,7,9,8)),
  "Tet16Face" : ( (1,2,0,3,6,7,8,9,4,5,11,12,10,14,15,13),(0,2,1,3,9,8,7,6,5,4,9,8,7,6,10,12,11,13,15,14)),
  "Hex8Face" : ( (1,2,3,0,5,6,7,4),(0,3,2,1,4,7,6,5)),
  "Hex20Face" : ( (1,2,3,0,5,6,7,4,9,10,11,8,13,14,15,12,17,18,19,16),(0,3,2,1,4,7,6,5,11,10,9,8,12,15,14,13,19,18,17,16)),
  "Hex27Face" : ( (1,2,3,0,5,6,7,4,9,10,11,8,13,14,15,12,17,18,19,16,20,22,23,24,22,25,26),(0,3,2,1,4,7,6,5,11,10,9,8,12,15,14,13,19,18,17,16,20,24,23,22,21,25,26)),
  "Hex32Face" : ( (1,2,3,0,5,6,7,4,10,11,12,13,14,15,8,9,17,18,19,16,21,22,23,20,26,27,28,29,30,31,34,25), (0,3,2,1,4,7,6,5,15,14,13,12,11,10,9,8,16,19,18,17,20,23,22,21,31,30,29,28,27,26,25,24))
}
def listToArrayStr(l):
   out="{"
   for s in l:
        if len(out)==1:
          out+=" %s"%s
        else:
          out+=", %s"%s
   return out+" }"
def LENLEN(l):
   if len(l) == 1:
      if l[0] == -1:
         return -1
      else:
         return 1
   else:
      return len(l)

outall="ReferenceElementInfo ReferenceElement_InfoList[]={\n"

for name in ["Point1", "Line2", "Line3", "Line4", "Tri3", "Tri6", "Tri9", "Tri10", "Rec4", "Rec8", "Rec9", "Rec12", "Rec16", "Tet4", "Tet10", "Tet16", "Hex8", "Hex20", "Hex27", "Hex32", "Line2Face", "Line3Face", "Line4Face", "Tri3Face", "Tri6Face", "Tri9Face", "Tri10Face", "Rec4Face", "Rec8Face", "Rec9Face", "Rec12Face", "Rec16Face", "Tet4Face", "Tet10Face", "Tet16Face", "Hex8Face", "Hex20Face", "Hex27Face", "Hex32Face", "Point1_Contact", "Line2_Contact", "Line3_Contact", "Line4_Contact", "Tri3_Contact", "Tri6_Contact", "Tri9_Contact", "Tri10_Contact", "Rec4_Contact", "Rec8_Contact", "Rec9_Contact", "Rec12_Contact", "Rec16_Contact", "Line2Face_Contact", "Line3Face_Contact", "Line4Face_Contact", "Tri3Face_Contact", "Tri6Face_Contact", "Tri9Face_Contact", "Tri10Face_Contact", "Rec4Face_Contact", "Rec8Face_Contact", "Rec9Face_Contact", "Rec12Face_Contact", "Rec16Face_Contact", "Tet4Face_Contact", "Tet10Face_Contact", "Tet16Face_Contact", "Hex8Face_Contact", "Hex20Face_Contact", "Hex27Face_Contact", "Hex32Face_Contact", "Line3Macro", "Tri6Macro", "Rec9Macro", "Tet10Macro", "Hex27Macro" ]:
        isFace=False
        isMacro=False
        isContact=False
        n=1
        z="Point"
        for z in [ "Point", "Line", "Tri", "Rec", "Tet", "Hex" ] :
            if name.startswith(z):
               s=name[len(z):]
               n=s
               if s.find("Macro")>=0:
                  isMacro=True
                  if s.find("Macro")>0: n=s[:s.find("Macro")]
                  s=s[s.find("Macro")+5:]
               if s.find("Face")>=0:
                  isFace=True
                  if s.find("Face")>0: n=s[:s.find("Face")]
                  s=s[s.find("Face")+4:]
               if s.find("_Contact")>=0:
                  isContact=True
                  if s.find("_Contact")>0: n=s[:s.find("_Contact")]
                  s=s[s.find("_Contact")+8:]
               n=int(n)
               break
        numLinearNodes=GEOBASE[z][1]
        if isFace:
           Quadrature=GEOBASE[z][2]
        else:
           Quadrature=z
        Parametrization="%s%s"%(z,n)
        if isContact:
            offsets=[0, n, 2*n]
            numSides=2
        else:
            numSides=1
            offsets=[0, n]
        numNodes=offsets[-1]
        if isMacro:
           numSubElements=GEOBASE[z][0]
           BasisFunctions="%s%s"%(z,numLinearNodes)
           subElementNodes=GEOBASE[z][3]
        else:
            numSubElements=1
            BasisFunctions=Parametrization
            subElementNodes=[ i for i in range(numNodes) ]

        linearTypeId="%s%s"%(z,numLinearNodes)
        linearNodes=[ i for i in range(numLinearNodes) ]
        if isFace: linearTypeId+="Face"
       
        if isContact: 
             linearTypeId+="_Contact"
             linearNodes+=[ n+i for i in range(numLinearNodes) ]
         
        if isFace:
           relevantGeoNodes=RELEVANTGEO["%s%s"%(z,n)]
        else:
           relevantGeoNodes=[ i for i in range(n) ]

        if isContact:
             faceNodes =[-1]
             shiftNodes = [-1]
             reverseNodes = [-1]
        else:
           if isFace:
                faceNodes=FACENODES["%s%s"%(z,n)]
                shiftNodes = SHIFTNODES["%s%sFace"%(z,n)][0]
                reverseNodes = SHIFTNODES["%s%sFace"%(z,n)][1]
           else:
                faceNodes=[  i for i in range(numNodes) ]
                shiftNodes = SHIFTNODES["%s%s"%(z,n)][0]
                reverseNodes = SHIFTNODES["%s%s"%(z,n)][1]


        out = '{ %s, "%s", %d, %d, %d, %s, %s,\n    %s, %sQuad, %sShape, %sShape,\n    %s,\n  %s, %s,\n  %s, %s,\n    %s,\n    %s },\n'% \
               (name, name, numNodes, numSubElements, numSides, listToArrayStr(offsets), linearTypeId, listToArrayStr(linearNodes),
                  Quadrature, Parametrization, BasisFunctions, listToArrayStr(subElementNodes),
                  len(relevantGeoNodes), listToArrayStr(relevantGeoNodes),
                  LENLEN(faceNodes), listToArrayStr(faceNodes), listToArrayStr(shiftNodes), listToArrayStr(reverseNodes) )
        outall+=out
out = '{ %s, "%s", %d, %d, %d, %s, %s,\n    %s, %s, %s, %s,\n    %s,\n  %s, %s,\n  %s, %s,\n    %s,\n    %s }\n'% \
        ("NoRef", "noElement", 0, 0, 0, listToArrayStr([0]), "NoRef", listToArrayStr([0]),
         "NoQuad", "NoShape", "NoShape", listToArrayStr([0]),
         -1, listToArrayStr([0]), -1, listToArrayStr([0]), listToArrayStr([0]), listToArrayStr([0]) )
outall+=out
print outall+"\n};"
