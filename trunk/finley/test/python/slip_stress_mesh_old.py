
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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""

generates   finley mesh simple vertical fault

THIS CODE CREATES RICH CONTACT ELEMENTS AND RICH FACE ELEMENTS
with fix for contact elements at FAULT ENDS

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Louise Kettle"

from esys.escript import *
from numpy import zeros,float64,array,size

#... generate domain ...
ne = 10
width  = 100000.
height =  30000.
fstart=array([width/2.,7.*width/16.,3.*height/8.])
fend=array([width/2.,9.*width/16.,5.*height/8.])

def faultL(l0,l1, l2,ne0, ne1, ne2,contact=False,xstart=zeros(3),xend=zeros(3)):
   meshfaultL=open('meshfault3D.fly','w')

   FaultError1="ERROR: fault defined on or too close to an outer surface"
   FaultError2="ERROR: the mesh is too coarse for fault"

   N0=ne0+1
   N1=ne1+1
   N2=ne2+1
   if (N0<=N1 and N0<=N2):
     if (N1 <= N2):
        M0=1
        M1=N0
        M2=N0*N1
        M0i=1
        M1i=ne0
        M2i=ne0*ne1
     else:
        M0=1
        M2=N0
        M1=N0*N2
        M0i=1
        M2i=ne0
        M1i=ne0*ne2
   elif (N1<=N2 and N1<=N0):
     if (N2 <= N0): 
        M1=1
        M2=N1
        M0=N2*N1
        M1i=1
        M2i=ne1
        M0i=ne2*ne1
     else:
        M1=1
        M0=N1
        M2=N1*N0
        M1i=1
        M0i=ne1
        M2i=ne0*ne1
   else:
     if (N0 <= N1):
        M2=1
        M0=N2
        M1=N2*N0
        M2i=1
        M0i=ne2
        M1i=ne0*ne2
     else: 
        M2=1
        M1=N2
        M0=N1*N2
        M2i=1
        M1i=ne2
        M0i=ne2*ne1
   
   dim=3
   Element_numNodes=8
   Element_Num=ne0*ne1*ne2
   if contact==False:
      numNodes=N0*N1*N2
   
   else:
      # define double (contact element) nodes on interior of fault 
      i0start=round(xstart[0]*ne0/l0)
      i1start=round(xstart[1]*ne1/l1)
      i2start=round(xstart[2]*ne2/l2)
      i0end=round(xend[0]*ne0/l0)
      i1end=round(xend[1]*ne1/l1)
      i2end=round(xend[2]*ne2/l2)
      n0double=int(i0end)-int(i0start)
      n1double=int(i1end)-int(i1start)
      n2double=int(i2end)-int(i2start)
      if (i0start == 0) or (i1start==0) or (i2start==0):
         raise FaultError1
         
      if (i0end == ne0) or (i1end==ne1) or (i2end==ne2):
         raise FaultError1

      if n0double==0:
         numNodes=N0*N1*N2+(n1double-1)*(n2double-1)

      elif n1double==0:
         numNodes=N0*N1*N2+(n0double-1)*(n2double-1) 
     
      elif n2double==0:
         numNodes=N0*N1*N2+(n0double-1)*(n1double-1)

   # define nodes for normal elements
   # there are N0*N1*N2 normal nodes
   
   Node=zeros([3,numNodes],float64)
   Node_ref=zeros(numNodes,float64)
   Node_DOF=zeros(numNodes,float64)
   Node_tag=zeros(numNodes,float64)

   meshfaultL.write("KettleFault\n")
   #print 'Nodes'
   meshfaultL.write("%dD-nodes %d\n"%(dim,numNodes))

   for i2 in range(N2):
      for i1 in range (N1):
         for i0 in range(N0):
            k=  i0 + N0*i1 + N0*N1*i2 # M0*i0+M1*i1+M2*i2;
            Node_ref[k]= i0 + N0*i1 + N0*N1*i2
            # no periodic boundary conditions
            Node_DOF[k]=Node_ref[k]
            Node_tag[k]=0
            Node[0][k]=(i0)*l0/(N0-1)
            Node[1][k]=(i1)*l1/(N1-1)
            Node[2][k]=(i2)*l2/(N2-1)

   # define double nodes on fault (will have same coordinates as some of nodes already defined)
   # only get double nodes on INTERIOR of fault

   if contact==True:
      Fault_NE=N0*N1*N2   
      if n0double==0:
        if(n1double<=n2double):       
             M1f=1
             M2f=n1double-1
        else:
             M2f=1
             M1f=n2double-1  
    
        for i2 in range(n2double-1):
           for i1 in range(n1double-1):
              # CHANGED:
              k=Fault_NE+i1+(n1double-1)*i2
              Node_ref[k]= k #Fault_NE + i1 + (n1double-1)*i2
              Node_DOF[k]=Node_ref[k]
              Node_tag[k]=1
              Node[0][k]=i0start*l0/ne0
              Node[1][k]=i1start*l1/ne1 + (i1+1)*l1/ne1
              Node[2][k]=i2start*l2/ne2 + (i2+1)*l2/ne2
      # elif n1double==0:
      # elif n2double==0:
      print("fstart = ",[i0start*l0/ne0, i1start*l1/ne1                  , i2start*l2/ne2])
      print("fend = ", [i0start*l0/ne0 , i1start*l1/ne1 + n1double*l1/ne1, i2start*l2/ne2 + n2double*l2/ne2])

   # write nodes to file
   for i in range(numNodes):
       meshfaultL.write("%d %d %d"%(Node_ref[i],Node_DOF[i],Node_tag[i]))
       for j in range(dim):
           meshfaultL.write(" %lf"%Node[j][i])
       meshfaultL.write("\n")




   # defining interior elements
   # there are ne0*ne1*ne2 interior elements

   Element_Nodes=zeros([8,ne0*ne1*ne2],float64)
   Element_ref=zeros(ne0*ne1*ne2,float64)
   Element_tag=zeros(ne0*ne1*ne2,float64)

   #print 'Interior elements'

   print("M0,M1,M2",M0,M1,M2)

   for i2 in range(ne2):
      for i1 in range (ne1):
         for i0 in range(ne0):
            k=i0 + ne0*i1 + ne0*ne1*i2;
            # define corner node (node0)
            node0=i0 + N0*i1 + N0*N1*i2;
            Element_ref[k]=k
            Element_tag[k]=0
            # for hex8 the interior elements are specified by 8 nodes
            #CHANGED:
            Element_Nodes[0][k]=node0;
            Element_Nodes[1][k]=node0+1;
            Element_Nodes[2][k]=node0+N0+1;
            Element_Nodes[3][k]=node0+N0;
            Element_Nodes[4][k]=node0+N0*N1;
            Element_Nodes[5][k]=node0+N0*N1+1;
            Element_Nodes[6][k]=node0+N0*N1+N0+1;
            Element_Nodes[7][k]=node0+N0*N1+N0;

   if contact==True:
      if n0double==0:
         x0s= i0start*l0/ne0
         x1s= i1start*l1/ne1
         x2s= i2start*l2/ne2
         x0e= i0end*l0/ne0
         x1e= i1end*l1/ne1
         x2e= i2end*l2/ne2
         #print "x0s,x1s,x2s,x0e,x1e,x2e",  x0s,x1s,x2s,x0e,x1e,x2e
         if (n1double==1) or (n2double==1):
            raise FaultError2
         for i2 in range(n2double):
            for i1 in range(n1double):
               # here the coordinates of kfault and kold are the same
               # Ref for fault node (only on interior nodes of fault):
               if (i1>0) and (i2>0):
                  kfault=Fault_NE+(i1-1.)+(n1double-1)*(i2-1.)
                  #print kfault , Node[0][int(kfault)],Node[1][int(kfault)],Node[2][int(kfault)]
               else:
                  kfault=0.
               # determine bottom corner node of each element

               # Ref for normal interior node:
               kold=int(i0start+N0*(i1start + i1) + (N0*N1)*(i2start+i2))
               #print kold, Node[0][kold],Node[1][kold],Node[2][kold]
               # Ref for interior element:
               kint=int(i0start + ne0*(i1start+i1) + (ne0*ne1)*(i2start+i2))          
               #print kint, Element_Nodes[0][kint]

               x0= (i0start)*l0/ne0
               x1= (i1start+i1)*l1/ne1
               x2= (i2start+i2)*l2/ne2
               
               # for x0 > xstart we need to overwrite old Nodes in interior element references 
               # with fault nodes:

               # for the interior elements with x1<x1s and x2<x2s the only nodes need changing 
               # are on the fault:
               if (i1==0) and (i2==0): 
                  # nearest fault node:
                  kfaultref=int(Fault_NE+i1+(n1double-1)*i2)
               elif (i1==0): 
                  # nearest fault node
                  kfaultref=int(Fault_NE+i1+(i2-1.)*(n1double-1))
               elif (i2==0): 
                  # nearest fault node
                  kfaultref=int(Fault_NE+(i1-1.) + i2*(n1double-1))
               else: 
                  # looking at element with fault node on bottom corner
                  kfaultref=int(kfault)

               #print x0,x1,x2
               #print kold, Node[0][kold],Node[1][kold],Node[2][kold]
               #print kfaultref, Node[0][kfaultref],Node[1][kfaultref],Node[2][kfaultref]

               # overwrite 4 outer corner elements of fault (only one node changed)
               if (i1==0 and i2==0):           
                  #nodecheck=int(Element_Nodes[7][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]

                  Element_Nodes[7][kint]=kfaultref

                  #print kfaultref, Node[0][kfaultref],Node[1][kfaultref],Node[2][kfaultref]


               elif (i1==0 and i2==n2double-1):

                  #nodecheck=int(Element_Nodes[3][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]

                  Element_Nodes[3][kint]=kfaultref

                  #print kfaultref, Node[0][kfaultref],Node[1][kfaultref],Node[2][kfaultref]

               elif (i1==n1double-1 and i2==0):
                  #nodecheck=int(Element_Nodes[4][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]

                  Element_Nodes[4][kint]=kfaultref
                  #print kfaultref, Node[0][kfaultref],Node[1][kfaultref],Node[2][kfaultref]

               elif (i1==n1double-1 and i2==n2double-1):
                  #nodecheck=int(Element_Nodes[0][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]

                  Element_Nodes[0][kint]=kfaultref
                  #print kfaultref, Node[0][kfaultref],Node[1][kfaultref],Node[2][kfaultref]

               # overwrite 4 sides of fault (only 2 nodes changed)
               elif (i1==0):

                  #nodecheck=int(Element_Nodes[3][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]
                  #nodecheck=int(Element_Nodes[7][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]

                  Element_Nodes[3][kint]=kfaultref
                  kfaultref1=int(kfaultref+(n1double-1))
                  Element_Nodes[7][kint]=kfaultref1

                  #print kfaultref, Node[0][kfaultref],Node[1][kfaultref],Node[2][kfaultref]
                  #print kfaultref1, Node[0][kfaultref1],Node[1][kfaultref1],Node[2][kfaultref1]



               elif (i1==n1double-1):

                  #nodecheck=int(Element_Nodes[0][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]
                  #nodecheck=int(Element_Nodes[4][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]


                  Element_Nodes[0][kint]=kfaultref                     
                  kfaultref1=kfaultref+(n1double-1)
                  Element_Nodes[4][kint]=kfaultref1

                  #print kfaultref, Node[0][kfaultref],Node[1][kfaultref],Node[2][kfaultref]
                  #print kfaultref1, Node[0][kfaultref1],Node[1][kfaultref1],Node[2][kfaultref1]

               elif (i2==0):

                  #nodecheck=int(Element_Nodes[4][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]
                  #nodecheck=int(Element_Nodes[7][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]

                  Element_Nodes[4][kint]=kfaultref
                  kfaultref1=kfaultref+1
                  Element_Nodes[7][kint]=kfaultref1

                  #print kfaultref, Node[0][kfaultref],Node[1][kfaultref],Node[2][kfaultref]
                  #print kfaultref1, Node[0][kfaultref1],Node[1][kfaultref1],Node[2][kfaultref1]

               elif (i2==n2double-1):

                  #nodecheck=int(Element_Nodes[0][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]
                  #nodecheck=int(Element_Nodes[3][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]

                  Element_Nodes[0][kint]=kfaultref
                  kfaultref1=kfaultref+1
                  Element_Nodes[3][kint]=kfaultref1

                  #print kfaultref, Node[0][kfaultref],Node[1][kfaultref],Node[2][kfaultref]
                  #print kfaultref1, Node[0][kfaultref1],Node[1][kfaultref1],Node[2][kfaultref1]

               # overwrite interior fault elements (4 nodes changed)
               else:
                  #nodecheck=int(Element_Nodes[0][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]
                  #print kfaultref, Node[0][kfaultref],Node[1][kfaultref],Node[2][kfaultref]

                  Element_Nodes[0][kint]=kfaultref
                  #if (x1<x1e and x2<x2e):
                  kfaultref1=kfaultref+1
                  kfaultref2=kfaultref+(n1double-1)
                  kfaultref3=kfaultref+1+(n1double-1)

                  #nodecheck=int(Element_Nodes[3][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]
                  #print kfaultref1, Node[0][kfaultref1],Node[1][kfaultref1],Node[2][kfaultref1]

                  #nodecheck=int(Element_Nodes[4][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]
                  #print kfaultref2, Node[0][kfaultref2],Node[1][kfaultref2],Node[2][kfaultref2]

                  #nodecheck=int(Element_Nodes[7][kint] )
                  #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]
                  #print kfaultref3, Node[0][kfaultref3],Node[1][kfaultref3],Node[2][kfaultref3]


                  Element_Nodes[3][kint]=kfaultref1
                  Element_Nodes[4][kint]=kfaultref2
                  Element_Nodes[7][kint]=kfaultref3
                   #elif x1<x1e:
                   #     kfaultref1=kfaultref+M1f
                   #     Element_Nodes[3][kint]=kfaultref1
                   #elif x2<x2e: 
                   #     kfaultref2=kfaultref+M2f
                   #     Element_Nodes[4][kint]=kfaultref2

  
               # print kint, kfaultref

   # write interior elements to file
   Element_Type='Hex8'
   meshfaultL.write("%s %d\n"%(Element_Type,Element_Num))

   for i in range(Element_Num):    
       meshfaultL.write("%d %d"%(Element_ref[i],Element_tag[i]))
       for j in range(Element_numNodes): 
           meshfaultL.write(" %d"%Element_Nodes[j][i])    
       meshfaultL.write("\n")

   # face elements
   FaceElement_Type='Hex8Face'
   FaceElement_Num= 2*(ne0*ne1 + ne0*ne2 + ne1*ne2)
   FaceElement_numNodes=8

   meshfaultL.write("%s %d\n"%(FaceElement_Type,FaceElement_Num))

   FaceElement_Nodes=zeros([FaceElement_numNodes,FaceElement_Num],float64)
   FaceElement_ref=zeros(FaceElement_Num,float64)
   FaceElement_tag=zeros(FaceElement_Num,float64)
   
   kcount=0

   # defining face elements on x2=0 face
   for i1 in range (ne1):
       for i0 in range(ne0):
          i2=0
          k=i0 + ne0*i1 + ne0*ne1*i2;
          # define corner node (node0)
          node0=i0 + N0*i1 + N0*N1*i2;
          FaceElement_ref[kcount]=kcount
          FaceElement_tag[kcount]=3
          # for hex8face the face elements are specified by 8 nodes
          FaceElement_Nodes[0][kcount]=node0;
          FaceElement_Nodes[1][kcount]=node0+1;
          FaceElement_Nodes[2][kcount]=node0+N0+1;
          FaceElement_Nodes[3][kcount]=node0+N0;
          FaceElement_Nodes[4][kcount]=node0+N0*N1;
          FaceElement_Nodes[5][kcount]=node0+N0*N1+1;
          FaceElement_Nodes[6][kcount]=node0+N0*N1+N0+1;
          FaceElement_Nodes[7][kcount]=node0+N0*N1+N0;
          kcount+=1

   # defining face elements on x2=L face
   for i1 in range (ne1):
       for i0 in range(ne0):
          i2=ne2-1
          k=i0 + ne0*i1 + ne0*ne1*i2;
          # define corner node (node0)
          node0=i0 + N0*i1 + N0*N1*i2;
          FaceElement_ref[kcount]=kcount
          FaceElement_tag[kcount]=3
          # for hex8face the face elements are specified by 8 nodes
          FaceElement_Nodes[0][kcount]=node0+N0*N1;
          FaceElement_Nodes[1][kcount]=node0+N0*N1+1;
          FaceElement_Nodes[2][kcount]=node0+N0*N1+N0+1;
          FaceElement_Nodes[3][kcount]=node0+N0*N1+N0;
          FaceElement_Nodes[4][kcount]=node0;
          FaceElement_Nodes[5][kcount]=node0+1;
          FaceElement_Nodes[6][kcount]=node0+N0+1;
          FaceElement_Nodes[7][kcount]=node0+N0;
          kcount+=1

   # defining face elements on x1=0 face
   for i2 in range (ne2):
       for i0 in range(ne0):
          i1=0
          k=i0 + ne0*i1 + ne0*ne1*i2;
          # define corner node (node0)
          node0=i0 + N0*i1 + N0*N1*i2;
          FaceElement_ref[kcount]=kcount
          FaceElement_tag[kcount]=3
          # for hex8face the face elements are specified by 8 nodes
          FaceElement_Nodes[0][kcount]=node0;
          FaceElement_Nodes[1][kcount]=node0+N0*N1;
          FaceElement_Nodes[2][kcount]=node0+N0*N1+1;
          FaceElement_Nodes[3][kcount]=node0+1;
          FaceElement_Nodes[4][kcount]=node0+N0;
          FaceElement_Nodes[5][kcount]=node0+N0*N1+N0;
          FaceElement_Nodes[6][kcount]=node0+N0*N1+N0+1;
          FaceElement_Nodes[7][kcount]=node0+N0+1;
          kcount+=1

   # defining face elements on x1=L face
   for i2 in range (ne2):
       for i0 in range(ne0):
          i1=ne1-1
          k=i0 + ne0*i1 + ne0*ne1*i2;
          # define corner node (node0)
          node0=i0 + N0*i1 + N0*N1*i2;
          FaceElement_ref[kcount]=kcount
          FaceElement_tag[kcount]=3
          # for hex8face the face elements are specified by 8 nodes
          FaceElement_Nodes[0][kcount]=node0+N0;
          FaceElement_Nodes[1][kcount]=node0+N0*N1+N0;
          FaceElement_Nodes[2][kcount]=node0+N0*N1+N0+1;
          FaceElement_Nodes[3][kcount]=node0+N0+1;
          FaceElement_Nodes[4][kcount]=node0;
          FaceElement_Nodes[5][kcount]=node0+N0*N1;
          FaceElement_Nodes[6][kcount]=node0+N0*N1+1;
          FaceElement_Nodes[7][kcount]=node0+1;
          kcount+=1

   # defining face elements on x0=0 face
   for i2 in range (ne2):
       for i1 in range(ne1):
          i0=0
          k=i0 + ne0*i1 + ne0*ne1*i2;
          # define corner node (node0)
          node0=i0 + N0*i1 + N0*N1*i2;
          FaceElement_ref[kcount]=kcount
          FaceElement_tag[kcount]=3
          # for hex8face the face elements are specified by 8 nodes
          FaceElement_Nodes[0][kcount]=node0;
          FaceElement_Nodes[1][kcount]=node0+N0;
          FaceElement_Nodes[2][kcount]=node0+N0*N1+N0;
          FaceElement_Nodes[3][kcount]=node0+N0*N1;
          FaceElement_Nodes[4][kcount]=node0+1;
          FaceElement_Nodes[5][kcount]=node0+N0+1;
          FaceElement_Nodes[6][kcount]=node0+N0*N1+N0+1;
          FaceElement_Nodes[7][kcount]=node0+N0*N1+1;
          kcount+=1

   # defining face elements on x0=L face
   for i2 in range (ne2):
       for i1 in range(ne1):
          i0=ne1-1
          k=i0 + ne0*i1 + ne0*ne1*i2;
          # define corner node (node0)
          node0=i0 + N0*i1 + N0*N1*i2;
          FaceElement_ref[kcount]=kcount
          FaceElement_tag[kcount]=3
          # for hex8face the face elements are specified by 8 nodes
          FaceElement_Nodes[0][kcount]=node0+1;
          FaceElement_Nodes[1][kcount]=node0+N0+1;
          FaceElement_Nodes[2][kcount]=node0+N0*N1+N0+1;
          FaceElement_Nodes[3][kcount]=node0+N0*N1+1;
          FaceElement_Nodes[4][kcount]=node0;
          FaceElement_Nodes[5][kcount]=node0+N0;
          FaceElement_Nodes[6][kcount]=node0+N0*N1+N0;
          FaceElement_Nodes[7][kcount]=node0+N0*N1;
          kcount+=1




   for i in range(FaceElement_Num):    
       meshfaultL.write("%d %d"%(FaceElement_ref[i],FaceElement_tag[i]))
       for j in range(FaceElement_numNodes): 
           meshfaultL.write(" %d"%FaceElement_Nodes[j][i])
       meshfaultL.write("\n")


   # contact elements   
   ContactElement_Type='Hex8Face_Contact'
   ContactElement_Num=0
   ContactElement_numNodes=16
   # print contact elements on fault
   if contact==True:
      if n0double==0:
         ContactElement_Num=(n1double)*(n2double)
         ContactElement_Nodes=zeros([ContactElement_numNodes,ContactElement_Num],float64)
         ContactElement_ref=zeros(ContactElement_Num,float64)
         ContactElement_tag=zeros(ContactElement_Num,float64)
         #print ContactElement_Num

         for i2 in range(n2double):
            for i1 in range(n1double):
               k=i1+(n1double)*i2
               #print k
               # define reference for interior elements with x0<=x0s
               # here the nodes are the old interior nodes
               kintold=int((i0start-1) + ne0*(i1start+i1) + ne0*ne1*(i2start+i2))

               # define reference for interior elements with x0>x0s
               # here the double nodes are the fault nodes

               kintfault=int(i0start + ne0*(i1start+i1) + ne0*ne1*(i2start+i2))

               #nodecheck=int(Element_Nodes[1][kintold] )
               #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]
               #nodecheck=int(Element_Nodes[0][kintfault] )
               #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]

               #nodecheck=int(Element_Nodes[2][kintold] )
               #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]
               #nodecheck=int(Element_Nodes[3][kintfault] )
               #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]

               #nodecheck=int(Element_Nodes[6][kintold] )
               #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]
               #nodecheck=int(Element_Nodes[7][kintfault] )
               #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]

               #nodecheck=int(Element_Nodes[5][kintold] )
               #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]
               #nodecheck=int(Element_Nodes[4][kintfault] )
               #print nodecheck, Node[0][nodecheck],Node[1][nodecheck],Node[2][nodecheck]

               ContactElement_ref[k]=k
               ContactElement_tag[k]=2

               ContactElement_Nodes[0][k]=Element_Nodes[1][kintold]  
               ContactElement_Nodes[1][k]=Element_Nodes[2][kintold]
               ContactElement_Nodes[2][k]=Element_Nodes[6][kintold]
               ContactElement_Nodes[3][k]=Element_Nodes[5][kintold]
               ContactElement_Nodes[4][k]=Element_Nodes[0][kintold]
               ContactElement_Nodes[5][k]=Element_Nodes[3][kintold]
               ContactElement_Nodes[6][k]=Element_Nodes[7][kintold]
               ContactElement_Nodes[7][k]=Element_Nodes[4][kintold]

               ContactElement_Nodes[8][k]=Element_Nodes[0][kintfault]
               ContactElement_Nodes[9][k]=Element_Nodes[3][kintfault]
               ContactElement_Nodes[10][k]=Element_Nodes[7][kintfault]
               ContactElement_Nodes[11][k]=Element_Nodes[4][kintfault]
               ContactElement_Nodes[12][k]=Element_Nodes[1][kintfault]
               ContactElement_Nodes[13][k]=Element_Nodes[2][kintfault]
               ContactElement_Nodes[14][k]=Element_Nodes[6][kintfault]
               ContactElement_Nodes[15][k]=Element_Nodes[5][kintfault]

   meshfaultL.write("%s %d\n"%(ContactElement_Type,ContactElement_Num))

   for i in range(ContactElement_Num):
        meshfaultL.write("%d %d"%(ContactElement_ref[i],ContactElement_tag[i]))
        for j in range(ContactElement_numNodes): 
            meshfaultL.write(" %d"%ContactElement_Nodes[j][i])
        meshfaultL.write("\n")

   # point sources (not supported yet)

   meshfaultL.write("Point1 0")



   meshfaultL.close() 


ne_w=int((ne/height)*width+0.5)
mydomainfile = faultL(width,width, height,ne, ne, ne_w,contact=True,xstart=fstart,xend=fend)
