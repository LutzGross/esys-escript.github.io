
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Finley: Reference elements */

/**************************************************************/

#include "ReferenceElements.h"

/**************************************************************/

Finley_RefElementInfo Finley_RefElement_InfoList[]={
  {Point1,     "Point1", 0,     0,  1,  1, 1, 1, 
    Point1, {0},
    Finley_Shape_Point1, Finley_Quad_getNodesPoint,       Finley_Quad_getNumNodesPoint,
     1,{0},
     1,{0}, {0},{-1}},
  {Line2,      "Line2", 1,      1,  2,  2, 1, 2, 
     Line2, {0,1},
     Finley_Shape_Line2,  Finley_Quad_getNodesLine,        Finley_Quad_getNumNodesLine,
     2,{0,1},
     2,{0,1}, {1,0},{-1}},
  {Line3, "Line3", 1,      1,  3,  3, 2, 2, 
     Line2, {0,1},
     Finley_Shape_Line3,  Finley_Quad_getNodesLine,        Finley_Quad_getNumNodesLine,
     3,{0,1,2},
     3,{0,1,2}, {1,0,2},{-1}},
  {Line4, "Line4", 1,      1,  4,  4, 3, 2, 
     Line2, {0,1},
     Finley_Shape_Line4,  Finley_Quad_getNodesLine,        Finley_Quad_getNumNodesLine,
     4,{0,1,2,3},
     4,{0,1,2,3}, {1,0,3,2},{-1}},
  {Tri3,       "Tri3", 2,       2,  3,  3, 1, 3, 
     Tri3, {0,1,2},
     Finley_Shape_Tri3,   Finley_Quad_getNodesTri,         Finley_Quad_getNumNodesTri,
     3,{0,1,2},
     3,{0,1,2}, {1,2,0},{0,2,1}},
  {Tri6,       "Tri6", 2,      2,  6,  6, 2, 3, 
     Tri3, {0,1,2},
     Finley_Shape_Tri6,   Finley_Quad_getNodesTri,         Finley_Quad_getNumNodesTri,
     6,{0,1,2,3,4,5}, 
     6,{0,1,2,3,4,5}, {1,2,0,4,5,3},{0,2,1,5,4,3}},
  {Tri9,       "Tri9", 2,      2,  9,  9, 3, 3, 
     Tri3, {0,1,2},
     Finley_Shape_Tri9,   Finley_Quad_getNodesTri,         Finley_Quad_getNumNodesTri,
     9,{0,1,2,3,4,5,6,7,8}, 
     9,{0,1,2,3,4,5,6,7,8}, {1,2,0,5,6,7,8,3,4},{0,2,1,8,7,6,5,4,3}},
  {Tri10, "Tri10", 2,      2, 10, 10, 3, 3, 
     Tri3, {0,1,2},
     Finley_Shape_Tri10,  Finley_Quad_getNodesTri,         Finley_Quad_getNumNodesTri,
     10,{0,1,2,3,4,5,6,7,8,9}, 
     10,{0,1,2,3,4,5,6,7,8,9}, {1,2,0,5,6,7,8,3,4,9},{0,2,1,8,7,6,5,4,3,9}},
  {Rec4, "Rec4", 2,       2,  4,  4, 1, 4, 
     Rec4, {0,1,2,3},
     Finley_Shape_Rec4,   Finley_Quad_getNodesRec,         Finley_Quad_getNumNodesRec,
     4,{0,1,2,3},
     4,{0,1,2,3}, {1,2,3,0},{0,3,2,1}},
  {Rec8,       "Rec8", 2,       2,  8,  8, 2, 4, 
     Rec4, {0,1,2,3},
     Finley_Shape_Rec8,   Finley_Quad_getNodesRec,         Finley_Quad_getNumNodesRec,
     8,{0,1,2,3,4,5,6,7},
     8,{0,1,2,3,4,5,6,7}, {1,2,3,0,5,6,7,4},{0,3,2,1,7,6,5,4}},
  {Rec9,       "Rec9", 2,       2,  9,  9, 2, 4, 
     Rec4, {0,1,2,3},
     Finley_Shape_Rec9,   Finley_Quad_getNodesRec,         Finley_Quad_getNumNodesRec, 
     9,{0,1,2,3,4,5,6,7,8}, 
     9,{0,1,2,3,4,5,6,7,8}, {1,2,3,0,5,6,7,4,8},{0,3,2,1,7,6,5,4,8}},
  {Rec12,      "Rec12", 2,      2, 12, 12, 3, 4, 
     Rec4, {0,1,2,3},
     Finley_Shape_Rec12,  Finley_Quad_getNodesRec,         Finley_Quad_getNumNodesRec,
     12,{0,1,2,3,4,5,6,7,8,9,10,11},
     12,{0,1,2,3,4,5,6,7,8,9,10,11}, {1,2,3,0,6,7,8,9,10,11,4,5},{0,3,2,1,11,10,9,8,7,6,5,4}},
  {Rec16,      "Rec16", 2,      2, 16, 16, 3, 4, 
     Rec4, {0,1,2,3},
     Finley_Shape_Rec16,  Finley_Quad_getNodesRec,         Finley_Quad_getNumNodesRec,
     16,{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15},
     16,{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}, {1,2,3,0,6,7,8,9,10,11,4,5,13,14,15,12},{0,3,2,1,11,10,9,8,7,6,5,4,12,15,14,13}},
  {Tet4,       "Tet4", 3,       3,  4,  4, 1, 4, 
     Tet4, {0,1,2,3},
     Finley_Shape_Tet4,   Finley_Quad_getNodesTet,         Finley_Quad_getNumNodesTet,
     4,{0,1,2,3},
     -1,{999}, {999},{999}},
  {Tet10,      "Tet10", 3,      3, 10, 10, 2, 4, 
     Tet4, {0,1,2,3},
     Finley_Shape_Tet10,  Finley_Quad_getNodesTet,         Finley_Quad_getNumNodesTet,
     10,{0,1,2,3,4,5,6,7,8,9},
     -1,{999}, {999},{999}},
  {Tet16,      "Tet16", 3,      3, 16, 16, 3, 4, 
     Tet4, {0,1,2,3},
     Finley_Shape_Tet16,  Finley_Quad_getNodesTet,         Finley_Quad_getNumNodesTet,
     16,{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15},
     -1,{999}, {999},{999}},
  {Hex8,       "Hex8", 3,       3,  8,  8, 1, 8, 
     Hex8, {0,1,2,3,4,5,6,7},
     Finley_Shape_Hex8,   Finley_Quad_getNodesHex,         Finley_Quad_getNumNodesHex,
     8,{0,1,2,3,4,5,6,7},
     -1,{999}, {999},{999}},
  {Hex20,      "Hex20", 3,      3, 20, 20, 2, 8, 
     Hex8, {0,1,2,3,4,5,6,7},
     Finley_Shape_Hex20,  Finley_Quad_getNodesHex,         Finley_Quad_getNumNodesHex,
     20,{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19},
     -1,{999}, {999},{999}},
  {Hex27,      "Hex27", 3,      3, 27, 27, 2, 8, 
     Hex8, {0,1,2,3,4,5,6,7},
     Finley_Shape_Hex27,  Finley_Quad_getNodesHex,         Finley_Quad_getNumNodesHex,
     27,{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26},
     -1,{999}, {999},{999}},
  {Hex32,      "Hex32", 3,      3, 32, 32, 3, 8, 
     Hex8, {0,1,2,3,4,5,6,7},
     Finley_Shape_Hex32,  Finley_Quad_getNodesHex,         Finley_Quad_getNumNodesHex,
     32,{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},
     -1,{999}, {999},{999}},
  {Line2Face,  "Line2Face", 0,  1,  2,  2, 1, 1, 
     Line2Face, {0,1},
     Finley_Shape_Line2,  Finley_Quad_getNodesPointOnFace, Finley_Quad_getNumNodesPoint,
     1,{0}, 
     1,{0}, {0,1,2},{-1}},
  {Line3Face,  "Line3Face", 0,  1,  3,  3, 2, 1, 
     Line2Face, {0,1},
     Finley_Shape_Line3,  Finley_Quad_getNodesPointOnFace, Finley_Quad_getNumNodesPoint,
     1,{0}, 
     1,{0}, {0,1,2},{-1}},
  {Line4Face,  "Line4Face", 0,  1,  4,  4, 3, 1, 
     Line2Face, {0,1},
     Finley_Shape_Line4,  Finley_Quad_getNodesPointOnFace, Finley_Quad_getNumNodesPoint,
     1,{0}, 
     1,{0}, {0,1,2},{-1}},
  {Tri3Face,   "Tri3Face", 1,   2,  3,  3, 1, 2, 
     Tri3Face, {0,1,2},
     Finley_Shape_Tri3,   Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     2,{0,1}, 
     2,{0,1}, {1,0,2},{-1}},
  {Tri6Face,   "Tri6Face", 1,   2,  6,  6, 2, 2, 
     Tri3Face, {0,1,2},
     Finley_Shape_Tri6,   Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     3,{0,1,3}, 
     3,{0,1,3}, {1,0,2,3,5,4},{-1}},
  {Tri9Face,   "Tri9Face", 1,   2,  9,  9, 3, 2, 
     Tri3Face, {0,1,2},
     Finley_Shape_Tri9,   Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     4,{0,1,3,4}, 
     4,{0,1,3,4}, {1,0,2,4,3,7,8,6,5},{-1}},
  {Tri10Face,  "Tri10Face", 1,  2, 10, 10, 3, 2, 
     Tri3Face, {0,1,2},
     Finley_Shape_Tri10,  Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     4,{0,1,3,4},
     4,{0,1,3,4}, {1,0,2,4,3,7,8,6,5,9},{-1}},
  {Rec4Face,   "Rec4Face", 1,   2,  4,  4, 1, 2,
     Rec4Face, {0,1,2,3},
     Finley_Shape_Rec4,   Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     2,{0,1},
     2,{0,1}, {1,0,3,2},{-1}},
  {Rec8Face,   "Rec8Face", 1,   2,  8,  8, 2, 2, 
     Rec4Face, {0,1,2,3},
     Finley_Shape_Rec8,   Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     3,{0,1,4},
     3,{0,1,4}, {1,0,3,2,4,7,6,5},{-1}},
  {Rec9Face,   "Rec9Face", 1,   2,  9,  9, 2, 2, 
     Rec4Face, {0,1,2,3},
     Finley_Shape_Rec9,   Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     3,{0,1,4},
     3,{0,1,4}, {1,0,3,2,4,7,6,5,8},{-1}},
  {Rec12Face,  "Rec12Face", 1,  2, 12, 12, 3, 2, 
     Rec4Face, {0,1,2,3},
     Finley_Shape_Rec12,  Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     4,{0,1,4,5}, 
     4,{0,1,4,5}, {1,0,3,2,5,4,11,10,9,8,7,6},{-1}},
  {Rec16Face,  "Rec16Face", 1,  2, 16, 16, 3, 2, 
     Rec4Face, {0,1,2,3},
     Finley_Shape_Rec16,  Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     4,{0,1,4,5}, 
     4,{0,1,4,5}, {1,0,3,2,5,4,11,10,9,8,7,6,13,12,15,14},{-1}},
  {Tet4Face,   "Tet4Face", 2,   3,  4,  4, 1, 3, 
     Tet4Face, {0,1,2,3},
     Finley_Shape_Tet4,   Finley_Quad_getNodesTriOnFace,   Finley_Quad_getNumNodesTri,
     3,{0,1,2}, 
     3,{0,1,2}, {1,2,0,3},{0,2,1,3}},
  {Tet10Face,  "Tet10Face", 2,  3, 10, 10, 2, 3, 
     Tet4Face, {0,1,2,3},
     Finley_Shape_Tet10,  Finley_Quad_getNodesTriOnFace,   Finley_Quad_getNumNodesTri,
     4,{0,1,2,4,5,6}, 
     4,{0,1,2,4,5,6}, {1,2,0,3,5,6,4,8,9,7},{0,2,1,3,6,7,9,8}},
  {Tet16Face,  "Tet16Face", 2,  3, 16, 16, 3, 3, 
     Tet4Face, {0,1,2,3},
     Finley_Shape_Tet16,  Finley_Quad_getNodesTriOnFace,   Finley_Quad_getNumNodesTri,
     5,{0,1,2,4,5,6,7,8,9}, 
     5,{0,1,2,4,5,6,7,8,9}, {1,2,0,3,6,7,8,9,4,5,11,12,10,14,15,13},{0,2,1,3,9,8,7,6,5,4,9,8,7,6,10,12,11,13,15,14}},
  {Hex8Face,   "Hex8Face", 2,   3,  8,  8, 1, 4, 
     Hex8Face, {0,1,2,3,4,5,6,7},
     Finley_Shape_Hex8,   Finley_Quad_getNodesRecOnFace,   Finley_Quad_getNumNodesRec,
     4,{0,1,2,3}, 
     4,{0,1,2,3}, {1,2,3,0,5,6,7,4},{0,3,2,1,4,7,6,5}},
  {Hex20Face,  "Hex20Face", 2,  3, 20, 20, 2, 4, 
     Hex8Face, {0,1,2,3,4,5,6,7},
     Finley_Shape_Hex20,  Finley_Quad_getNodesRecOnFace,   Finley_Quad_getNumNodesRec,
     8,{0,1,2,3,8,9,10,11},
     8,{0,1,2,3,8,9,10,11}, {1,2,3,0,5,6,7,4,9,10,11,8,13,14,15,12,17,18,19,16},{0,3,2,1,4,7,6,5,11,10,9,8,12,15,14,13,19,18,17,16}},
  {Hex27Face,  "Hex27Face", 2,  3, 27, 27, 2, 4, 
     Hex8Face, {0,1,2,3,4,5,6,7},
     Finley_Shape_Hex27,  Finley_Quad_getNodesRecOnFace,   Finley_Quad_getNumNodesRec,
     9,{0,1,2,3,8,9,10,11,20},
     9,{0,1,2,3,8,9,10,11,20}, {1,2,3,0,5,6,7,4,9,10,11,8,13,14,15,12,17,18,19,16,20,22,23,24,22,25,26},{0,3,2,1,4,7,6,5,11,10,9,8,12,15,14,13,19,18,17,16,20,24,23,22,21,25,26}},
  {Hex32Face,  "Hex32Face", 2,  3, 32, 32, 3, 4, 
     Hex8Face, {0,1,2,3,4,5,6,7},
     Finley_Shape_Hex32,  Finley_Quad_getNodesRecOnFace,   Finley_Quad_getNumNodesRec,
     12,{0,1,2,3,8,9,10,11,12,13,14,15},
     12,{0,1,2,3,8,9,10,11,12,13,14,15},
     {1,2,3,0,5,6,7,4,10,11,12,13,14,15,8,9,17,18,19,16,21,22,23,20,26,27,28,29,30,31,34,25},
                   {0,3,2,1,4,7,6,5,15,14,13,12,11,10,9,8,16,19,18,17,20,23,22,21,31,30,29,28,27,26,25,24}},

  {Point1_Contact, "Point1_Contact", 0, 0,  2,  1, 1, 1,
     Point1_Contact, {0,1},
     Finley_Shape_Point1, Finley_Quad_getNodesPoint, Finley_Quad_getNumNodesPoint,
     1,{0},
     -1,{999}, {999},{999}},
  {Line2_Contact,  "Line2_Contact", 1,  1,  4,  2, 1, 2, 
     Line2_Contact, {0,1,2,3},
     Finley_Shape_Line2,   Finley_Quad_getNodesLine,  Finley_Quad_getNumNodesLine,
     2,{0,1},
     -1,{999}, {999},{999}},
  {Line3_Contact,  "Line3_Contact", 1,  1,  6,  3, 2, 2, 
     Line2_Contact, {0,1,3,4},
     Finley_Shape_Line3,   Finley_Quad_getNodesLine,  Finley_Quad_getNumNodesLine,
     3,{0,1,2},
     -1,{999}, {999},{999}},
  {Line4_Contact,  "Line4_Contact", 1,  1,  8,  4, 3, 2, 
     Line2_Contact, {0,1,4,5},
     Finley_Shape_Line4,   Finley_Quad_getNodesLine,  Finley_Quad_getNumNodesLine,
     4,{0,1,2,3},
     -1,{999}, {999},{999}},
  {Tri3_Contact,   "Tri3_Contact",  2,  2,  6,  3, 1, 3, 
     Tri3_Contact, {0,1,2,3,4,5},
     Finley_Shape_Tri3,     Finley_Quad_getNodesTri,   Finley_Quad_getNumNodesTri,
     3,{0,1,2},
     -1,{999}, {999},{999}},
  {Tri6_Contact,   "Tri6_Contact", 2,   2, 12,  6, 2, 3, 
     Tri3_Contact, {0,1,2,6,7,8},
     Finley_Shape_Tri6,     Finley_Quad_getNodesTri,   Finley_Quad_getNumNodesTri,
     6,{0,1,2,3,4,5},
     -1,{999}, {999},{999}},
  {Tri9_Contact,   "Tri9_Contact", 2,   2, 18,  9, 3, 3, 
     Tri3_Contact, {0,1,2,9,10,11},
     Finley_Shape_Tri9,     Finley_Quad_getNodesTri,   Finley_Quad_getNumNodesTri,
     9,{0,1,2,3,4,5,6,7,8},
     -1,{999}, {999},{999}},
  {Tri10_Contact,  "Tri10_Contact", 2,  2, 20, 10, 3, 3, 
     Tri3_Contact, {0,1,2,10,11,12},
     Finley_Shape_Tri10,    Finley_Quad_getNodesTri,   Finley_Quad_getNumNodesTri,
     10,{0,1,2,3,4,5,6,7,8,9},
     -1,{999}, {999},{999}},
  {Rec4_Contact,   "Rec4_Contact", 2,   2,  8,  4, 1, 4, 
     Rec4_Contact, {0,1,2,3,4,5,6,7},
     Finley_Shape_Rec4,     Finley_Quad_getNodesRec,   Finley_Quad_getNumNodesRec,
     4,{0,1,2,3},
     -1,{999}, {999},{999}},
  {Rec8_Contact,   "Rec8_Contact", 2,   2, 16,  8, 2, 4, 
     Rec4_Contact, {0,1,2,3,8,9,10,11},
     Finley_Shape_Rec8,     Finley_Quad_getNodesRec,   Finley_Quad_getNumNodesRec,
     8,{0,1,2,3,4,5,6,7},
     -1,{999}, {999},{999}},
  {Rec9_Contact,   "Rec9_Contact", 2,   2, 18,  9, 2, 4, 
     Rec4_Contact, {0,1,2,3,9,10,11,12},
     Finley_Shape_Rec9,     Finley_Quad_getNodesRec,   Finley_Quad_getNumNodesRec,
     9,{0,1,2,3,4,5,6,7,8},
     -1,{999}, {999},{999}},
  {Rec12_Contact,  "Rec12_Contact", 2,  2, 24, 12, 3, 4, 
     Rec4_Contact, {0,1,2,3,12,13,14,15},
     Finley_Shape_Rec12,    Finley_Quad_getNodesRec,   Finley_Quad_getNumNodesRec,
     12,{0,1,2,3,4,5,6,7,8,9,10,11},
     -1,{999}, {999},{999}},
  {Rec16_Contact,  "Rec16_Contact", 2,  2, 32, 16, 3, 4, 
     Rec4_Contact, {0,1,2,3,16,17,18,19},
     Finley_Shape_Rec16,    Finley_Quad_getNodesRec,   Finley_Quad_getNumNodesRec,
     16,{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15},
     -1,{999}, {999},{999}},
  {Line2Face_Contact,  "Line2Face_Contact", 0,  1,  4,  2, 1, 1, 
     Line2Face_Contact,  {0,1,2,3},
     Finley_Shape_Line2,  Finley_Quad_getNodesPointOnFace, Finley_Quad_getNumNodesPoint,
     1,{0},
     -1,{999}, {999},{999}},
  {Line3Face_Contact,  "Line3Face_Contact", 0,  1,  6,  3, 2, 1, 
     Line2Face_Contact, {0,1,3,4},
     Finley_Shape_Line3,  Finley_Quad_getNodesPointOnFace, Finley_Quad_getNumNodesPoint,
     1,{0},
     -1,{999}, {999},{999}},
  {Line4Face_Contact,  "Line4Face_Contact", 0,  1,  8,  4, 3, 1, 
     Line2Face_Contact, {0,1,4,5},
     Finley_Shape_Line4,  Finley_Quad_getNodesPointOnFace, Finley_Quad_getNumNodesPoint,
     1,{0},
     -1,{999}, {999},{999}},
  {Tri3Face_Contact,   "Tri3Face_Contact", 1,   2,  6,  3, 1, 2, 
     Tri3Face_Contact, {0,1,2,3,4,5},
     Finley_Shape_Tri3,   Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     2,{0,1},
     -1,{999}, {999},{999}},
  {Tri6Face_Contact,   "Tri6Face_Contact", 1,   2,  12,  6, 2, 2, 
     Tri3Face_Contact, {0,1,2,6,7,8},
     Finley_Shape_Tri6,   Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     3,{0,1,3},
     -1,{999}, {999},{999}},
  {Tri9Face_Contact,   "Tri9Face_Contact", 1,   2,  18,  9, 3, 2, 
     Tri3Face_Contact, {0,1,2,9,10,11},
     Finley_Shape_Tri9,   Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     4,{0,1,3,4},
     -1,{999}, {999},{999}},
  {Tri10Face_Contact,  "Tri10Face_Contact", 1,  2, 20, 10, 3, 2, 
     Tri3Face_Contact, {0,1,2,10,11,12},
     Finley_Shape_Tri10,  Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     4,{0,1,3,4},
     -1,{999}, {999},{999}},
  {Rec4Face_Contact,   "Rec4Face_Contact",  1,  2,  8,  4, 1, 2, 
     Rec4Face_Contact, {0,1,2,3,4,5,6,7},
     Finley_Shape_Rec4,   Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     2,{0,1},
     -1,{999}, {999},{999}},
  {Rec8Face_Contact,   "Rec8Face_Contact", 1,   2, 16,  8, 2, 2, 
     Rec4Face_Contact, {0,1,2,3,8,9,10,11},
     Finley_Shape_Rec8,   Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     3,{0,1,4},
     -1,{999}, {999},{999}},
  {Rec9Face_Contact,   "Rec9Face_Contact", 1,   2, 18,  9, 2, 2, 
     Rec4Face_Contact, {0,1,2,3,9,10,11,12},
     Finley_Shape_Rec9,   Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     3,{0,1,4},
     -1,{999}, {999},{999}},
  {Rec12Face_Contact,  "Rec12Face_Contact", 1,  2, 24, 12, 3, 2, 
     Rec4Face_Contact, {0,1,2,3,12,13,14,15},
     Finley_Shape_Rec12,  Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     4,{0,1,4,5},
     -1,{999}, {999},{999}},
  {Rec16Face_Contact,  "Rec16Face_Contact", 1,  2, 32, 16, 3, 2, 
     Rec4Face_Contact, {0,1,2,3,16,17,18,19},
     Finley_Shape_Rec16,  Finley_Quad_getNodesLineOnFace,  Finley_Quad_getNumNodesLine,
     4,{0,1,4,5},
     -1,{999}, {999},{999}},
  {Tet4Face_Contact,   "Tet4Face_Contact",  2,  3,  8,  4, 1, 3, 
     Tet4Face_Contact, {0,1,2,3,4,5,6,7},
     Finley_Shape_Tet4,   Finley_Quad_getNodesTriOnFace,   Finley_Quad_getNumNodesTri,
     3,{0,1,2},
     -1,{999}, {999},{999}},
  {Tet10Face_Contact,  "Tet10Face_Contact", 2,  3, 20, 10, 2, 3, 
     Tet4Face_Contact, {0,1,2,3,10,11,12,13},
     Finley_Shape_Tet10,  Finley_Quad_getNodesTriOnFace,   Finley_Quad_getNumNodesTri,
     4,{0,1,2,4,5,6},
     -1,{999}, {999},{999}},
  {Tet16Face_Contact,  "Tet16Face_Contact", 2,  3, 32, 16, 3, 3, 
     Tet4Face_Contact, {0,1,2,3,16,17,18,19},
     Finley_Shape_Tet16,  Finley_Quad_getNodesTriOnFace,   Finley_Quad_getNumNodesTri,
     5,{0,1,2,4,5,6,7,8,9},
     -1,{999}, {999},{999}},
  {Hex8Face_Contact,   "Hex8Face_Contact", 2,   3, 16,  8, 1, 4, 
     Hex8Face_Contact, {0,1,2,3,4,5,6,7, 8,9,10,11,12,13,14,15},
     Finley_Shape_Hex8,   Finley_Quad_getNodesRecOnFace,   Finley_Quad_getNumNodesRec,
     4,{0,1,2,3},
     -1,{999}, {999},{999}},
  {Hex20Face_Contact,  "Hex20Face_Contact", 2,  3, 40, 20, 2, 4, 
     Hex8Face_Contact, {0,1,2,3,4,5,6,7, 20,21,22,23,24,25,26,27},
     Finley_Shape_Hex20,  Finley_Quad_getNodesRecOnFace,   Finley_Quad_getNumNodesRec,
     8,{0,1,2,3,8,9,10,11},
     -1,{999}, {999},{999}},
  {Hex27Face_Contact,  "Hex27Face_Contact", 2,  3, 54, 27, 2, 4, 
     Hex8Face_Contact, {0,1,2,3,4,5,6,7, 27,28,29,30,31,32,33,34},
     Finley_Shape_Hex27,  Finley_Quad_getNodesRecOnFace,   Finley_Quad_getNumNodesRec,
     9,{0,1,2,3,8,9,10,11,20},
     -1,{999}, {999},{999}},
  {Hex32Face_Contact,  "Hex32Face_Contact", 2,  3, 64, 32, 3, 4, 
     Hex8Face_Contact, {0,1,2,3,4,5,6,7, 32,33,34,35,36,37,38,39},
     Finley_Shape_Hex32,  Finley_Quad_getNodesRecOnFace,   Finley_Quad_getNumNodesRec,
     12,{0,1,2,3,8,9,10,11,12,13,14,15},
     -1,{999}, {999},{999}},

  {NoType, "noElement", 0, 0,  0,  0, 0, 0, NoType, {999}, Finley_Shape_Point1, Finley_Quad_getNodesPoint, Finley_Quad_getNumNodesPoint,0,{999},-1,{999},{999},{999}}  /* marks end of list */
};

/**************************************************************/

/*   get a quadrature scheme with NumQuadNodes quadrature nodes for the tri  */
/*   as a queezed scheme on a quad [0,1]^2 */

Finley_RefElement* Finley_RefElement_alloc(ElementTypeId id,int numQuadNodes) {
  Finley_RefElement *out=NULL;
  int Ndim, NS;
  
  /*  allocate the Finley_RefElement to be returned: */
  
  out=MEMALLOC(1,Finley_RefElement);
  if (Finley_checkPtr(out)) return NULL;
  out->Type=&(Finley_RefElement_InfoList[id]);
  out->numQuadNodes=numQuadNodes;
  out->QuadNodes=NULL;
  out->QuadWeights=NULL;
  out->S=NULL;
  out->dSdv=NULL;
  
  /*  allocate memory: */
  
  Ndim=Finley_RefElement_InfoList[id].numDim;
  NS=Finley_RefElement_InfoList[id].numShapes;
  out->QuadNodes=MEMALLOC(numQuadNodes*Ndim,double);
  out->QuadWeights=MEMALLOC(numQuadNodes,double);
  out->S=MEMALLOC(NS*numQuadNodes,double);
  out->dSdv=MEMALLOC(NS*Ndim*numQuadNodes,double);
  if ( Finley_checkPtr(out->QuadNodes) || Finley_checkPtr(out->QuadWeights) || Finley_checkPtr(out->S) || Finley_checkPtr(out->dSdv) ) {
         Finley_RefElement_dealloc(out);
         return NULL;
  }
  
  /*  set the quadrature nodes: */
  
  Finley_RefElement_InfoList[id].getQuadNodes(numQuadNodes,out->QuadNodes,out->QuadWeights);
  if (! Finley_noError()) {
         Finley_RefElement_dealloc(out);
         return NULL;
  } 
  
  /*  eval shape functions on quadrature node: */
  
  Finley_RefElement_InfoList[id].getValues(numQuadNodes,out->QuadNodes,out->S,out->dSdv);
  if (! Finley_noError()) {
         Finley_RefElement_dealloc(out);
         return NULL;
  } 
  
  /*  all done: */
  
  #ifdef Finley_TRACE
  printf("reference element %s with %d quadrature nodes allocated.\n",Finley_RefElement_InfoList[id].Name,numQuadNodes);
  #endif
  return out;
}

/**************************************************************/

void Finley_RefElement_dealloc(Finley_RefElement* in) {
  if (in!=NULL) {
     #ifdef Finley_TRACE
     printf("reference element %s is deallocated.\n",in->Type->Name);
     #endif
     MEMFREE(in->QuadNodes);
     MEMFREE(in->QuadWeights);
     MEMFREE(in->S);
     MEMFREE(in->dSdv);
     MEMFREE(in);
  }
}

/**************************************************************/

ElementTypeId Finley_RefElement_getTypeId(char* element_type) {
    int ptr=0;
    ElementTypeId out=NoType;
    while (Finley_RefElement_InfoList[ptr].TypeId!=NoType && out==NoType) {
       if (strcmp(element_type,Finley_RefElement_InfoList[ptr].Name)==0) out=Finley_RefElement_InfoList[ptr].TypeId;
       ptr++;
    }
    return out;
}
