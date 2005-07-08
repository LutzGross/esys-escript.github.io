/* $Id$ */

#ifndef INC_FINLEY_REFERENCEELEMENTS
#define INC_FINLEY_REFERENCEELEMENTS

/**************************************************************/

/*   Finley: Reference elements */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"
#include "ShapeFunctions.h"
#include "Quadrature.h"

/**************************************************************/

/*     The ids of the allowed reference ellements: */

#define MAX_numNodes 64

typedef enum {
  Point1,
  Line2,
  Line3,
  Line4,
  Tri3,
  Tri6,
  Tri9,
  Tri10,
  Rec4,
  Rec8,
  Rec9,
  Rec12,
  Rec16,
  Tet4,
  Tet10,
  Tet16,
  Hex8,
  Hex20,
  Hex32,
  Line2Face,
  Line3Face,
  Line4Face,
  Tri3Face,
  Tri6Face,
  Tri9Face,
  Tri10Face,
  Rec4Face,
  Rec8Face,
  Rec9Face,
  Rec12Face,
  Rec16Face,
  Tet4Face, 
  Tet10Face, 
  Tet16Face,
  Hex8Face,
  Hex20Face, 
  Hex32Face, 
  Point1_Contact,
  Line2_Contact,
  Line3_Contact,
  Line4_Contact,
  Tri3_Contact,
  Tri6_Contact,
  Tri9_Contact,
  Tri10_Contact,
  Rec4_Contact,
  Rec8_Contact,
  Rec9_Contact,
  Rec12_Contact,
  Rec16_Contact,
  Line2Face_Contact,
  Line3Face_Contact,
  Line4Face_Contact,
  Tri3Face_Contact,
  Tri6Face_Contact,
  Tri9Face_Contact,
  Tri10Face_Contact,
  Rec4Face_Contact,
  Rec8Face_Contact,
  Rec9Face_Contact,
  Rec12Face_Contact,
  Rec16Face_Contact,
  Tet4Face_Contact, 
  Tet10Face_Contact, 
  Tet16Face_Contact,
  Hex8Face_Contact,
  Hex20Face_Contact, 
  Hex32Face_Contact, 
  NoType   /* marks end of list */
} ElementTypeId;

/**************************************************************/

/*  this struct holds the definition of the reference element: */

typedef struct Finley_RefElementInfo {
  ElementTypeId TypeId;                      /* the id */
  char* Name;                                /* the name in text form e.g. Line1,Rec12,... */
  dim_t numDim;                              /* dimension of the element */
  dim_t numNodes;                            /* number of nodes defining the element*/
  dim_t numShapes;                           /* number of shape functions, typically = numNodes*/
  dim_t numOrder;                            /* order of the shape functions */
  dim_t numVertices;                         /* number of vertices of the element */
  ElementTypeId LinearTypeId;                /* id of the linear version of the element */
  index_t linearNodes[MAX_numNodes];         /* gives the list of nodes defining the linear element, typically it is linearNodes[i]=i */
  Finley_Shape_Function* getValues;          /* function to evaluate the shape functions at a set of points */
  Finley_Quad_getNodes* getQuadNodes;        /* function to set the quadrature points */
  Finley_Quad_getNumNodes* getNumQuadNodes;  /* function selects the number of quadrature nodes for a given accuracy order */
  dim_t numGeoNodes;                         /* nuber of nodes used to describe the geometry of the geometrically relevant part of the element */
                                             /* typically this is numNodes but for volumenic elements used to descrbe faces this is the number of */
                                             /* nodes on the particular face */
  index_t geoNodes[MAX_numNodes];            /* list to gather the geometrically relevant nodes */
  dim_t numNodesOnFace;                      /* if the element is allowed as a face element, numNodesOnFace defines the number of nodes */
                                             /* defining the face */
                /* the following lists are only used for face elements defined by numNodesOnFace>0 */
  index_t faceNode[MAX_numNodes];             /* list of the nodes defining the face */
  index_t shiftNodes[MAX_numNodes];           /* defines a permutation of the nodes which rotates the nodes on the face */
  index_t reverseNodes[MAX_numNodes];         /* reverses the order of the nodes on a face. teh permutation has keep 0 fixed. */
                                              /* shiftNodes={-1} or reverseNodes={-1} are ignored. */
}  Finley_RefElementInfo;

/**************************************************************/

/*  this struct holds the realization of a reference element */

typedef struct Finley_RefElement {
  Finley_RefElementInfo* Type;     /* type of the reference element */
  int numQuadNodes;                /* number of quadrature points */
  double *QuadNodes;               /* coordinates of quadrature nodes */
  double *QuadWeights;             /* weights of the quadrature scheme */
  double *S;                       /* shape functions at quadrature nodes */
  double *dSdv;                    /* derivative of the shape functions at quadrature nodes */
}  Finley_RefElement;

/**************************************************************/

/*    interfaces: */

Finley_RefElement* Finley_RefElement_alloc(ElementTypeId,int);
void Finley_RefElement_dealloc(Finley_RefElement*);
ElementTypeId Finley_RefElement_getTypeId(char*);

#endif /* #ifndef INC_FINLEY_REFERENCEELEMENTS */

/*
 * $Log$
 * Revision 1.2  2005/07/08 04:07:56  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.1  2005/06/29 02:34:55  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
 */
