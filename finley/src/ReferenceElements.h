
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
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

#ifndef INC_FINLEY_REFERENCEELEMENTS
#define INC_FINLEY_REFERENCEELEMENTS


/**************************************************************/

#include "Finley.h"
#include "ShapeFunctions.h"
#include "Quadrature.h"

/**************************************************************/

/*     The ids of the allowed reference elements: */

#define MAX_numNodes 64
#define MAX_numSubElements 8
#define MAX_numSides 2

typedef enum {
  Finley_Point1,
  Finley_Line2,
  Finley_Line3,
  Finley_Line4,
  Finley_Tri3,
  Finley_Tri6,
  Finley_Tri9,
  Finley_Tri10,
  Finley_Rec4,
  Finley_Rec8,
  Finley_Rec9,
  Finley_Rec12,
  Finley_Rec16,
  Finley_Tet4,
  Finley_Tet10,
  Finley_Tet16,
  Finley_Hex8,
  Finley_Hex20,
  Finley_Hex27,
  Finley_Hex32,
  Finley_Line2Face,
  Finley_Line3Face,
  Finley_Line4Face,
  Finley_Tri3Face,
  Finley_Tri6Face,
  Finley_Tri9Face,
  Finley_Tri10Face,
  Finley_Rec4Face,
  Finley_Rec8Face,
  Finley_Rec9Face,
  Finley_Rec12Face,
  Finley_Rec16Face,
  Finley_Tet4Face, 
  Finley_Tet10Face, 
  Finley_Tet16Face,
  Finley_Hex8Face,
  Finley_Hex20Face, 
  Finley_Hex27Face, 
  Finley_Hex32Face, 
  Finley_Point1_Contact,
  Finley_Line2_Contact,
  Finley_Line3_Contact,
  Finley_Line4_Contact,
  Finley_Tri3_Contact,
  Finley_Tri6_Contact,
  Finley_Tri9_Contact,
  Finley_Tri10_Contact,
  Finley_Rec4_Contact,
  Finley_Rec8_Contact,
  Finley_Rec9_Contact,
  Finley_Rec12_Contact,
  Finley_Rec16_Contact,
  Finley_Line2Face_Contact,
  Finley_Line3Face_Contact,
  Finley_Line4Face_Contact,
  Finley_Tri3Face_Contact,
  Finley_Tri6Face_Contact,
  Finley_Tri9Face_Contact,
  Finley_Tri10Face_Contact,
  Finley_Rec4Face_Contact,
  Finley_Rec8Face_Contact,
  Finley_Rec9Face_Contact,
  Finley_Rec12Face_Contact,
  Finley_Rec16Face_Contact,
  Finley_Tet4Face_Contact, 
  Finley_Tet10Face_Contact, 
  Finley_Tet16Face_Contact,
  Finley_Hex8Face_Contact,
  Finley_Hex20Face_Contact,
  Finley_Hex27Face_Contact, 
  Finley_Hex32Face_Contact, 
  Finley_Line3Macro, 
  Finley_Tri6Macro, 
  Finley_Rec9Macro,
  Finley_Tet10Macro,
  Finley_Hex27Macro,
  Finley_NoRef   /* marks end of list */
} Finley_ElementTypeId;

/**************************************************************/

/*  this struct holds the definition of the reference element: */

typedef struct Finley_ReferenceElementInfo {
  Finley_ElementTypeId TypeId;               /* the id */
  char* Name;                                /* the name in text form e.g. Line1,Rec12,... */
  dim_t numNodes;                            /* number of nodes defining the element*/
  dim_t numSubElements;                      /* number of subelements. >1 if macro elements are used. */
  dim_t numSides;							 /* specifies the number of sides the element supports. This =2 if contact elements are used
                                                otherwise =1. */
                                                

  index_t offsets[MAX_numSides+1];			 /* offset to the side nodes: offsets[s]...offset[s+1]-1 refers to the nodes to be used for side s*/								

  
  Finley_ElementTypeId LinearTypeId;         /* id of the linear version of the element */
  
  index_t linearNodes[MAX_numNodes*MAX_numSides];  /* gives the list of nodes defining the linear or macro element */
  
  Finley_QuadTypeId Quadrature;                /* quadrature scheme */
  Finley_ShapeFunctionTypeId Parametrization;  /* shape function for parametrization of the element */
  Finley_ShapeFunctionTypeId BasisFunctions;   /* shape function for the basis functions */ 

  index_t subElementNodes[MAX_numNodes*MAX_numSides*MAX_numSubElements];         /* gives the list of nodes defining the subelements:
																		subElementNodes[INDEX2(i,s,BasisFunctions->numShape*numSides)] is the i-th node in the s-th subelement.*/ 
/*********************************************************************************************************************************** */  
  dim_t numRelevantGeoNodes;                 /* number of nodes used to describe the geometry of the geometrically relevant part of the element
                                                typically this is numNodes but for 'Face' elements where the quadrature points are defined on face of the element 
						this is the number of nodes on the particular face. */
  index_t relevantGeoNodes[MAX_numNodes];    /* list to gather the geometrically relevant nodes (length used is numRelevantGeoNodes)
                                                this list is used for the VTK interface */
  
  dim_t numNodesOnFace;                       /* if the element is allowed as a face element, numNodesOnFace defines the number of nodes defining the face */
                                              /* the following lists are only used for face elements defined by numNodesOnFace>0 */
  index_t faceNodes[MAX_numNodes];             /* list of the nodes defining the face */
  index_t shiftNodes[MAX_numNodes];           /* defines a permutation of the nodes which rotates the nodes on the face */
  index_t reverseNodes[MAX_numNodes];         /* reverses the order of the nodes on a face. The permutation has to keep 0 fixed. */
                                              /* shiftNodes={-1} or reverseNodes={-1} are ignored. */
}  Finley_ReferenceElementInfo;


/**************************************************************/

/*  this struct holds the realization of a reference element */

typedef struct Finley_ReferenceElement {
	Finley_ReferenceElementInfo* Type;     /* type of the reference element */
	Finley_ReferenceElementInfo* LinearType;     /* type of the linear reference element */
	index_t reference_counter;	       /* reference counter */
        dim_t integrationOrder;                /* used integration order */
	dim_t numNodes;
        dim_t numLocalDim;
	dim_t numLinearNodes;
	Finley_ShapeFunction* Parametrization;
	Finley_ShapeFunction* BasisFunctions;
	Finley_ShapeFunction* LinearBasisFunctions;
        double* DBasisFunctionDv;                              /* pointer to derivatives to basis function corresponding to the Parametrization quad points */
        bool_t DBasisFunctionDvShared;                /* TRUE to indicate that DBasisFunctionDv is shared with another object which is managing it */

}  Finley_ReferenceElement;

/**************************************************************/

/*    interfaces: */

Finley_ReferenceElement* Finley_ReferenceElement_alloc(Finley_ElementTypeId,int);
void Finley_ReferenceElement_dealloc(Finley_ReferenceElement*);
Finley_ElementTypeId Finley_ReferenceElement_getTypeId(char*);
Finley_ReferenceElement* Finley_ReferenceElement_reference(Finley_ReferenceElement* in);
Finley_ReferenceElementInfo* Finley_ReferenceElement_getInfo(Finley_ElementTypeId id);


#define Finley_ReferenceElement_getNumNodes(__in__) (__in__)->Type->numNodes

#endif /* #ifndef INC_FINLEY_REFERENCEELEMENTS */

