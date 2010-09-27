
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Dudley: Reference elements */

/**************************************************************/

#ifndef INC_DUDLEY_REFERENCEELEMENTS
#define INC_DUDLEY_REFERENCEELEMENTS


/**************************************************************/

#include "Dudley.h"
#include "ShapeFunctions.h"
#include "Quadrature.h"

/**************************************************************/

/*     The ids of the allowed reference ellements: */

#define MAX_numNodes 64

typedef enum {
  Point1=0,
  Line2=1,
  Tri3=2,
  Tet4=3,
  Line2Face=4,
  Tri3Face=5,
  Tet4Face=6, 
  NoRef=7   /* marks end of list */
} ElementTypeId;

/**************************************************************/

/*  this struct holds the definition of the reference element: */

typedef struct Dudley_ReferenceElementInfo {
  ElementTypeId TypeId;                      /* the id */
  char* Name;                                /* the name in text form e.g. Line1,Rec12,... */
  dim_t numNodes;                            /* number of nodes defining the element*/
  
  Dudley_QuadTypeId Quadrature;                /* quadrature scheme */
  Dudley_ShapeFunctionTypeId BasisFunctions;   /* shape function for the basis functions */ 

/*********************************************************************************************************************************** */  
 
  dim_t numNodesOnFace;                       /* if the element is allowed as a face element, numNodesOnFace defines the number of nodes defining the face */
                                              /* the following lists are only used for face elements defined by numNodesOnFace>0 */

  index_t shiftNodes[MAX_numNodes];           /* defines a permutation of the nodes which rotates the nodes on the face */
  index_t reverseNodes[MAX_numNodes];         /* reverses the order of the nodes on a face. the permutation has keep 0 fixed. */
                                              /* shiftNodes={-1} or reverseNodes={-1} are ignored. */
}  Dudley_ReferenceElementInfo;


/**************************************************************/

/*  this struct holds the realization of a reference element */

typedef struct Dudley_ReferenceElement {
	Dudley_ReferenceElementInfo* Type;     /* type of the reference element  - dudley only supports linear elements*/
	index_t reference_counter;	       /* reference counter */
        dim_t integrationOrder;                /* used integration order */
	dim_t numNodes;
        dim_t numLocalDim;
	Dudley_ShapeFunction* BasisFunctions;
        double* DBasisFunctionDv;                              /* pointer to derivatives to basis function corresponding to the Parametrization quad points */
        bool_t DBasisFunctionDvShared;                /* TRUE to indicate that DBasisFunctionDv is shared with another object which is managing it */

}  Dudley_ReferenceElement;

/**************************************************************/

/*    interfaces: */

Dudley_ReferenceElement* Dudley_ReferenceElement_alloc(ElementTypeId,int);
void Dudley_ReferenceElement_dealloc(Dudley_ReferenceElement*);
ElementTypeId Dudley_ReferenceElement_getTypeId(char*);
Dudley_ReferenceElement* Dudley_ReferenceElement_reference(Dudley_ReferenceElement* in);
Dudley_ReferenceElementInfo* Dudley_ReferenceElement_getInfo(ElementTypeId id);


#define Dudley_ReferenceElement_getNumNodes(__in__) (__in__)->Type->numNodes

#endif /* #ifndef INC_DUDLEY_REFERENCEELEMENTS */
