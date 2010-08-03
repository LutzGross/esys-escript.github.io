
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

/*   Finley: Shape functions header file */

/**************************************************************/

#ifndef INC_DUDLEY_SHAPEFUNCTIONS
#define INC_DUDLEY_SHAPEFUNCTIONS

/**************************************************************/

#include "Dudley.h"

/**************************************************************/

#define S_INDEX(_J_,_I_,_NUMNODES_)                INDEX2(_J_,_I_,_NUMNODES_)
#define DSDV_INDEX(_J_,_K_,_I_,_NUMNODES_,_DIM_)   INDEX3(_J_,_K_,_I_,_NUMNODES_,_DIM_)


typedef enum {
  Point1Shape,
  Line2Shape,
  Line3Shape,
  Line4Shape,
  Tri3Shape,
  Tri6Shape,
  Tri9Shape,
  Tri10Shape,
  Rec4Shape,
  Rec8Shape,
  Rec9Shape,
  Rec12Shape,
  Rec16Shape,
  Tet4Shape,
  Tet10Shape,
  Tet16Shape,
  Hex8Shape,
  Hex20Shape,
  Hex27Shape,
  Hex32Shape,
  NoShape   /* marks end of list */
} Finley_ShapeFunctionTypeId;

/**************************************************************/

/*  this struct holds the definition of the shape functions on element: */

typedef void (Finley_ShapeFunction_Evaluation) (dim_t,double*,double*,double*);

typedef struct Finley_ShapeFunctionInfo {
 
  Finley_ShapeFunctionTypeId TypeId;                        /* the id */
  char* Name;                                /* the name in text form e.g. Line1,Rec12,... */
  dim_t numDim;                              /* spacial dimension */
  dim_t numShapes;                           /* number of shape functions */
  dim_t numOrder;                            /* order of the shape functions */
  dim_t numVertices;                         /* number of vertices of the element */
  Finley_ShapeFunction_Evaluation* getValues;           /* function to evaluate the shape functions at a set of points */
}  Finley_ShapeFunctionInfo;


/**************************************************************/

/*  this struct holds the evaluation of a shape function on a quadrature scheme: */

typedef struct Finley_ShapeFunction {
  Finley_ShapeFunctionInfo* Type;     /* type of the reference element */
  int numQuadNodes;                /* number of quadrature points */
  double *QuadNodes;               /* coordinates of quadrature nodes */
  double *QuadWeights;             /* weights of the quadrature scheme */
  double *S;                       /* shape functions at quadrature nodes */
  double *dSdv;                    /* derivative of the shape functions at quadrature nodes */
  index_t reference_counter;	   /* reference counter */
}  Finley_ShapeFunction;

/**************************************************************/
/*   Interfaces: */

Finley_ShapeFunction_Evaluation Finley_Shape_Point1;
Finley_ShapeFunction_Evaluation Finley_Shape_Line2;
Finley_ShapeFunction_Evaluation Finley_Shape_Line3;
Finley_ShapeFunction_Evaluation Finley_Shape_Line4;
Finley_ShapeFunction_Evaluation Finley_Shape_Tri3;
Finley_ShapeFunction_Evaluation Finley_Shape_Tri6;
Finley_ShapeFunction_Evaluation Finley_Shape_Tri9;
Finley_ShapeFunction_Evaluation Finley_Shape_Tri10;
Finley_ShapeFunction_Evaluation Finley_Shape_Rec4;
Finley_ShapeFunction_Evaluation Finley_Shape_Rec8;
Finley_ShapeFunction_Evaluation Finley_Shape_Rec9;
Finley_ShapeFunction_Evaluation Finley_Shape_Rec12;
Finley_ShapeFunction_Evaluation Finley_Shape_Rec16;
Finley_ShapeFunction_Evaluation Finley_Shape_Tet4;
Finley_ShapeFunction_Evaluation Finley_Shape_Tet10;
Finley_ShapeFunction_Evaluation Finley_Shape_Tet16;
Finley_ShapeFunction_Evaluation Finley_Shape_Hex8;
Finley_ShapeFunction_Evaluation Finley_Shape_Hex20;
Finley_ShapeFunction_Evaluation Finley_Shape_Hex27;
Finley_ShapeFunction_Evaluation Finley_Shape_Hex32;

Finley_ShapeFunction* Finley_ShapeFunction_alloc(Finley_ShapeFunctionTypeId id,int numQuadDim, int numQuadNodes, double *QuadNodes, double *QuadWeights);
void Finley_ShapeFunction_dealloc(Finley_ShapeFunction*);
Finley_ShapeFunctionTypeId Finley_ShapeFunction_getTypeId(char*);
Finley_ShapeFunction* Finley_ShapeFunction_reference(Finley_ShapeFunction* in);
Finley_ShapeFunctionInfo* Finley_ShapeFunction_getInfo(Finley_ShapeFunctionTypeId id);
#endif /* #ifndef INC_FINLEY_SHAPEFUNCTIONS */
